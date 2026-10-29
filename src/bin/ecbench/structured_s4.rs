//! Product-space S4 relaxation in the native IC pipeline. Every returned
//! relation is lifted and verified by the counted group. Native solver words
//! are retained for the common band; raw S has a legacy wall-priced solver
//! term. Incomplete construction/split accounting prevents a full-cost claim.
use crypto_lib::cryptanalysis::ic_boundary::{
    binary_subspace_factor_base, lift_abscissae, BinaryGroup, BinaryInstance, FactorBase, GroupOps,
    OracleCounters,
};
use crypto_lib::cryptanalysis::ic_framework::stages::{
    BooleanSystem, DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, SolverCost,
    SolverTotals, SystemShape,
};
use crypto_lib::cryptanalysis::inherited_f4::ReducedBasis;
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    f4_word_ops_thread, solve_boolean_system, sym_semaev_s4, FieldStructure, SolveOptions,
    SolverEngine, SplitRule, SymElement,
};
use std::collections::BTreeMap;
use std::time::Instant;

#[allow(dead_code)]
#[path = "structured_s4_math.rs"]
mod math;

pub struct StructuredBase<'i>(&'i BinaryInstance);
impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for StructuredBase<'_> {
    fn name(&self) -> &str {
        "structured-s4-base-v1"
    }
    fn describe(&self, _: &Params) -> String {
        "Explicit seeded binary subspace; prefix, geometric, trace-zero geometric or random".into()
    }
    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("dimension", "basis dimension"),
            ("family", "prefix|geometric|geomtraceu|random"),
            ("seed", "public base seed"),
        ]
    }
    fn build(
        &self,
        _: &InstanceCtx<BinaryGroup<'a>>,
        p: &Params,
        _: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        let ell = p.u64("dimension")? as usize;
        let seed = p.u64("seed")?;
        let family = p.get("family").ok_or("missing family")?;
        if ell == 0 || ell > 10 || ell >= self.0.n as usize || self.0.n > 62 {
            return Err("structured base needs 1<=dimension<=10 and dimension<n<=62".into());
        }
        if !["prefix", "geometric", "geomtraceu", "random"].contains(&family) {
            return Err("unknown structured family".into());
        }
        let started = Instant::now();
        let basis = math::make_basis(&self.0.gf, ell, family, seed);
        let mut fb = binary_subspace_factor_base(self.0, &basis);
        fb.cost.wall_ns = started.elapsed().as_nanos() as u64;
        Ok(fb)
    }
}

pub struct StructuredOracle<'i> {
    instance: &'i BinaryInstance,
    summands: u32,
    form: String,
    degree: u32,
    preclose: bool,
    node_budget: usize,
    max_solutions: usize,
    max_cols: u128,
    bs: Vec<u64>,
    ps: Vec<Vec<u64>>,
    st: Option<FieldStructure>,
    totals: SolverTotals,
}
impl<'i> StructuredOracle<'i> {
    fn new(summands: u32, instance: &'i BinaryInstance) -> Self {
        Self {
            instance,
            summands,
            form: String::new(),
            degree: 0,
            preclose: false,
            node_budget: 0,
            max_solutions: 0,
            max_cols: 0,
            bs: vec![],
            ps: vec![],
            st: None,
            totals: SolverTotals::default(),
        }
    }
    fn system(
        &self,
        xr: u64,
    ) -> (
        Vec<crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolPoly>,
        usize,
    ) {
        let n = self.instance.n;
        let st = self.st.as_ref().unwrap();
        if self.form == "symmetric" {
            let (eqs, ps) = math::symmetric_system(&self.instance.gf, &self.bs, xr, st).unwrap();
            (eqs, ps.iter().map(Vec::len).sum())
        } else {
            let base: Vec<_> = self.bs.iter().map(|&x| math::element(x, n)).collect();
            let nv = 3 * self.bs.len();
            let xs: Vec<_> = (0..3)
                .map(|i| SymElement::from_subspace_vars(&base, i * self.bs.len(), n, nv))
                .collect();
            (
                sym_semaev_s4(&xs[0], &xs[1], &xs[2], &math::element(xr, n), st),
                nv,
            )
        }
    }
}
impl<'a> DecompositionOracle<BinaryGroup<'a>> for StructuredOracle<'_> {
    fn name(&self) -> &str {
        "structured-s4-oracle-v1"
    }
    fn summands(&self) -> u32 {
        self.summands
    }
    fn describe(&self, _: &Params) -> String {
        "S4 direct or product-space symmetric relaxation; counted sign lift; solver words and construction unpriced".into()
    }
    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("m", "must be 3"),
            ("form", "direct|symmetric"),
            (
                "degree",
                "requested hybrid degree; input-degree floor applies",
            ),
            ("preclose", "0|1: fixed-degree full closure before hybrid"),
            ("node_budget", "deterministic node limit"),
            ("max_solutions", "root cap"),
            ("max_cols", "closure column cap"),
        ]
    }
    fn prepare(
        &mut self,
        _: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        p: &Params,
        _: &mut GroupOps,
    ) -> Result<(), String> {
        if std::env::vars_os().any(|(k, _)| {
            ["KIC_", "F4_", "SOLVER_"]
                .iter()
                .any(|s| k.to_string_lossy().starts_with(s))
        }) {
            return Err("structured S4 refuses solver environment overrides".into());
        }
        if self.summands != 3 || self.instance.b != 1 {
            return Err("structured S4 requires m=3 and curve b=1".into());
        }
        self.form = p.get("form").ok_or("missing form")?.into();
        if self.form != "direct" && self.form != "symmetric" {
            return Err("form must be direct or symmetric".into());
        }
        self.degree = p
            .u64("degree")?
            .try_into()
            .map_err(|_| "degree too large")?;
        self.preclose = match p.u64("preclose")? {
            0 => false,
            1 => true,
            _ => return Err("preclose must be 0 or 1".into()),
        };
        self.node_budget = p
            .u64("node_budget")?
            .try_into()
            .map_err(|_| "node_budget too large")?;
        self.max_solutions = p
            .u64("max_solutions")?
            .try_into()
            .map_err(|_| "max_solutions too large")?;
        self.max_cols = p.u64("max_cols")? as u128;
        if !(2..=4).contains(&self.degree) || self.max_solutions == 0 || self.node_budget == 0 {
            return Err("invalid deterministic limits".into());
        }
        self.bs = fb
            .subspace_basis
            .clone()
            .ok_or("factor base must record an exact subspace basis")?;
        self.ps = math::products(&self.instance.gf, &self.bs, 3);
        let nv = if self.form == "symmetric" {
            self.ps.iter().map(Vec::len).sum()
        } else {
            3 * self.bs.len()
        };
        if nv > 21 {
            return Err(format!(
                "pilot root extraction supports at most 21 variables, got {nv}"
            ));
        }
        self.st = Some(FieldStructure::new(
            self.instance.n,
            &self.instance.irreducible,
        ));
        self.totals.solver = "structured-s4-v1: inherited F4; native solver words only".into();
        Ok(())
    }
    fn decompose(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
    ) -> Option<Vec<usize>> {
        if point.infinity {
            return None;
        }
        let started = Instant::now();
        let (mut eqs, nv) = self.system(point.x);
        let shape = BooleanSystem {
            equations: eqs.clone(),
            n_vars: nv,
        }
        .shape();
        let terms = eqs.iter().map(|p| p.terms.len() as u64).sum();
        let mut closure_words = 0;
        let mut refuted = false;
        let mut capped = false;
        if self.preclose {
            let cols = (0..=self.degree)
                .fold((0u128, 1u128), |(sum, c), k| {
                    let next = if k < (nv as u32) {
                        c * (nv as u32 - k) as u128 / (k + 1) as u128
                    } else {
                        0
                    };
                    (sum + c, next)
                })
                .0;
            if cols > self.max_cols {
                capped = true;
            } else {
                match ReducedBasis::from_system_closed_full(&eqs, nv, self.degree, 1_000_000) {
                    Some((mut basis, mut cost)) => {
                        refuted = basis.is_refuted();
                        if !refuted {
                            eqs.extend(basis.decisive_rows(&mut cost));
                        }
                        closure_words = cost.word_ops();
                    }
                    None => capped = true,
                }
            }
        }
        let mut roots = vec![];
        let before = f4_word_ops_thread();
        let mut max_degree = self.degree;
        let mut exhausted = capped;
        if !refuted && !capped {
            let opts = SolveOptions {
                engine: SolverEngine::InheritedF4 {
                    max_degree: self.degree,
                },
                max_solutions: self.max_solutions,
                node_budget: self.node_budget,
                split_rule: SplitRule::HighestFree,
            };
            let (found, stats) = solve_boolean_system(&eqs, nv, &opts);
            max_degree = stats.max_degree_built;
            exhausted = stats.exhausted || found.len() >= self.max_solutions;
            roots = found;
        }
        let words = closure_words + f4_word_ops_thread().wrapping_sub(before);
        let (mut split_tests, mut nonsplit) = (0u64, 0u64);
        let mut relation = None;
        let algebraic_roots = roots.len() as u64;
        for root in roots {
            assert!(eqs.iter().all(|p| p.eval(root) == 0));
            let xs = if self.form == "symmetric" {
                math::split(
                    &self.instance.gf,
                    &fb.abscissa_list,
                    math::decode(&self.ps, root),
                    &mut split_tests,
                )
            } else {
                Some(std::array::from_fn(|i| {
                    self.bs
                        .iter()
                        .enumerate()
                        .filter(|(j, _)| root >> (i * self.bs.len() + j) & 1 != 0)
                        .fold(0, |a, (_, b)| a ^ b)
                }))
            };
            if let Some(xs) = xs {
                if let Some(indices) = lift_abscissae(ctx.group, fb, ops, &xs, point) {
                    relation = Some(indices);
                    break;
                }
                counters.lift_failures += 1;
            } else {
                nonsplit += 1;
            }
        }
        let extra = BTreeMap::from([
            ("symbolic_input_terms".into(), terms),
            ("algebraic_roots".into(), algebraic_roots),
            ("cubic_split_tests".into(), split_tests),
            ("nonsplitting_roots".into(), nonsplit),
            ("closure_refutations".into(), u64::from(refuted)),
            ("incomplete_queries".into(), u64::from(exhausted)),
        ]);
        let cost = SolverCost {
            ops: words,
            op_unit: super::pdp3_koblitz::WORD_XOR_UNIT.into(),
            wall_ns: started.elapsed().as_nanos() as u64,
            peak_bytes: 0,
            degree_reached: Some(max_degree),
            solving_degree: None,
            timed_out: exhausted,
            extra,
        };
        self.totals.absorb(&shape, Some(&cost), exhausted);
        relation
    }
    fn last_system(&self) -> Option<SystemShape> {
        self.totals.shape.clone()
    }
    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}
pub fn base<'i>(
    name: &str,
    instance: &'i BinaryInstance,
) -> Option<Box<dyn FactorBaseBuilder<BinaryGroup<'i>> + 'i>> {
    (name == "structured-s4-base-v1").then(|| Box::new(StructuredBase(instance)) as _)
}
pub fn oracle<'i>(
    name: &str,
    summands: u32,
    instance: &'i BinaryInstance,
) -> Option<Box<dyn DecompositionOracle<BinaryGroup<'i>> + 'i>> {
    (name == "structured-s4-oracle-v1")
        .then(|| Box::new(StructuredOracle::new(summands, instance)) as _)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_lib::cryptanalysis::ic_boundary::{koblitz_instance, CountedGroup};
    #[test]
    fn every_family_and_form_recovers_a_counted_repeated_root_witness() {
        let inst = koblitz_instance(0, 13).unwrap();
        let group = BinaryGroup(&inst.fast);
        let ctx = InstanceCtx {
            group: &group,
            generator: inst.generator,
            target: inst.generator,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: "public repeated-root test".into(),
            field_degree: Some(13),
        };
        for family in ["prefix", "geometric", "geomtraceu", "random"] {
            let params = Params(BTreeMap::from([
                ("dimension".into(), "3".into()),
                ("family".into(), family.into()),
                ("seed".into(), "17".into()),
            ]));
            let mut setup = GroupOps::default();
            let fb = StructuredBase(&inst)
                .build(&ctx, &params, &mut setup)
                .unwrap();
            assert!(
                !fb.points.is_empty(),
                "positive control has no points for {family}"
            );
            let p = fb.points[0];
            let twice = group.double(&mut setup, p);
            let point = group.add(&mut setup, p, twice);
            for form in ["direct", "symmetric"] {
                let params = Params(BTreeMap::from([
                    ("form".into(), form.into()),
                    ("degree".into(), "2".into()),
                    (
                        "preclose".into(),
                        if form == "symmetric" { "1" } else { "0" }.into(),
                    ),
                    ("node_budget".into(), "20000".into()),
                    ("max_solutions".into(), "4096".into()),
                    ("max_cols".into(), "30000".into()),
                ]));
                let mut oracle = StructuredOracle::new(3, &inst);
                oracle.prepare(&ctx, &fb, &params, &mut setup).unwrap();
                let mut ops = GroupOps::default();
                let mut counters = OracleCounters::default();
                let found = oracle
                    .decompose(&ctx, &fb, &mut ops, &mut counters, point)
                    .expect("known public witness must lift");
                assert_eq!(found.len(), 3);
                let mut check = GroupOps::default();
                let mut sum = group.identity();
                for i in found {
                    sum = group.add(&mut check, sum, fb.points[i]);
                }
                assert_eq!(sum, point);
                assert!(ops.adds + ops.doubles > 0);
                assert!(oracle.totals.ops > 0);
                assert_eq!(oracle.totals.calls, 1);
            }
        }
    }
    #[test]
    fn oracle_refuses_incomplete_search_parameters() {
        let inst = koblitz_instance(0, 7).unwrap();
        let group = BinaryGroup(&inst.fast);
        let ctx = InstanceCtx {
            group: &group,
            generator: inst.generator,
            target: inst.generator,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: "parameter control".into(),
            field_degree: Some(7),
        };
        let mut ops = GroupOps::default();
        let fb = binary_subspace_factor_base(&inst, &[1, 2]);
        assert!(StructuredOracle::new(3, &inst)
            .prepare(&ctx, &fb, &Params::default(), &mut ops)
            .is_err());
    }
}
