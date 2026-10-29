// Ported from the tested native pursuit instrument at d67f6085; finite controls in that study.
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_groebner::{FieldStructure, SymElement};
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolPoly, Mask};
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use num_bigint::BigUint;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde_json::{json, Value};
pub(super) fn element(x: u64, n: u32) -> F2mElement {
    F2mElement::from_biguint(&BigUint::from(x), n)
}
pub(super) fn word(x: &F2mElement) -> u64 {
    x.raw_bits().first().copied().unwrap_or(0)
}
pub(super) fn point_json(p: &BinaryPoint) -> Value {
    match p {
        BinaryPoint::Infinity => json!({"infinity":true}),
        BinaryPoint::Affine { x, y } => json!({"x":word(x),"y":word(y)}),
    }
}

/// Deterministic binary echelon basis. The span, rather than its ordering, matters.
pub(super) fn basis(xs: impl IntoIterator<Item = u64>) -> Vec<u64> {
    let mut piv = [0u64; 64];
    for mut x in xs {
        while x != 0 {
            let i = 63 - x.leading_zeros() as usize;
            if piv[i] == 0 {
                piv[i] = x;
                break;
            }
            x ^= piv[i];
        }
    }
    piv.into_iter().filter(|&x| x != 0).collect()
}
pub(super) fn span(bs: &[u64]) -> Vec<u64> {
    let mut xs = vec![0];
    for &b in bs {
        let len = xs.len();
        for i in 0..len {
            xs.push(xs[i] ^ b);
        }
    }
    xs.sort_unstable();
    xs
}
pub(super) fn trace(gf: &Gf2, x: u64) -> u64 {
    let (mut y, mut t) = (x, 0);
    for _ in 0..gf.n {
        t ^= y;
        y = gf.sqr(y);
    }
    assert!(t <= 1);
    t
}
pub(super) fn kernel(rows: &[u64], n: u32) -> Vec<u64> {
    let mut rows = rows.to_vec();
    let mut pivots = Vec::new();
    let mut r = 0;
    for col in 0..n as usize {
        let Some(p) = (r..rows.len()).find(|&i| rows[i] >> col & 1 != 0) else {
            continue;
        };
        rows.swap(r, p);
        for i in 0..rows.len() {
            if i != r && rows[i] >> col & 1 != 0 {
                rows[i] ^= rows[r];
            }
        }
        pivots.push(col);
        r += 1;
    }
    (0..n as usize)
        .filter(|j| !pivots.contains(j))
        .map(|j| {
            let mut v = 1u64 << j;
            for (i, &p) in pivots.iter().enumerate() {
                if rows[i] >> j & 1 != 0 {
                    v |= 1u64 << p;
                }
            }
            v
        })
        .collect()
}
pub(super) fn make_basis(gf: &Gf2, ell: usize, family: &str, seed: u64) -> Vec<u64> {
    assert!(ell > 0 && ell < gf.n as usize && gf.n < 63);
    let mut rng = StdRng::seed_from_u64(seed);
    if family == "prefix" {
        return (0..ell).map(|j| 1 << j).collect();
    }
    if family == "random" {
        let mut bs = Vec::new();
        while bs.len() < ell {
            bs = basis(bs.into_iter().chain([rng.gen_range(1..1u64 << gf.n)]));
        }
        return bs;
    }
    assert!(
        family == "geometric" || family == "geomtraceu",
        "unknown family"
    );
    loop {
        let g = rng.gen_range(2..1u64 << gf.n);
        let mut powers = vec![1];
        for j in 1..ell {
            powers.push(gf.mul(powers[j - 1], g));
        }
        if basis(powers.iter().copied()).len() != ell {
            continue;
        }
        let c = if family == "geometric" {
            rng.gen_range(1..1u64 << gf.n)
        } else {
            let rows: Vec<u64> = powers
                .iter()
                .map(|&p| (0..gf.n).fold(0, |a, j| a | (trace(gf, gf.mul(p, 1 << j)) << j)))
                .collect();
            let ker = kernel(&rows, gf.n);
            assert!(!ker.is_empty());
            let mut c = 0;
            while c == 0 {
                // XOR, not integer sum: uniform nonzero vector of the kernel.
                for &v in &ker {
                    if rng.gen::<bool>() {
                        c ^= v;
                    }
                }
            }
            c
        };
        let bs: Vec<_> = powers.into_iter().map(|p| gf.mul(c, p)).collect();
        assert_eq!(basis(bs.iter().copied()).len(), ell);
        if family == "geomtraceu" {
            assert!(bs.iter().all(|&b| trace(gf, b) == 0));
        }
        return bs;
    }
}
pub(super) fn products(gf: &Gf2, bs: &[u64], m: usize) -> Vec<Vec<u64>> {
    let mut out = vec![basis(bs.iter().copied())];
    for k in 1..m {
        out.push(basis(
            out[k - 1]
                .iter()
                .flat_map(|&x| bs.iter().map(move |&y| gf.mul(x, y))),
        ));
    }
    out
}

/// FPPR-style symmetric relaxation on the actual product spaces V^(k).
pub(super) fn symmetric_system(
    gf: &Gf2,
    bs: &[u64],
    xr: u64,
    st: &FieldStructure,
) -> Option<(Vec<F2BoolPoly>, Vec<Vec<u64>>)> {
    let ps = products(gf, bs, 3);
    let nv = ps.iter().map(Vec::len).sum::<usize>();
    if nv > 64 {
        return None;
    }
    let mut offset = 0;
    let es: Vec<_> = ps
        .iter()
        .map(|p| {
            let xs: Vec<_> = p.iter().map(|&x| element(x, gf.n)).collect();
            let e = SymElement::from_subspace_vars(&xs, offset, gf.n, nv);
            offset += p.len();
            e
        })
        .collect();
    let (e1, e2, e3) = (&es[0], &es[1], &es[2]);
    let (e12, e22, e32) = (e1.square(st), e2.square(st), e3.square(st));
    let r1 = SymElement::constant(&element(xr, gf.n), gf.n, nv);
    let r2 = r1.square(st);
    let r3 = r2.mul(&r1, st);
    let r4 = r2.square(st);
    let out = r4
        .add(&e12.square(st))
        .add(&e32.square(st))
        .add(&e22.square(st).mul(&r4, st))
        .add(&e32.mul(e3, st).mul(&r1, st))
        .add(&e3.mul(&e22, st).mul(&r3, st))
        .add(&e3.mul(&e12, st).mul(&r1, st))
        .add(&e3.mul(&r3, st))
        .add(&e12.mul(&e32, st).mul(&r2, st))
        .add(&e32.mul(&r4, st))
        .add(&e32)
        .add(&e22.mul(&r2, st));
    assert!(out
        .coords
        .iter()
        .flat_map(|p| &p.terms)
        .all(|t| t.mask.count_ones() <= 2));
    Some((out.coords, ps))
}
pub(super) fn symmetric_values(gf: &Gf2, xs: [u64; 3]) -> [u64; 3] {
    [
        xs[0] ^ xs[1] ^ xs[2],
        gf.mul(xs[0], xs[1]) ^ gf.mul(xs[0], xs[2]) ^ gf.mul(xs[1], xs[2]),
        gf.mul(gf.mul(xs[0], xs[1]), xs[2]),
    ]
}
pub(super) fn decode(ps: &[Vec<u64>], assignment: Mask) -> [u64; 3] {
    let mut offset = 0;
    std::array::from_fn(|i| {
        let mut x = 0;
        for (j, &b) in ps[i].iter().enumerate() {
            if assignment >> (offset + j) & 1 != 0 {
                x ^= b;
            }
        }
        offset += ps[i].len();
        x
    })
}
pub(super) fn encode(ps: &[Vec<u64>], es: [u64; 3]) -> u64 {
    let mut offset = 0;
    let mut assignment = 0;
    for i in 0..3 {
        let j = span(&ps[i]).into_iter().position(|x| x == es[i]);
        assert!(j.is_some()); // membership only; recover coefficients separately.
        let mut found = None;
        for bits in 0..1u64 << ps[i].len() {
            let x = ps[i]
                .iter()
                .enumerate()
                .filter(|(j, _)| bits >> j & 1 != 0)
                .fold(0, |a, (_, b)| a ^ b);
            if x == es[i] {
                found = Some(bits);
                break;
            }
        }
        assignment |= found.expect("product-space membership") << offset;
        offset += ps[i].len();
    }
    assignment
}
/// Exact synthetic division, including repeated roots. No distinct-root shortcut.
pub(super) fn split(gf: &Gf2, domain: &[u64], es: [u64; 3], tests: &mut u64) -> Option<[u64; 3]> {
    let mut coeff = vec![1, es[0], es[1], es[2]];
    let mut roots = [0; 3];
    for dst in &mut roots {
        let root = domain.iter().copied().find(|&x| {
            *tests += 1;
            coeff.iter().fold(0, |a, &c| gf.mul(a, x) ^ c) == 0
        })?;
        *dst = root;
        let mut q = vec![coeff[0]];
        for i in 1..coeff.len() - 1 {
            q.push(coeff[i] ^ gf.mul(root, q[i - 1]));
        }
        assert_eq!(coeff.last().unwrap() ^ gf.mul(root, *q.last().unwrap()), 0);
        coeff = q;
    }
    roots.sort_unstable();
    assert_eq!(symmetric_values(gf, roots), es);
    Some(roots)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode_preserves_the_top_assignment_bit() {
        let bits = Mask::BITS as usize;
        let ps = vec![vec![0; bits - 1], vec![7], vec![]];
        assert_eq!(decode(&ps, (1 as Mask) << (bits - 1)), [0, 7, 0]);
        assert_eq!(decode(&ps, 0), [0, 0, 0]);
    }
}
