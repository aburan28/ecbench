#!/usr/bin/env python3
"""Build the measured ecbench PDF addendum and append it to the field report.

Run from any directory. All numerical tables are read from the release files;
the PDF layout and metadata are deterministic for a given input snapshot.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab import rl_config
from reportlab.graphics.shapes import Drawing, Rect, String
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate, Frame, PageBreak, PageTemplate, Paragraph, Spacer, Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[1]
INK = colors.HexColor("#18334a")
TEAL = colors.HexColor("#207c85")
PALE = colors.HexColor("#edf4f6")
GRAY = colors.HexColor("#52616d")
AMBER = colors.HexColor("#d79527")


def read(path: Path):
    return json.loads(path.read_text())


def register_fonts():
    candidates = [Path("/usr/share/fonts/truetype/dejavu"),
                  Path("/usr/local/share/fonts/dejavu")]
    for font_root in candidates:
        if (font_root / "DejaVuSans.ttf").exists():
            for name, file in [("Body", "DejaVuSans.ttf"),
                               ("BodyBold", "DejaVuSans-Bold.ttf"),
                               ("Mono", "DejaVuSansMono.ttf")]:
                pdfmetrics.registerFont(TTFont(name, str(font_root / file)))
            pdfmetrics.registerFontFamily("Body", normal="Body", bold="BodyBold")
            return
    raise RuntimeError("DejaVu Sans fonts are required to render mathematical text")


def build(addendum: Path, original: Path, combined: Path):
    rl_config.invariant = 1
    register_fonts()
    data = ROOT / "ecbench" / "data"
    work = ROOT / "build"
    manifest = read(data / "manifest.json")
    cases = read(work / "algebra_cases.json")
    coverage_path = work / "algebra_coverage_cases.json"
    coverage = read(coverage_path) if coverage_path.exists() else None
    endo = read(work / "endomorphism_cases.json")
    glv = read(work / "glv4.json")
    results = [read(work / f"algebra_{method}.json")
               for method in ["buchberger", "f5b"]]
    ev = read(work / "endomorphism_verification.json")
    eb = read(work / "endomorphism_benchmark.json")
    core = manifest["counts"]
    ac = cases["counts"]
    total_cases = ac["cases"] + (len(coverage["cases"]) if coverage else 0)
    result_by_id = {r["case_id"]: r for r in results[0]["results"]}
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle("ReportBody", fontName="Body", fontSize=9.8,
                              leading=14, spaceAfter=9, textColor=INK))
    styles.add(ParagraphStyle("ReportSmall", parent=styles["ReportBody"],
                              fontSize=8.3, leading=11.7, spaceAfter=7))
    styles.add(ParagraphStyle("ReportTitle", fontName="BodyBold", fontSize=22,
                              leading=28, textColor=INK, spaceAfter=14))
    styles.add(ParagraphStyle("ReportHeading", fontName="BodyBold", fontSize=12.2,
                              leading=17, textColor=TEAL, spaceBefore=9, spaceAfter=7))
    styles.add(ParagraphStyle("ReportTable", fontName="Body", fontSize=8.1,
                              leading=11.1, textColor=INK, alignment=TA_LEFT))
    styles.add(ParagraphStyle("ReportTableHead", parent=styles["ReportTable"],
                              fontName="BodyBold", textColor=colors.white))
    styles.add(ParagraphStyle("ReportCode", fontName="Mono", fontSize=8.1,
                              leading=12, backColor=PALE, borderPadding=8,
                              spaceBefore=3, spaceAfter=11, textColor=INK))
    story = []

    def p(text, kind="ReportBody"):
        story.append(Paragraph(text, styles[kind]))

    def h(text):
        p(text, "ReportHeading")

    def title(number, text, subtitle):
        p(f"ECBENCH / ADDENDUM {number}", "ReportSmall")
        p(text, "ReportTitle")
        p(subtitle)

    def table(rows, widths, small=False):
        cells = [[Paragraph(escape(str(x)), styles["ReportTableHead" if i == 0
                                                   else "ReportTable"])
                  for x in row] for i, row in enumerate(rows)]
        t = Table(cells, colWidths=widths, repeatRows=1, hAlign="LEFT")
        t.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), INK),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 8),
            ("RIGHTPADDING", (0, 0), (-1, -1), 8),
            ("TOPPADDING", (0, 0), (-1, -1), 6 if small else 8),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 6 if small else 8),
            ("LINEBELOW", (0, 0), (-1, 0), 0.5, INK),
        ]))
        story.extend([t, Spacer(1, 10)])

    # Page A1: the frozen arithmetic corpus and the newly measured research layer.
    title(1, "A golden EC benchmark battery",
          "Updated experimental scope • release 1.0.0 • 12 October 2026")
    p("The original field-and-curve report is extended here with executable "
      "polynomial relation-search fixtures, measured solver results, operational "
      "isogenies, two-term endomorphism decompositions, and a separate genuine "
      "GLV/GLS four-term example. This addendum gives the current artifact counts "
      "and measurement contract; these counts supersede any earlier vector inventory.")
    table([
        ["Layer", "Current content", "Purpose"],
        ["Arithmetic corpus", f"{core['fields']} fields; {core['curves']} curves; "
         f"{core['vectors']:,} known answers", "Arithmetic, encoding, rejection, subgroup membership, maps"],
        ["Polynomial battery", f"{total_cases} fixtures; main battery "
         f"{ac['cases']} ({ac['binary_cases']} binary, {ac['prime_cases']} prime)",
         "Two-factor Semaev relation solving with exact finite domains"],
        ["GLV2 / Frobenius", f"{ev['supported_cases']} supported of {ev['curves_classified']}; "
         f"{ev['scalar_vectors_checked']} scalar vectors", "Map stability, eigenvalues, lattices, signed reconstruction"],
        ["GLV + GLS4 supplement", "82-bit quadratic field; 80-bit subgroup", "Verified four-term lattice with 20-bit component bound"],
    ], [135, 166, 198])
    h("Which curves can be grouped?")
    p(f"Over each fixed field, equal Frobenius trace characterizes an isogeny class; "
      f"the {core['curves']} represented models occupy {core['represented_isogeny_classes']} "
      "such classes. Equal j-invariant groups geometric isomorphism classes and "
      "twists, and is a different partition. An ordinary endomorphism-order "
      "discriminant records further structure inside an isogeny class. None of "
      "these labels determines the performance of a particular polynomial presentation.")
    h("Fields retained from the original report")
    p("Prime fields: 2<super>61</super>−1, 2<super>64</super>−2<super>32</super>+1, "
      "2<super>83</super>−2<super>59</super>−1, and 2<super>83</super>−55. "
      "Binary fields: degrees 53, 79, and 83 with the explicit irreducible "
      "polynomial bases in the dataset. The binary examples are custom small-field "
      "models analogous in shape to NIST binary curves, rather than NIST standardized parameters.")
    if coverage:
        p(f"A separate {len(coverage['cases'])}-fixture coverage extension adds "
          "degree-53 and degree-79 binary cases plus the other three prime fields. "
          f"The combined {total_cases} polynomial fixtures cover all seven core fields; "
          "the main complexity battery retains its original 55 cases and measurements.",
          "ReportSmall")
    p("These small fields are experimental benchmarks. An approximately 80-bit "
      "prime subgroup admits generic attacks of roughly 2<super>40</super> group "
      "operations; this corpus makes no production-security recommendation.", "ReportSmall")
    story.append(PageBreak())

    # Page A2: exact polynomial systems and honest factor-base populations.
    title(2, "The polynomial experiments",
          "Index-calculus relation finding, with fixed two-summand targets")
    p("Each case seeks actual selected-subgroup points P and Q satisfying "
      "P+Q=T. The target is fixed before solving. Its x-coordinate is substituted "
      "into f₃(x₁,x₂,x₃), and all coefficients, variables, monomial order, finite "
      "domains, and preprocessing choices are recorded as canonical sparse input.")
    h("Semaev equations used")
    p("For y²=x³+Ax+B in odd characteristic:", "ReportSmall")
    p("f₃ = (x₁−x₂)² x₃²<br/>"
      "     − 2[(x₁+x₂)(x₁x₂+A)+2B] x₃<br/>"
      "     + (x₁x₂−A)² − 4B(x₁+x₂)", "ReportCode")
    p("For y²+xy=x³+ax²+b in characteristic two:", "ReportSmall")
    p("f₃ = (x₁x₂+x₁x₃+x₂x₃)² + x₁x₂x₃ + b", "ReportCode")
    p("The binary system expands every extension-field coordinate over F₂. "
      "A factor-base subspace of dimension d produces 2d Boolean variables. "
      "Semaev-derived equations are Boolean-normalized, while v²−v=0 field "
      "equations remain explicit. With a fixed target the normalized f₃ "
      "coordinates have degree at most two; subgroup-membership indicators can "
      "have higher degree and must be retained in the input-degree statistics.")
    h(f"Main {ac['cases']}-fixture complexity grid")
    table([
        ["Cases", "Factor bases and target cohorts"],
        ["49 binary", "K83-a0, S83-a1, B83-a0; low monomials d=2–6 and recorded subgroup-seeded linear bases d=2–5"],
        ["6 prime", "p83-S-0; 4, 8, or 16 factor-base x values from small generator multiples; exact split membership polynomials"],
        ["25 planted", "Constructed from eligible factor-base points; exercise recovery, not random-target success probability"],
        ["30 seeded subgroup", "Independently sampled reproducible nonzero subgroup targets; includes 5 empty-factor-base diagnostics"],
    ], [104, 395], small=True)
    p("Low-monomial B83-a0 bases are empty after the selected-subgroup filter "
      "at every tested d=2–6. They remain useful unit-ideal diagnostics. "
      "Subgroup-seeded bases explicitly include eligible generator x values; "
      "their construction is biased and has its own result stratum.", "ReportSmall")
    p("A polynomial root is only a candidate. The runner lifts it to the "
      "recorded rational points, verifies membership in the chosen cyclic subgroup, "
      "and checks the exact signed sum T. It retains ordered signed decompositions "
      "and their swap orbits, allowing repeated factors. Characteristic-two "
      "polynomial coefficients do not turn −P into P in the odd-order subgroup.", "ReportSmall")
    if coverage:
        p("The 14 added coverage fixtures use K53-a1 and K79-a1 with seeded "
          "dimension-2/3 bases, and p61-G-0, p64-G-0, and p83s-G-0 with "
          "four-x-value bases. Each configuration has a planted and a sampled "
          "target. They extend field coverage rather than the complexity grid.", "ReportSmall")
    story.append(PageBreak())

    # Page A3: metric definitions, including the omitted generation interval.
    title(3, "What is measured",
          "Separate yield, presentation degree, solver work, and elapsed time")
    table([
        ["Quantity", "Definition in this battery"],
        ["Target success", "Targets with at least one verified relation / all attempted targets; a lower bound when solves are censored"],
        ["Relations per attempt", "Verified distinct unordered signed decompositions / all attempted targets; may exceed one"],
        ["Final basis degree", "Maximum total degree among polynomials in the computed final Gröbner basis"],
        ["Degree of regularity", "First d where the degree-d component of S/J is zero; J is generated by the highest-degree homogeneous parts of the exact recorded affine inputs"],
        ["Solving degree", "Maximum degree reached during the solver computation; unmeasured (null), because SymPy does not expose it"],
        ["Solve time", "Prebuilt-equation parsing/setup + Gröbner basis + finite-domain root extraction + independent subgroup/point verification"],
    ], [135, 364], small=True)
    h("Timing boundary")
    p("The runs solve prebuilt fixtures. Factor-base and target construction, "
      "Semaev construction, Weil restriction, and golden-truth enumeration occur "
      "before the benchmark. Their generation time was not recorded. "
      "Both fixture_generation_seconds and generation_inclusive_total_pipeline_seconds "
      "are null; the reported solve time is not a generation-inclusive attack time.")
    p("A separate observed subprocess wall interval includes Python/SymPy startup "
      "and the optional homogeneous-ideal regularity diagnostic. Each case has a "
      "2.5-second subprocess budget. A completed solve remains completed if only "
      "its later diagnostic times out. An unfinished solve is censored, with "
      "unknown relation yield, rather than an empty mathematical relation set.")
    h("Why the degree labels matter")
    p("Regularity here is an invariant of the recorded generator presentation. "
      "Membership filters, field equations, and preprocessing affect J. An "
      "affine unit ideal can therefore have final basis degree zero and a "
      "positive homogeneous-input regularity. Neither number measures the "
      "maximum intermediate degree reached by the solver.")
    degree_rows = [["Buchberger example", "Input degree", "Final GB degree", "Regularity"]]
    for cid, label in [
        ("K83-a0-f3-low_monomials-d2-planted", "K83-a0, low d=2, planted"),
        ("S83-a1-f3-low_monomials-d5-planted", "S83-a1, low d=5, planted"),
        ("p83-S-0-f3-fb4-planted", "p83-S-0, base 4, planted"),
        ("p83-S-0-f3-fb4-seeded_subgroup", "p83-S-0, base 4, sampled"),
    ]:
        result = result_by_id[cid]
        basis = result["groebner_basis"]
        degree_rows.append([label, result["input"]["max_total_degree"],
                            str(basis["max_total_degree"]) + (" (unit ideal)" if basis["unit_ideal"] else ""),
                            result["degree_of_regularity"]["value"]])
    table(degree_rows, [242, 75, 102, 80], small=True)
    p("Roots are extracted by evaluating the computed basis over the declared "
      "finite grid. The baseline uses actual Gröbner computation, followed by "
      "an explicit enumeration stage; it makes no fast triangular/FGLM-root-solver claim. "
      "Golden answers are compared only after solver roots have been obtained.", "ReportSmall")
    story.append(PageBreak())

    # Page A4: measured overall and separate target cohorts.
    title(4, "Measured solver baseline",
          "SymPy 1.14.0 • grevlex • 2.5 seconds per case • Linux x86_64, Python 3.12.14")
    rows = [["Method", "Complete / 55", "Solve timeouts", "Diagnostic-only timeouts",
             "Suite wall (s)", "Σ solved (s)"]]
    for d in results:
        s = d["summary"]
        rows.append([d["configuration"]["method"], s["completed"], s["timeouts_censored"],
                     s["diagnostic_timeouts"], f"{d['suite_wall_seconds']:.2f}",
                     f"{s['total_solve_seconds_completed']:.2f}"])
    table(rows, [94, 82, 77, 100, 73, 73], small=True)
    d = Drawing(499, 94)
    for i, result in enumerate(results):
        y = 60 - i * 31
        s = result["summary"]
        d.add(String(1, y + 6, result["configuration"]["method"], fontName="Body", fontSize=9, fillColor=INK))
        x = 91
        scale = 6.0
        d.add(Rect(x, y, s["completed"] * scale, 17, fillColor=TEAL, strokeColor=None))
        d.add(Rect(x + s["completed"] * scale, y, s["timeouts_censored"] * scale,
                   17, fillColor=AMBER, strokeColor=None))
        d.add(String(428, y + 4, f"{s['completed']} / {s['timeouts_censored']}",
                     fontName="Body", fontSize=9, fillColor=INK))
    d.add(Rect(91, 4, 9, 9, fillColor=TEAL, strokeColor=None))
    d.add(String(105, 4, "completed solves", fontName="Body", fontSize=8, fillColor=INK))
    d.add(Rect(247, 4, 9, 9, fillColor=AMBER, strokeColor=None))
    d.add(String(261, 4, "censored solve timeouts", fontName="Body", fontSize=8, fillColor=INK))
    story.extend([d, Spacer(1, 7)])
    h("Yield by target population")
    rows = [["Method / cohort", "Done", "Success / attempts", "Success bounds",
             "Unordered relations", "Relations / attempt"]]
    for result in results:
        for name in ["planted", "seeded_subgroup"]:
            s = result["cohorts"][name]
            label = f"{result['configuration']['method']} / " + ("planted" if name == "planted" else "sampled")
            rows.append([label, f"{s['completed']}/{s['cases']}",
                         f"{s['successful_targets_completed']}/{s['cases']}",
                         f"{100*s['all_attempt_success_fraction']:.1f}–{100*s['all_attempt_success_fraction_upper_bound']:.1f}%",
                         s["unordered_relations_completed"], f"{s['relations_per_attempt']:.3f}"])
    table(rows, [135, 57, 80, 88, 72, 67], small=True)
    p("Success bounds allow all unfinished targets to succeed for the upper "
      "endpoint. Relation counts and relations per attempt are observed lower "
      "bounds when timeouts exist. No finite upper bound on the unknown relation "
      "counts is claimed in this table.", "ReportSmall")
    s0, s1 = [r["summary"] for r in results]
    p(f"Buchberger recovered {s0['ordered_relations_completed']} ordered / "
      f"{s0['unordered_relations_completed']} unordered signed relations; f5b "
      f"recovered {s1['ordered_relations_completed']} / {s1['unordered_relations_completed']}. "
      "All completed cases matched the independently enumerated "
      "golden roots and signed point decompositions, with zero solver errors. "
      "No completed sampled target yielded a relation at these small factor-base "
      "sizes; the unfinished sampled target in each run remains unknown.")
    p("The method comparison describes this fixed fixture set, implementation, "
      "machine, and deadline. It does not establish a general ranking of Gröbner "
      "algorithms. The sampled cohort also includes five empty-base diagnostics; "
      "basis construction and model variant remain separate strata in the JSON "
      "rather than being interpreted as an unbiased cryptanalytic population.", "ReportSmall")
    coverage_runs = [work / f"algebra_coverage_{method}.json"
                     for method in ["buchberger", "f5b"]]
    if coverage and all(path.exists() for path in coverage_runs):
        h("Separate 14-fixture field-coverage batch")
        rows = [["Method", "Done / 14", "Solve timeouts", "Unordered relations", "Σ solved (s)"]]
        for path in coverage_runs:
            result = read(path)
            summary = result["summary"]
            rows.append([result["configuration"]["method"], summary["completed"],
                         summary["timeouts_censored"], summary["unordered_relations_completed"],
                         f"{summary['total_solve_seconds_completed']:.3f}"])
        table(rows, [119, 80, 93, 122, 85], small=True)
        p("Each method recovered all seven planted targets; all seven sampled "
          "targets completed without a relation. Both coverage batches had zero "
          "errors or diagnostic timeouts. These are separate runs of the added "
          "fixtures, preserving the original 55-case timings above.", "ReportSmall")
    story.append(PageBreak())

    # Page A5: operational maps, selected-subgroup stability, and timings.
    title(5, "Isogenies and two-term maps",
          "Explicit known answers are distinguished from class-level metadata")
    table([
        ["Map battery", "Known answers", "Scope"],
        ["Separable degree-2 isogeny", "20", "E: y²=x³−x to E′: y²=x³+4x; kernel {O,(0,0)}"],
        ["Gaussian unit", "40", "φ(x,y)=(−x, i y), where i lies in the base field"],
        ["Eisenstein unit", "78", "φ(x,y)=(ζ₃x,y), where ζ₃ lies in the base field"],
        ["Relative binary Frobenius", "70", "(x,y) ↦ (x²,y²), with the coefficient-squared target model"],
    ], [142, 77, 280], small=True)
    p("The explicit 2-isogeny has x′=x−1/x and y′=y(1+1/x²) away from "
      "its kernel. The kernel and infinity map to infinity; exceptional inputs "
      "are tested. A class identified by equal trace can lack an explicit map in "
      "this release. The dataset records this distinction and does not supply a "
      "complete isogeny graph or volcano.")
    h("Selected-subgroup endomorphism tests")
    p(f"The classifier supports {ev['supported_cases']} of the {ev['curves_classified']} "
      f"core curves, including {ev['binary_frobenius_cases']} binary Koblitz Frobenius "
      f"cases. All {ev['scalar_vectors_checked']} scalar decompositions verify a "
      "kernel-lattice determinant, the exact congruence k≡k₁+k₂λ (mod r), "
      "nearest-plane bounds, and the signed point reconstruction "
      "[k]G=[k₁]G+[k₂]φ(G). Unsupported cases carry explicit reasons.")
    p("For binary Koblitz coefficients in F₂, τ(x,y)=(x²,y²) is an "
      "endomorphism of the same curve, with τ²−(2a−1)τ+[2]=0. For binary "
      "coefficients outside F₂, this squaring map ordinarily lands on a different "
      "coefficient-squared curve. It must not be used as a same-curve GLV map.")
    p("A map must preserve the chosen cyclic subgroup before a subgroup "
      "eigenvalue or scalar decomposition is valid. The rank-two examples test "
      "this directly: [r]P=O alone cannot establish P∈〈G〉. Powers of an "
      "endomorphism with a quadratic minimal polynomial do not by themselves "
      "supply a balanced four-term GLV/GLS construction.", "ReportSmall")
    h("Reference timing samples")
    timing = {(r["curve_id"], r["operation"]): r for r in eb["results"]}
    selected = [cid for cid in ["p61-E3-5", "p64-E4-2", "K83-a0"]
                if (cid, "scalar_mul_subgroup_baseline") in timing]
    rows = [["Curve", "Scalar baseline (µs)", "Full two-term path (µs)"]]
    for cid in selected:
        rows.append([cid, f"{timing[cid, 'scalar_mul_subgroup_baseline']['median']/1000:.1f}",
                     f"{timing[cid, 'glv2_scalar_mul_full']['median']/1000:.1f}"])
    table(rows, [180, 160, 159], small=True)
    p("Medians of three variable-time Python samples, one operation per sample. "
      "The full path includes decomposition, the endomorphism, and joint signed "
      "multiplication; precomputed joint-multiplication timings are separately "
      "named in JSON. These samples illustrate the timing interface and are "
      "separate from polynomial relation-search performance.", "ReportSmall")
    story.append(PageBreak())

    # Page A6: genuine supplemental GLV/GLS and next measurement boundaries.
    title(6, "A genuine GLV + GLS4 example",
          "A separate quadratic-extension supplement; the 56-curve corpus is preserved")
    f = glv["field"]
    c = glv["curve"]
    table([
        ["Parameter", "Verified value"],
        ["Base prime / field", f"p=2⁴¹−853={f['base_prime']}; Fp²=Fp[u]/(u²−2)"],
        ["Field cardinality", f"q={f['cardinality']} (82 bits)"],
        ["Curve", "E: y²=x³+(7+5u); ordinary, j=0"],
        ["Rational group", "C₂ × C₂r; N=4r"],
        ["Prime subgroup", f"r={c['subgroup_order']} (80 bits)"],
        ["Geometric End ring", "Z[ζ₃], maximal quadratic order of discriminant −3"],
    ], [135, 364], small=True)
    p("Let η=1+u. The commuting efficient maps are "
      "φ(x,y)=(ζ₃x,y) and "
      "ψ(x,y)=(η<super>1−p</super>x<super>p</super>, "
      "η<super>3(1−p)/2</super>y<super>p</super>). The subgroup eigenvalues "
      "are λφ=174560118547837207691337 and λψ=896867993705181372003749. "
      "On geometric endomorphisms ψ²=−πq; on Fp²-rational points πq acts as "
      "identity, so ψ² acts as [−1]. These are different statements.")
    p("The four-term congruence uses eigenweights 1, λφ, λψ, λφλψ. "
      "The supplied reduced rank-four lattice has determinant −r. Exact rational "
      "Gram–Schmidt and Babai nearest-plane decomposition establish universal "
      "absolute component ceilings:")
    bounds = glv["lattices"]["Babai_universal_bounds"]
    p("[" + ", ".join(bounds["glv4"]["coordinate_absolute_bound_ceilings"]) + "]", "ReportCode")
    p("Every component is therefore at most 20 bits in absolute magnitude. "
      "The accompanying two-term GLV lattice guarantees 40-bit components. "
      "Sixteen scalar vectors verify both congruences and actual point "
      "reconstructions; the supplement also records field/map answers, certified "
      "prime factors, and 51 independent PARI checks.")
    p("Four-dimensional refers to a four-term congruence lattice using these "
      "GLV and GLS maps. The geometric endomorphism algebra remains quadratic. "
      "This supplemental curve has its own scalar/map benchmark family; it is "
      "not one of the current 55 polynomial fixtures.", "ReportSmall")
    h("What remains before a full index-calculus attack benchmark")
    p("This release measures two-summand relation solving on finite factor "
      "bases. A full ECDLP pipeline additionally needs larger-arity construction, "
      "generation-inclusive timing, relation-matrix rank, linear algebra modulo "
      "r, and individual-log recovery. An instrumented solver must provide "
      "maximum working degree before a solving-degree comparison is possible. "
      "The Sage/Singular adapter exports actual equations, but no unavailable "
      "external-backend timing is fabricated.", "ReportSmall")
    p("References: I. Semaev, <i>Summation polynomials and the discrete logarithm "
      "problem on elliptic curves</i>, ePrint 2004/031, "
      '<link href="https://eprint.iacr.org/2004/031">eprint.iacr.org/2004/031</link>. '
      "P. Longa and F. Sica, <i>Four-dimensional Gallant–Lambert–Vanstone scalar "
      "multiplication</i>, ePrint 2011/608, "
      '<link href="https://eprint.iacr.org/2011/608">eprint.iacr.org/2011/608</link>. '
      "Detailed contracts: docs/ALGEBRA_METRICS.md, docs/BACKEND_PROTOCOL.md, "
      "and the canonical JSON artifacts in this release.", "ReportSmall")

    def decorate(canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(TEAL)
        canvas.setLineWidth(0.7)
        canvas.line(48, 36, A4[0] - 48, 36)
        canvas.setFont("Body", 7.2)
        canvas.setFillColor(GRAY)
        canvas.drawString(48, 23, "ecbench-golden 1.0.0 • measured-fixture addendum")
        canvas.drawRightString(A4[0] - 48, 23, f"A{doc.page}")
        canvas.restoreState()

    addendum.parent.mkdir(parents=True, exist_ok=True)
    doc = BaseDocTemplate(str(addendum), pagesize=A4, leftMargin=48,
                          rightMargin=48, topMargin=44, bottomMargin=49,
                          title="ecbench: measured polynomial and endomorphism benchmark addendum",
                          author="ecbench reproducible calculation project")
    doc.addPageTemplates(PageTemplate(id="Addendum", frames=[Frame(
        48, 49, A4[0] - 96, A4[1] - 93, leftPadding=0, rightPadding=0,
        topPadding=0, bottomPadding=0)], onPage=decorate))
    doc.build(story)

    from pypdf import PdfReader, PdfWriter
    annex = PdfReader(addendum)
    if len(annex.pages) != 6:
        raise RuntimeError(f"Addendum overflow: expected 6 pages, got {len(annex.pages)}")
    writer = PdfWriter()
    writer.append(str(original), outline_item="Original finite-field and curve report")
    writer.append(str(addendum), outline_item="Measured ecbench benchmark addendum")
    writer.add_metadata({"/Title": "Elliptic-curve fields and ecbench benchmark battery",
                         "/Author": "ecbench reproducible calculation project"})
    with combined.open("wb") as stream:
        writer.write(stream)
    print(json.dumps({"addendum": str(addendum), "addendum_pages": 6,
                      "complete_report": str(combined),
                      "complete_report_pages": len(PdfReader(combined).pages)}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path,
        default=ROOT / "docs" / "elliptic_curve_fields_report.pdf")
    parser.add_argument("--addendum", type=Path,
                        default=ROOT / "docs" / "ecbench_benchmark_addendum.pdf")
    parser.add_argument("--combined", type=Path,
                        default=ROOT / "docs" / "ecbench_complete_report.pdf")
    args = parser.parse_args()
    build(args.addendum, args.original, args.combined)


if __name__ == "__main__":
    main()
