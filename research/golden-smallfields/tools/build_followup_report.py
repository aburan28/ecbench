#!/usr/bin/env python3
"""Build a reproducible six-page follow-up and append it to the existing PDF."""
import hashlib
import json
from pathlib import Path
import sys
from xml.sax.saxutils import escape

from reportlab import rl_config
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from build_benchmark_report import register_fonts

OUT = ROOT / 'followups' / '2026-10-15'
INK = colors.HexColor('#18334a')
TEAL = colors.HexColor('#207c85')


def load(name): return json.loads((OUT / name).read_text())


def build():
    rl_config.invariant = 1
    register_fonts()
    proof = load('row_equivalence_proofs.json')
    recovery = load('timeout_recovery.json')
    bench = load('preprocessing_benchmark.json')
    reg = load('regularity_preprocessing_comparison.json')
    density = load('relation_density.json')
    iso = load('represented_isogeny.json')
    complete = load('complete_mathematical_metrics.json')
    if len(bench.get('summary', [])) != 8:
        raise ValueError('The paired benchmark must finish before building the report')
    styles = {
        'title': ParagraphStyle('title', fontName='BodyBold', fontSize=21, leading=27, textColor=INK, spaceAfter=13),
        'body': ParagraphStyle('body', fontName='Body', fontSize=9.5, leading=14, textColor=INK, spaceAfter=9),
        'small': ParagraphStyle('small', fontName='Body', fontSize=8, leading=11, textColor=INK, spaceAfter=6),
        'cell': ParagraphStyle('cell', fontName='Body', fontSize=8, leading=10, textColor=INK),
        'head': ParagraphStyle('head', fontName='BodyBold', fontSize=8, leading=10, textColor=colors.white),
        'mono': ParagraphStyle('mono', fontName='Mono', fontSize=8.2, leading=12, textColor=INK, spaceAfter=10),
    }
    story = []
    def para(text, style='body'):
        story.append(Paragraph(escape(str(text)).replace('\n', '<br/>'), styles[style]))
    def page(title):
        if story: story.append(PageBreak())
        para('ecbench • Research follow-up • 15 October 2026', 'small')
        para(title, 'title')
    def table(headers, rows, widths):
        data = [[Paragraph(escape(str(cell)), styles['head']) for cell in headers]]
        data += [[Paragraph(escape(str(cell)), styles['cell']) for cell in row] for row in rows]
        value = Table(data, colWidths=widths, repeatRows=1, hAlign='LEFT')
        value.setStyle(TableStyle([('BACKGROUND', (0,0), (-1,0), TEAL),
            ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.HexColor('#edf4f6'), colors.white]),
            ('VALIGN', (0,0), (-1,-1), 'TOP'), ('LEFTPADDING', (0,0), (-1,-1), 7),
            ('RIGHTPADDING', (0,0), (-1,-1), 7), ('TOPPADDING', (0,0), (-1,-1), 6),
            ('BOTTOMPADDING', (0,0), (-1,-1), 6)]))
        story.append(value); story.append(Spacer(1, 11))
    def seconds(value): return '—' if value is None else f'{value:.3f}'
    def footer(canvas, doc):
        canvas.setFont('Body', 8); canvas.setFillColor(INK)
        canvas.drawString(36, 23, 'ecbench | Exact arithmetic, measured experiments, explicit proof boundaries')
        canvas.drawRightString(A4[0]-36, 23, f'Follow-up {doc.page}')

    page('Open leads: what the follow-up found')
    raw = [r for r in bench['rows'] if r['variant'] == 'original']
    reduced = [r for r in bench['rows'] if r['variant'] == 'row_reduced']
    raw_ok = sum(r['result']['status'] == 'ok' for r in raw)
    reduced_ok = sum(r['result']['status'] == 'ok' for r in reduced)
    para('The original 56 curves, 69 polynomial fixtures, historical measurements, and golden release remain unchanged. This follow-up investigates unresolved solver timeouts, hidden linear constraints, random-target yield, and a concrete isogeny path.')
    table(['Lead', 'New evidence', 'Status'], [
        ['GF(2) preprocessing', f'{len(proof["transformations"])} exact transformations, two-way row-span proofs, and exhaustive original-domain checks', 'Verified'],
        ['Difficult relation solves', f'Paired reruns: {reduced_ok}/{len(reduced)} transformed and {raw_ok}/{len(raw)} original runs completed at a 20 s budget', 'Measured'],
        ['Regularity as a grouping feature', 'Equal affine ideals can have different highest-input-part regularities after row reduction', 'Generator dependent'],
        ['Sampled-target zero yield', f'Exact pair-sum supports for {density["distinct_factor_bases"]} factor bases explain the extremely low uniform-target density', 'Quantified'],
        ['Same order, different groups', 'An explicit dual pair of degree-two isogenies connects two existing p61 models', 'Constructive'],
    ], [114, 319, 90])
    para('The practical improvement is ordinary Gaussian elimination applied before Gröbner computation. Every output equation has an exact combination witness, and every original equation is checked to lie in the output row span. These checks establish equality of polynomial ideals, rather than only agreement on sampled roots.')
    para('The timing experiment covers four selected difficult binary cases and two actual SymPy algorithms, with three paired repetitions. It does not establish a universal backend ranking. The original solver’s maximum working degree remains unavailable; no final-basis degree is substituted for it.')

    page('A certified preprocessing step')
    para('Collect the monomials occurring in the exact recorded GF(2) equations. Interpret each polynomial as a bit row over those columns, row-reduce with XOR, retain nonzero rows, and put any discovered constant unit first. No monomial multiples or guessed relations are introduced.')
    para('For original generators F and reduced generators R, the certificate proves both inclusions of their constant-coefficient row spans:\nR ⊆ span(F), and F ⊆ span(R).\nConsequently ⟨F⟩ = ⟨R⟩ as affine polynomial ideals. Variables, coordinate domains, factor bases, signed targets, and expected point decompositions are unchanged.', 'mono')
    selected = {}
    for item in proof['transformations']:
        p = item['proof']
        if p['source_case_id'] in {row['source_case_id'] for row in bench['summary']}:
            selected[p['source_case_id']] = p
    table(['Case', 'Input rows', 'Output rows', 'Linear output rows', 'Unit'], [
        [cid, p['input_generator_count'], p['output_generator_count'], p['new_linear_generators'], p['unit_generator_found']]
        for cid, p in selected.items()], [267, 57, 62, 78, 59])
    para(f'All {len(proof["transformations"])} binary fixtures admit this exact transformation. Their derived systems pass exhaustive evaluation over {proof["original_system_truth_rechecked"]["domain_assignments_checked"]:,} declared assignments, followed by independent subgroup and signed-point checks. Prime-characteristic cases are explicitly rejected by this GF(2)-specific implementation.')
    para('Proof masks identify the original equations XORed to obtain each output equation. The verifier also checks reverse inclusion: checking only that outputs are valid combinations would permit dropping necessary equations. Tests reject forged units, false combination witnesses, omitted rows, changed targets, and an unsupported coefficient characteristic.')

    page('Paired solver measurements')
    para('Each row below summarizes three new paired repetitions. Both variants use fresh isolated workers, identical targets and finite domains, grevlex, and a 20-second worker budget. Diagnostic regularity work is disabled during these performance pairs. Variant order alternates across pairs.')
    table(['Case / method', 'Complete: raw / RREF', 'Median raw wall (s)', 'Median RREF wall (s)', 'Paired wall speedup'], [
        [row['source_case_id'] + '\n' + row['method'],
         f'{row["variants"]["original"]["completed"]}/3 · {row["variants"]["row_reduced"]["completed"]}/3',
         seconds(row['variants']['original']['median_pipeline_wall_seconds_completed']),
         seconds(row['variants']['row_reduced']['median_pipeline_wall_seconds_completed']),
         '—' if row['median_paired_wall_speedup'] is None else f'{row["median_paired_wall_speedup"]:.2f}×']
        for row in bench['summary']], [236, 76, 69, 73, 69])
    para('Wall time starts before preprocessing and exact proof verification and ends after worker completion. It includes interpreter/SymPy startup, serialization, root extraction, and signed relation checks. The paired core metric separately adds preprocessing and proof verification to the worker’s core solve timer. Neither includes factor-base construction or fixture generation.')
    para('A dash in a median or speedup means no matched completed observations. Timeouts remain censored; the 20-second budget is not treated as an exact solve duration. Completed-only medians are shown with completion counts. Three repetitions provide local measurements, not confidence intervals or production implementation results.')

    page('Regularity and relation-yield power')
    by_case = {}
    for row in reg['rows']:
        by_case.setdefault(row['source_case_id'], {})[row['variant']] = row['result'].get('degree_of_regularity', {}).get('value')
    table(['Equivalent affine system', 'Original top-part regularity', 'RREF top-part regularity'], [
        [cid, values.get('original'), values.get('row_reduced')] for cid, values in by_case.items()], [335, 94, 94])
    para('Regularity here is the first zero Hilbert-function value of J = ⟨top(f) : f is a recorded input generator⟩. Row-equivalent affine generator lists can define different J. Therefore this number is meaningful only with a fixed generator list and preprocessing policy; it is not an intrinsic curve invariant.')
    maximum = {}
    for row in density['factor_bases']:
        if row['field_id'] not in maximum or row['uniform_nonzero_target_success_probability']['approximation'] > maximum[row['field_id']]['uniform_nonzero_target_success_probability']['approximation']:
            maximum[row['field_id']] = row
    table(['Field', 'Highest-density recorded factor base', 'Points', 'Exact-support probability ≈'], [
        [fid, row['curve_id'], row['factor_base_points'], f'{row["uniform_nonzero_target_success_probability"]["approximation"]:.3g}']
        for fid, row in sorted(maximum.items())], [54, 257, 54, 158])
    expected = density['cohort']['uniform_model_expected_successes']['approximation']
    para(f'For each pool, enumerate every unordered signed pair, combine repeated factors correctly, and exclude identity targets. If S distinct nonzero sums occur in an r-order subgroup, an ideal uniform nonzero target has success probability S/(r−1). Across the 37 seeded-case contexts, the uniform-model expected successes total only {expected:.6g}; the same sum bounds the chance of any success.')
    para('These are exact support calculations under a hypothetical uniform-target model. The existing hash-seeded targets are deterministic and modulo-biased, so the figure is not a p-value for those recorded trials. The zero-yield observations have very little power to compare curve families. Exact fractions and every supported target are retained in the JSON.')

    page('An explicit isogeny between represented curves')
    para('Over p = 2^61 − 1, the two p61 E4 records have the same trace and order 2^61, yet their rational groups differ:\np61-E4-0: Z/2 × Z/2^60\np61-E4-1: Z/2^61.', 'mono')
    para('The previously recorded degree-two quotient can be connected all the way to the second stored representative. This closes the gap between a trace-based class-membership statement and a concrete map for this pair.')
    para('E: y² = x³ − x\nE′: Y² = X³ + 4X\nφ(x,y) = (x − 1/x, y(1 + 1/x²)).\nIts kernel is {O,(0,0)}.', 'mono')
    para(f'The stored target is y² = x³ − 37x. Choose u = {iso["representative_scale_u"]} modulo p, with 4u⁴ = −37. Scale the quotient output by (u²,u³) to reach that exact model.')
    para('For the dual, undo the scale to recover (X,Y), then apply\nψ(X,Y) = ((X + 4/X)/4, Y(1 − 4/X²)/8).\nAt infinity or X=0, the appropriate map returns infinity.', 'mono')
    para(f'Symbolic rational-function checks verify both curve equations and both dual compositions, ψφ=[2] and φψ=[2], after the appropriate scaling. The dataset adds {len(iso["vectors"])} explicit point/dual-composition vectors and {iso["homomorphism_pair_checks"]} homomorphism checks, including exceptional points.')
    para('This example illustrates the distinction between isogeny and rational group isomorphism. It supplies one verified path between existing records; it does not construct every path or a complete isogeny volcano. The primary curve and class counts are unchanged.')

    page('Recovered work and remaining leads')
    s = recovery['summary']
    para(f'The timeout follow-up retried {s["attempts"]} historical incomplete solves or diagnostics at a 20-second budget. It recovered {s["recovered_solve_timeouts"]} of {s["previous_solve_timeouts"]} prior solve timeouts; {s["completed_regularity_diagnostics"]} diagnostic runs finished. Every completed relation solve matched the golden roots and signed decompositions. Historical 2.5-second observations remain unchanged.')
    para(f'The consolidated metric catalog now contains exact highest-input-part Hilbert diagnostics and verified final-basis metrics for all {complete["cases"]} original systems. Every value points to an actual completed computation. An algorithm-independent regularity value can be shared across methods; the timing of a censored F5b diagnostic is never filled from a Buchberger run.', 'small')
    table(['Priority', 'Open lead', 'Concrete next experiment'], [
        ['1', 'Positive random-target controls', 'Use small-field controls with enough pair-sum density to observe successes; hold subgroup size, factor-base policy, and sampling model fixed.'],
        ['2', 'Actual working-degree instrumentation', 'Add a precisely defined intermediate-degree trace or instrumented Macaulay/F4 backend; retain unavailable values until measured.'],
        ['3', 'Larger polynomial systems', 'Extend relation arity and vary subspace dimensions after profiling certified preprocessing; record new case identities and rank/yield.'],
        ['4', 'Independent external backend', 'Execute the existing Sage/Singular adapter on supported coefficient fields; compare exact roots and actual stage timings.'],
        ['5', 'Isogeny / endomorphism-aware factor bases', 'Compare matched factor-base pools and targets under explicit maps, including kernel and image restrictions.'],
        ['6', 'Complete discrete-log pipeline', 'Collect independent relation rows, measure rank and linear algebra modulo r, and verify individual-log recovery before claiming an attack runtime.'],
    ], [46, 151, 326])
    para('Reproduce from the project checkout:', 'small')
    para('python tools/follow_up_leads.py --phase all --timeout 20 --repetitions 3\npython -m unittest discover -s tests -v\npython tools/build_followup_report.py', 'mono')
    para('Sources and method context: the bundled SymPy 1.14.0 runner; Semaev, ePrint 2004/031; the original report’s Tate/isogeny references. New observations are supported by the exact transformations, pair-sum census, symbolic map checks, and raw measurements in followups/2026-10-15/.')

    addendum = ROOT / 'docs' / 'ecbench_open_leads_followup.pdf'
    document = SimpleDocTemplate(str(addendum), pagesize=A4, rightMargin=36, leftMargin=36, topMargin=32, bottomMargin=42,
        title='ecbench: open-lead follow-up, 15 October 2026', author='ecbench research')
    document.build(story, onFirstPage=footer, onLaterPages=footer)
    if len(PdfReader(addendum).pages) != 6:
        raise RuntimeError('Follow-up layout must fit six pages; inspect overflow before release')
    writer = PdfWriter()
    for pdf in [ROOT / 'docs' / 'ecbench_complete_report.pdf', addendum]:
        for pdf_page in PdfReader(pdf).pages: writer.add_page(pdf_page)
    writer.add_metadata({'/Title': 'ecbench: complete report with open-lead follow-up', '/Author': 'ecbench research'})
    combined = ROOT / 'docs' / 'ecbench_updated_report.pdf'
    with combined.open('wb') as stream: writer.write(stream)
    print(json.dumps({'followup_pages': len(PdfReader(addendum).pages), 'combined_pages': len(PdfReader(combined).pages),
        'pdf': str(combined), 'sha256': hashlib.sha256(combined.read_bytes()).hexdigest()}, indent=2))


if __name__ == '__main__': build()
