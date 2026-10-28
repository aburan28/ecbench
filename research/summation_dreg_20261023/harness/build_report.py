"""Build REPORT.pdf from REPORT.md in the study directory (repo layout)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import md2pdf
md = open(os.path.join(STUDY, 'REPORT.md')).read()
out = os.path.join(STUDY, 'REPORT.pdf')
md2pdf.render(
    md, out, base_dir=STUDY,
    title='Controlling the Degree of Regularity in the Summation-Polynomial Step',
    subtitle='An independent, twice-reproduced measurement study of elliptic-curve index calculus on binary-Koblitz and prime-Koblitz curves',
    meta_lines=['Independent cross-check of the aburan28/crypto research program · v2, 2026-10-23',
                'Engines: msolve 0.6.5 · Singular 4.3.3.0 · PARI/GP 2.15.4  —  every result run twice'])
print('PDF bytes:', os.path.getsize(out))
