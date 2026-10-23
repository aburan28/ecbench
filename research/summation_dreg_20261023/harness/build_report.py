import md2pdf, os
md = open('REPORT.md').read()
md2pdf.render(
    md, 'out/summation_dreg_report.pdf',
    base_dir=os.path.abspath('.'),
    title='Controlling the Degree of Regularity in the Summation-Polynomial Step',
    subtitle='An independent, twice-reproduced measurement study of elliptic-curve index calculus on binary-Koblitz and prime-Koblitz curves',
    meta_lines=['Independent cross-check of the aburan28/crypto research program · 2026-10-23',
                'Engines: msolve 0.6.5 · Singular 4.3.3.0 · PARI/GP 2.15.4  —  all results run twice'])
print('PDF bytes:', os.path.getsize('out/summation_dreg_report.pdf'))
