"""Minimal, dependency-light Markdown -> PDF renderer (reportlab) for the report.
Supports: # ## ### #### headings; paragraphs; - and * bullets (one level);
1. ordered; | pipe | tables | with --- header rule; ![alt](path) images on their
own line; **bold**, *italic*, `code`; --- horizontal rule; > blockquote;
```fenced``` code blocks; [PAGEBREAK] marker. Deterministic."""
import re, os
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.fonts import addMapping

def _register_fonts():
    import matplotlib
    d = os.path.join(os.path.dirname(matplotlib.__file__), 'mpl-data', 'fonts', 'ttf')
    reg = [('DejaVu', 'DejaVuSans.ttf'), ('DejaVu-Bold', 'DejaVuSans-Bold.ttf'),
           ('DejaVu-Oblique', 'DejaVuSans-Oblique.ttf'), ('DejaVu-BoldOblique', 'DejaVuSans-BoldOblique.ttf'),
           ('DejaVuMono', 'DejaVuSansMono.ttf')]
    for name, fn in reg:
        p = os.path.join(d, fn)
        if os.path.exists(p):
            try: pdfmetrics.registerFont(TTFont(name, p))
            except Exception: pass
    # map so <b>/<i> inside paragraphs resolve to the right TTF
    addMapping('DejaVu', 0, 0, 'DejaVu'); addMapping('DejaVu', 1, 0, 'DejaVu-Bold')
    addMapping('DejaVu', 0, 1, 'DejaVu-Oblique'); addMapping('DejaVu', 1, 1, 'DejaVu-BoldOblique')
_register_fonts()
HAVE_DJV = 'DejaVu' in pdfmetrics.getRegisteredFontNames()
FSANS = 'DejaVu' if HAVE_DJV else FSANS
FBOLD = 'DejaVu-Bold' if HAVE_DJV else FBOLD
FOBL = 'DejaVu-Oblique' if HAVE_DJV else FOBL
FBOLDOBL = 'DejaVu-BoldOblique' if HAVE_DJV else FBOLDOBL
FMONO = 'DejaVuMono' if HAVE_DJV else FMONO
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
                                Image, HRFlowable, PageBreak, KeepTogether, ListFlowable, ListItem)

INK, INK2, MUT = colors.HexColor('#0b0b0b'), colors.HexColor('#52514e'), colors.HexColor('#8a8982')
BLUE, AQUA, ORANGE, RED, GREEN = (colors.HexColor(c) for c in ('#2a78d6', '#1baf7a', '#eb6834', '#e34948', '#008300'))
SURF, RULE, BG2 = colors.HexColor('#fcfcfb'), colors.HexColor('#dcdbd6'), colors.HexColor('#f2f1 ee'.replace(' ', ''))

def _styles():
    ss = getSampleStyleSheet()
    base = dict(fontName=FSANS, textColor=INK, leading=14.5, fontSize=10.3)
    S = {}
    S['body'] = ParagraphStyle('body', **base, alignment=TA_JUSTIFY, spaceAfter=6)
    S['h1'] = ParagraphStyle('h1', fontName=FBOLD, fontSize=19, leading=23, textColor=INK, spaceBefore=6, spaceAfter=10)
    S['h2'] = ParagraphStyle('h2', fontName=FBOLD, fontSize=14.5, leading=18, textColor=BLUE, spaceBefore=15, spaceAfter=6)
    S['h3'] = ParagraphStyle('h3', fontName=FBOLD, fontSize=11.8, leading=15, textColor=INK, spaceBefore=11, spaceAfter=4)
    S['h4'] = ParagraphStyle('h4', fontName=FBOLDOBL, fontSize=10.5, leading=13.5, textColor=INK2, spaceBefore=8, spaceAfter=3)
    S['bullet'] = ParagraphStyle('bullet', **base, alignment=TA_LEFT, spaceAfter=3, leftIndent=6)
    S['cap'] = ParagraphStyle('cap', fontName=FOBL, fontSize=8.8, leading=11.5, textColor=INK2, alignment=TA_CENTER, spaceBefore=3, spaceAfter=10)
    S['quote'] = ParagraphStyle('quote', fontName=FSANS, fontSize=9.8, leading=13.5, alignment=TA_LEFT, leftIndent=10, borderPadding=(2, 2, 2, 8), textColor=INK2, backColor=colors.HexColor('#f4f3ef'), spaceAfter=7, borderColor=BLUE)
    S['code'] = ParagraphStyle('code', fontName=FMONO, fontSize=8.3, leading=11, textColor=INK, backColor=colors.HexColor('#f4f3ef'), borderPadding=5, spaceAfter=7, leftIndent=2)
    S['tcell'] = ParagraphStyle('tcell', fontName=FSANS, fontSize=8.0, leading=9.8, textColor=INK)
    S['thead'] = ParagraphStyle('thead', fontName=FBOLD, fontSize=8.1, leading=10, textColor=colors.white)
    S['title'] = ParagraphStyle('title', fontName=FBOLD, fontSize=24, leading=28, textColor=INK, spaceAfter=4, alignment=TA_LEFT)
    S['sub'] = ParagraphStyle('sub', fontName=FSANS, fontSize=12, leading=16, textColor=INK2, spaceAfter=3, alignment=TA_LEFT)
    return S

def _inline(text):
    text = text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
    # protect code spans FIRST so * _ { } inside them are not treated as markup
    spans = []
    def _stash(m):
        spans.append(m.group(1))
        return '\x00%d\x00' % (len(spans) - 1)
    text = re.sub(r'`([^`]+)`', _stash, text)
    # links before emphasis (URLs may contain special chars but not * here)
    text = re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)', r'<link href="\2" color="#2a78d6">\1</link>', text)
    text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'(?<!\*)\*(?!\*)([^*]+?)\*(?!\*)', r'<i>\1</i>', text)
    # restore code spans as monospace font, escaping their own content
    def _restore(m):
        s = spans[int(m.group(1))]
        return '<font face="%s" size="8.8">%s</font>' % (FMONO, s)
    text = re.sub(r'\x00(\d+)\x00', _restore, text)
    return text

def render(md, out_path, base_dir='', title=None, subtitle=None, meta_lines=None):
    S = _styles()
    story = []
    if title:
        story += [Paragraph(_inline(title), S['title'])]
        if subtitle: story += [Paragraph(_inline(subtitle), S['sub'])]
        if meta_lines:
            for ml in meta_lines: story += [Paragraph(_inline(ml), S['sub'])]
        story += [Spacer(1, 4), HRFlowable(width='100%', color=BLUE, thickness=1.4), Spacer(1, 10)]
    lines = md.split('\n')
    i = 0
    def flush_table(tbl):
        if not tbl: return
        # tbl: list of rows (list of cell strs); second row is --- separator
        header = tbl[0]
        body = tbl[2:] if len(tbl) > 1 and set(''.join(tbl[1]).replace('|', '').replace(' ', '')) <= set('-:') else tbl[1:]
        data = [[Paragraph(_inline(c), S['thead']) for c in header]]
        for r in body:
            data.append([Paragraph(_inline(c), S['tcell']) for c in r])
        ncol = len(header)
        avail = 170 * mm
        # width heuristic: first column a bit wider
        w = [avail / ncol] * ncol
        t = Table(data, colWidths=w, repeatRows=1)
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), BLUE),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [SURF, colors.HexColor('#f3f2ee')]),
            ('GRID', (0, 0), (-1, -1), 0.4, RULE),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('LEFTPADDING', (0, 0), (-1, -1), 4), ('RIGHTPADDING', (0, 0), (-1, -1), 4),
            ('TOPPADDING', (0, 0), (-1, -1), 3), ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
        ]))
        story.append(t); story.append(Spacer(1, 8))
    while i < len(lines):
        ln = lines[i]
        if ln.strip() == '[PAGEBREAK]':
            story.append(PageBreak()); i += 1; continue
        if ln.startswith('```'):
            buf = []; i += 1
            while i < len(lines) and not lines[i].startswith('```'):
                buf.append(lines[i]); i += 1
            i += 1
            code = '<br/>'.join(l.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace(' ', '&nbsp;') for l in buf)
            story.append(Paragraph(code, S['code'])); continue
        if re.match(r'^!\[.*\]\(.+\)\s*$', ln):
            m = re.match(r'^!\[(.*)\]\((.+)\)\s*$', ln)
            alt, path = m.group(1), m.group(2)
            fp = path if os.path.isabs(path) else os.path.join(base_dir, path)
            if os.path.exists(fp):
                from reportlab.lib.utils import ImageReader
                iw, ih = ImageReader(fp).getSize()
                maxw = 165 * mm
                w = min(maxw, iw * 0.75)
                h = w * ih / iw
                story.append(Image(fp, width=w, height=h))
                if alt: story.append(Paragraph(_inline(alt), S['cap']))
            i += 1; continue
        if ln.startswith('| '):
            tbl = []
            while i < len(lines) and lines[i].startswith('|'):
                cells = [c.strip() for c in lines[i].strip().strip('|').split('|')]
                tbl.append(cells); i += 1
            flush_table(tbl); continue
        if ln.startswith('#### '):
            story.append(Paragraph(_inline(ln[5:]), S['h4'])); i += 1; continue
        if ln.startswith('### '):
            story.append(Paragraph(_inline(ln[4:]), S['h3'])); i += 1; continue
        if ln.startswith('## '):
            story.append(Paragraph(_inline(ln[3:]), S['h2'])); i += 1; continue
        if ln.startswith('# '):
            story.append(Paragraph(_inline(ln[2:]), S['h1'])); i += 1; continue
        if ln.strip() == '---':
            story.append(Spacer(1, 3)); story.append(HRFlowable(width='100%', color=RULE, thickness=0.7)); story.append(Spacer(1, 3)); i += 1; continue
        if ln.startswith('> '):
            story.append(Paragraph(_inline(ln[2:]), S['quote'])); i += 1; continue
        if re.match(r'^\s*[-*] ', ln):
            items = []
            while i < len(lines) and re.match(r'^\s*[-*] ', lines[i]):
                items.append(Paragraph(_inline(re.sub(r'^\s*[-*] ', '', lines[i])), S['bullet'])); i += 1
            story.append(ListFlowable([ListItem(it, leftIndent=12, value='•') for it in items], bulletType='bullet', start='•', leftIndent=10))
            story.append(Spacer(1, 4)); continue
        if re.match(r'^\s*\d+\. ', ln):
            items = []
            while i < len(lines) and re.match(r'^\s*\d+\. ', lines[i]):
                items.append(Paragraph(_inline(re.sub(r'^\s*\d+\. ', '', lines[i])), S['bullet'])); i += 1
            story.append(ListFlowable([ListItem(it, leftIndent=12) for it in items], bulletType='1', leftIndent=12))
            story.append(Spacer(1, 4)); continue
        if ln.strip() == '':
            i += 1; continue
        # paragraph (gather until blank)
        buf = [ln]; i += 1
        while i < len(lines) and lines[i].strip() != '' and not re.match(r'^(#|\||!\[|> |```|\s*[-*] |\s*\d+\. |---)', lines[i]):
            buf.append(lines[i]); i += 1
        story.append(Paragraph(_inline(' '.join(buf)), S['body']))
    def _footer(canvas, doc):
        canvas.saveState()
        canvas.setFont(FSANS, 7.5); canvas.setFillColor(MUT)
        canvas.drawString(20 * mm, 12 * mm, 'Summation-polynomial degree-of-regularity study — independent measurement report')
        canvas.drawRightString(190 * mm, 12 * mm, 'p. %d' % doc.page)
        canvas.restoreState()
    doc = SimpleDocTemplate(out_path, pagesize=A4, leftMargin=20 * mm, rightMargin=20 * mm,
                            topMargin=18 * mm, bottomMargin=18 * mm, title=title or 'Report')
    doc.build(story, onFirstPage=_footer, onLaterPages=_footer)
    return out_path

if __name__ == '__main__':
    import sys
    md = open(sys.argv[1]).read()
    render(md, sys.argv[2], base_dir=os.path.dirname(os.path.abspath(sys.argv[1])))
    print('wrote', sys.argv[2])
