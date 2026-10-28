"""Build the SVG artwork for the Steinbeck Lab organisation README.

Run from this folder:
    uv run --with rdkit --with fonttools --with brotli python build.py

Inputs: logos in src/, fonts in fonts/. Output: *.svg next to this file.
Every image owns its background, so one file serves light and dark GitHub.
All structures are drawn by RDKit from SMILES; all numbers are real
(the descriptor values come from a live call to api.naturalproducts.net).
"""
import base64
import io
import math
import re
import string
from pathlib import Path

from fontTools import subset
from fontTools.ttLib import TTFont
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

HERE = Path(__file__).parent
SRC = HERE / "src"
FONTS = HERE / "fonts"

SANS = "'Schibsted Grotesk', 'Segoe UI', Helvetica, Arial, sans-serif"
MONO = "'IBM Plex Mono', ui-monospace, SFMono-Regular, Menlo, monospace"

rdDepictor.SetPreferCoordGen(True)


# ---------------------------------------------------------------- fonts

def font_face(family, file, weight):
    font = TTFont(FONTS / file)
    opts = subset.Options()
    opts.flavor = "woff2"
    sub = subset.Subsetter(opts)
    sub.populate(text=string.printable + "·")
    sub.subset(font)
    buf = io.BytesIO()
    font.flavor = "woff2"
    font.save(buf)
    data = base64.b64encode(buf.getvalue()).decode()
    return (f"@font-face{{font-family:'{family}';font-weight:{weight};"
            f"src:url(data:font/woff2;base64,{data}) format('woff2')}}")


FONT_CSS = "".join([
    font_face("Schibsted Grotesk", "schibsted-grotesk-latin-800-normal.woff2", 800),
    font_face("Schibsted Grotesk", "schibsted-grotesk-latin-400-normal.woff2", 400),
    font_face("IBM Plex Mono", "ibm-plex-mono-latin-500-normal.woff2", 500),
])
MOTION_OFF = "@media (prefers-reduced-motion:reduce){*{animation:none!important}}"


# ---------------------------------------------------------------- helpers

def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def png_uri(name):
    return "data:image/png;base64," + base64.b64encode((SRC / name).read_bytes()).decode()


def png_size(name):
    b = (SRC / name).read_bytes()
    return int.from_bytes(b[16:20], "big"), int.from_bytes(b[20:24], "big")


def svg(w, h, body, css="", title=""):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
            f'viewBox="0 0 {w} {h}" role="img"><title>{esc(title)}</title>'
            f"<style>{FONT_CSS}{css}{MOTION_OFF}</style>{body}</svg>")


def text(x, y, s, size, fill, font=SANS, weight=400, anchor="start", extra=""):
    return (f'<text x="{x}" y="{y}" font-family="{font}" font-weight="{weight}" font-size="{size}" '
            f'fill="{fill}" text-anchor="{anchor}"{extra}>{s}</text>')


def arrow_ne(x, y, color, size=14):
    s = size
    return (f'<path d="M{x} {y + s}L{x + s} {y}M{x + 4} {y}H{x + s}V{y + s - 4}" fill="none" '
            f'stroke="{color}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>')


def arrow(x1, y1, x2, y2, color, width=2):
    """A straight arrow with an open head, pointing at (x2, y2)."""
    a = math.atan2(y2 - y1, x2 - x1)
    h = 10
    p1 = (x2 - h * math.cos(a - 0.5), y2 - h * math.sin(a - 0.5))
    p2 = (x2 - h * math.cos(a + 0.5), y2 - h * math.sin(a + 0.5))
    return (f'<path d="M{x1:.1f} {y1:.1f}L{x2:.1f} {y2:.1f}M{p1[0]:.1f} {p1[1]:.1f}L{x2:.1f} {y2:.1f}L{p2[0]:.1f} {p2[1]:.1f}" '
            f'fill="none" stroke="{color}" stroke-width="{width}" stroke-linecap="round" stroke-linejoin="round"/>')


def check(x, y, color):
    return (f'<path d="M{x} {y + 7}l6 6l12 -14" fill="none" stroke="{color}" stroke-width="3" '
            f'stroke-linecap="round" stroke-linejoin="round"/>')


def logo(name, x, y, box_w, box_h, multiply=False):
    """Left-align a logo inside a box, keeping its aspect ratio."""
    w, h = png_size(name)
    s = min(box_w / w, box_h / h)
    blend = ' style="mix-blend-mode:multiply"' if multiply else ""
    return (f'<image href="{png_uri(name)}" x="{x}" y="{y + (box_h - h * s) / 2:.1f}" '
            f'width="{w * s:.1f}" height="{h * s:.1f}"{blend}/>')


def rgb(c):
    return tuple(int(c[i:i + 2], 16) / 255 for i in (1, 3, 5))


def mol(smiles, w, h, ink, n, o, s=None, bond=22, lw=2, font=0.6, attach=False, cls=None):
    """RDKit drawing of one molecule as an SVG fragment, transparent background."""
    m = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DSVG(w, h)
    op = d.drawOptions()
    op.clearBackground = False
    op.bondLineWidth = lw
    op.scaleBondWidth = False
    op.fixedBondLength = bond
    op.padding = 0.06
    op.baseFontSize = font
    op.dummiesAreAttachments = attach
    op.setAtomPalette({-1: rgb(ink), 0: rgb(ink), 6: rgb(ink), 7: rgb(n), 8: rgb(o), 16: rgb(s or o)})
    d.DrawMolecule(m)
    d.FinishDrawing()
    body = d.GetDrawingText()
    body = re.sub(r"<\?xml.*?\?>|<svg[^>]*>|</svg>|<!-- .*? -->|<rect[^>]*/>", "", body, flags=re.S)
    body = re.sub(r" class='[^']*'", "", body)
    if cls:
        body = cls(body)
    return f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" overflow="visible">{body}</svg>'


def place(x, y, inner):
    return f'<g transform="translate({x} {y})">{inner}</g>'


INVERT = ('<filter id="inv"><feColorMatrix type="matrix" values="-1 0 0 0 1 0 -1 0 0 1 0 0 -1 0 1 0 0 0 1 0"/>'
          '<feColorMatrix type="hueRotate" values="180"/></filter>')


# ---------------------------------------------------------------- hero

HERO_MOLS = [
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("camptothecin", "CC[C@@]1(O)C(=O)OCc2c1cc1-c3nc4ccccc4cc3Cn1c2=O"),
    ("capsaicin", "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O"),
    ("quercetin", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    ("nicotine", "CN1CCC[C@H]1c1cccnc1"),
    ("penicillin G", "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O"),
    ("piperine", "O=C(/C=C/C=C/c1ccc2c(c1)OCO2)N1CCCCC1"),
    ("menthol", "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O"),
    ("resveratrol", "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1"),
    ("salicin", "OC[C@H]1O[C@@H](Oc2ccccc2CO)[C@H](O)[C@@H](O)[C@@H]1O"),
    ("vanillin", "COc1cc(C=O)ccc1O"),
    ("theobromine", "Cn1cnc2c1c(=O)[nH]c(=O)n2C"),
]


def hero():
    W, H = 1280, 600
    ink, dim, n, o, s = "#efe9dd", "#a39b8d", "#8fb0ff", "#ff8a73", "#e8c15a"
    # One row of real natural products, drawn twice so the loop is seamless.
    cell_h, row, x = 190, [], 0
    for name, smi in HERO_MOLS:
        m = Chem.MolFromSmiles(smi)
        rdDepictor.Compute2DCoords(m)
        pts = m.GetConformer().GetPositions()
        cw = int(max(160, (pts[:, 0].max() - pts[:, 0].min()) * 27 + 80))
        row.append(place(x, 0, mol(smi, cw, cell_h, "#d9d2c4", n, o, s, bond=27, lw=2.2, font=0.55)))
        row.append(text(x + cw / 2, cell_h + 26, esc(name), 14, dim, MONO, 500, "middle"))
        x += cw + 40
    loop = x
    strip = "".join(row)
    css = (f".mq{{animation:mq 70s linear infinite}}@keyframes mq{{to{{transform:translateX(-{loop}px)}}}}"
           )
    fade = ('<linearGradient id="edge"><stop offset="0" stop-color="#fff" stop-opacity="0"/>'
            '<stop offset=".12" stop-color="#fff"/><stop offset=".88" stop-color="#fff"/>'
            '<stop offset="1" stop-color="#fff" stop-opacity="0"/></linearGradient>'
            f'<mask id="m"><rect x="0" y="0" width="{W}" height="{H}" fill="url(#edge)"/></mask>')
    body = f"""<defs>{INVERT}{fade}<clipPath id="r"><rect width="{W}" height="{H}" rx="20"/></clipPath></defs>
<g clip-path="url(#r)">
<rect width="{W}" height="{H}" fill="#0f0f0e"/>
<image href="{png_uri('labmark.png')}" x="1124" y="52" width="100" height="100" filter="url(#inv)"/>
<g>
{text(60, 92, "NATURAL PRODUCTS CHEMINFORMATICS · FRIEDRICH SCHILLER UNIVERSITY JENA", 15, dim, MONO, 500)}
{text(52, 236, "Steinbeck Lab", 150, ink, SANS, 800, extra=' letter-spacing="-4.5"')}
{text(60, 296, "Open software for chemistry and natural products.", 30, ink)}
</g>
<g mask="url(#m)"><g transform="translate(0 330)"><g class="mq"><g id="row">{strip}</g><use href="#row" x="{loop}"/></g></g></g>
</g>"""
    return svg(W, H, body, css, "Steinbeck Lab: open software for chemistry and natural products, Friedrich Schiller University Jena")


# ---------------------------------------------------------------- posters

POSTER_W, POSTER_H = 1280, 460
PX, PY, PW, PH = 640, 40, 600, 380  # the right-hand panel


def poster(p, visual, css=""):
    c = p["c"]
    lines = "".join(f'<tspan x="64" dy="{0 if i == 0 else 42}">{esc(l)}</tspan>' for i, l in enumerate(p["lines"]))
    ref = text(64, 372, esc(p["ref"]), 19, c["dim"]) if p.get("ref") else ""
    body = f"""<defs><clipPath id="r"><rect width="{POSTER_W}" height="{POSTER_H}" rx="20"/></clipPath></defs>
<g clip-path="url(#r)">
<rect width="{POSTER_W}" height="{POSTER_H}" fill="{c['field']}"/>
{logo(p['logo'], 64, 56, 470, p.get('logo_h', 104), p.get('multiply', False))}
<text x="64" y="238" font-family="{SANS}" font-size="34" fill="{c['ink']}" letter-spacing="-0.5">{lines}</text>
{ref}
{text(64, 410, esc(p['domain']), 20, c['accent'], MONO, 500)}
{arrow_ne(64 + len(p['domain']) * 12 + 14, 395, c['accent'])}
{visual}
</g>"""
    return svg(POSTER_W, POSTER_H, body, css, p["name"] + ": " + " ".join(p["lines"]))


def v_coconut(c):
    picks = [("quercetin", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
             ("nicotine", "CN1CCC[C@H]1c1cccnc1"),
             ("menthol", "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O"),
             ("resveratrol", "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1"),
             ("vanillin", "COc1cc(C=O)ccc1O"),
             ("salicin", "OC[C@H]1O[C@@H](Oc2ccccc2CO)[C@H](O)[C@@H](O)[C@@H]1O")]
    out, cw, ch = [], 200, 180
    for i, (name, smi) in enumerate(picks):
        x, y = PX + (i % 3) * cw, PY + 10 + (i // 3) * ch
        out.append(f'<rect x="{x + 6}" y="{y}" width="{cw - 12}" height="{ch - 14}" rx="12" fill="{c["tile"]}"/>')
        out.append(place(x + 10, y + 4, mol(smi, cw - 20, 132, c["ink"], c["n"], c["o"], bond=20, lw=1.7, font=0.5)))
        out.append(text(x + cw / 2, y + 150, name, 14, c["dim"], MONO, 500, "middle"))
    return "".join(out)


def v_decimer(c):
    smi = "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O"
    canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    px, py, pw, ph = PX + 150, PY + 18, 300, 220
    page = mol(smi, pw - 20, ph - 20, "#3b3833", "#3b3833", "#3b3833", bond=20, lw=1.8, font=0.55)
    sw = len(canon) * 12.6
    css = (f".scan{{animation:scan 3.6s cubic-bezier(.45,0,.55,1) infinite alternate}}"
           f"@keyframes scan{{from{{transform:translateY(0)}}to{{transform:translateY({ph - 6}px)}}}}"
           )
    body = f"""<filter id="soft"><feGaussianBlur stdDeviation="0.55"/></filter>
<g transform="rotate(-3 {px + pw / 2} {py + ph / 2})">
<rect x="{px}" y="{py}" width="{pw}" height="{ph}" rx="4" fill="#efe6d2"/>
<g filter="url(#soft)">{place(px + 10, py + 10, page)}</g>
<g class="scan"><rect x="{px - 8}" y="{py + 2}" width="{pw + 16}" height="3" fill="{c['accent']}"/>
<rect x="{px}" y="{py - 26}" width="{pw}" height="28" fill="{c['accent']}" opacity=".14"/></g>
</g>
{arrow(PX + 300, PY + 256, PX + 300, PY + 290, c['accent'], 2.4)}
{text(round(PX + 300 - sw / 2, 1), PY + 334, esc(canon), 21, c['ink'], MONO, 500)}
{text(PX + 300, PY + 370, "SMILES", 13, c['dim'], MONO, 500, "middle", ' letter-spacing="2"')}"""
    return body, css


def v_chemaudit(c):
    salt = "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]"
    parent = "CC(=O)Oc1ccccc1C(=O)O"
    rows = ["validated", "standardized", "quality scored"]
    out = [place(PX, PY + 10, mol(salt, 270, 200, c["ink"], c["n"], c["o"], bond=27, lw=2.2, font=0.6)),
           arrow(PX + 280, PY + 110, PX + 330, PY + 110, c["accent"], 2.4),
           place(PX + 340, PY + 10, mol(parent, 250, 200, c["ink"], c["n"], c["o"], bond=27, lw=2.2, font=0.6))]
    for i, r in enumerate(rows):
        y = PY + 268 + i * 38
        out.append(check(PX + 40, y, c["accent"]))
        out.append(text(PX + 78, y + 18, r, 22, c["ink"]))
    out.append(text(PX + 135, PY + 230, "salt form in", 14, c["dim"], MONO, 500, "middle"))
    out.append(text(PX + 465, PY + 230, "parent out", 14, c["dim"], MONO, 500, "middle"))
    return "".join(out)


def v_cms(c):
    # Response values from GET api.naturalproducts.net/latest/chem/descriptors (caffeine, rdkit).
    k, v, d = c["accent"], c["ink"], c["dim"]
    lines = [
        [("GET ", d), ("/latest/chem/descriptors", v)],
        [("    ?smiles=", d), ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", k)],
        [("    &toolkit=", d), ("rdkit", k)],
        [],
        [("{", d)],
        [('  "molecular_weight"', k), (": ", d), ("194.19", v), (",", d)],
        [('  "alogp"', k), (": ", d), ("-1.03", v), (",", d)],
        [('  "topological_polar_surface_area"', k), (": ", d), ("61.82", v), (",", d)],
        [('  "qed_drug_likeliness"', k), (": ", d), ("0.54", v), (",", d)],
        [("  ...", d)],
        [("}", d)],
    ]
    out = [f'<rect x="{PX}" y="{PY}" width="{PW - 20}" height="{PH}" rx="14" fill="{c["tile"]}"/>',
           text(PX + 28, PY + 40, "api.naturalproducts.net", 15, d, MONO, 500)]
    for i, parts in enumerate(lines):
        spans = "".join(f'<tspan fill="{col}">{esc(t)}</tspan>' for t, col in parts)
        out.append(f'<text x="{PX + 28}" y="{PY + 88 + i * 27}" font-family="{MONO}" font-weight="500" '
                   f'font-size="17" xml:space="preserve">{spans}</text>')
    return "".join(out)


def v_mortar(c):
    whole = "O=C(/C=C/C=C/c1ccc2c(c1)OCO2)N1CCCCC1"
    m = Chem.MolFromSmiles(whole)
    # Cut every bond between a ring atom and a chain atom: ring systems fall away from the linker.
    cuts = [b.GetIdx() for b in m.GetBonds()
            if not b.IsInRing() and (b.GetBeginAtom().IsInRing() != b.GetEndAtom().IsInRing())]
    frags = Chem.GetMolFrags(Chem.FragmentOnBonds(m, cuts), asMols=True)
    out = [place(PX + 20, PY - 10, mol(whole, 560, 170, c["ink"], c["n"], c["o"], bond=22, lw=2, font=0.6)),
           f'<line x1="{PX + 30}" y1="{PY + 178}" x2="{PX + 570}" y2="{PY + 178}" stroke="{c["accent"]}" '
           f'stroke-width="2" stroke-dasharray="2 8" stroke-linecap="round"/>']
    widths = [200, 190, 170]
    x = PX
    for f, w in zip(frags, widths):
        smi = Chem.MolToSmiles(f)
        out.append(f'<rect x="{x + 6}" y="{PY + 200}" width="{w - 12}" height="170" rx="12" fill="{c["tile"]}"/>')
        out.append(place(x + 10, PY + 210, mol(smi, w - 20, 150, c["ink"], c["n"], c["o"], bond=19, lw=1.8, font=0.55, attach=True)))
        x += w
    return "".join(out)


def v_orthonym(c):
    smi = "COc1cc(C=O)ccc1O"
    out = [place(PX, PY + 40, mol(smi, 250, 250, c["ink"], c["n"], c["o"], bond=34, lw=2.6, font=0.7)),
           text(PX + 262, PY + 150, "4-hydroxy-3-", 31, c["ink"], SANS, 800, extra=' letter-spacing="-0.8"'),
           text(PX + 262, PY + 190, "methoxybenzaldehyde", 31, c["ink"], SANS, 800, extra=' letter-spacing="-0.8"')]
    # The round trip: structure to name, and OPSIN reading the name back.
    out.append(f'<path d="M{PX + 150} {PY + 34}C{PX + 230} {PY - 6} {PX + 330} {PY + 10} {PX + 380} {PY + 100}" '
               f'fill="none" stroke="{c["accent"]}" stroke-width="2"/>')
    out.append(arrow(PX + 372, PY + 86, PX + 381, PY + 102, c["accent"], 2))
    out.append(text(PX + 290, PY + 24, "Orthonym", 15, c["accent"], MONO, 500))
    out.append(f'<path d="M{PX + 420} {PY + 222}C{PX + 400} {PY + 330} {PX + 260} {PY + 350} {PX + 150} {PY + 300}" '
               f'fill="none" stroke="{c["dim"]}" stroke-width="2" stroke-dasharray="1 7" stroke-linecap="round"/>')
    out.append(arrow(PX + 166, PY + 306, PX + 150, PY + 300, c["dim"], 2))
    out.append(text(PX + 330, PY + 352, "OPSIN reads it back", 15, c["dim"], MONO, 500))
    out.append(check(PX + 516, PY + 334, c["accent"]))
    return "".join(out)


def v_surge(c):
    isomers = ["CCCCO", "CCC(C)O", "CC(C)CO", "CC(C)(C)O", "CCOCC", "CCCOC", "CC(C)OC"]
    out = [f'<rect x="{PX}" y="{PY}" width="{PW - 20}" height="56" rx="12" fill="{c["tile"]}"/>',
           f'<text x="{PX + 24}" y="{PY + 36}" font-family="{MONO}" font-weight="500" font-size="21" xml:space="preserve">'
           f'<tspan fill="{c["dim"]}">$ </tspan><tspan fill="{c["ink"]}">surge -S C4H10O</tspan></text>',
           text(PX + PW - 44, PY + 36, "7 isomers", 17, c["accent"], MONO, 500, "end")]
    cw, ch = 145, 150
    for i, smi in enumerate(isomers):
        col, rw = (i, 0) if i < 4 else (i - 4 + 0.5, 1)
        x, y = PX + col * cw, PY + 76 + rw * (ch + 4)
        out.append(place(x, y, mol(smi, cw - 10, ch, c["ink"], c["n"], c["o"], bond=31, lw=2.4, font=0.65)))
    return "".join(out)


def v_sherlock(c):
    # 13C NMR of caffeine (CDCl3), ppm; Sherlock works from spectra like this to candidate structures.
    peaks = [(155.4, .78), (151.7, .7), (148.7, .74), (141.4, .9), (107.6, .62), (33.6, 1), (29.7, .95), (27.9, .92)]
    ax0, ax1, base, top = PX + 10, PX + 400, PY + 300, PY + 60
    sx = lambda ppm: ax1 - (ppm / 180) * (ax1 - ax0)
    out = [f'<line x1="{ax0}" y1="{base}" x2="{ax1}" y2="{base}" stroke="{c["dim"]}" stroke-width="1.5"/>']
    for t in (0, 50, 100, 150):
        out.append(f'<line x1="{sx(t):.1f}" y1="{base}" x2="{sx(t):.1f}" y2="{base + 8}" stroke="{c["dim"]}" stroke-width="1.5"/>')
        out.append(text(round(sx(t), 1), base + 30, str(t), 14, c["dim"], MONO, 500, "middle"))
    out.append(text(ax1, base + 56, "13C  ppm", 14, c["dim"], MONO, 500, "end"))
    for i, (ppm, hgt) in enumerate(peaks):
        x = sx(ppm)
        y = base - hgt * (base - top)
        out.append(f'<line x1="{x:.1f}" y1="{base}" x2="{x:.1f}" '
                   f'y2="{y:.1f}" stroke="{c["accent"]}" stroke-width="2.6" stroke-linecap="round"/>')
    out.append(arrow(PX + 418, PY + 190, PX + 450, PY + 190, c["accent"], 2.4))
    out.append(place(PX + 450, PY + 90, mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C", 150, 200, c["ink"], c["n"], c["o"],
                                            bond=20, lw=2, font=0.55)))
    return "".join(out)


POSTERS = [
    dict(slug="coconut", name="COCONUT", logo="coconut.png", logo_h=96,
         lines=["The COlleCtion of Open", "Natural prodUcTs: deposit,", "curate and reuse NP data."],
         ref="Nucleic Acids Research, 2025", domain="coconut.naturalproducts.net",
         c=dict(field="#f1e8da", tile="#e9dccb", ink="#2e2118", dim="#6b5848", accent="#8a4f26", n="#2f55b5", o="#b4432f"),
         v=v_coconut),
    dict(slug="decimer", name="DECIMER.ai", logo="decimer_w.png", logo_h=96,
         lines=["Extract chemical structures", "from PDFs and images, and", "turn them into SMILES."],
         ref="Nature Communications, 2023", domain="decimer.ai",
         c=dict(field="#0f1d17", tile="#16281f", ink="#f1ece2", dim="#a7b0a5", accent="#f2d27a", n="#8fb0ff", o="#ff8a73"),
         v=v_decimer),
    dict(slug="chemaudit", name="ChemAudit", logo="chemaudit.png", logo_h=104,
         lines=["Validate, standardize and", "score structures before they", "feed a model or a database."],
         ref="Software on Zenodo", domain="chemaudit.naturalproducts.net",
         c=dict(field="#fbf2ec", tile="#f3e3d8", ink="#26140f", dim="#6d5850", accent="#c8411b", n="#2f55b5", o="#c8411b"),
         v=v_chemaudit),
    dict(slug="cms", name="Cheminformatics Microservice", logo="cms_w.png", logo_h=120,
         lines=["CDK, RDKit and Open Babel,", "unified behind one REST API", "and web app."],
         ref="Journal of Cheminformatics, 2023", domain="app.naturalproducts.net",
         c=dict(field="#0b1d3a", tile="#07152c", ink="#f4ead8", dim="#a9b4c8", accent="#f5c99b", n="#8fb0ff", o="#ff8a73"),
         v=v_cms),
    dict(slug="mortar", name="MORTAR", logo="mortar.png", logo_h=96, multiply=True,
         lines=["The MOlecule fRagmenTAtion", "fRamework: in silico", "fragmentation on the desktop."],
         ref="Journal of Cheminformatics, 2022", domain="felixbaensch.github.io/MORTAR",
         c=dict(field="#f5eedf", tile="#ece1c9", ink="#2f2618", dim="#6e604a", accent="#a86c0c", n="#2f55b5", o="#b4432f"),
         v=v_mortar),
    dict(slug="orthonym", name="Orthonym", logo="orthonym_w.png", logo_h=88,
         lines=["Deterministic structure-to-", "IUPAC names, each checked", "by an OPSIN round trip."],
         domain="orthonym.decimer.ai",
         c=dict(field="#0c0c0d", tile="#18181a", ink="#f2f2f0", dim="#a3a3a3", accent="#d8cfbd", n="#8fb0ff", o="#ff8a73"),
         v=v_orthonym),
    dict(slug="surge", name="SURGE", logo="surge.png", logo_h=96,
         lines=["A fast open-source chemical", "graph generator: every", "isomer of a formula."],
         ref="Journal of Cheminformatics, 2022", domain="github.com/StructureGenerator/surge",
         c=dict(field="#edf1fb", tile="#dfe6f8", ink="#14166b", dim="#4d5378", accent="#2a63c9", n="#2a63c9", o="#c23b2e"),
         v=v_surge),
    dict(slug="sherlock", name="Sherlock", logo="sherlock.png", logo_h=150,
         lines=["Computer-assisted structure", "elucidation: from NMR", "spectra to structures."],
         ref="Molecules, 2023", domain="github.com/Steinbeck-Lab/sherlock",
         c=dict(field="#14161c", tile="#1c1f27", ink="#ece8e0", dim="#a4a9b5", accent="#e6b65c", n="#8fb0ff", o="#ff8a73"),
         v=v_sherlock),
]


def footer():
    W, H = 1280, 140
    ink, dim = "#efe9dd", "#a39b8d"
    body = f"""<defs>{INVERT}</defs><rect width="{W}" height="{H}" rx="20" fill="#0f0f0e"/>
<image href="{png_uri('labmark.png')}" x="48" y="26" width="88" height="88" filter="url(#inv)"/>
{text(156, 78, "Brewed in Jena since 1997", 34, ink, SANS, 800, extra=' letter-spacing="-0.8"')}
{text(158, 106, "the logo is caffeine, for a reason", 15, dim, MONO, 500)}
{text(1232, 78, "cheminf.uni-jena.de", 18, dim, MONO, 500, "end")}
{text(1232, 106, "naturalproducts.net", 18, dim, MONO, 500, "end")}"""
    return svg(W, H, body, title="Brewed in Jena since 1997")


if __name__ == "__main__":
    for old in HERE.glob("*.svg"):
        old.unlink()
    out = {"hero.svg": hero(), "footer.svg": footer()}
    for p in POSTERS:
        v = p["v"](p["c"])
        visual, css = v if isinstance(v, tuple) else (v, "")
        out[f"poster-{p['slug']}.svg"] = poster(p, visual, css)
    for name, data in out.items():
        (HERE / name).write_text(data)
        print(f"{name:26} {len(data) // 1024:4} KB")
