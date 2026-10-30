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
    sub.populate(text=string.printable + "·→−Δδ")
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
# One image per tool: logo, tagline, paper, link and an evidence strip on the left,
# the tool's Anatomy machine (machines/how-<slug>.svg) on the right, on one field.

POSTER_W, POSTER_H = 1280, 640
SX, SY, SW, SH = 56, 392, 528, 212  # the evidence strip under the link
MX, MY, MW, MH = 600, 16, 664, 608  # the machine


def machine(slug):
    """The machine SVG without its own card, nested on the right; returns (svg, field colour, label)."""
    s = (HERE / "machines" / f"how-{slug}.svg").read_text()
    root = re.search(r"<svg\b[^>]*>", s).group(0)
    w, h = map(float, re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', root).group(1, 2))
    label = re.search(r'aria-label="([^"]*)"', root)
    field = re.findall(r"--anatomy-card:\s*(#[0-9a-fA-F]{6})", s)[-1]
    inner = s[s.index(root) + len(root):s.rindex("</svg>")]
    # The first rect is the machine's own card; the poster field replaces it.
    inner = re.sub(r'<rect\b[^>]*class="[a-z]+-(?:card|frame)"[^>]*/>', "", inner, count=1)
    k = min(MW / w, MH / h)
    x, y = MX + (MW - w * k) / 2, MY + (MH - h * k) / 2
    nested = (f'<svg x="{x:.1f}" y="{y:.1f}" width="{w * k:.1f}" height="{h * k:.1f}" viewBox="0 0 {w:g} {h:g}" '
              f'overflow="visible"><g class="iso-svg iso" data-theme="light">{inner}</g></svg>')
    return nested, field, (label.group(1) if label else "")


def poster(p, strip):
    c = p["c"]
    mach, field, label = machine(p["slug"])
    lines = "".join(f'<tspan x="{SX}" dy="{0 if i == 0 else 38}">{esc(l)}</tspan>' for i, l in enumerate(p["lines"]))
    ref = text(SX, 318, esc(p["ref"]), 18, c["dim"]) if p.get("ref") else ""
    body = f"""<defs><clipPath id="pp-r"><rect width="{POSTER_W}" height="{POSTER_H}" rx="20"/></clipPath></defs>
<g clip-path="url(#pp-r)">
<rect width="{POSTER_W}" height="{POSTER_H}" fill="{field}"/>
{logo(p['logo'], SX, 44, 440, p.get('logo_h', 84), p.get('multiply', False))}
<text x="{SX}" y="196" font-family="{SANS}" font-size="30" fill="{c['ink']}" letter-spacing="-0.4">{lines}</text>
{ref}
{text(SX, 354, esc(p['domain']), 19, c['accent'], MONO, 500)}
{arrow_ne(SX + len(p['domain']) * 11.4 + 12, 340, c['accent'], 13)}
{strip}
{mach}
</g>"""
    return svg(POSTER_W, POSTER_H, body, "", p["name"] + ": " + " ".join(p["lines"]) + " " + label)


def tile(x, y, w, h, fill):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="{fill}"/>'


def s_coconut(c):
    picks = [("quercetin", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
             ("nicotine", "CN1CCC[C@H]1c1cccnc1"),
             ("menthol", "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O"),
             ("salicin", "OC[C@H]1O[C@@H](Oc2ccccc2CO)[C@H](O)[C@@H](O)[C@@H]1O")]
    out, cw = [], SW / 4
    for i, (name, smi) in enumerate(picks):
        x = SX + i * cw
        out.append(tile(x + 4, SY, cw - 8, SH, c["tile"]))
        out.append(place(round(x + 8, 1), SY + 8, mol(smi, int(cw - 16), 160, c["ink"], c["n"], c["o"], bond=15, lw=1.5, font=0.45)))
        out.append(text(round(x + cw / 2, 1), SY + SH - 18, name, 14, c["dim"], MONO, 500, "middle"))
    return "".join(out)


def step_tiles(c, steps):
    out, cw = [], SW / len(steps)
    for i, (head, l1, l2) in enumerate(steps):
        x = SX + i * cw
        out.append(tile(round(x + 4, 1), SY, round(cw - 8, 1), SH, c["tile"]))
        out.append(text(round(x + 22, 1), SY + 56, str(i + 1), 40, c["accent"], SANS, 800))
        out.append(text(round(x + 22, 1), SY + 104, esc(head), 21, c["ink"], SANS, 800, extra=' letter-spacing="-0.3"'))
        out.append(text(round(x + 22, 1), SY + 140, esc(l1), 15, c["dim"]))
        out.append(text(round(x + 22, 1), SY + 162, esc(l2), 15, c["dim"]))
    return "".join(out)


def s_decimer(c):
    return step_tiles(c, [("Segment", "finds each structure", "drawn on the page"),
                          ("Classify", "flags crops that are", "chemical structures"),
                          ("Recognise", "the Image Transformer", "writes the SMILES")])


def s_chemaudit(c):
    # Score rule from ValidationEngine._calculate_score; bands from the ChemAudit frontend.
    x0, x1, by = SX + 24, SX + SW - 24, SY + 132
    sx = lambda v: x0 + (x1 - x0) * v / 100
    out = [tile(SX, SY, SW, SH, c["tile"]),
           text(SX + 24, SY + 40, "score = 100 − 50·critical − 20·error − 5·warning", 16, c["ink"], MONO, 500),
           text(SX + 24, SY + 70, "27 checks per structure", 15, c["dim"])]
    for a, b, col in ((0, 50, "#d63a2f"), (50, 80, "#e0a020"), (80, 100, "#2e9e6b")):
        out.append(f'<rect x="{sx(a):.1f}" y="{by}" width="{sx(b) - sx(a) - 2:.1f}" height="14" rx="3" fill="{col}"/>')
    for t in (0, 50, 80, 100):
        out.append(text(round(sx(t), 1), by + 38, str(t), 13, c["dim"], MONO, 500, "middle"))
    for v, name, dy in ((85, "salt 85", 20), (100, "parent 100", 40)):
        out.append(f'<path d="M{sx(v):.1f} {by - 4}l-6 -10h12z" fill="{c["ink"]}"/>')
        out.append(f'<line x1="{sx(v):.1f}" y1="{by - 14}" x2="{sx(v):.1f}" y2="{by - dy + 4}" stroke="{c["ink"]}" stroke-width="1"/>')
        out.append(text(round(sx(v) - 6, 1), by - dy, name, 13, c["ink"], MONO, 500, "end"))
    return "".join(out)


def s_cms(c):
    # Live response from GET api.naturalproducts.net/latest/chem/descriptors (caffeine, rdkit).
    k, v, d = "#f5c99b", "#f4ead8", "#a9b4c8"
    lines = [
        [("GET ", d), ("/latest/chem/descriptors", v), ("?toolkit=", d), ("rdkit", k)],
        [("    &smiles=", d), ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", k)],
        [('{ "molecular_weight"', k), (": ", d), ("194.19", v), (",", d)],
        [('  "alogp"', k), (": ", d), ("-1.03", v), (",", d)],
        [('  "topological_polar_surface_area"', k), (": ", d), ("61.82", v), (",", d)],
        [('  "qed_drug_likeliness"', k), (": ", d), ("0.54", v), (", ... }", d)],
    ]
    out = [tile(SX, SY, SW, SH, "#16294d")]
    for i, parts in enumerate(lines):
        spans = "".join(f'<tspan fill="{col}">{esc(t)}</tspan>' for t, col in parts)
        out.append(f'<text x="{SX + 22}" y="{SY + 40 + i * 30 + (12 if i > 1 else 0)}" font-family="{MONO}" '
                   f'font-weight="500" font-size="14.5" xml:space="preserve">{spans}</text>')
    return "".join(out)


def s_mortar(c):
    # Algorithms and their default outputs, from MORTAR's *Fragmenter.java defaults.
    algos = [("Ertl", "algorithm", "functional groups", "and alkanes", "#d9534f"),
             ("Sugar", "Removal Utility", "aglycone", "and sugars", "#2e9e6b"),
             ("Scaffold", "Generator", "scaffolds", "", "#3e63dd"),
             ("MolWURCS", "", "sugar", "fragments", "#9b59b6")]
    out, cw = [], SW / 4
    for i, (a, b, o1, o2, col) in enumerate(algos):
        x = round(SX + i * cw + 4, 1)
        out.append(tile(x, SY, round(cw - 8, 1), SH, c["tile"]))
        out.append(f'<rect x="{x}" y="{SY}" width="{cw - 8:.1f}" height="10" rx="5" fill="{col}"/>')
        out.append(text(x + 14, SY + 50, esc(a), 17, c["ink"], SANS, 800, extra=' letter-spacing="-0.3"'))
        out.append(text(x + 14, SY + 74, esc(b), 15, c["ink"]))
        out.append(text(x + 14, SY + 146, esc(o1), 13, c["dim"]))
        out.append(text(x + 14, SY + 166, esc(o2), 13, c["dim"]))
        out.append(text(x + 14, SY + 120, "returns", 13, col, MONO, 500))
    return "".join(out)


def s_orthonym(c):
    smi = "COc1cc(C=O)ccc1O"
    return "".join([tile(SX, SY, SW, SH, c["tile"]),
                    place(SX + 8, SY + 6, mol(smi, 200, 200, c["ink"], c["n"], c["o"], bond=26, lw=2.2, font=0.6)),
                    text(SX + 220, SY + 84, "4-hydroxy-3-", 26, c["ink"], SANS, 800, extra=' letter-spacing="-0.5"'),
                    text(SX + 220, SY + 118, "methoxybenzaldehyde", 26, c["ink"], SANS, 800, extra=' letter-spacing="-0.5"'),
                    text(SX + 220, SY + 154, "vanillin: structure in, name out", 15, c["dim"])])


def s_surge(c):
    # From the Surge README: `surge -u C10H16O` makes the 452458 isomers of C10H16O in 0.1 s (2021 PC).
    return "".join([tile(SX, SY, SW, SH, c["tile"]),
                    f'<text x="{SX + 24}" y="{SY + 42}" font-family="{MONO}" font-weight="500" font-size="18" xml:space="preserve">'
                    f'<tspan fill="{c["dim"]}">$ </tspan><tspan fill="{c["ink"]}">surge -u C10H16O</tspan></text>',
                    text(SX + 22, SY + 122, "452,458", 64, c["accent"], SANS, 800, extra=' letter-spacing="-2"'),
                    text(SX + 300, SY + 100, "isomers in 0.1 s", 19, c["ink"]),
                    text(SX + 300, SY + 124, "on a 2021 PC", 15, c["dim"]),
                    text(SX + 24, SY + 178, "nauty · canonical augmentation · no duplicates", 15, c["dim"])])


def s_sherlock(c):
    # Ranking rule from casekit FilterAndRank.rank, which Sherlock calls; PyLSD from the README.
    return "".join([tile(SX, SY, SW, SH, c["tile"]),
                    place(SX + 10, SY + 6, mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C", 190, 200, c["ink"], c["n"], c["o"], bond=22, lw=2, font=0.55)),
                    text(SX + 214, SY + 62, "PyLSD builds candidates,", 19, c["ink"]),
                    text(SX + 214, SY + 92, "ranked by matched signals,", 19, c["ink"]),
                    text(SX + 214, SY + 122, "then by mean |Δδ|", 19, c["ink"]),
                    text(SX + 214, SY + 162, "caffeine, 8 carbon peaks", 15, c["dim"])])


POSTERS = [
    dict(slug="coconut", name="COCONUT", logo="coconut.png", logo_h=80,
         lines=["The COlleCtion of Open Natural", "prodUcTs: deposit, curate", "and reuse natural product data."],
         ref="Nucleic Acids Research, 2025", domain="coconut.naturalproducts.net",
         c=dict(tile="#ece0cf", ink="#2e2118", dim="#6b5848", accent="#8a4f26", n="#2f55b5", o="#b4432f"), s=s_coconut),
    dict(slug="decimer", name="DECIMER.ai", logo="decimer.png", logo_h=80,
         lines=["Extract chemical structures", "from PDFs and images, and", "turn them into SMILES."],
         ref="Nature Communications, 2023", domain="decimer.ai",
         c=dict(tile="#dde9e0", ink="#17302a", dim="#4d6358", accent="#2f5d4a", n="#2f55b5", o="#b4432f"), s=s_decimer),
    dict(slug="chemaudit", name="ChemAudit", logo="chemaudit.png", logo_h=88,
         lines=["Validate, standardize and score", "structures before they feed", "a model or a database."],
         ref="Software on Zenodo", domain="chemaudit.naturalproducts.net",
         c=dict(tile="#f3e3d8", ink="#26140f", dim="#6d5850", accent="#c8411b", n="#2f55b5", o="#c8411b"), s=s_chemaudit),
    dict(slug="cms", name="Cheminformatics Microservice", logo="cms.png", logo_h=100,
         lines=["CDK, RDKit and Open Babel,", "unified behind one REST API", "and web app."],
         ref="Journal of Cheminformatics, 2023", domain="app.naturalproducts.net",
         c=dict(tile="#dde5f3", ink="#0f1e3a", dim="#4a5a78", accent="#b6531a", n="#2f55b5", o="#b4432f"), s=s_cms),
    dict(slug="mortar", name="MORTAR", logo="mortar.png", logo_h=80, multiply=True,
         lines=["The MOlecule fRagmenTAtion", "fRamework: in silico fragmentation", "of whole compound sets."],
         ref="Journal of Cheminformatics, 2022", domain="felixbaensch.github.io/MORTAR",
         c=dict(tile="#ece1c9", ink="#2f2618", dim="#6e604a", accent="#8f5c0a", n="#2f55b5", o="#b4432f"), s=s_mortar),
    dict(slug="orthonym", name="Orthonym", logo="orthonym.png", logo_h=70,
         lines=["Deterministic structure-to-IUPAC", "names, each checked by an", "OPSIN round trip."],
         domain="orthonym.decimer.ai",
         c=dict(tile="#e4e4ea", ink="#18181b", dim="#55565e", accent="#5a3ddb", n="#2f55b5", o="#b4432f"), s=s_orthonym),
    dict(slug="surge", name="SURGE", logo="surge.png", logo_h=80,
         lines=["A fast open-source chemical", "graph generator: every", "isomer of a formula."],
         ref="Journal of Cheminformatics, 2022", domain="github.com/StructureGenerator/surge",
         c=dict(tile="#dfe6f8", ink="#14166b", dim="#4d5378", accent="#2a63c9", n="#2a63c9", o="#c23b2e"), s=s_surge),
    dict(slug="sherlock", name="Sherlock", logo="sherlock.png", logo_h=120,
         lines=["Computer-assisted structure", "elucidation: from NMR", "spectra to structures."],
         ref="Molecules, 2023", domain="github.com/Steinbeck-Lab/sherlock",
         c=dict(tile="#efe4cf", ink="#1b2830", dim="#4c5a64", accent="#94600f", n="#2f55b5", o="#b4432f"), s=s_sherlock),
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
    out = {"hero.svg": hero(), "footer.svg": footer()}
    for p in POSTERS:
        out[f"poster-{p['slug']}.svg"] = poster(p, p["s"](p["c"]))
    for name, data in out.items():
        (HERE / name).write_text(data)
        print(f"{name:26} {len(data) // 1024:4} KB")
