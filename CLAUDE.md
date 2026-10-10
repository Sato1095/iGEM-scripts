# iGEM Thessaloniki 2026 — AETh3R

Context for Claude Code sessions on this repository. Read this first.

---

## What the project is

**AETh3R** is a peptide-PROTAC (p-PROTAC) therapeutic for non-small cell lung cancer
(NSCLC), targeting two proteins for degradation rather than inhibition:

- **MDM2** — the negative regulator of p53. Degrading it releases p53-mediated tumour
  suppression. This is the arm the dry lab designed binders for.
- **KRAS** — mutated in ~30% of NSCLC. Dry-lab work here is exploratory (see below).

The degrader recruits the **VHL** E3 ligase via a literature-derived peptide (pHIF1),
joined to the warhead by a glycine–serine linker.

This repo holds the **wiki** (Flask + Frozen-Flask, deployed by GitLab CI) and, under
`analysis/`, the **dry lab scripts**.

---

## Hard rules

### 1. The analysis scripts are a record, not a product

They were written for our own pipeline. They are **not** a software tool for others to
use, and the wiki must never claim otherwise. We are not entering the Software & AI
Village and not nominating for Best Software Tool.

### 2. Never change script behaviour

Every number on the Model page came out of these scripts. A refactor that alters output
silently invalidates the wiki.

When cleaning up a script, you may change:
- comments and docstrings
- variable and function names
- formatting, import order
- dead code removal

You may **not** change:
- any logic, threshold, filter or ordering
- default parameter values
- output format or file names

After editing, re-run the script where inputs still exist and confirm the output is
unchanged. If it cannot be re-run, say so in that script's README rather than assuming.

### 3. The three selection scripts stay separate

Stage 1, Stage 2 and Stage 3 could be merged into one cleaner script. **Do not.** The
fact that the funnel evolved through three attempts — each added because the previous
one missed something — is the team's main modelling contribution. Merging them erases
the evidence for the argument the wiki makes.

### 4. Figures and images

Images, icons and fonts go on `static.igem.wiki` via the uploads tool, referenced as
`{{ 'https://static.igem.wiki/teams/6118/wiki/<page>/<name>.avif' }}`. Raster images are
served as `.avif`; SVGs stay `.svg`. No external CDNs — everything must load from iGEM
infrastructure.

---

## The dry lab pipeline, in order

### Warhead: de novo MDM2 binders

1. **Interface analysis** — three hydrophobic sub-pockets mapped from 1YCR (the p53
   anchors Phe19, Trp23, Leu26).
2. **Generation** — RFdiffusion (500 backbones) → ProteinMPNN (4 sequences each) →
   AlphaFold Multimer. 2,000 candidates.
3. **Three-stage selection funnel** (see `analysis/`):
   - *Stage 1* — composite Z-score of AlphaFold confidence metrics. Picked **Z1 (42.3)**.
   - *Stage 2* — confidence filter + pharmacophore flags. Picked **Z2 (267.3)**.
   - *Stage 3* — interface geometry and chemistry from atomic coordinates. Picked
     **Z3 (131.2)**.

**The point of the sequence:** Stage 2's top pick, Z2, carries net charge −8 — chemically
implausible for a hydrophobic groove, but the confidence metrics rank it highly because
they score predicted geometry and are blind to chemistry. Stage 3 was written to catch
exactly that.

### Recruiter

**pHIF1** (LA-Hyp-Y-Hle-P), from Chen et al. Contains two non-canonical residues —
4-hydroxyproline and homoleucine — which constrain every downstream tool.

### Assembly and validation

Four constructs: PMI-PROTAC, PZ1, PZ2, PZ3. Simulated as ternary complexes with MDM2 and
**VBC** (VHL + ElonginB + ElonginC) for 100 ns.

**The headline result:** PZ3 carries the best MDM2 warhead and makes the most contacts
with the target, but forms the *weakest* induced MDM2–VHL interface. PZ1 is the opposite
and is the leading degrader candidate. Target binding and degradation efficiency are
separable.

---

## Canonical naming — use these, consistently

| Use | Not |
|---|---|
| AETh3R | AETHER, AETH3R, Aether |
| Z1 (42.3) · Z2 (267.3) · Z3 (131.2) | 42.3, 267.3, 131.2 bare; "Z4"; 130.3 |
| PZ1 · PZ2 · PZ3 · PMI-PROTAC | — |
| VBC | "the VHL complex" when VBC is meant |
| Dr. Alessio Ciulli | Alessandro Ciulli |

---

## Numbers that appear on the wiki

If a cleanup changes any of these, something has gone wrong.

- 2,000 candidates → 1,275 fill 3/3 pockets → 409 pass all gates → **18** meet or exceed PMI
- Z3: 42.0 apolar Å²/res, 49 contacts, 3/3 pockets, net charge −1, pLDDT 89.5
- Z2: net charge −8 (rejected at the chemistry gate)
- PMI benchmark: 40.0 apolar Å²/res, 35 contacts, 3/3 pockets
- Induced MDM2–VHL contacts: PZ1 ~450–600, PMI-PROTAC ~300–380, PZ2 ~200–300, PZ3 ~130–200
- Pathways: 117 edges after filtering → **29 nodes, 99 edges** after curation; 29 rules
  (12 default, 17 overrides); 1 fixed point per condition from 6.7 × 10⁷ states;
  32 rule variants, 4 satisfy all criteria, **4 of 32** show true synergy (28 kill, but
  24 of those are KRAS-alone)

---

## Open issues — do not paper over these

1. **KRAS arm.** The wiki says AETh3R degrades KRAS *and* MDM2, but the dry lab only
   designed MDM2 binders. The KRAS work is an unstapled SOS1 helix that gave a clear
   negative result (ipTM 0.27; no RMSD plateau) and is framed as future work. The Design
   page overview may still overstate this — check before claiming a KRAS degrader exists.

2. **MDMX PROTAC status.** Model page says "under construction for simulation"; an earlier
   Design draft said "successfully assembled". Settle which is true.

3. **Ternary scope.** We simulated MDM2·PROTAC·**VBC**. We did *not* model the Cullin-2 /
   Rbx1 scaffold. Any text claiming the full CRL2^VHL assembly is wrong.

4. **MM/PBSA** showed negative cooperativity for the constructs tested — reported honestly
   on the wiki. Do not quietly drop it.

5. **Single trajectories.** One 100 ns run per system, no replicates. Stated as a
   limitation; keep it stated.

---

## Wiki conventions

- Pages live in `wiki/pages/*.html`, extend `layout.html`, use
  `<div class="content-layout">` with a `<nav class="toc">` sidebar that `toc.js` fills
  from the headings.
- Collapsible blocks: `.accordion > .acc-item > .acc-header` + `.acc-body`, toggled by
  `description_images/accordion.js` via an `is-open` class.
- References: one numbered list per page, **in order of first appearance**, deduplicated,
  each citation a link (`[<a href="#ref-N" class="cite">N</a>]`), list wrapped in a
  `.ref-accordion`. `references.js` opens the accordion when a citation is clicked.
- Per-page CSS: `safety.css`, `sustainability.css`, `inclusivity.css`, `design.css`.
  Shared vocabulary across them: `.*-lede` (section standfirst), `.*-note` (a result or
  caveat lifted out of the prose), `.*-todo` (unwritten, deliberately loud).

### Known CSS trap

Full-bleed bands use `width: 100vw; margin-left: calc(-50vw + 50%)`. That centres on the
*containing block*, not the viewport — so inside `.content-body` (which sits right of the
18% TOC sidebar) every band is displaced. If a band looks offset, this is why. The fix is
structural: move the band outside `.content-layout`.

---

## Still to do

- Fill the Design page placeholders: integrin αvβ6 expression, antibody selection, LNP
  carrier choice.
- Replace figure placeholders marked `TODO` in the HTML.
- Resolve the open issues above.
