# Pathways — a Boolean model of apoptotic commitment

**The question:** does degrading MDM2 and KRAS *together* kill an A549 cell
that survives either one alone?

A549 is the right background for asking it — TP53 wild-type, so there is
intact p53 to release, KRAS G12S, so the other arm is independently
druggable, and CDKN2A deleted. This is a separate strand from the structural
work; it shares the project's targets, not its methods.

Run the whole thing with `run_pipeline.sh` (or `run_pipeline.ps1` on
Windows). Python 3.8+, standard library only.

---

## How the model was built

### 1 · Build the network — `build_network.py`

Every edge comes from OmniPath, with its sign, its supporting resources and
its PMIDs carried along. We declared the 31 nodes we wanted
(`model_scope_v1.txt`) and fetched every qualified edge among them, rather
than growing the network outward from KRAS and MDM2 — that way the boundary
of the model is a decision someone can argue with, not a side effect of how
far a search happened to run.

An edge had to be directed, reported by at least 2 independent resources,
with at least 1 reference and curation effort of 5 or more. Everything
rejected is written out with the reason.

**203,001 interactions → 117 edges among 31 nodes.**

### 2 · Curate it — `curate_network.py`

Evidence thresholds tell you how well an edge is attested, not whether it
points the right way. Forty edges were fixed by hand, and they fall into
five recurring problems with any aggregated interactome:

- **reverse records** — `APAF1→CYCS`, when cytochrome c is released *from*
  mitochondria, not made by APAF1
- **signs that contradict their own annotations** — `AKT1→FOXO3` recorded as
  activation, when the phosphorylation it cites excludes FOXO3 from the nucleus
- **reciprocal sequestration** — one binding event recorded in both
  directions, which in Boolean form becomes a mutual-inhibition toggle
- **text-mining artefacts** — records resting only on automated readers,
  which routinely swap enzyme and substrate
- **transcriptional bypasses** — `TP53→CASP3`, letting p53 fire the
  executioner caspase without ever permeabilising the mitochondrion

PIK3CA and PIK3R1 became a single PI3K node (recorded as mutually
activating, which in Boolean form latches on forever). PARP1 came out — once
`CASP3→PARP1` is removed as a readout out-edge, nothing is left for it to
read.

**One edge was put back against our own thresholds.** `KRAS→PI3K` has a
single supporting resource, so step 1 rejected it — which left the whole
PI3K/AKT arm with no route from the perturbed node, and the model unable to
ask its own question. It is tagged `MANUAL_ADD` in the edge list, and that
tag should stay visible wherever the network is shown.

**117 edges → 29 nodes, 99 edges.** Every change is in
`network_changelog.csv` with its reason.

### 3 · Write the rules — `make_rules.py`

An edge list cannot tell you whether two activators combine with AND or OR.
Group a node's regulators and the rule comes down to two questions: is one
activator enough to switch it on, and is one inhibitor enough to switch it
off? Saying yes to both is the default, and it assumes nothing the edges do
not already claim.

So the only things worth justifying are the departures. Of 29 rules, 3 are
clamped inputs, 15 are left at the default, and 11 depart from it for
reasons of mechanism — the apoptosome needs both cytochrome c *and* Apaf-1,
BAK is already sitting at the membrane rather than waiting to be made, and
so on. `rules_justification.md` goes through them one at a time.

This is the only step with human judgement in the middle of it:
`rules_structural.txt` was written by hand from the draft, so it is an input
to the pipeline rather than something the pipeline regenerates. The script's
check mode is what keeps that honest — it re-derives each rule's behaviour
and compares it against the edge signs, and reports every edge no rule uses,
so a deliberate omission stays distinguishable from an accident.

### 4 · Find the attractors — `attractor_analysis.py`

Clamp the perturbed inputs, then repeatedly fix any node whose rule has
become constant no matter what the rest of the network does. In all four
conditions this determines all 29 nodes within six rounds.

That is a proof rather than a sample: **each condition has one fixed point,
it is the only attractor, and every one of the 67,108,864 starting states
reaches it.** No cycles, and the outcome does not depend on where you start.
Random-start sampling was run alongside and agreed every time.

### 5 · Screen the alternatives — `run_screen.py`

Five rules are defensible in two forms each, so all 32 combinations are
scored against criteria that were **frozen before any variant was run**
(`screen_criteria.md`). That ordering is the whole point of the exercise.

| | |
|---|---|
| satisfy every criterion | **4 of 32** |
| show true synergy — the combination kills, neither single agent does | **4 of 32** |
| apoptosis under MDM2 degradation alone | 0 of 32 |
| under KRAS degradation alone | 24 of 32 |
| under the combination | 28 of 32 |

The two fours are **different, overlapping sets** sharing two members — not
one result quoted twice.

---

## What came out, and what it is worth

Every passing variant needs R1 in form B and R5 in form B. Two rules stay
open, and both are reported rather than resolved:

- **R2, the PUMA rule.** The synergy argument turns on whether PUMA is
  p53-only or also FOXO3-driven, and this data cannot say: FOXO3 is off in
  every condition the criteria constrain, so both forms behave identically.
- **R1.** The screen needs form B, a coincidence detector requiring input
  from both arms. The binding literature leans the other way, since BAD is
  the canonical high-affinity BCL-XL antagonist. The tension is stated, not
  settled.

**This is model selection, not validation.** The criteria come from the same
proteomics that made us look at these rules in the first place, so passing
is a constraint on the model, not evidence it is right. Specifically:

- Nothing tests the KRAS arm — no dataset in hand perturbs KRAS, so the
  KRAS-degraded and both-degraded outcomes are predictions.
- Nutlin is not degradation. It blocks MDM2–p53 binding; the model clamps
  MDM2 to zero. MDM2 protein actually *rises* 2.84 log2 under nutlin,
  because p53 transcribes it — something a clamped model cannot show.
- BCL2L11, PMAIP1, MCL1 and MDM4 were never detected, so no criterion
  touches them.
- Bortezomib and seven other drugs in the panel are deliberately held back,
  as the only genuine out-of-sample test we have.

### One caveat that changes how to read R4

`BAX` (+0.70) and `BAK1` (+0.71) both rise clearly in the proteomics, and
neither was made a pass criterion. The reason is recorded in
`screen_criteria.md` before the screen ran: these nodes mean "active", not
"present", and the data pulls both ways — BAX and BAK1 rise while cytochrome
c does not move.

That exclusion decides R4. Either one as a criterion would cut the passing
set from four to two and fix R4 as form B. So "R4 is undetermined" means
*undetermined because we would not treat abundance as activity* — not that
the data is silent. We did not amend the criteria or re-run, because
changing a pass criterion after seeing results is exactly what freezing them
is meant to prevent.

---

## Files

| | |
|---|---|
| `model_scope_v1.txt` | the 31 declared nodes — the model's stated boundary |
| `decisions_final.csv` | 25 reviewed sign calls, each with a rationale |
| `run01/raw/` | the archived OmniPath responses |
| `final/` | step 1 output: the 117-edge network and its run manifest |
| `network_{edges,nodes,changelog}.csv` | the curated network and every change to it |
| `rules_structural.txt` | the 29 Boolean rules |
| `rules_justification.md` | why each rule is what it is |
| `screen_criteria.md` | the pre-registered criteria |
| `validation_nutlin_A549.csv` | the proteomics, mapped onto model nodes |
| `attractor_report.md` · `screen_report.md` · `screen_results.json` | results |
| `model.{bnd,cfg,bnet}` | MaBoSS and BoolNet exports |
| `CHECKSUMS.md` · `FINAL_REPORT_MANIFEST.md` · `checksums.txt` | checksums for every artefact |

**The OmniPath snapshot is pinned.** `run01/raw/` holds the raw API
responses, each with its URL, timestamp and checksum. OmniPath changes
continuously, so re-fetching would quietly give a different network and
invalidate everything downstream — the model is only reproducible against
this archive.

If you open the model in MaBoSS, note that `model.cfg` starts every node at
0, which is just what the exporter writes. Set the initial states to the
control attractor in `rules_justification.md` instead.
