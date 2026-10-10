# How the Boolean rules were written

Model: the curated network — 29 nodes, 99 signed edges — in A549.
Rule set: `rules_structural.txt`.

This document covers everything up to the point where the rules meet the
proteomics. It ends by listing the places they disagree. What happens to
those disagreements is `screen_criteria.md` and the screen itself; **no
result from the screen appears here**, because every judgement below was
made before it.

---

## 1. What a signed edge list cannot tell you

An edge asserts that the source can change the target. It does not say how
the inputs to a node combine. Group a node's regulators into activators and
inhibitors and the rule is fully determined by exactly two questions:

> **Is one activator enough to switch the node on?**
> **Is one inhibitor enough to switch it off?**

Answer yes to both and you get `(any activator) AND NOT (any inhibitor)`.
That is the default convention, and it is the weakest assumption consistent
with a signed edge list — it adds nothing the edges do not already claim.

So the default is the starting point everywhere, and the only thing worth
documenting is where we answered *no*, and why. That is what the tiers
below are.

## 2. The draft, and where it flagged itself

`make_rules.py` writes the draft under that convention:

```
python3 make_rules.py --edges network_edges.csv --nodes network_nodes.csv \
    --fixed KRAS MDM2 MDM4 --out rules_draft.txt
```

It marks the nodes where the default is most likely to be wrong — where
several inputs have to be combined, or where an activator's mechanism names
a physical complex. Eight came back flagged: **BAX, BCL2L1, BCL2L11, CASP9,
CDKN1A, MCL1, PI3K, TP53**.

The flags are a prompt to check, not a verdict. Checking them changed four
(BCL2L1, CASP9, PI3K, TP53) and left four alone (BAX, BCL2L11, CDKN1A,
MCL1). Seven further nodes were changed that the draft had not flagged —
those came from reading the edges themselves rather than from the input
count.

Of the 29 rules: **3 are clamped inputs, 15 are the unmodified default, and
11 depart from it.**

## 3. Tier 1 — the default, unchanged

Fifteen rules. Both questions answered yes; nothing added.

> AKT1, APAF1, BAD, BAX, BCL2L11, CASP3, CASP7, CDKN1A, CYCS, DIABLO,
> MAPK1, MAPK3, MCL1, PTEN, XIAP

Three are worth pointing at, because they look like decisions and are not:

- `BAD = not (AKT1 or RAF1)` — two kinases, either sufficient. Worked out
  by hand as a truth table first; the default produces the same thing.
- `BAX = (TP53 or BCL2L11 or GSK3B or BBC3) and not (BCL2L1 or MCL1 or AKT1)`
  — the guards enter under OR inside the NOT, so **both** must be absent for
  BAX to fire. That is the stringent reading, and it is what the default
  gives. No intervention was needed to get it.
- `MCL1 = (MAPK1 or MAPK3) and not (BCL2L11 or GSK3B or PMAIP1 or TP53 or BBC3)`
  — ERK stabilises; any one BH3-only protein neutralises. MCL1 is therefore
  lost under either single agent, which is why it cannot be the node that
  distinguishes them.

`CDKN1A` is also here: `(TP53 or FOXO3) and not (GSK3B or MDM2 or MAPK1)` is
exactly the default. It was flagged by the draft, checked, and kept.

## 4. Tier 2 — structural and semantic overrides

Eleven rules where one of the two questions is answered *no*, for reasons
that come from what the node means or how the complex works. Every one of
these would be made whatever the proteomics said.

| rule | what changed | why |
|---|---|---|
| `CASP9 = CYCS and APAF1 and not (…)` | AND between activators | The apoptosome is a stoichiometric complex. Under OR, cytochrome c activates caspase-9 with no Apaf-1 present and the execution module becomes decorative. A composition constraint, not a preference. |
| `FOXO3 = not (AKT1 or MAPK1)` | activator requirement dropped | FOXO3 is constitutively expressed; regulation is nuclear exclusion by AKT. Its only activator edge is a weak p53 transcriptional record (effort 19), which under the default would make FOXO3 p53-dependent. |
| `GSK3B = not (AKT1 or MAPK1 or MAPK3)` | activator requirement dropped | Constitutively active kinase held off by AKT (effort 128). Its only activator edge is TP53 at effort 6. |
| `BAK1 = not (BCL2L1 or MCL1)` | activator requirement dropped | BAK is constitutively resident at the outer membrane; control is guard sequestration, not synthesis. |
| `PMAIP1 = TP53` | inhibitor dropped | `BCL2L1 ⊣ PMAIP1` is the reverse record of a sequestration event already represented as antagonist→guard. Keeping both directions creates a mutual-inhibition toggle. |
| `BBC3 = TP53 or FOXO3` | inhibitor dropped | `MCL1 ⊣ BBC3` is the same reciprocal-record problem as PMAIP1 above. **Both activators are kept.** |
| `RAF1 = KRAS`, `MAP2K1 = RAF1` | ERK feedback and AKT ⊣ RAF1 dropped | Modulatory phosphorylations that tune amplitude and adaptation. Read as Boolean switches they shut the MAPK arm off at baseline, since AKT is constitutively active in A549, and create a sustained oscillator. Neither matches a cell with constitutive ERK activity. |
| `PI3K = KRAS and not (TP53 or PTEN)` | FOXO3 → PI3K dropped | Creates an AKT/FOXO3 negative-feedback oscillator that made the KRAS-degraded condition non-deterministic. |
| `TP53 = not (MDM2 and MDM4)` | activators dropped; AND between suppressors | ERK-dependent p53 would confound the KRAS comparison. Under OR with MDM4 clamped ON, p53 could never activate and MDM2 degradation would be inert — a degenerate model. |

**A correction worth recording.** An earlier version of this document wrote
the BBC3 rule as `BBC3 = TP53`, grouped with PMAIP1 under "reciprocal
records removed". That was wrong, and in a way that mattered. The default
for BBC3 is `(TP53 or FOXO3) and not MCL1`; removing the reciprocal `MCL1`
record gets you to `TP53 or FOXO3` and no further. Dropping `FOXO3 → BBC3`
is a second, independent change that the reciprocal-record argument does not
license — and it is one of the alternatives the screen exists to test. The
rule set has always carried both activators; it was the description that was
wrong.

### One network change made on the same grounds

`KRAS → PI3K` was restored as a documented manual addition. It is
single-resource in OmniPath and was excluded by the evidence threshold,
which left the entire PI3K/AKT arm disconnected from the perturbed node.
Without it, KRAS degradation cannot reach AKT, BAD, FOXO3 or GSK3B and the
model cannot address the question it was built for. It is flagged
`provenance=MANUAL_ADD` in the edge list.

## 5. Checking the rules against the proteomics

Up to here nothing has looked at data. The rules now get compared with the
MDM2-inhibition proteomics in `validation_nutlin_A549.csv`, node by node:
the model
"changes" a node if its state differs between the control and MDM2-degraded
attractors, and that direction is compared with the measured call.

Fourteen nodes agree. Four were not detected (PMAIP1, BCL2L11, MCL1, MDM4).
**Ten disagree:**

| node | model | observed | |
|---|---|---|---|
| `BCL2L1` | down | **up** +0.78 | the guard survives MDM2 inhibition; the model clears it |
| `CDKN1A` | unchanged | **up** +2.13 | p21 cannot rise while ERK is active |
| `GSK3B` | unchanged | **up** +0.36 | held off by MAPK1, which C1 and C2 both require to be on |
| `CASP9` | unchanged | **down** −0.42 | already off in the model; a fall from off is not representable |
| `CYCS` | up | unchanged | |
| `DIABLO` | up | unchanged | |
| `PTEN` | up | unchanged | |
| `AKT1` | down | unchanged | |
| `XIAP` | down | unchanged | |
| `MDM2` | down | **up** +2.84 | nutlin blocks MDM2–p53 binding; the model clamps MDM2 to 0, so the rise cannot be represented at all |

### What was done with each

**Two became screened rules.** `BCL2L1` and `CDKN1A` are both places where a
defensible alternative answer to one of the two questions would remove the
disagreement, so both alternatives go into the screen rather than being
chosen here — R1 and R3.

**Two more rules go in alongside them.** `TP53` (R5) had previously been set
to the AND form by hand with the Nutlin data in view; it is returned to the
screen as a free variant so the choice is made by the stated procedure
instead. The guard logic at BAX and BAK1 (R4) goes in because it is coupled
to R1: if the alternative BCL2L1 rule keeps BCL-XL present, whether the
guards still block the effectors is no longer something the default can be
assumed to settle. `BBC3` (R2) goes in as the remaining open question from
tier 2 — whether PUMA is p53-only or also FOXO3-driven decides whether the
KRAS arm can reach the mitochondrion at all.

**The rest are not screened, and are open.** Five of them —
`CYCS`, `DIABLO`, `PTEN`, `AKT1` and `XIAP` — are cases where the model moves
a node and the measurement does not. `GSK3B` is the reverse. These are not
measurement failures and the screen does not address them; they are recorded
here as standing disagreements. Two have a structural cause worth stating:

- `GSK3B = not (AKT1 or MAPK1 or MAPK3)` cannot be on while MAPK1 is on, and
  a viable control cell requires MAPK1 on. The observed rise is
  unreachable for any rule set built on this network.
- The same holds for `APAF1 = (TP53 or CYCS) and not BCL2L1` once BCL-XL is
  required to survive: APAF1 is forced off. APAF1 is measured up (+0.82).
  It does not appear in the mismatch table above only because the structural
  rules clear BCL-XL, which the data says they should not.

`MDM2` is a category of its own: nutlin and degradation are different
perturbations, and a clamped node cannot show the rise that p53-driven
transcription produces. No rule set can fix this.

**Abundance is not activity.** Several of these disagreements sit on that
line. `BAX` and `BAK1` are measured up while `CYCS` is unchanged — more
effector protein, no reported cytochrome c release. A Boolean node that
means "active" cannot be scored against a measurement that means "present",
and where the two readings conflict the model is being asked a question the
data cannot answer.

## 6. Initial conditions

**The outcome does not depend on them.** Each condition was run from 120
random initial states across all 29 nodes and converged to a single
attractor in 120/120. The perturbation, not the starting point, determines
the result, and no cyclic attractors were found.

Control attractor, ON: KRAS, PI3K, AKT1, RAF1, MAP2K1, MAPK1, MAPK3, MCL1,
BCL2L1, XIAP, MDM2, MDM4. Everything else OFF.

This is itself a check on the rules. It is a proliferating A549 cell with
constitutive KRAS/ERK signalling, active AKT, both mitochondrial guards
intact, XIAP restraining the caspases, p53 held off by MDM2, and no caspase
activity. It was not imposed — it emerges, and a rule set that failed to
produce it would be wrong before any perturbation was tested.

**Fixed nodes.** KRAS and MDM2 are clamped (1 = present, 0 = degraded)
because they are the perturbations. MDM4 is clamped ON because A549
expresses it and MDM2-directed degraders do not target it; its in-edges were
removed for the same reason.

## 7. What a reader should take from this

The control attractor and the single-attractor result are properties of the
network and the convention, not of any choice made to fit data. The tier 2
overrides are mechanism, defensible on their own terms, and were settled
before the proteomics was opened.

What the proteomics then did was identify ten disagreements. Two of them,
plus three related rules, were turned into the alternatives the screen
tests. The rest remain open and are listed above rather than resolved.

That the rules themselves were not fitted does not make the screen blind:
the five rules it varies were chosen with the disagreements in view. The
screen is a consistency check against data already seen, and
`screen_criteria.md` states the same thing before any variant is run.
