# Reproducing a Published Boolean Model of Blood Stem Cell Fate

## What This Project Is

Herrera et al. (2024) built a 21-node gene regulatory network describing how hematopoietic stem cells (HSCs) choose between becoming megakaryocyte-erythroid progenitors (MEPs), granulocyte-monocyte progenitors (GMPs) or common lymphoid progenitors (CLPs). Their analysis was done in R using the BoolNet package.

My goal was to rebuild their Boolean model independently in Python, from the rules published in the paper and without using their analysis code, and check whether my version reproduces their results. I used the authors' rules file and BoolNet only as references to check my work against.

## The Model

Each of the 21 nodes (transcription factors, signalling molecules, ROS and metabolic states) is either ON (1) or OFF (0). Twenty nodes have a logical rule that decides their next state from the current state of their regulators. The twenty-first, environmental oxygen, is the network's input and stays fixed.

One naming point that matters: `oxygen` is environmental oxygen, the input. `o2` is superoxide (O₂·⁻), a reactive oxygen species. They are separate nodes.

## Getting the Rules Right

I transcribed the rules from Table 1 into Python. Copying from the PDF scrambled the logic symbols, so I cross-checked every rule against the authors' own BoolNet rules file in their repository. This caught two transcription errors: a missing NOT in the HIF1 rule, and incorrect logic in the PU.1 rule. The final rules below match the authors' file exactly.

```python
runx1  = (pu1 or gata2 or runx1) and not ikzf1
meis1  = (meis1 or runx1 or pu1) and not gfi1
hif1   = (not oxygen and o2 and not (p53 or foxo3 or ampk)) or (not oxygen and meis1)
foxo3  = (ampk and (hif1 or p53)) and not (akt or cebpa)
p53    = hif1 and not (akt or gata1 or pu1)
gata2  = (gata2 or p53) and not (gata1 or pu1 or gfi1)
gata1  = (gata1 or gata2 or akt or runx1) and not (pu1 or p53 or ikzf1)
pu1    = (pu1 or runx1 or (cebpa and ikzf1)) and not (gata1 or gata2 or gfi1)
cebpa  = (pu1 and runx1) and not mef2c
ikzf1  = (mef2c or ikzf1 or runx1) and not (cebpa or pu1 or gata1)
gfi1   = (ikzf1 or cebpa) and not (pu1 or p53)
mef2c  = pu1 and not cebpa
mtor   = akt and not ampk and not p53
ampk   = (not akt or not oxphos or not hif1) and p53
akt    = h2o2 or not (p53 or foxo3)
h2o2   = (gata1 or hif1 or pu1 or oxphos or (sod and o2)) and not antiox
o2     = not sod or oxphos
sod    = p53 or foxo3
antiox = p53 or foxo3
oxphos = (mtor or akt) and not (foxo3 or hif1)
```

In the code, these are packaged into a single `update_state()` function that computes every node's next value from the current state, so no rule ever reads a value that has already been updated in the same step.

## Validation 1: Checking the Published Cell States

Figure 3 of the paper shows 13 stable cell states (fixed-point attractors): 1 HSC, 6 MEP, 4 GMP and 2 CLP. A fixed point, by definition, comes back unchanged after one update.

I encoded all 13 columns from the figure and ran each through `update_state()` once. Eleven came back unchanged. Two, GMP1 and GMP3, did not: in both, Runx1 flipped from OFF to ON.

## Validation 2: Finding the Fixed Points from Scratch

Checking known states can't reveal missing or extra ones, so I searched all 2²¹ = 2,097,152 possible states and kept every state that maps to itself in one update.

**Result**: exactly 13 fixed points.** Eleven match Figure 3 exactly. The other two match GMP1 and GMP3 except at a single node, Runx1.

This search finds only fixed points, not cycles, since a cycle never maps to itself in a single step. Cycles are identified in Validation 3.

## A Discrepancy Between the Rules and Figure 3

The published Runx1 rule is:

```
Runx1 = (PU.1 or GATA2 or Runx1) and not Ikzf1
```

Every GMP state has PU.1 ON and Ikzf1 OFF, so under this rule Runx1 must be ON in all GMP attractors. Figure 3 shows it OFF in two of them.

The authors' analysis notebook writes its results to PDF reports that aren't included in the repository, so there was no saved attractor list to compare against. To settle the question, I ran the authors' own rules file through the same R package they used (BoolNet, `getAttractors()` with synchronous updating; output files are in the `herrera_check/` folder). Their software returns Runx1 active in all four GMP attractors, matching my independent Python implementation exactly.

Figure 3 shows the attractors from asynchronous updating, while this check used synchronous updating. That doesn't affect the comparison: a fixed point is a state that every rule maps to itself, so it's the same whichever order the genes update in. The two methods differ only in cycles. Synchronous updating also produces five two-state cycles that don't appear under asynchronous updating, and these aren't part of Figure 3.

So the discrepancy is in Figure 3, not in the model: two GMP columns appear to have been drawn with Runx1 inactive when the model itself places it active. The number and identity of all 13 stable cell states are otherwise fully reproduced.

## Validation 3: Basin Sizes

The basin of an attractor is the set of all starting states that eventually end up in it. For each of the 2,097,152 possible starting states, I applied the update rules repeatedly until a state repeated. The attractor is the repeated state and everything after its first appearance, which handles both fixed points and cycles. Each attractor is named by its states written as 0/1 strings, sorted and joined, so a cycle gets the same name whichever state a trajectory enters it at. Each starting state then adds 1 to its attractor's count.

**Result: 18 attractors (13 fixed points and 5 two-state cycles), with counts adding up to exactly 2,097,152.** I ran the authors' rules file through BoolNet to get its per-attractor basin counts, and all 18 of my counts match BoolNet's exactly.

Grouping fixed points by marker gene (GATA2 → HSC, GATA1 → MEP, PU.1 → GMP, Gfi1 → CLP) and excluding cycles, as the authors' notebook does:

| Cell type | Starting states | My result | Paper |
|---|---|---|---|
| HSC | 3,496 | 0.18% | 0.18% |
| MEP | 1,617,312 | 84.37% | 84.37% |
| GMP | 81,920 | 4.27% | 4.27% |
| CLP | 214,200 | 11.17% | 11.17% |
| Cycles (excluded) | 180,224 | | |

All four published percentages are reproduced.

Each of the five cycles alternates between two states that differ in mutually repressing genes, such as Cebpa and Mef2c. Under synchronous updating, both genes respond at the same instant to each other's previous value and flip back and forth indefinitely. Under asynchronous updating, these cycles resolve into fixed points, which is why the authors treat them as timing artifacts rather than cell states.

## Mutant Simulations

To simulate a mutation, I force one gene permanently OFF (loss of function, LOF) or ON (gain of function, GOF): the normal rules run, and then the mutated gene is overwritten with its fixed value. Only starting states where the mutated gene already has its fixed value are used, giving 2²⁰ = 1,048,576 starting states per mutant.

Following the authors' notebook, mutant percentages are calculated out of all basins, cycles included.

**Quantitative results reported in the paper's text:**

| Mutation | Attractors | Effect checked | My result | Paper |
|---|---|---|---|---|
| FOXO3 LOF | 18 | HSC basin | 0.067% | 0.067% |
| Antioxidant GOF | 18 | HSC basin | 0.331% | 0.3311% |
| Mef2c LOF | 14 | CLP basin | 6.134% | 6.13% |
| AKT GOF | 16 | HSC attractor | absent | absent |

All four are reproduced. Not every mutant has 18 attractors: removing or forcing a gene can merge or eliminate attractors.

**Qualitative claims from the paper's gain-of-function table (Table 2):**

| Mutation | Attractors | Paper's claim | My result (paper's method) |
|---|---|---|---|
| Meis1 GOF | 18 | MEP basin reduced | MEP 65.12%, reduced ✓ |
| Runx1 GOF | 16 | HSC and MEP basins reduced | HSC 0.049%, MEP 8.01%, both reduced ✓ |

Using the paper's own method, these claims are reproduced too. The next two sections show that how robust they are depends on how the comparison is made.

What the mutants show biologically:

- **FOXO3 loss shrinks the HSC basin** while leaving an HSC state in place, just much harder to reach. FOXO3 drives the antioxidant response that protects stem cells.
- **Antioxidant gain roughly doubles the HSC basin.** Suppressing ROS makes the stem-cell state easier to reach, the paper's central oxygen/ROS argument showing up directly in the numbers.
- **AKT gain removes the HSC state entirely.** AKT pushes cells out of quiescence; with it permanently on, no stable stem state exists.
- **Mef2c loss shifts cells from lymphoid to myeloid:** the CLP basin falls while the GMP basin rises from 3.9% to 6.25%.

## How the Paper Compares Wild Type and Mutants

The paper calculates wild-type percentages with cycles excluded, but mutant percentages with cycles included, and compares the two directly. The paper reports both sets of numbers without noting the difference in denominators.

Because the wild type's cycles-excluded percentages are each divided by 0.914 (the fraction of wild-type starting states that reach a cell fate rather than a cycle: 1,916,928 ÷ 2,097,152), every fold change calculated the paper's way is exactly 0.914 times the fold change calculated with consistent denominators:

```
paper's fold change = consistent fold change × 0.914
```

I confirmed this holds for every mutation and cell type. The effect is systematic: decreases look about 9% stronger than they are, increases look about 9% weaker, and no change looks like a 9% decrease. Any real increase smaller than 9.4% (1 ÷ 0.914) appears as a decrease.

| Mutation | Cell type | Paper's method | Consistent (cycles included) |
|---|---|---|---|
| FOXO3 LOF | HSC | 0.368 | 0.403 |
| Antioxidant GOF | HSC | 1.816 | 1.986 |
| Mef2c LOF | CLP | 0.549 | 0.601 |
| AKT GOF | HSC | eliminated | eliminated |
| FOXO3 LOF | GMP | 0.914 | 1.000 |
| FOXO3 LOF | MEP | 0.928 | 1.015 |
| AKT GOF | MEP | 0.961 | 1.051 |

*(Fold change = mutant percentage ÷ wild-type percentage. Below 1 is a decrease, above 1 an increase.)*

What stands out:

- **The paper's headline mutant results keep their direction.** Every effect it reports in the text points the same way under both methods; only the size changes.
- **An apparent effect that isn't real:** FOXO3, Antioxidant and AKT mutants appear to reduce the GMP basin by about 9% under the paper's method. With consistent denominators, GMP is unchanged at exactly 3.906% in all three.
- **Two small effects reverse direction:** under FOXO3 LOF and AKT GOF, the MEP basin actually increases slightly (by 1.5% and 5.1%), but the paper's method shows both as decreases. The paper's results for these two mutants concern only the HSC basin, so this doesn't overturn a published claim. It does show the mixed comparison can reverse small effects for anyone using it to interpret other cell types.

## Sensitivity to How Cycles Are Treated

A consistent comparison can either include cycles for every condition or exclude them for every condition. The two agree only when a mutation leaves the share of cycles unchanged. When it doesn't, every fold change for that mutant shifts by the same factor:

```
fold change (cycles excluded) = fold change (cycles included) × (0.914 ÷ mutant's fate fraction)
```

where a mutant's fate fraction is the share of its starting states that reach a cell type rather than a cycle.

| Mutant | Cycles | Fate fraction | Multiplier |
|---|---|---|---|
| FOXO3 LOF | 8.6% | 0.914 | 1.000 |
| Antioxidant GOF | 8.6% | 0.914 | 1.000 |
| AKT GOF | 8.6% | 0.914 | 1.000 |
| Mef2c LOF | 11.7% | 0.883 | 1.035 |
| Meis1 GOF | 23.2% | 0.768 | 1.191 |
| Runx1 GOF | 85.7% | 0.143 | 6.40 |

For FOXO3, Antioxidant and AKT the choice makes no difference. For Meis1 and Runx1 GOF it changes the answer:

| Mutation | Cell type | Paper's claim | Paper's method | Consistent, cycles included | Consistent, cycles excluded |
|---|---|---|---|---|---|
| Meis1 GOF | MEP | reduced | 0.772 | 0.844 | **1.005** |
| Runx1 GOF | MEP | reduced | 0.095 | 0.104 | 0.664 |
| Runx1 GOF | HSC | reduced | 0.268 | 0.293 | **1.873** |

- **Meis1 GOF, MEP:** among states that reach a cell fate, MEP makes up 84.83%, essentially the wild type's 84.37%. The reported reduction comes from states being diverted into cycles, not from fewer MEPs among real fates.
- **Runx1 GOF, HSC:** with cycles excluded, the HSC share nearly doubles. The reported reduction reverses.
- **Runx1 GOF, MEP:** reduced under every method. This claim is robust.

The authors exclude cycles from the wild type because they consider them artifacts of synchronous updating. Applying the same reasoning to the mutants removes the reported Meis1 MEP reduction and reverses the Runx1 HSC reduction.

Neither choice is clearly correct. Excluding cycles assumes their states would resolve into cell types in the same proportions as all other states; including them treats "cycle" as a separate outcome. Neither assumption is tested here. For mutants that strongly change the share of cycles, the conclusion depends on that assumption, and synchronous basin sizes alone can't settle it. With 86% of Runx1 GOF's starting states ending in cycles, any basin percentage for that mutant should be read with caution.

## Limitations

- Only the Boolean model was reproduced. The paper's continuous (ODE) and stochastic models were not.
- All basin calculations use synchronous updating, as in the paper's basin analysis. Resolving the cycle question above would need an asynchronous analysis, where the cycles disappear.
- Mutant results were checked against the paper's reported percentages and qualitative claims, not against per-attractor BoolNet counts for each mutant.
- Cell types are assigned by a single marker gene, following the paper.

## Repository Contents

- `model.py` / notebook: the rules, all three validations, the mutant simulations and the comparison tables
- `boolnet_attractors_active.txt`, `boolnet_attractors_full.txt`: BoolNet output from the authors' rules file, used as reference data
- `basin_counts_*.json`: saved basin counts for the wild type and each mutant

AI tools were used for code review and debugging.

## References

Herrera, J., Bensussen, A., García-Gómez, M. L., Garay-Arroyo, A., & Álvarez-Buylla, E. R. (2024). A system-level model reveals that transcriptional stochasticity is required for hematopoietic stem cell differentiation. *npj Systems Biology and Applications*, 10, 145. https://doi.org/10.1038/s41540-024-00469-8

Authors' code: https://github.com/joelhema/Herrera-et-al-2024