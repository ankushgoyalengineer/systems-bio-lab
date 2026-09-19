# Gene Expression with Negative Feedback: When Does Repression Cause Overshoot?

## Why I Built This

After working through enzyme kinetics, I wanted to move from a single reaction into an actual gene regulatory circuit — the mRNA → protein cascade, with a protein that represses its own transcription. This motif shows up everywhere in biology (circadian clocks, stress responses, developmental switches), and "negative feedback stabilizes the system" is one of those phrases that gets repeated a lot without anyone showing *why*, or what "stabilizes" actually looks like dynamically. I wanted to build the system myself and find out.

The core question I set out to answer: "Does adding self-repression to a simple gene expression cascade just shift the system to a new steady state, or does it change the *behavior* of how the system gets there — and if so, under what conditions?"

## Grounding the Parameters

Rather than picking round numbers, I parameterized the model using values reported for *E. coli* gene expression:

- **mRNA and protein decay rates** (`ym = 0.14 min⁻¹`, `yp = 0.02 min⁻¹`) come from standard *E. coli* half-lives: bacterial mRNA typically degrades with a half-life around 5 minutes (`ln(2)/5 ≈ 0.14`), while stable proteins are diluted mainly by cell division, roughly a 35-minute cycle (`ln(2)/35 ≈ 0.02`) [Alon, 2019].
- **Transcription and translation rates** (`alpha = 10.0`, `beta = 5.0`) represent canonical strong-promoter and ribosome-binding-site efficiencies used as baseline modeling values in the same text [Alon, 2019].
- **Repressor mechanics** (`Kd = 100.0` molecules/cell, `n = 2–4`) are grounded in the theoretical framework Elowitz and Leibler used to design the repressilator circuit computationally before building it *in vivo* [Elowitz & Leibler, 2000].

Units throughout: time in minutes, mRNA and protein in molecules per cell.

## Stage 1: Baseline — No Feedback

Before touching any repression, I built the plain two-variable cascade: constant transcription produces mRNA, mRNA is translated into protein, and both decay linearly.

Before trusting the simulation, I solved for the steady state by hand:

dM/dt = 0  →  M = alpha / ym
dP/dt = 0  →  P = beta * M / yp

With the grounded parameters above (alpha=10.0, ym=0.14, beta=5.0, yp=0.02):

M = 10.0 / 0.14        = 71.43 molecules/cell
P = 5.0 * 71.43 / 0.02 = 17,857.14 molecules/cell

The simulation converged to M ≈ 71.43 and P ≈ 17,857.14 (printed final value: 17857.142857132494) — agreement to roughly 8 significant figures, with the tiny residual explained by Euler discretization over many steps. This checkpoint mattered before adding any complexity: if the baseline hadn't matched the hand calculation, nothing built on top of it would have been trustworthy.

The unregulated system rises smoothly to steady state with no overshoot — there's nothing in the equations that would cause it to overcorrect, since both species just decay toward equilibrium. That's the baseline everything else is compared against.

NOTE: mRNA plateaus at ~71 molecules/cell, too small to see on this scale next to protein's ~17,857 — see printed values above for the exact number

![Baseline: no feedback](03-gene-expression/no_feedback_baseline.png)

## Stage 2: Adding Repression

I added a Hill function so transcription decreases as protein accumulates:

transcription_rate = alpha / (1 + (P / Kd)^n)

`n` controls how sharply repression switches on as protein crosses the threshold `Kd`. With biologically realistic repressor cooperativity (n=4, consistent with the repressilator's design), the system no longer approaches steady state smoothly. Protein peaks at 588.31 molecules/cell before settling to a final steady state of 281.18 — a peak-to-final ratio of 2.09, meaning the transient overshoot is roughly double the eventual steady-state level. That overshoot is the first sign that feedback doesn't just relocate the equilibrium, it changes the *path* the system takes to get there.

![Feedback with n=4: overshoot and settle](03-gene-expression/feedback_n4.png)

## Stage 3: A Numerical Artifact I Almost Missed

Curious how far this went, I initially pushed the Hill coefficient to a mathematical extreme (n=1000) to see whether the system could be made to oscillate indefinitely rather than settle. I want to be explicit that **n=1000 has no biological meaning** — real transcription factor cooperativity, including the repressilator's, sits in the range of roughly 1–4; there's no physical repressor complex of a thousand subunits. This was a stress test of the equation, not a biological scenario.

At my original step size (dt=0.01), that stress test produced what looked like small, persistent wiggling that never fully died out. I didn't trust it outright: a two-variable feedback loop like this one, with no third species and no explicit delay, is mathematically incapable of sustaining a genuine oscillation, since both variables only ever decay — this keeps the system net-damping regardless of how steep the repression is. Real sustained genetic oscillators (like the repressilator itself, which needs three interlocking repressors) require at least three species or a built-in delay. A "sustained oscillation" showing up here was a red flag rather than a discovery.

I checked it against a much finer step size (dt=0.001), extending the comparison to the full 1440-minute (24-hour) window used consistently elsewhere in this analysis, rather than the shorter window I'd initially tested. That correction mattered: at a shorter window, the system hadn't yet reached true steady state, so an earlier version of this check was inadvertently comparing two mid-transient values rather than two converged endpoints. Over the full window, both step sizes converge to the same steady state — 160.0864741377557 (dt=0.01) versus 160.0864741377811 (dt=0.001) — agreeing to roughly 10 significant figures.

![n=10: dt sensitivity check](03-gene-expression/dt_sensitivity_n10.png)

That confirmed dt=0.01 is trustworthy in the range I actually use for the rest of the analysis (n up to 10), once compared on equal, fully-converged footing. The earlier appearance of sustained wiggling at n=1000 reflects the integrator struggling to resolve an artificially sharp, non-physiological switch — not a property of any biologically realistic version of this circuit.

This ended up being one of the more useful things to come out of the project — not because I found a bug in the biologically relevant range, but because a surprising result is worth checking against known theory and hard numbers, on a properly matched comparison, before it's worth reporting. Catching that my first version of this check hadn't actually reached steady state was as important as the original artifact investigation itself.

## Stage 4: How Overshoot Scales with Feedback Steepness

With n=1000 set aside as a non-biological stress test, I swept the Hill coefficient across the biologically realistic range — n = 1, 2, 4, 8, 10, consistent with real repressor cooperativity — and measured max(protein) / final(protein) for each, using the grounded parameters throughout.

My prediction going in: steeper repression should produce a sharper, later-triggering brake on transcription, so overshoot should increase with n.

![Hill coefficient vs overshoot ratio](03-gene-expression/hill_sweep_overshoot.png)

n	Overshoot ratio
1	1.0000
2	1.3163
4	2.0923
8	3.1336
10	3.4580

(Exact values from overshoot_pairs. The n=4 ratio, 2.0923, matches the standalone n=4 run in Stage 2 to every printed digit — an independent cross-check that the sweep and the illustrative single run agree.)

The result confirmed the prediction's direction: overshoot increases with n across the entire range I tested, with no sign of leveling off even at n=10. This differs from what I initially assumed while exploring further into non-physiological territory (n well past 100), where the curve does eventually bend and flatten mathematically — but that flattening happens far outside any biologically plausible repressor cooperativity, so it isn't the headline result here. Within the range that actually corresponds to real gene circuits, the honest finding is simpler: steeper, more switch-like repression produces progressively larger overshoot, monotonically, across the full biologically realistic range.

What This Project Actually Showed

Negative feedback in this system doesn't just relocate the steady state — it introduces a transient overshoot-and-settle behavior that gets more pronounced as repression sharpens, and within realistic bounds on cooperativity, that relationship keeps climbing rather than saturating. Separately: a surprising simulation result is worth doubting and checking against known theory and hard numbers before it's worth reporting. That check itself needs the same scrutiny — my first attempt to validate the n=1000 finding compared two values that turned out not to be true steady states at all, and catching that was as important as the original question. The corrected check confirmed the apparent instability was a resolution artifact tied to a non-physiological parameter value, not a property of the biology — but only once the verification itself was verified.

References

Alon, U. (2019). An Introduction to Systems Biology: Design Principles of Biological Circuits (2nd ed.). CRC Press.

Elowitz, M. B., & Leibler, S. (2000). A synthetic oscillatory network of transcriptional regulators. Nature, 403(6767), 335–338.

(Full implementation, including the GeneNetwork class and the sweep script, is available in model.py and the accompanying notebook.)