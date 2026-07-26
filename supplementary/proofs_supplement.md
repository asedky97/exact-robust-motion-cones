# Supplementary Proofs

Companion to *Exact Robust Motion Cones: Guaranteed Sticking Manipulation under
Parametric and Geometric Contact Uncertainty* (IEEE Access submission). This
document holds the proof of Lemma 1, deferred from the main text for space.

## Lemma 1 (Branch alignment)

**Statement.** Let {C_i}_{i=1}^n be intervals on the unit circle, each of angular
width < π. If ⋂_i C_i ≠ ∅, then there exists a reference direction r such that
(i) r ∈ ⋂_i C_i, and (ii) for every i, the angular representation of C_i on the
branch (∠r − π, ∠r + π] is a single contiguous interval.

**Proof.** Since the intersection is nonempty, it is itself an interval (possibly
degenerate). Pick any r in the intersection. For each i, r ∈ C_i and
width(C_i) < π, so the angular distance from ∠r to any point in C_i is strictly
less than π in either direction. Hence C_i is contained within the 2π-wide branch
centred at ∠r, and its representation on that branch is a single interval
[a_i, b_i] with a_i < ∠r < b_i. ∎

## Corollary 1 and Remark 1

Corollary 1 (one-recenter fallback) and Remark 1 (implementation) appear in
the main text immediately after the lemma statement. Corollary 1 is not
proved here or anywhere in the package: Lemma 1 only guarantees that *some*
common reference exists when the true intersection is nonempty; it says
nothing about whether the paper's specific mechanical procedure (try
r0 = ∠(M0n0), and on failure retry once from r1 = (Φ⁻+Φ⁺)/2) actually finds
that reference in at most one retry. The claim rests on a continuity
argument over the connected uncertainty box that makes it plausible for the
full model class, not a proof of it.

It is also empirically untested: `branch_alignment_check.py` samples 10,000
realistic geometries and 10,000 adversarial ones (full-circle normal angle,
lever arms down to 3 cm) against the paper's uncertainty box, and the
recenter step is triggered zero times in both sweeps: r0 alone is always
sufficient. That means the experiments confirm r0 is sufficient in the
regimes tested; they say nothing about whether the one-recenter fallback
behaves correctly, since it has never actually fired.
