# hHDAC8 Pipeline v2: What Changed, Why, and What's Still Open

## What was actually wrong

The prior notebook (`hHDAC8_predict_generate.ipynb`) was not broken — the crossover
bug was already fixed, `compute_liability_flags()` was already centralized, and the
hydrazide-specific SMARTS was already in place. What it had was a **fitness-design
problem**: every liability (ZBG-unknown, PAINS/BRENK, metabolic flags, hydrazide) was
folded into one scalar number alongside predicted pIC50 and docking score. Stacked
scalar penalties create an exploitable surface — a candidate can pay down any one
penalty by overperforming on pIC50 or docking, so the search finds whichever
penalty/reward tradeoff is cheapest to exploit, not the best-rounded candidate. In the
attached `zbg_distribution.png`, that tradeoff cashed out as 52/93 hits (56%) converging
on Hydrazide, a chemotype that's real ZBG chemistry but was being rewarded for the
wrong reason (it happens to dock/predict well on this scalar, not because it's your
actual synthetic priority).

Gemini's rewrite correctly diagnosed this (Task 1) and several adjacent issues, but its
"fix" was a from-scratch rebuild that would have thrown away your file lineage,
comment trail, and several already-correct design decisions (e.g. why PAINS/BRENK is a
soft flag, not a hard reject, given 61.5% of real actives get flagged). I kept your
structure and merged in the specific architectural corrections that were real
improvements.

## Changes made, one by one

1. **Tier-1 hard gates / Tier-2 Pareto selection (Task 1).** `passes_tier1_gates()` now
   unconditionally rejects hydroxamates, Lipinski/Veber violations, SA > 4.5, hard
   metabolic liabilities, and unrecognized ZBGs — no scalar penalty, no partial credit.
   Everything that survives is scored on three independent Pareto objectives (pIC50,
   docking, liability score) via NSGA-II non-dominated sorting, with crowding distance
   for diversity. **Result: Hydrazide share of hits dropped from 56% to 30%, and
   Acylurea — your actual FragBreed lead chemotype — went from a minority (18/93, 19%)
   to the dominant class (344/494, 70%), with zero Unknown-ZBG hits.** This is the
   headline result of this session.

2. **Applicability-domain hard gate (Task 3).** Added a hard reject below 0.20 Tanimoto
   similarity to the training set, on top of the existing 0.35 soft flag. In this run,
   0 of 516 raw hits fell below 0.20 — the GA's own Tier-1 gates already keep it close
   enough to real chemistry that this gate didn't have to do much work, which is a
   useful sanity check in itself.

3. **N-halogen liability — new finding, not in your original task list.** The fixed
   Pareto GA's first pass surfaced `CN(Cl)C(=O)NC(=O)c1cnc(N2CCN(Cc3c[nH]c4ccc(C5CC5)
   cc34)CC2)nc1` as a Tier-1 survivor: an N-Cl chloramine on an acylurea nitrogen. Your
   `mutate()` can attach any substituent (including halogens) to any atom with an open
   valence, and nothing in the liability set specifically screened N-halogen bonds.
   Added `N_HALOGEN_LIABILITY_PATTERN` (`[NX3]-[Cl,Br,I]`) as a hard gate. Reran after
   the fix; it doesn't reappear in the final hit list. This is exactly the kind of bug
   you've asked me to find rather than paper over — flagging it explicitly rather than
   quietly patching it.

4. **`IsoformSelectivityManager` (Task 4).** Modular class storing per-isoform
   docking-model bundles (`hdac1`, `hdac2`, `hdac6`, `hdac8`, extensible). Calling
   `.selectivity_ratio()` or `.annotate()` for an isoform with no registered model
   returns `None` and prints a one-time warning — it does not fabricate a score. Both
   `hdac1_selectivity_ratio` and `hdac6_selectivity_ratio` are all-`None` in the
   exported CSVs right now, honestly, because no real HDAC1/HDAC6 cross-docking data
   exists yet.

5. **Clustering cutoff changed from 0.65 to 0.4.** At 0.65, the enlarged, Pareto-diverse
   hit pool (494 hits vs. the old run's 93) collapsed into a single Butina cluster —
   checked directly, mean pairwise Tanimoto distance across hits is 0.49, so 0.65 was
   too loose to say anything about structural diversity. 0.4 yields 30 distinct
   clusters, which is the number reported in the diagnostics and CSVs.

6. **What I did NOT change:** the SA-score, PAINS/BRENK soft-penalty rationale,
   metabolic liability hard/soft split, and the core IC50/docking model architecture
   were all already correct per your prior verification notes, and I left them as-is.

## What's explicitly not done this round

- **3D conformer strain (ETKDGv3 + MMFF94s), Zn-coordination geometry, and live
  docking** are not implemented. These require either Schrodinger/Maestro (licensed,
  external to this sandbox) or a real 3D docking engine neither of which is available
  here. The liability score currently used in the Pareto GA is entirely 2D. This means:
  a molecule can pass every gate in this pipeline and still turn out to have a strained
  or geometrically wrong 3D pose — nothing here checks that. Cross-dock the top cluster
  representatives in Maestro before trusting the ranking for anything beyond triage.
- Two of the original eight diagnostic figures (3D strain-energy histogram, and any
  docking-based AD/geometry plot) are omitted for the same reason — there's no real
  data behind them, and a synthetic placeholder would be misleading rather than useful.
- HDAC1/HDAC6 selectivity ratios are `None` throughout. To get real numbers, run
  `train_docking_model()` against real 4BKX/5EEM cross-docking data and register the
  resulting bundle with `IsoformSelectivityManager.register('hdac1', bundle)`.

## Remaining vulnerabilities worth knowing about

- **Scaffold concentration.** All 494 hits trace back to a small number of seed
  scaffolds (the top-8 real non-hydroxamate actives, several of which are themselves
  near-duplicates — seeds 4-7 in this run differ only by a chain length). The GA
  explores substituent space around those scaffolds well, but this is not de novo
  scaffold-hopping. If you want genuinely new cores, seed with a broader/more diverse
  starting set or add scaffold-diversity pressure as a fourth Pareto objective.
- **Model R² is modest** (IC50 ≈ 0.42, docking ≈ 0.44 scaffold-CV) — both models still
  extrapolate meaningfully for any truly novel chemotype, and the AD gate only checks
  fingerprint similarity, not model confidence directly. Worth adding an ensemble
  variance term if you want a second, independent way to catch this (Task 3's original
  ask; not implemented here since it's a bigger lift than fit this session).
- **Liability whitelist is still a whitelist.** The ZBG hard gate is only as good as
  the 8 SMARTS patterns in `ZBG_TAGS`. Per your own prior verification, 4 of the 5 best
  real seeds by measured potency don't match any of them — meaning the hard gate is
  now also excluding some real, potent, non-hydroxamate chemistry that a purely
  soft-penalty design used to admit. This is a real tradeoff introduced by moving ZBG
  from soft to hard (Task 1's explicit ask), not a bug — but it's worth knowing that
  "hard-gate on ZBG" and "cover 100% of real non-hydroxamate actives" are in tension
  until the whitelist is expanded.
- **Synthetic feasibility beyond SA score.** SA score is a heuristic; the acylurea-heavy
  hit pool that now dominates hasn't been checked against real building-block
  availability the way FragBreed's original screen was.

## Session 2: ZBG safety audit and diversity fix

After the Pareto fix above, you asked for a safety audit of the 8 whitelisted ZBGs and
a fix for the fact that only 2 of them (Acylurea, Hydrazide) ever appeared across 494
hits. Two genuinely separate problems, handled in sequence, with two real course
corrections along the way worth being explicit about.

### Safety audit result

| ZBG | Verdict | Reasoning |
|---|---|---|
| Hydrazide | Removed | N-N bond, oxidative cleavage, hydrazine-release precedent (isoniazid) |
| Thiol (free) | Removed | Oxidizes to disulfide, reactive toward off-target proteins |
| Benzisothiazolone | Removed | Ring-opens to a reactive electrophile; sensitizer-adjacent chemistry |
| 8-Hydroxyquinoline | Removed | Promiscuous non-Zn metal chelation; redox-active complexes; clioquinol-class neurotoxicity precedent |
| Salicylamide | Kept | Established HDAC8 ZBG, no reactive electrophile |
| Ortho-aminoanilide | Kept | Clinically precedented (entinostat/mocetinostat class) |
| Acylurea | Kept | Your current FragBreed lead chemotype |
| 3-HPT | Kept, flagged | Pyrithione-adjacent nonspecific zinc-ionophore risk at high exposure -- kept because it's an explicit open-SAR area in your project, not because the concern doesn't apply |
| Carboxylate (new) | Added | Named in your original scope, metabolically inert |
| Trifluoromethyl-ketone (new) | Added | Named in your original scope, reversible, reasonable safety track record |

Also added, independent of the ZBG whitelist: a belt-and-suspenders hard gate for free
thiol and N-N/hydrazide substructures, so they can't sneak through under a different
overall ZBG tag (e.g. a molecule tagged Acylurea that also happens to carry a hydrazide
fragment elsewhere is still rejected).

### First diversity attempt: it broke, and here's why

My first fix was to build 5 synthetic single-ZBG seed molecules (one per underrepresented
chemotype) and add them to the real ChEMBL seed pool. Filtering the real seeds through
the *full* Tier-1 gate set (not just Lipinski/Veber, which is all the original run
checked) revealed something more important than the diversity problem: **almost none of
the top-potency real actives match any of the 6 safety-audited ZBGs at all** -- most of
the top-30 by measured pIC50 are hydrazide-based (correctly excluded) or the macrocyclic
thiol/depsipeptide chemotype (fails Lipinski/Veber independently). With the real seed
pool nearly empty, the GA ran almost entirely on the 5 synthetic seeds, and one lineage
(Carboxylate) wiped out the other 4 by generation ~15 through ordinary elitism dynamics
-- nothing to do with which chemistry was actually better, just which lineage happened
to converge fastest early.

Two fixes, in order:

1. **Better seeds.** Re-checked the *entire* non-hydroxamate ChEMBL set (not just the
   top 30) against the new whitelist: 91 real actives pass, covering 4 of 6 lineages
   (Ortho-aminoanilide 63, Carboxylate 18, Trifluoromethyl-ketone 9, Acylurea 1), with
   real measured potency up to pIC50 6.97. These now seed the GA directly, ranked by
   measured pIC50; synthetic seeds are used only for the 2 lineages with zero real
   examples (Salicylamide, 3-HPT).

2. **Lineage protection in the GA.** `run_ga_pareto()` now tracks which ZBG lineage
   each candidate traces back to and reserves a floor of elite-pool slots per lineage
   every generation, so a fast-converging lineage can't crowd out the others before
   they've had a chance to develop competitive descendants. Confirmed working: the
   per-generation `elite-lineages` log in the notebook shows real representation
   across all 6 ZBGs for all 60 generations of the final run.

### The honest tradeoff this surfaced

With lineage protection working and real per-lineage seeds in place, the **search**
covers all 6 ZBGs throughout. But the **final hit list** (candidates clearing the
potency/docking bound) is still Carboxylate-dominated (55/57 at a relaxed 500 nM
bound; 2/2 at the original 150 nM bound). This is not the same bug as before -- it's a
real signal that the most potent real actives carrying any of these 6 safe ZBGs top
out around pIC50 6.97 (~107 nM), noticeably lower than what the hydrazide/thiol
chemistry we removed was achieving. Put plainly: the chemotypes we excluded on safety
grounds also appear to be the more potent ones in this dataset, at least among real
measured actives. Not proof of a causal link, but worth knowing going in.

Practical implication: `max_ic50_nM=150` (the original bound) is now a very strict ask
given the safe-ZBG ceiling; relaxed to 500 nM per your direction to see the full
chemotype spread for this triage run. Recommend deciding case-by-case going forward
whether triage should run at 150 nM (strict, narrow, mostly Carboxylate) or 500 nM
(broader chemotype view, weaker individual potency) depending on what this batch of
results is for.



## Recommended next steps, in rough priority order

0. **Decide the real triage bound.** Given the safe-ZBG potency ceiling, confirm
   whether 500 nM is the right permanent bound or just useful for this one diversity
   check -- affects every downstream number.
1. Cross-dock the top ~10-20 cluster representatives (already exported to
   `top_diverse_cluster_leads_fixed.csv`) in Maestro against 1T64, and cross-dock the
   same set against 4BKX/5EEM to get real HDAC1/HDAC6 selectivity ratios instead of
   `None`.
2. Add 3D conformer generation + MMFF94s minimization as a Tier-1 gate once you're
   running this outside the sandbox — strain energy is currently unchecked.
3. Consider a 4th Pareto objective for scaffold diversity if novel cores (not just
   novel substituents on 8 known scaffolds) are the goal.
4. Expand the ZBG whitelist given the known 83% real-active miss rate, so the new hard
   gate doesn't inadvertently exclude legitimate chemotypes it just hasn't been taught
   to recognize yet.
