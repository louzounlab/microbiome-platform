# LOCATE reference for the product tests

Evidence from the approved V2 validation (validation/locate/PROTOCOL.md, RESULTS-V2.md):
the product (locate_lightning2, Python 3.13) vs the original LOCATE-model 0.0.3 as
published (Python 3.10 reference environment, never shipped), on example E2.

The rule (V2): d = product − original, per seed, in mean test Spearman over metabolites;
pass if the 90% confidence interval of the mean of d lies within ±0.02.

## What the product test checks (since Amendment 3, 2026-11-17)

`test_locate_product_passes_the_approved_rule_on_this_machine` trains the product on seeds
70-109 of E2 on whatever machine runs the tests, and applies the rule against
`original_scores_seeds70-109.csv` (seed, A: the original's score on each seed, recorded on
the laptop). A neural network's training rounds differently on another processor, operating
system or thread count, so the test does not ask for identical numbers, only for results
equivalent to the original's within the approved margin.
(validation/locate/SERVER-peptibase10-PLAN.md, Amendment 3; RESULTS-peptibase10-amend3.md.)

## Recorded evidence

- `v2_scores.csv`: seeds 10-29, the original (A) and the product on the laptop (C).
  Result: PASS, mean −0.0013, interval [−0.0193, +0.0168], no systematic direction.
- `platforms/linux-x86_64/scores.csv`: the lab server peptibase10 under AlmaLinux 9.8,
  Amendment 2, seeds 30-69 (seed, A, S), PASS (RESULTS-peptibase10-amend2.md).
- `platforms/linux-x86_64/scores-amend3.csv`: the same server after its upgrade to
  AlmaLinux 10.3, Amendment 3, seeds 70-109 (RESULTS-peptibase10-amend3.md).

## For reference only (no longer checked)

- `C_pred_val_seed{10,11,12}.csv`, `C_run_seed*.json`: the product's exact test-set
  predictions on the laptop (Windows, torch 2.14.1, Lightning 2.6.6).
- `platforms/linux-x86_64/pred_val_seed{10,11,12}.csv`, `seed*.json`: the same on the lab
  server under AlmaLinux 9.8. After the upgrade to 10.3 they were no longer reproduced
  (seed 10: same stopping epoch, predictions moved by up to 0.014; seeds 11 and 12 stopped at
  other epochs), which is why the test became the tolerance rule.
