"""Amendment 3 (SERVER-peptibase10-PLAN.md, recorded before running): the server after its OS upgrade
(AlmaLinux 10.3) (S) vs the original (A) on 40 fresh seeds (70-109), V2's rule with t(0.95, 39).
Writes RESULTS-peptibase10-amend3.md and the scores the tests use:
tests/data/locate_reference/original_scores_seeds70-109.csv (seed, A: the product test's reference) and
tests/data/locate_reference/platforms/linux-x86_64/scores-amend3.csv (seed, A, S).

    python validation/locate/compare_amend3.py [<locate-reference folder>]

A: <reference>/runs-amend3/A/E2/seed<n>/ (run_all.ps1, the Python 3.10 reference setup).
S: <reference>/server-amend3/seed<n>.json (server_check.py run on peptibase10; copied from there).
Nothing here is a parameter.
"""

import json
import math
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import t as student_t

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from compare_v2 import mean_test_spearman  # noqa: E402  (V2's score, unchanged)

REFERENCE = Path(sys.argv[1] if len(sys.argv) > 1 else r"C:\Users\noapa\code\locate-reference")
SEEDS = range(70, 110)
TESTS = HERE.parents[1] / "tests" / "data" / "locate_reference"
T_95_39 = 1.685  # t(0.95, df = 39), as recorded in the amendment
MARGIN = 0.02


def main() -> None:
    assert len(SEEDS) == 40 and abs(student_t.ppf(0.95, 39) - T_95_39) < 5e-4  # the recorded constant is right
    measured = pd.read_csv(HERE / "data" / "E2" / "Y_val.csv", index_col=0, float_precision="round_trip")
    rows, machines = [], set()
    for seed in SEEDS:
        folder = REFERENCE / "runs-amend3" / "A" / "E2" / f"seed{seed}"
        info = json.loads((folder / "run.json").read_text())
        if (info["impl"], info["torch"], info.get("lstsq")) != ("original", "1.13.0+cpu", "default"):
            raise SystemExit(f"A seed {seed} is not the setup the amendment requires: {info}")
        a = mean_test_spearman(pd.read_csv(folder / "pred_val.csv", index_col=0, float_precision="round_trip"),
                               measured)
        server = json.loads((REFERENCE / "server-amend3" / f"seed{seed}.json").read_text())
        if server["torch"] != "2.14.1+cpu" or server["seed"] != seed:
            raise SystemExit(f"S seed {seed} is not the deployed product: {server}")
        machines.add(server["machine"])
        rows.append({"seed": seed, "epochs_A": info["epochs"], "epochs_S": server["epochs"], "A": a,
                     "S": server["score"], "d": server["score"] - a})
    table = pd.DataFrame(rows)
    mean, sd = float(table["d"].mean()), float(table["d"].std(ddof=1))
    half = T_95_39 * sd / math.sqrt(len(table))
    low, high = mean - half, mean + half
    passed = -MARGIN <= low and high <= MARGIN
    direction = ("systematically lower" if high < 0 else "systematically higher" if low > 0
                 else "no systematic direction")
    lines = [
        "# LOCATE on peptibase10 after the OS upgrade: Amendment 3 results (40 fresh seeds, V2's rule)",
        "",
        "Rule (SERVER-peptibase10-PLAN.md, Amendment 3, recorded before running): E2, seeds 70-109; "
        "d = mean test Spearman(S) − (A); pass if the 90% confidence interval of the mean of d, "
        f"mean ± {T_95_39} × SD / √40, lies within ±{MARGIN}.",
        "",
        f"A: the original LOCATE-model 0.0.3 (Python 3.10, torch 1.13.0+cpu) on the laptop. "
        f"S: the product on {', '.join(sorted(machines))} (AlmaLinux 10.3, torch 2.14.1+cpu).",
        "",
        f"## Verdict: **{'PASS' if passed else 'FAIL'}**",
        "",
        f"- Mean difference d (S − A): **{mean:+.4f}**; 90% confidence interval **[{low:+.4f}, {high:+.4f}]** "
        f"(margin ±{MARGIN}).",
        f"- Direction: **{direction}**. S was lower in {(table['d'] < 0).sum()} of {len(table)} seeds.",
        f"- Mean test Spearman: A {table['A'].mean():.4f}, S {table['S'].mean():.4f}. SD of d: {sd:.4f}.",
        "",
        "| Seed | Epochs A | Epochs S | Spearman A | Spearman S | d = S − A |",
        "|---|---|---|---|---|---|",
    ]
    lines += [f"| {r.seed} | {r.epochs_A} | {r.epochs_S} | {r.A:.4f} | {r.S:.4f} | {r.d:+.4f} |"
              for r in table.itertuples()]
    (HERE / "RESULTS-peptibase10-amend3.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    table[["seed", "A"]].to_csv(TESTS / "original_scores_seeds70-109.csv", index=False, float_format="%.17g")
    table[["seed", "A", "S"]].to_csv(TESTS / "platforms" / "linux-x86_64" / "scores-amend3.csv", index=False,
                                     float_format="%.17g")
    sys.stdout.reconfigure(encoding="utf-8")  # "−" and "±" on a Windows console
    print("\n".join(lines[6:11]))


if __name__ == "__main__":
    main()
