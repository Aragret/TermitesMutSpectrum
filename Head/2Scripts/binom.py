import argparse
import sys
from typing import Dict, Iterable, List, Tuple, Union

import pandas as pd
import scipy.stats
import statsmodels.api as sm

PATH_TO_MUTSPEC = "/home/glebo/Documents/practice/python/etc/mutspec.csv"
COLS = ["NucSubst", "ObsToExp"]
sub_type = [
    "A>C",
    "A>G",
    "A>T",
    "C>A",
    "C>G",
    "C>T",
    "G>A",
    "G>C",
    "G>T",
    "T>A",
    "T>C",
    "T>G"
]
directional_pairs = [
    ("A>C", "C>A"),
    ("A>G", "G>A"),
    ("A>T", "T>A"),
    ("C>G", "G>C"),
    ("C>T", "T>C"),
    ("G>T", "T>G"),
]
reciprocal_pairs = [
    ("A>C", "T>G"),
    ("A>G", "T>C"),
    ("A>T", "T>A"),
    ("C>G", "G>C"),
    ("C>T", "G>A"),
    ("G>T", "C>A"),
]


def asterics_for_vector(pvals: Iterable) -> List[str]:
    asterics = []
    for val in pvals:
        if val == 0.0:
            asterics.append("")
        elif val < 0.001:
            asterics.append("***")
        elif val < 0.01:
            asterics.append("**")
        elif val < 0.05:
            asterics.append("*")
        else:
            asterics.append("")
    return asterics


def binom_testing(
        subs: List[Tuple[str]],
        pairs: List[Tuple[str]],
        mut_num: Dict[str, Union[int, float]],
        label: str = None):
    data = []
    min_mut_num = min([num for num in mut_num.values() if num != 0])
    for sub in subs:
        mut_num[sub] = mut_num[sub] / min_mut_num
    for mut1, mut2 in pairs:
        n1, n2 = round(mut_num[mut1]), round(mut_num[mut2])
        if n2 == 0 or n1 == 0:
            row = (mut1, mut2, 0, 0) if label is None else (
                label, mut1, mut2, 0, 0)
            data.append(row)
        else:
            ratio = n1 / n2
            res = scipy.stats.binom_test(n1, n1 + n2, p=0.5)
            pval = res
            row = (mut1, mut2, ratio, pval) if label is None else (
                label, mut1, mut2, ratio, pval)
            data.append(row)
    base_cols = ["mut1", "mut2", "ratio", "pval"]
    cols = base_cols if label is None else ["label"] + base_cols
    data = pd.DataFrame(data, columns=cols)

    _, qval, _, _ = sm.stats.multipletests(
        data["pval"].values, method="fdr_bh")  # adjust pval
    data["pval_adj"] = qval
    data["asterics"] = asterics_for_vector(qval)
    return data


def main():
    parser = argparse.ArgumentParser(
        description="Do statistics (binomial test) for 12-comp. mutspec distinct values.\n"
        "For each pair of reciprocal and directional substitutions count significance of its difference."
    )
    parser.add_argument("inp", type=argparse.FileType("r"),
                        help="Path to input 16-comp MutSpec table. Required columns [{}, {}]"
                        .format(*COLS))
    parser.add_argument(
        "out", default=None, nargs="?", type=argparse.FileType("w"),
        help="Path to output csv with statistics, (by default write file in the same dir: `inp`.compared.csv)"
    )
    args = parser.parse_args()

    inp = args.inp
    out = args.out or inp.name.replace(".csv", ".compared.csv")
    print(f"input file: {inp.name}", file=sys.stderr)
    print(
        f"output file: {out if isinstance(out, str) else out.name}", file=sys.stderr)

    mutspec = pd.read_csv(inp, usecols=COLS)
    for c in COLS:
        assert c in mutspec.columns, f"Column {c} required in MutSpec table"

    mut_num = dict(zip(mutspec.NucSubst, mutspec.ObsToExp))

    dirp = binom_testing(sub_type, directional_pairs, mut_num, "directional")
    recp = binom_testing(sub_type, reciprocal_pairs, mut_num, "reciprocal")

    compared = pd.concat([dirp, recp], axis=0).reset_index(drop=True)
    compared.to_csv(out, index=None)
    print("\nDone", file=sys.stderr)


if __name__ == "__main__":
    main()

