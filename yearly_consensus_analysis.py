#!/usr/bin/env python3

import argparse
import os
import re
from collections import Counter

import numpy as np
import pandas as pd
from Bio import SeqIO


UNKNOWN_AA = set(["X", "B", "Z", "J", "U", "O", "*", "?"])


def parse_year_from_header(header):
    """
    从 fasta header 中解析年份。

    适用于:
        >ADV76673|1968
        >APO19006|2016

    年份只取 | 后面的最后一个字段，避免 APO19006 被误识别成 1900。
    """

    main_id = header.split()[0]

    if "|" not in main_id:
        return None

    year_str = main_id.split("|")[-1].strip()

    if re.fullmatch(r"\d{4}", year_str):
        year = int(year_str)
        if 1800 <= year <= 2200:
            return year

    return None


def read_aligned_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))

    if len(records) == 0:
        raise ValueError("No sequences found in fasta file.")

    lengths = [len(r.seq) for r in records]

    if len(set(lengths)) != 1:
        raise ValueError(
            "Input fasta is not aligned: sequence lengths are not identical."
        )

    data = []

    for r in records:
        seq_id = r.id
        desc = r.description
        year = parse_year_from_header(desc)

        data.append(
            {
                "seq_id": seq_id,
                "description": desc,
                "year": year,
                "seq": str(r.seq).upper(),
            }
        )

    df = pd.DataFrame(data)

    n_missing_year = df["year"].isna().sum()

    if n_missing_year > 0:
        print(f"[Warning] {n_missing_year} sequences have no valid year and will be removed.")

    df = df.dropna(subset=["year"]).copy()
    df["year"] = df["year"].astype(int)

    return df, lengths[0]


def choose_reference_sequence(df, ref_id=None):
    """
    用于建立 position label。
    默认选择最早年份的第一条序列。
    """

    if ref_id is not None:
        hit = df[df["seq_id"] == ref_id]

        if len(hit) == 0:
            hit = df[df["description"].str.contains(re.escape(ref_id), regex=True)]

        if len(hit) == 0:
            raise ValueError(f"Reference sequence not found: {ref_id}")

        ref_row = hit.iloc[0]

    else:
        earliest_year = df["year"].min()
        ref_row = df[df["year"] == earliest_year].iloc[0]

    return ref_row["seq_id"], ref_row["year"], ref_row["seq"], ref_row["description"]


def build_position_map(ref_seq):
    """
    建立 alignment position 到 reference ungapped position 的映射。

    aln_pos: 比对坐标，包含 gap
    ref_pos: 参考序列去 gap 后的位置
    """

    mapping = []

    ref_pos = 0

    for i, aa in enumerate(ref_seq):
        aln_pos = i + 1

        if aa != "-":
            ref_pos += 1
            ref_pos_value = ref_pos
        else:
            ref_pos_value = np.nan

        mapping.append(
            {
                "aln_pos": aln_pos,
                "ref_pos": ref_pos_value,
                "ref_aa": aa,
            }
        )

    return pd.DataFrame(mapping)


def consensus_at_position(aas, ignore_gap=True, ignore_unknown=True):
    """
    对某一年某个位点的 amino acids 计算 consensus AA 和频率。
    """

    n_total = len(aas)

    valid = []

    for aa in aas:
        if ignore_gap and aa == "-":
            continue

        if ignore_unknown and aa in UNKNOWN_AA:
            continue

        valid.append(aa)

    n_valid = len(valid)

    if n_valid == 0:
        return {
            "consensus_aa": "X",
            "consensus_count": 0,
            "consensus_fraction": np.nan,
            "consensus_percent": np.nan,
            "n_valid": 0,
            "n_total": n_total,
            "aa_counts": "",
        }

    c = Counter(valid)
    consensus_aa, consensus_count = c.most_common(1)[0]

    consensus_fraction = consensus_count / n_valid
    consensus_percent = consensus_fraction * 100

    aa_counts = ";".join([f"{aa}:{count}" for aa, count in c.most_common()])

    return {
        "consensus_aa": consensus_aa,
        "consensus_count": consensus_count,
        "consensus_fraction": consensus_fraction,
        "consensus_percent": consensus_percent,
        "n_valid": n_valid,
        "n_total": n_total,
        "aa_counts": aa_counts,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Calculate yearly consensus amino acids and frequencies from aligned HA fasta."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input aligned fasta file, e.g. HA_1968_2026.aln.fasta",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        default="yearly_consensus_output",
        help="Output directory",
    )

    parser.add_argument(
        "--ref-id",
        default=None,
        help="Reference sequence ID. If not provided, use first sequence from earliest year.",
    )

    parser.add_argument(
        "--include-ref-gaps",
        action="store_true",
        help="Include alignment columns where reference sequence has gap.",
    )

    parser.add_argument(
        "--keep-gaps-in-consensus",
        action="store_true",
        help="Allow '-' to be used in consensus calculation.",
    )

    parser.add_argument(
        "--keep-unknown",
        action="store_true",
        help="Allow X/B/Z/* etc. to be used in consensus calculation.",
    )

    parser.add_argument(
        "--percent-decimals",
        type=int,
        default=1,
        help="Number of decimals for percentage output.",
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("[1] Reading aligned fasta...")
    df, aln_len = read_aligned_fasta(args.input)

    print(f"Number of sequences with valid year: {len(df)}")
    print(f"Year range: {df['year'].min()} - {df['year'].max()}")
    print(f"Alignment length: {aln_len}")

    print("[2] Choosing reference sequence for position mapping...")
    ref_id, ref_year, ref_seq, ref_desc = choose_reference_sequence(df, args.ref_id)

    print(f"Reference ID: {ref_id}")
    print(f"Reference year: {ref_year}")
    print(f"Reference description: {ref_desc}")

    pos_map = build_position_map(ref_seq)

    pos_map.to_csv(
        os.path.join(args.outdir, "alignment_to_reference_position_map.tsv"),
        sep="\t",
        index=False,
    )

    if args.include_ref_gaps:
        used_pos_map = pos_map.copy()
    else:
        used_pos_map = pos_map[pos_map["ref_pos"].notna()].copy()

    used_pos_map["position_label"] = used_pos_map.apply(
        lambda row: f"pos_{int(row['ref_pos'])}"
        if pd.notna(row["ref_pos"])
        else f"aln_{int(row['aln_pos'])}",
        axis=1,
    )

    print(f"Number of analyzed positions: {len(used_pos_map)}")

    ignore_gap = not args.keep_gaps_in_consensus
    ignore_unknown = not args.keep_unknown

    print(f"Ignore gaps in consensus: {ignore_gap}")
    print(f"Ignore unknown amino acids in consensus: {ignore_unknown}")

    print("[3] Calculating yearly consensus...")

    long_records = []

    years = sorted(df["year"].unique())

    for year in years:
        sub = df[df["year"] == year]
        seqs = sub["seq"].tolist()
        n_seq = len(seqs)

        for _, prow in used_pos_map.iterrows():
            aln_pos = int(prow["aln_pos"])
            idx = aln_pos - 1

            aas = [s[idx] for s in seqs]

            result = consensus_at_position(
                aas,
                ignore_gap=ignore_gap,
                ignore_unknown=ignore_unknown,
            )

            long_records.append(
                {
                    "year": year,
                    "aln_pos": aln_pos,
                    "ref_pos": prow["ref_pos"],
                    "position_label": prow["position_label"],
                    "ref_aa": prow["ref_aa"],
                    "n_sequences_this_year": n_seq,
                    **result,
                }
            )

    long_df = pd.DataFrame(long_records)

    # ref_pos 整理
    long_df["ref_pos"] = long_df["ref_pos"].apply(
        lambda x: int(x) if pd.notna(x) else np.nan
    )

    long_out = os.path.join(args.outdir, "yearly_consensus_long.tsv")
    long_df.to_csv(long_out, sep="\t", index=False)

    print(f"Long-format table saved: {long_out}")

    print("[4] Creating matrix tables...")

    aa_mat = long_df.pivot(
        index="year",
        columns="position_label",
        values="consensus_aa",
    )

    percent_mat = long_df.pivot(
        index="year",
        columns="position_label",
        values="consensus_percent",
    )

    fraction_mat = long_df.pivot(
        index="year",
        columns="position_label",
        values="consensus_fraction",
    )

    # 保持 position 顺序
    position_order = used_pos_map["position_label"].tolist()

    aa_mat = aa_mat[position_order]
    percent_mat = percent_mat[position_order]
    fraction_mat = fraction_mat[position_order]

    aa_out = os.path.join(args.outdir, "yearly_consensus_AA_matrix.tsv")
    percent_out = os.path.join(args.outdir, "yearly_consensus_percent_matrix.tsv")
    fraction_out = os.path.join(args.outdir, "yearly_consensus_fraction_matrix.tsv")

    aa_mat.to_csv(aa_out, sep="\t")
    percent_mat.to_csv(percent_out, sep="\t", float_format=f"%.{args.percent_decimals}f")
    fraction_mat.to_csv(fraction_out, sep="\t", float_format="%.6f")

    print(f"Consensus AA matrix saved: {aa_out}")
    print(f"Consensus percent matrix saved: {percent_out}")
    print(f"Consensus fraction matrix saved: {fraction_out}")

    print("[5] Creating combined AA(percent) matrix...")

    combined_mat = aa_mat.copy()

    for year in aa_mat.index:
        for pos in aa_mat.columns:
            aa = aa_mat.loc[year, pos]
            pct = percent_mat.loc[year, pos]

            if pd.isna(pct):
                combined_mat.loc[year, pos] = f"{aa}(NA)"
            else:
                combined_mat.loc[year, pos] = f"{aa}({pct:.{args.percent_decimals}f}%)"

    combined_out = os.path.join(args.outdir, "yearly_consensus_AA_percent_combined_matrix.tsv")
    combined_mat.to_csv(combined_out, sep="\t")

    print(f"Combined AA(percent) matrix saved: {combined_out}")

    print("[6] Exporting yearly consensus fasta...")

    fasta_out = os.path.join(args.outdir, "yearly_consensus_sequences.fasta")

    with open(fasta_out, "w") as f:
        for year in years:
            sub_long = long_df[long_df["year"] == year].copy()
            sub_long = sub_long.set_index("position_label").loc[position_order].reset_index()

            seq = "".join(sub_long["consensus_aa"].tolist())

            n_seq = df[df["year"] == year].shape[0]

            f.write(f">consensus|{year}|n={n_seq}\n")

            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")

    print(f"Yearly consensus fasta saved: {fasta_out}")

    print("[7] Exporting yearly sequence counts...")

    year_counts = (
        df.groupby("year")
        .size()
        .reset_index(name="n_sequences")
    )

    year_counts_out = os.path.join(args.outdir, "yearly_sequence_counts.tsv")
    year_counts.to_csv(year_counts_out, sep="\t", index=False)

    print(f"Yearly sequence counts saved: {year_counts_out}")

    print("\nDone.")


if __name__ == "__main__":
    main()