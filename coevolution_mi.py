#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import numpy as np


AA = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {a: i for i, a in enumerate(AA)}


def read_fasta(path):
    names = []
    seqs = []

    name = None
    seq_chunks = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    names.append(name)
                    seqs.append("".join(seq_chunks).upper())

                name = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if name is not None:
            names.append(name)
            seqs.append("".join(seq_chunks).upper())

    if len(seqs) == 0:
        raise ValueError("没有读到任何 FASTA 序列")

    lengths = [len(s) for s in seqs]
    if len(set(lengths)) != 1:
        raise ValueError("输入 FASTA 必须是已经比对好的序列，所有序列长度应一致")

    return names, seqs


def seqs_to_index_array(seqs):
    """
    将氨基酸序列转换成整数矩阵：
    A,C,D... -> 0,1,2...
    gap、X、B、Z 等非标准字符 -> -1
    """
    arr = np.array([list(s) for s in seqs])
    idx = np.full(arr.shape, -1, dtype=np.int16)

    for aa, k in AA_TO_IDX.items():
        idx[arr == aa] = k

    return idx


def compute_sequence_weights(idx, theta=0.8):
    """
    根据序列相似度计算权重，降低冗余序列影响。

    如果两条序列的相似度 >= theta，则认为它们属于同一相似序列簇。
    每条序列权重 = 1 / 与其相似的序列数量。
    """
    n_seq, length = idx.shape
    neighbors = np.ones(n_seq, dtype=np.float64)

    for i in range(n_seq):
        for j in range(i + 1, n_seq):
            valid = (idx[i] >= 0) & (idx[j] >= 0)

            if valid.sum() == 0:
                sim = 0.0
            else:
                sim = np.mean(idx[i, valid] == idx[j, valid])

            if sim >= theta:
                neighbors[i] += 1
                neighbors[j] += 1

    weights = 1.0 / neighbors
    return weights


def consensus_residue(col, weights):
    counts = np.zeros(len(AA), dtype=np.float64)
    valid = col >= 0

    if valid.sum() == 0:
        return "X"

    np.add.at(counts, col[valid], weights[valid])
    return AA[int(np.argmax(counts))]


def mutual_information(col_i, col_j, weights, pseudo=0.01, min_eff_pair=5.0):
    """
    计算两个比对列之间的加权 MI。

    对于含 gap 或非标准氨基酸的位置，该序列在该位点对中被忽略。
    """
    q = len(AA)
    valid = (col_i >= 0) & (col_j >= 0)

    if valid.sum() == 0:
        return np.nan, 0.0

    w = weights[valid]
    eff_n = w.sum()

    if eff_n < min_eff_pair:
        return np.nan, eff_n

    counts = np.zeros((q, q), dtype=np.float64)

    ai = col_i[valid]
    aj = col_j[valid]

    np.add.at(counts, (ai, aj), w)

    pxy_obs = counts / eff_n

    # 加入少量伪计数，避免低样本下概率为 0
    if pseudo > 0:
        pxy = (1.0 - pseudo) * pxy_obs + pseudo / (q * q)
    else:
        pxy = pxy_obs

    px = pxy.sum(axis=1)
    py = pxy.sum(axis=0)

    denom = np.outer(px, py)

    mask = (pxy > 0) & (denom > 0)

    mi = np.sum(pxy[mask] * np.log2(pxy[mask] / denom[mask]))

    return mi, eff_n


def apc_correction(mi_matrix):
    """
    Average Product Correction, APC 校正。
    用于降低 MI 中由于保守性、系统发育背景等带来的假阳性。
    """
    m = mi_matrix.copy()

    row_mean = np.nanmean(m, axis=1)
    global_mean = np.nanmean(m)

    if not np.isfinite(global_mean) or global_mean == 0:
        return m

    apc = m - np.outer(row_mean, row_mean) / global_mean
    return apc


def zscore_matrix(score_matrix):
    vals = score_matrix[np.triu_indices_from(score_matrix, k=1)]
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return np.full_like(score_matrix, np.nan)

    mean = np.mean(vals)
    std = np.std(vals)

    if std == 0:
        return np.full_like(score_matrix, np.nan)

    return (score_matrix - mean) / std


def main():
    parser = argparse.ArgumentParser(
        description="Detect coevolving amino-acid positions from an aligned FASTA file using MI/APC."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="已经比对好的 FASTA 文件"
    )

    parser.add_argument(
        "-o", "--output",
        default="coevolution_pairs.csv",
        help="输出 CSV 文件，默认：coevolution_pairs.csv"
    )

    parser.add_argument(
        "--max-gap",
        type=float,
        default=0.5,
        help="过滤 gap 比例高于该值的列，默认：0.5"
    )

    parser.add_argument(
        "--theta",
        type=float,
        default=0.8,
        help="序列权重计算中的相似度阈值，默认：0.8"
    )

    parser.add_argument(
        "--no-weights",
        action="store_true",
        help="不进行序列权重校正"
    )

    parser.add_argument(
        "--pseudo",
        type=float,
        default=0.01,
        help="MI 计算中的伪计数比例，默认：0.01"
    )

    parser.add_argument(
        "--min-eff-pair",
        type=float,
        default=5.0,
        help="某个位点对至少需要的有效序列数，默认：5.0"
    )

    parser.add_argument(
        "--z",
        type=float,
        default=3.0,
        help="输出的最小 Z-score 阈值，默认：3.0；若想输出所有位点对可设为 -999"
    )

    parser.add_argument(
        "--top",
        type=int,
        default=200,
        help="最多输出前多少个位点对，默认：200；设为 0 表示全部输出"
    )

    args = parser.parse_args()

    names, seqs = read_fasta(args.input)
    idx = seqs_to_index_array(seqs)

    n_seq, aln_len = idx.shape

    print(f"读取序列数：{n_seq}")
    print(f"比对长度：{aln_len}")

    # 根据 gap 比例过滤列
    valid_fraction = np.mean(idx >= 0, axis=0)
    keep_cols = valid_fraction >= (1.0 - args.max_gap)

    if keep_cols.sum() < 2:
        raise ValueError("过滤后剩余列数少于 2，请调大 --max-gap")

    original_positions = np.arange(1, aln_len + 1)[keep_cols]
    idx = idx[:, keep_cols]

    n_col = idx.shape[1]

    print(f"过滤后用于分析的位点数：{n_col}")

    # 序列权重
    if args.no_weights:
        weights = np.ones(n_seq, dtype=np.float64)
    else:
        print("计算序列权重...")
        weights = compute_sequence_weights(idx, theta=args.theta)

    neff = weights.sum()
    print(f"有效序列数 Neff：{neff:.2f}")

    # 每个位点的 consensus 氨基酸
    consensus = [
        consensus_residue(idx[:, k], weights)
        for k in range(n_col)
    ]

    # 计算 MI 矩阵
    mi_matrix = np.full((n_col, n_col), np.nan, dtype=np.float64)
    eff_pair_matrix = np.zeros((n_col, n_col), dtype=np.float64)

    print("计算位点对 MI...")

    for i in range(n_col):
        if (i + 1) % 20 == 0 or i == n_col - 1:
            print(f"  进度：{i + 1}/{n_col}")

        for j in range(i + 1, n_col):
            mi, eff_n = mutual_information(
                idx[:, i],
                idx[:, j],
                weights,
                pseudo=args.pseudo,
                min_eff_pair=args.min_eff_pair
            )

            mi_matrix[i, j] = mi
            mi_matrix[j, i] = mi

            eff_pair_matrix[i, j] = eff_n
            eff_pair_matrix[j, i] = eff_n

    # APC 校正
    mi_apc = apc_correction(mi_matrix)

    # Z-score
    zmat = zscore_matrix(mi_apc)

    # 收集结果
    results = []

    for i in range(n_col):
        for j in range(i + 1, n_col):
            mi = mi_matrix[i, j]
            score = mi_apc[i, j]
            z = zmat[i, j]
            eff_n = eff_pair_matrix[i, j]

            if not np.isfinite(mi) or not np.isfinite(score) or not np.isfinite(z):
                continue

            if z < args.z:
                continue

            results.append({
                "pos_i": int(original_positions[i]),
                "pos_j": int(original_positions[j]),
                "aa_i": consensus[i],
                "aa_j": consensus[j],
                "MI_bits": mi,
                "MI_APC": score,
                "Z_score": z,
                "eff_pair_n": eff_n
            })

    # 按 Z-score 从高到低排序
    results.sort(key=lambda x: x["Z_score"], reverse=True)

    if args.top > 0:
        results = results[:args.top]

    # 输出 CSV
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "pos_i",
                "pos_j",
                "aa_i",
                "aa_j",
                "MI_bits",
                "MI_APC",
                "Z_score",
                "eff_pair_n"
            ]
        )

        writer.writeheader()

        for r in results:
            writer.writerow(r)

    print(f"完成。共输出 {len(results)} 个候选共进化位点对。")
    print(f"结果文件：{args.output}")


if __name__ == "__main__":
    main()