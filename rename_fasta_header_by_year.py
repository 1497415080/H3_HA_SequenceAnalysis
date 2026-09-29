import csv
import re

csv_file = "HA_20260920-Nvirus.csv"
fasta_file = "HA_20260920-Nvirus.fasta"
out_file = "HA_20260920-Nvirus.year.1968-2026.fasta"

MIN_YEAR = 1968
MAX_YEAR = 2026


def extract_year(date_text):
    """
    从 Collection_Date 中灵活提取年份。

    支持：
    1968
    5/31/2023
    2023-05-31
    2025-12
    2025/12
    2025.12
    Dec-2025
    2025 Dec

    不保留：
    NA
    N/A
    unknown
    空值
    无法提取年份
    """

    if date_text is None:
        return None

    date_text = str(date_text).strip()

    if not date_text:
        return None

    # 常见缺失值
    na_values = {
        "NA", "N/A", "na", "n/a",
        "None", "none",
        "NULL", "null",
        "Unknown", "unknown",
        "-", "--"
    }

    if date_text in na_values:
        return None

    # 提取 4 位年份，支持 1968、5/31/2023、2025-12 等
    match = re.search(r"\b(19\d{2}|20\d{2})\b", date_text)

    if not match:
        return None

    year = int(match.group(1))

    # 只保留 1968-2026
    if MIN_YEAR <= year <= MAX_YEAR:
        return str(year)

    return None


# 读取 CSV，建立 Accession -> year 映射
acc_to_year = {}

total_csv = 0
valid_csv = 0
invalid_csv = 0

with open(csv_file, "r", encoding="utf-8-sig", newline="") as f:
    reader = csv.DictReader(f)

    for row in reader:
        total_csv += 1

        acc = row.get("Accession", "").strip()
        collection_date = row.get("Collection_Date", "").strip()

        year = extract_year(collection_date)

        if acc and year:
            valid_csv += 1

            # CSV 中一般是不带版本号的 accession，例如 YFY77345
            acc_to_year[acc] = year

            # 兼容 CSV 中可能带版本号的情况，例如 YFY77345.1
            acc_to_year[acc.split(".")[0]] = year
        else:
            invalid_csv += 1


def get_accession_from_header(header_line):
    """
    从 FASTA 表头提取 accession。

    示例：
    >YFY77345.1 |hemagglutinin [Influenza A virus]

    返回：
    YFY77345.1
    """
    header = header_line[1:].strip()
    acc_full = header.split()[0]
    return acc_full


kept = 0
removed = 0
removed_no_year_or_out_range = 0

current_header = None
current_seq_lines = []

with open(fasta_file, "r", encoding="utf-8") as fin, \
     open(out_file, "w", encoding="utf-8") as fout:

    def process_record(header, seq_lines):
        global kept, removed, removed_no_year_or_out_range

        if header is None:
            return

        acc_full = get_accession_from_header(header)
        acc_no_version = acc_full.split(".")[0]

        year = acc_to_year.get(acc_full) or acc_to_year.get(acc_no_version)

        if year:
            fout.write(f">{acc_full} | {year}\n")
            for seq_line in seq_lines:
                fout.write(seq_line)
            kept += 1
        else:
            # 不在 CSV、Collection_Date 为 NA、无年份、年份不在 1968-2026 的序列都不保留
            removed += 1
            removed_no_year_or_out_range += 1

    for line in fin:
        if line.startswith(">"):
            # 处理上一条记录
            process_record(current_header, current_seq_lines)

            # 开始新记录
            current_header = line.rstrip("\n")
            current_seq_lines = []
        else:
            current_seq_lines.append(line)

    # 处理最后一条记录
    process_record(current_header, current_seq_lines)


print("完成")
print(f"输出文件：{out_file}")
print()
print("CSV 统计：")
print(f"CSV 总记录数：{total_csv}")
print(f"CSV 有效年份且在 {MIN_YEAR}-{MAX_YEAR} 范围内：{valid_csv}")
print(f"CSV 无效年份、NA、空值或超出范围：{invalid_csv}")
print()
print("FASTA 统计：")
print(f"保留序列数：{kept}")
print(f"删除序列数：{removed}")