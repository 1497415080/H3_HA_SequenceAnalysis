# H3 HA Sequence Analysis

This repository contains H3 hemagglutinin (HA) sequence data and code associated with the study *A broadly neutralizing H3N2 antibody defines a conserved hemagglutinin head epitope breached by recent stepwise antigenic drift*.

> **Manuscript status:** In preparation

## Repository contents

| File or directory | Description |
| --- | --- |
| `HA_20260920-Nvirus.csv` | H3 sequence information downloaded from NCBI Virus. |
| `HA_20260920-Nvirus.fasta` | H3 sequences downloaded from NCBI Virus. |
| `HA_20260920-Nvirus.year.fasta` | H3 sequences with year annotations in their FASTA headers. |
| `HA_20260920-Nvirus.year.1968-2026.fasta` | H3 sequence file described as deduplicated and standardized for the 1968–2026 analysis. |
| `year_consensus/` | Directory used for analysis outputs. |
| `rename_fasta_header_by_year.py` | Script for FASTA header processing. |
| `yearly_consensus_analysis.py` | Calculates yearly consensus amino acids and frequencies from an aligned, year-annotated FASTA file. |
| `coevolution_mi.py` | Scores pairs of alignment positions using mutual information (MI), average product correction (APC), and Z-scores. |
| `FinalCode.ipynb` | Analysis notebook. |

## Data source

The H3 sequence-information CSV and FASTA file were downloaded from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/) as of September 2026.

The sequence files are derived from a public sequence database. Before reusing or redistributing them, consult the current requirements of NCBI Virus and the applicable requirements for the underlying records.

## Software requirements

Install the Python packages listed in `requirements.txt`:

```bash
python -m pip install -r requirements.txt
```

The dependency list covers the third-party imports in the scripts and import statements supplied for this repository. Package versions are not pinned.

## Sequence-analysis inputs

Both `yearly_consensus_analysis.py` and `coevolution_mi.py` require an **aligned FASTA file whose sequences all have the same length**. Equal length alone does not establish that sequences have been correctly aligned.

For yearly consensus analysis, each sequence header must also contain a four-digit year as the final pipe-delimited field of its first token, for example:

```text
>ADV76673|1968
>APO19006|2016
```

`yearly_consensus_analysis.py` excludes records whose headers do not contain a valid year in this format.

## Yearly consensus analysis

Run the script from the repository root, substituting an aligned, year-annotated amino-acid FASTA file for the input path:

```bash
python yearly_consensus_analysis.py \
  -i PATH_TO_ALIGNED_YEAR_ANNOTATED_FASTA \
  -o year_consensus
```

By default, the script uses the first sequence from the earliest year as its reference. Its position labels count residues in that reference after removing gaps. Alignment columns where the reference has a gap are excluded from the consensus outputs by default. Gaps and specified unknown amino-acid characters are ignored when calculating consensus residues.

The script writes the following files to the output directory:

| File | Contents |
| --- | --- |
| `alignment_to_reference_position_map.tsv` | Mapping between alignment positions and ungapped reference positions. |
| `yearly_consensus_long.tsv` | Consensus residues, frequencies, counts, and position information by year and position. |
| `yearly_consensus_AA_matrix.tsv` | Year-by-position consensus amino-acid matrix. |
| `yearly_consensus_percent_matrix.tsv` | Year-by-position consensus percentage matrix. |
| `yearly_consensus_fraction_matrix.tsv` | Year-by-position consensus fraction matrix. |
| `yearly_consensus_AA_percent_combined_matrix.tsv` | Year-by-position matrix combining consensus amino acids and percentages. |
| `yearly_consensus_sequences.fasta` | One consensus sequence for each year represented in the input. |
| `yearly_sequence_counts.tsv` | Number of retained sequences by year. |

The script generates tables and a FASTA file; it does not generate figures. Run `python yearly_consensus_analysis.py --help` for its available options.

## Coevolution analysis

Run the script on an aligned amino-acid FASTA file:

```bash
python coevolution_mi.py \
  -i PATH_TO_ALIGNED_AMINO_ACID_FASTA \
  -o coevolution_pairs.csv
```

The script calculates MI for pairs of alignment columns. By default, it:

- excludes columns with more than 50% gaps or nonstandard amino-acid characters;
- weights sequences using a similarity threshold of 0.8;
- uses a pseudocount proportion of 0.01;
- does not score position pairs with an effective sample size below 5;
- applies APC and calculates Z-scores; and
- exports up to 200 position pairs with a Z-score of at least 3, sorted by descending Z-score.

The output CSV contains `pos_i`, `pos_j`, `aa_i`, `aa_j`, `MI_bits`, `MI_APC`, `Z_score`, and `eff_pair_n`. `pos_i` and `pos_j` are **1-based positions in the input alignment**, not automatically positions in a reference sequence or an H3 numbering scheme.

Run `python coevolution_mi.py --help` for all available parameters.

## Citation and data use

The associated manuscript is in preparation:

> *A broadly neutralizing H3N2 antibody defines a conserved hemagglutinin head epitope breached by recent stepwise antigenic drift.*

Sequence data should be reused and redistributed in accordance with the terms applicable to their original source.
