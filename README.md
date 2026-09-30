# H3 HA Sequence Analysis

This repository contains data and code associated with a study of influenza A(H3N2) hemagglutinin (HA) sequences. It is intended to document the input data, sequence-processing steps, and analyses used in the study.

> **Manuscript status:** [In preparation]  
> **Associated publication:** [Add the paper title, citation, and DOI when available]

## Repository contents

| File or directory | Description |
| --- | --- |
| `HA_20260920-Nvirus.csv` | [Describe the columns and whether this is an input or an output.] |
| `HA_20260920-Nvirus.fasta` | [Describe which sequences this file contains and how it was obtained.] |
| `HA_20260920-Nvirus.year.fasta` | [Explain how this differs from the FASTA file above.] |
| `HA_20260920-Nvirus.year.1968-2026.fasta` | [Explain the year range and how this file was generated.] |
| `year_consensus/` | [Describe the files in this directory and how they were generated.] |
| `rename_fasta_header_by_year.py` | [Describe its input and output.] |
| `yearly_consensus_analysis.py` | [Describe its input and output.] |
| `coevolution_mi.py` | [Describe the analysis it performs and its output.] |
| `FinalCode.ipynb` | [Describe the analyses or figures produced by the notebook.] |

## Data source

The HA sequences were obtained from [database name and URL] on [access date].

- Search criteria: [Describe the query or download criteria.]
- Inclusion criteria: [Describe which sequences were retained.]
- Exclusion criteria: [Describe which sequences were removed and why.]
- Sequence identifiers: [Explain where accession numbers or other identifiers can be found.]
- Data-use terms: [Add the database's citation and redistribution requirements.]

**Note:** Before redistributing sequence data, users should check the terms of the original data source.

## Analysis workflow

The intended order of analysis is:

1. Start with `[input filename]`.
2. Run `rename_fasta_header_by_year.py` to [describe what it does].
3. Run `yearly_consensus_analysis.py` to [describe what it does].
4. Run `coevolution_mi.py` to [describe what it does].
5. Open `FinalCode.ipynb` to [describe the remaining analyses and outputs].

[Correct the order above to match the workflow you actually used. Add example commands once the scripts' arguments are confirmed.]

## Software requirements

- Python: [version]
- Required Python packages: [to be listed in `requirements.txt`]

To install the dependencies once `requirements.txt` has been added:

```bash
python -m pip install -r requirements.txt
```

## Reproducing the results

1. Download or prepare the input data as described in **Data source**.
2. Install the software requirements.
3. Run the scripts in the order described in **Analysis workflow**.
4. Run `FinalCode.ipynb` from start to finish.
5. Compare the generated files with [specify the expected outputs or manuscript figures/tables].

[Document any manual processing steps, random seeds, required directory paths, or external software here.]

## Citation

If you use this repository, please cite:

> [Authors]. [Article title]. [Journal or preprint server], [year]. [DOI or URL]

If the associated article has not yet been published, this section will be updated when citation information becomes available.

## License and data use

Code license: [Specify a license, if you choose to provide one.]

Sequence data may be subject to the terms of the original data provider. See **Data source** before reusing or redistributing the data.

## Contact

For questions about this repository, contact [name or email, optional].
