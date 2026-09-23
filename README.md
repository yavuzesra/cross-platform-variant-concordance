# Illumina–ONT Small-Variant Concordance

[![tests](https://github.com/yavuzesra/cross-platform-variant-concordance/actions/workflows/tests.yml/badge.svg)](https://github.com/yavuzesra/cross-platform-variant-concordance/actions/workflows/tests.yml)

A Bash/Python workflow for comparing normalized small-variant calls from matched Illumina and Oxford Nanopore Technologies (ONT) sequencing data within a defined target region.

The workflow was developed as part of the BSc thesis **“Evaluation of a Hybrid-Capture Long-Read Sequencing Workflow for Autoinflammatory Disorders.”** It was used to evaluate agreement between an existing Illumina workflow and targeted ONT sequencing for five clinical workflow-application samples.

## What the workflow measures

For each matched sample pair, the workflow:

1. validates the reference, target BED, VCFs and ONT alignment;
2. restricts both VCFs to `FILTER=PASS`, primary chromosomes and the target BED;
3. normalizes records with `bcftools norm` against GRCh38;
4. matches variants by `CHROM:POS:REF:ALT`;
5. compares genotypes for shared records while ignoring phasing order;
6. identifies Illumina-only and ONT-only records;
7. measures ONT alignment depth at Illumina-only positions;
8. produces per-pair and combined summary tables.

This is a **VCF-to-VCF technical concordance analysis**, not a benchmark against an independent truth set.

## Study-level results

Across the five matched clinical sample pairs used in the thesis evaluation:

| Metric | Result |
|---|---:|
| Illumina PASS records within the comparison territory | 2,772 |
| Records also present in ONT | 2,594 |
| Illumina-to-ONT record overlap | 93.6% |
| Shared records with concordant genotype | 2,491 |
| Shared records with discordant genotype | 103 |
| Genotype concordance among shared records | 96.0% |

Patient-level VCFs, alignments, sample identifiers and variant-level result tables are not distributed in this public repository.

## Comparison logic

### Comparison territory

The included `resources/autoinflammatory.bed` defines the targeted autoinflammatory panel regions used in the project. Only records within these regions are compared.

### Variant preparation

Both platforms are processed using the same preparation rules:

- `FILTER=PASS` only;
- primary chromosomes (`chr1`–`chr22`, `chrX`, `chrY`);
- restriction to the target BED;
- multiallelic splitting;
- left normalization and reference checking with `bcftools norm`.

A normalized variant is identified by:

```text
CHROM + POS + REF + ALT
```

### Genotype concordance

Phasing is ignored for this comparison. For example, `0/1`, `0|1` and `1|0` are treated as the same heterozygous genotype.

### Interpretation of platform-only records

Illumina-only records are additionally annotated with ONT alignment depth at the same genomic position, using a project threshold of **20×**.

ONT-only records are reported but are not interpreted as Illumina false negatives because Illumina BAM/CRAM files, gVCFs and complete callable-region information were not available for this analysis.

## Repository structure

```text
.
├── config/
│   └── config.example.sh
├── examples/
│   ├── illumina.example.vcf
│   └── ont.example.vcf
├── metadata/
│   └── sample_pairs.example.tsv
├── resources/
│   └── autoinflammatory.bed
├── scripts/
│   ├── combine_summaries.py
│   ├── compare_small_variants.py
│   ├── prepare_vcfs.sh
│   ├── run_all.sh
│   └── validate_inputs.py
├── tests/
├── environment.yml
└── README.md
```

## Requirements

The original analysis used:

- Python 3.13.13
- bcftools 1.19
- samtools 1.19.2
- bgzip/tabix from htslib 1.19
- Bash

A Conda environment specification is provided in `environment.yml`.

```bash
conda env create -f environment.yml
conda activate illumina-ont-concordance
```

## Configuration

Copy the example configuration and metadata files:

```bash
cp config/config.example.sh config/config.sh
cp metadata/sample_pairs.example.tsv metadata/sample_pairs.tsv
```

Edit `config/config.sh` to point to the GRCh38 no-alt reference FASTA. Edit `metadata/sample_pairs.tsv` to provide the matched Illumina VCF, ONT VCF and indexed ONT BAM/CRAM paths for each pair.

The expected metadata columns are:

```text
sample_key
comparison_id
illumina_sample_id
ont_sample_id
illumina_vcf
ont_vcf
ont_alignment
```

The files in `examples/` are synthetic format examples only; they are not a complete end-to-end validation dataset.

## Run the workflow

```bash
bash scripts/run_all.sh
```

The pipeline writes generated files to `work/`, `logs/` and `results/`. These directories are intentionally excluded from version control.

## Main outputs

For each comparison the workflow generates variant-level tables including shared records, genotype-discordant shared records, Illumina-only records and ONT-only records, plus a summary table. `combine_summaries.py` creates a cohort-level concordance summary across all comparisons.

## Testing

Portable unit tests cover genotype normalization, variant classification, missing-value handling and summary aggregation:

```bash
pytest -q
```

GitHub Actions runs these tests and syntax checks on every push and pull request. End-to-end execution additionally requires the external bioinformatics tools and user-provided sequencing inputs.

## Limitations

- The analysis compares VCF records, not raw-read evidence across both platforms.
- ONT-only records cannot be classified as Illumina false negatives without Illumina alignment/callability data.
- `FILTER=PASS` indicates that a call passed the originating workflow filters; it does not establish clinical relevance or independent truth.
- Variant concordance depends on consistent genome build, normalization and comparison territory.

## Data availability

Clinical sequencing files and patient-level result tables are not included. Synthetic metadata and VCF examples are provided to document the expected input format.
