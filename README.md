# QC Pipeline (`qcpipe_integrated.py`)

## Summary
`qcpipe_integrated.py` is an integrated quality control (QC) pipeline for targeted next-generation sequencing (tNGS) data. It performs end-to-end QC across FASTQ reads, BAM alignments, enrichment regions, insert size distributions, taxonomic classification, and optional [Castanet](https://github.com/MultipathogenGenomics/castanet) outputs.

It generates a **single per-sample QC report (`*_qc.csv`)** and optional diagnostic plots.

### Key features

- FASTQ read length QC (raw + trimmed)
- Alignment QC using `samtools`
- Insert size distribution analysis from BAM files
- Enrichment analysis using BED-defined regions, usually in mitochondrial regions
- Kraken2 taxonomic summarisation
- [Castanet](https://github.com/MultipathogenGenomics/castanet) BAM + depth integration
- Insert size visualisation plots

---

## Dependencies

* **`bwa-mem2`**
* **`samtools` ($\\ge$ v1.10)**
* **`seqkit`**
* ** Python 3.x Environment**: Requires `pandas`, `numpy`, and `matplotlib`.

---

## Command Line Usage Reference

```bash
python qcpipe_integrated.py --sample <SAMPLE_ID> --batch <BATCH_ID> --outdir <OUTPUT_DIR> [OPTIONS]
```

### Argument Parameters Definition

| Flag | Long Flag | Requirement | Description                                                                                                            |
| :--- | :--- | :---: |:-----------------------------------------------------------------------------------------------------------------------|
| `-s` | `--sample` | **Required** | Unique identifier name string for the processing sample.                                                               |
| `-b` | `--batch` | **Required** | Running batch code identity identifier.                                                                                |
| `-o` | `--outdir` | **Required** | Targeted target directory route path to dump generated files.                                                          |
| `-r` | `--raw` | Optional | Path to raw forward ($R_1$) and reverse ($R_2$) files (space separated).                                               |
| `-t` | `--trimmed` | Optional | Path to trimmed forward ($R_1$) and reverse ($R_2$) files (space separated).                                           |
| `-k` | `--kraken` | Optional | Input data file path pointing to a standard kraken2 report.                                                            |
| `-m` | `--mttarget` | Optional | Fasta file of enrichment target used in panel (select from enrichment_targets folder). (triggers enrichment analysis)  |
| `--bedfile` | `--bedfile` | Optional | BED file corresponding to --mttarget fasta file that contains enriched regions (select from enrichment_targets folder) |
| `--castanetbam` | `--castanetbam` | Optional | Path to castanet BAM file.                                                                                             |
| `--castanetdepth`| `--castanetdepth`| Optional | Path to castanet _depth.csv file.                                                                                      |
| | `--keeptmp` | Optional Flag| Prevents cleanup deletion routines targeting intermediate alignments.                                                  |

---

## Execution Examples

### Example 1: Basic Processing (Read Geometry Evaluation Only)
Evaluates basic spatial track sizes across processed paired data files without parsing deeper coordinate indexes.

```bash
python qcpipe_integrated.py \
  --sample SRR1234567 \
  --batch Run_2026_06_A \
  --trimmed data/trimmed_R1.fastq data/trimmed_R2.fastq \
  --outdir ./qc_results
```

### Example 2: Full Enrichment Pipeline with Background Classification
A complete pipeline run involving target alignments to an engine index (e.g., mitochondrial genome), sorting coordinate enrichment via a structural BED template, parsing taxonomy distributions, and extracting deduplication metadata tables.

```bash
python qcpipe_integrated.py \
  --sample Patient_Sample_04 \
  --batch Clinical_Batch_42 \
  --raw data/raw_R1.fastq.gz data/raw_R2.fastq.gz \
  --trimmed data/trimmed_R1.fastq.gz data/trimmed_R2.fastq.gz \
  --mttarget enrichment_targets/HumanMt.fasta \
  --bedfile enrichment_targets/HumanMt_regions.bed \
  --kraken kraken2_report.txt \
  --castanetbam sample.bam \
  --castanetdepth sample_depth.csv \
  --outdir ./qc_reports
  ```
  ---

## Output Statistics Reference Guide

The script generates a single comma-delimited output file (`[sample_id]_qc.csv`) containing headers mapped dynamically based on input flags. The full data schema is broken down below:

### 1. Base Sample & Sequencing Statistics
| CSV Header Block | Data Type | Description                                                                       |
| :--- | :---: |:----------------------------------------------------------------------------------|
| `batch` | String | User-specified batch run identifier.                                              |
| `sampleid` | String | Extracted or provided specimen id                                                 |
| `R1` / `R2` | String | path locations for forward and reverse raw reads.                                 |
| `trimmedR1` / `trimmedR2`| String | path locations for filtered forward and reverse trimmed reads.                    |
| `rawreads` | Integer | Total absolute read count combined across both raw paired-end fastq inputs.       |
| `trimmedreads` | Integer | Total absolute read count combined across both processed paired-end fastq inputs. |
| `avglen` / `medlen` | Float | mean and median read lengths for all trimmed reads.                               |
| `r1count` / `r2count` | Integer | Individual absolute processed read counts split across $R_1$ and $R_2$ files.     |
| `r1avglen` / `r1medlen` | Float | average and median length for $R_1$ reads.                                        |
| `r2avglen` / `r2medlen` | Float | average and median length for $R_2$ reads.                                        |

### 2. Enrichment & Target Insert Dynamics (Mitochondrial / Target Profiling)
*Generated when supplying `--mttarget` and `--bedfile` arguments.*

| CSV Header Block | Data Type | Description                                                                                                                                     |
| :--- | :---: |:------------------------------------------------------------------------------------------------------------------------------------------------|
| `enrichmentloci_avginsert`| Float | Mean fragment size from enrichment target mapped reads.                                                                                         |
| `enrichmentloci_stdinsert`| Float | Standard deviation fragment size from enrichment target mapped reads.                                                                           |
| `enrichmentloci_insert25` | Float | 25th Percentile fragment size from enrichment target mapped reads.                                                                              |
| `enrichmentloci_insert50` | Float | Median fragment size from enrichment target mapped reads.                                                                                       |
| `enrichmentloci_insert75` | Float | 75th Percentile fragment size from enrichment target mapped reads.                                                                              |
| `enrichedMedian` | Float | Computed median mapping depth across coordinates tagged as `enriched`.                                                                          |
| `unenrichedMedian` | Float | Computed median mapping depth across coordinates tagged as `unenriched`.                                                                        |
| `enrichmentRatio` | Float | Enrichment ratio $\\frac{\\text{enrichedMedian}}{\\text{unenrichedMedian}}$ (utilizes smooth scaling padding $+1$ if background drops below 1). |
| `enrichmentloci` | String | enrichment locus used (selected by highest read depth).                                                                                         |

### 3. Comprehensive Metagenomic Taxonomy Metrics
*Generated when supplying a `--kraken` classification file.*

| CSV Header Block | Data Type | Description                                               |
| :--- | :---: |:----------------------------------------------------------|
| `kraken:Eukaryota` | Integer | Read count for Eukaryota from krakenreport.               |
| `kraken:Bacteria` | Integer | Read count for Bacteria from krakenreport.                |
| `kraken:Archaea` | Integer | Read count for Archaea from krakenreport.                 |
| `kraken:Viruses` | Integer | Read count for Viruses from krakenreport.                 |
| `kraken:Fungi` | Integer | Read count for Fungi from krakenreport.                   |
| `kraken:Caudoviricetes`| Integer | Read count for Caudoviricetes (phages) from krakenreport. |
| `kraken:Homo sapiens` | Integer | Read count for Homo sapiens from krakenreport             |

### 4. Global Structural Mapping & Deduplication Profile
*Generated when supplying `--castanetbam` and `--castanetdepth` metrics.*

| CSV Header Block | Data Type | Description                                                                                 |
| :--- | :---: |:--------------------------------------------------------------------------------------------|
| `all_mapped_avginsert` | Float | Raw mean insert size from castanet bam mapped reads.                                        |
| `all_mapped_stdinsert` | Float | Raw standard deviation insert size from castanet bam mapped reads.                          |
| `all_mapped_insert[25/50/75]`| Float | Raw percentile metrics ($Q_1, Q_2, Q_3$) from castanet bam mapped reads.                    |
| `filtered_mapped_avginsert` | Float | High-confidence mean insert length (Filtered to MapQ $>20$ and fragments $>50\\text{ bp}$). |
| `filtered_mapped_stdinsert` | Float | High-confidence standard deviation of insert size.                                          |
| `filtered_mapped_insert[25/50/75]`| Float | Percentile metrics ($Q_1, Q_2, Q_3$) for filtered insert size.                              |
| `castanet_total_mapped_reads`| Integer | Sum of all reads mapped to all targets from castanet depths file.                           |
| `castanet_dedup_reads` | Integer | Sum of all deduplicated reads mapped to all targets from castanet depths file.              |
