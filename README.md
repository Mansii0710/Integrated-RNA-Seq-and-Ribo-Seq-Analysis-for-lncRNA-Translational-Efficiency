# Estimating the Translational Efficiency of lncRNAs in Thyroid Tissue

This project integrates **RNA-seq** and **Ribo-seq** data from normal human thyroid tissue to estimate the **Translational Efficiency (TE)** of long non-coding RNAs (lncRNAs). Long non-coding RNAs were traditionally considered non-coding, but growing evidence shows that a subset of short open reading frames are capable of producing functional micropeptides. TE — defined as the ratio of ribosome occupancy to transcript abundance — is used here to identify lncRNAs with evidence of active translation.

The workflow is divided into three stages: RNA-seq analysis, Ribo-seq analysis, and TE calculation.
1. **RNA-seq Analysis** — quality control, alignment, transcript assembly, and quantification
2. **Ribo-seq Analysis** — quality control, rRNA removal, genome alignment, and P-site analysis
3. **Translational Efficiency Calculation** — integration of RNA-seq and Ribo-seq TPMs, filtering, and annotation
--------------------------------------------------

## Pipeline Workflow

<img width="1280" height="274" alt="manuscript figures2" src="https://github.com/user-attachments/assets/b65af23c-6469-436b-94a3-d832fc4ef46f" />

The diagram summarizes the parallel RNA‑seq and Ribo‑seq pipelines from setup and raw data acquisition to translational efficiency estimation.

--------------------------------------------------

## Repository Structure

```
.
├── rna_seq/
│   ├── 1.RNA_setup_env.sh        # Setup conda environment and install RNA-Seq tools
│   ├── 2.RNAprepro.sh           # Quality control (FastQC), trimming (Trim Galore), MultiQC report
│   ├── 3.RNApostpro.sh          # Alignment (STAR) and transcript assembly (StringTie)
│   ├── 4.salmonRNAseq.sh        # Quantification of transcript expression using Salmon (TPM values)
│   ├── 5.mergesalmon.py         # Merge TPM values from all samples into a single matrix
│   └── 6.removeintron2.py       # Filter transcripts to retain relevant lncRNA entries
│
├── ribo_seq/
│   ├── 1.Ribo_setup_env.sh      # Setup conda environment and install Ribo-Seq tools
│   ├── 2.ribopreqc.sh           # Initial quality control using FastQC and MultiQC
│   ├── 3.trimming.sh            # Adapter trimming and read length filtering using Cutadapt
│   ├── 4.sortmernaindex.txt     # Instructions for building SortMeRNA rRNA reference index
│   ├── 5_r_rnaremoval.sh        # Removal of rRNA contamination using SortMeRNA
│   ├── 6.staralmn.sh            # Alignment of reads using STAR and BAM processing with SAMtools
│   └── 7.ribowaltzfinal.R       # P-site detection and ribosome occupancy calculation using Ribowaltz
│
├── te_calculation/
│   ├── 1.TEnew.py               # Calculate translational efficiency (TE) from RNA and Ribo TPM values
│   └── 2.TEannotate.py          # Annotate TE results with gene metadata (GTF annotation)
│
└──README.md  
```
--------------------------------------------------

## **Datasets**

### RNA-seq
8 paired-end human normal thyroid tissue RNA-seq samples (Total RNA) were provided by professor at the University.

### Ribo-seq
4 single-end Ribo-seq samples from normal human thyroid tissue were obtained from NCBI GEO under accession **[GSE212031](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212031)**.

> Although the raw data was generated as paired-end, only single-end reads were used, as Ribo-seq analysis is optimized for single-end footprint data.

| Sample | Download |
|--------|----------|
| SRR21198995 | (ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR211/095/SRR21198995/SRR21198995_1.fastq.gz) |
| SRR21198996 | (ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR211/096/SRR21198996/SRR21198996_1.fastq.gz) |
| SRR21198997 | (ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR211/097/SRR21198997/SRR21198997_1.fastq.gz) |
| SRR21198998 | (ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR211/098/SRR21198998/SRR21198998_1.fastq.gz) |

--------------------------------------------------

## **Environment Setup**

### **RNA-seq Environment**
```
bash rna_seq/1.RNA_setup_env.sh 
```

| Tool | Purpose |
|------|---------|
| STAR 2.7.10b | Splice-aware genome alignment |
| StringTie | Transcript assembly and abundance estimation |
| FastQC | Per-read quality assessment |
| Trim Galore + Cutadapt 5.1 | Adapter and quality trimming |
| GFFread / GFFcompare | GTF-to-FASTA conversion and annotation |
| MultiQC | Aggregated QC reporting |
| SAMtools | BAM sorting and indexing |
| Salmon | Quasi-mapping transcript quantification |

### Ribo-seq Environment
```
bash ribo_seq/1.Ribo_setup_env.sh 
```

| Tool | Purpose |
|------|---------|
| FastQC 0.12.1 | Quality assessment of raw reads |
| Cutadapt 4.5 | 3′-adapter trimming for ribosome footprints |
| STAR 2.7.10b | Genome alignment |
| SAMtools 1.17 | BAM sorting and indexing |
| SortMeRNA 4.3.4 | rRNA contamination removal |
| UMI-tools 1.1.4 | UMI deduplication |
| MultiQC 1.15 | Aggregated QC reporting |
| RiboWaltz (R) | P-site offset estimation and footprint TPM |
| Ribotricer / Ribo-TISH | Ribosome periodicity and translation detection |

--------------------------------------------------

## Reference Files

| File | Source |
|------|--------|
| GRCh38 genome FASTA | [GENCODE](https://www.gencodegenes.org/human/) |
| Full annotation GTF (v48) | [GENCODE](https://www.gencodegenes.org/human/) |
| lncRNA annotation GTF (v48) | [GENCODE](https://www.gencodegenes.org/human/) |
| 18S + 28S rRNA sequences | [SILVA Release 138.2](https://www.arb-silva.de/) |
| 5S rRNA | [Rfam RF00001](https://rfam.xfam.org/family/RF00001) |
| 5.8S rRNA | [Rfam RF00002](https://rfam.xfam.org/family/RF00002) |

--------------------------------------------------

## Workflow Summary

### STAR Index Generation

Before alignment, the genome index must be generated using STAR. This step is performed only once.
The same STAR genome index generated for RNA-seq is reused for Ribo-seq alignment to ensure consistency between datasets.

```bash
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles GRCh38.fa \
     --sjdbGTFfile gencode.v48.annotation.gtf \
     --sjdbOverhang 100
```
--------------------------------------------------


**RNA-Seq Workflow**

1. **Raw Data Handling**
   - FASTQ files obtained and organized

2. **Quality Control**
   - FastQC used to evaluate base quality, GC content, and adapter contamination
   - MultiQC aggregates all QC reports into a single summary

3. **Trimming**
   - Trim Galore (Cutadapt-based) removes adapters and low-quality bases
   - Improves alignment accuracy

4. **Alignment**
   - STAR aligner maps reads to GRCh38 reference genome
   - Produces sorted BAM files

5. **Transcript Assembly**
   - StringTie reconstructs transcript structures
   - Distinguishes known vs novel transcripts

6. **Transcript Extraction**
   - GFFread converts GTF → FASTA
   - Required for downstream quantification

7. **Quantification**
   - Salmon calculates TPM values for each transcript
   - Provides transcript-level expression estimates

8. **Merging TPM Data**
   - Python script merges TPM values from all samples
   - Creates a unified transcript expression matrix

--------------------------------------------------

**Ribo-Seq Workflow**

1. **Raw Data Handling**
   - FASTQ files obtained and organized

2. **Quality Control**
   - FastQC identifies quality issues and adapter contamination
   - MultiQC provides summary reports across samples

3. **Adapter Trimming**
   - Cutadapt removes 3' adapters
   - Filters reads to ribosome footprint length (18–32 nt)

4. **rRNA Removal**
   - SortMeRNA removes contaminating rRNA reads
   - Uses SILVA (18S/28S) and Rfam (5S/5.8S) databases (explained in file - ribo_seq/4.sortmernaindex.txt)
   - Ensures only translation-relevant reads are retained

5. **Alignment**
   - STAR maps reads to genome with Ribo-Seq-specific parameters:
     - Unique mapping only
     - No intron gaps
     - End-to-end alignment

6. **Post-processing**
   - SAMtools sorts and indexes BAM files
   - Ensures compatibility for downstream analysis

7. **Ribosome Profiling Analysis**
   - Ribowaltz identifies P-site positions
   - Generates ribosome occupancy data per transcript

8. **Quantification**
   - Ribosome-protected fragments converted to TPM values
   - Reflects translation activity

--------------------------------------------------

**Translational Efficiency (TE) Calculation**

- RNA-Seq TPM and Ribo-Seq TPM datasets are merged
- Only common transcripts retained
- Low expression filtering applied:
  - RNA TPM ≥ 1
  - Ribo TPM ≥ 0.5

Formula:
TE = (Ribo TPM + 1) / (RNA TPM + 1)

- Pseudo-count avoids division errors
- Output ranked by TE

--------------------------------------------------

**Annotation**

- lncRNA transcripts annotated using GENCODE reference
- Information added:
  - Transcript ID
  - Gene ID
  - Gene Name
  - Gene Type

- Annotation step improves biological interpretation
- Enables linking TE values to functional genes

--------------------------------------------------

## Input / Output Summary

| Stage | Input | Output |
|-------|-------|--------|
| RNA-seq pre-processing | Raw paired-end FASTQ | Trimmed FASTQ + QC reports |
| RNA-seq alignment & assembly | Trimmed FASTQ | Sorted BAMs, assembled GTFs, transcript FASTAs |
| Salmon quantification | Trimmed FASTQ + lncRNA index | Per-sample TPM files → merged TPM matrix |
| Ribo-seq trimming | Raw FASTQ | Trimmed FASTQ (18–32 nt) |
| rRNA removal | Trimmed FASTQ | Non-rRNA FASTQ |
| Ribo-seq alignment | Non-rRNA FASTQ | Sorted genome + transcriptome BAMs |
| RiboWaltz | Transcriptome BAM + lncRNA GTF | P-site tables, TPM matrix, QC plots |
| TE calculation | RNA TPM + Ribo TPM | Filtered TE table |
| TE annotation | TE table + lncRNA GTF | Final annotated TE table |

--------------------------------------------------

## **Requirements**

- Linux / HPC environment
- Conda / Miniconda
- Python 3.9
- R (for Ribowaltz)

--------------------------------------------------

## Learning Outcomes

Working through this project develops practical skills in:

- Setting up reproducible bioinformatics environments using Conda and Mamba on HPC systems
- Running end-to-end RNA-seq and Ribo-seq pipelines from raw data to quantification
- Understanding the key differences between RNA-seq and Ribo-seq analysis requirements
- Performing quality control and adapter trimming with Ribo-seq-specific parameters
- Removing rRNA contamination using curated databases (SILVA, Rfam) and SortMeRNA
- Aligning reads with STAR using both general and ribosome-profiling-specific parameters
- Estimating P-site offsets and validating ribosome footprint periodicity using RiboWaltz
- Quantifying lncRNA expression with Salmon using a filtered reference
- Integrating multi-omics data (RNA-seq + Ribo-seq) to compute translational efficiency
- Interpreting TE values in the context of lncRNA biology and potential micropeptide coding
