## eDNA Illumina Metabarcoding Workflow
Author: Rodgee Mae Guden

This repository contains a modular Bash workflow for analysing Illumina eDNA metabarcoding data from African freshwater systems. The workflow is designed for multi-marker metabarcoding datasets targeting COI, 12S rRNA, 16S rRNA and 18S rRNA.
The pipeline follows the same step-based structure as the eDNA metagenomics workflow, where each step can be run independently using the `-s` argument.

---
Pipeline script:
```bash
edna_metabarcoding_pipeline.sh
```
Make the script executable:
```bash
chmod +x edna_metabarcoding_pipeline.sh
```
---
## Workflow overview

The workflow includes the following steps:

| Step | Description |
|---|---|
| `check` | Check input files and dependencies. |
| `concat` | Concatenate lane-level paired-end reads per sample. |
| `qc_raw` | Run FastQC on concatenated raw reads. |
| `fastp` | Perform quality filtering and adapter trimming. |
| `qc_filtered` | Run FastQC on filtered reads. |
| `demux` | Demultiplex reads by marker using cutadapt. |
| `asv` | Infer ASVs using QIIME2/DADA2. |
| `otu` | Cluster OTUs at 97% identity using VSEARCH. |
| `taxonomy_asv` | Assign taxonomy to ASVs using BLAST and optionally TaxonKit. |
| `taxonomy_otu` | Assign taxonomy to OTUs using BLAST and optionally TaxonKit. |
| `multiqc` | Summarise QC and log files with MultiQC. |
| `all` | Run the full pipeline. |

---
## Required input files

### 1. Raw reads
Raw Illumina paired-end reads should be stored in `READS_DIR`.
The pipeline assumes that each indexed sample may have multiple lane-level FASTQ files. The `concat` step combines them into one R1 and one R2 file per sample.
Example raw files:
```text
S16_L001_R1.fastq.gz
S16_L001_R2.fastq.gz
S16_L002_R1.fastq.gz
S16_L002_R2.fastq.gz
```
### 2. Sample metadata file

Default:
```bash
/staging/leuven/stg_00100/Rodgee/metabarcoding_data/metadata/samples.tsv
```
Recommended format:
```text
sample_id	r1_glob	r2_glob
S16	S16*_R1*.fastq.gz	S16*_R2*.fastq.gz
S17	S17*_R1*.fastq.gz	S17*_R2*.fastq.gz
S18	S18*_R1*.fastq.gz	S18*_R2*.fastq.gz
```
The explicit glob columns are recommended because they avoid accidental matching of similar sample IDs, for example `S1` and `S16`.
A minimal format is also accepted:
```text
sample_id
S16
S17
S18
```
### 3. Primer files
Primer FASTA files should be stored in `PRIMER_DIR`.
Expected files:
```text
COI_F.fasta
COI_R.fasta
12S_F.fasta
12S_R.fasta
16S_F.fasta
16S_R.fasta
18S_F.fasta
18S_R.fasta
```
Example:
```text
>COI_forward
GGDRCWGGWTGAACWGTWTAYCCNCC
```
---
## Software requirements

The following tools are required:

- FastQC
- MultiQC
- fastp
- cutadapt
- VSEARCH
- BLAST+
- QIIME2, for ASV analysis
- TaxonKit, optional but recommended

The ASV workflow requires **QIIME2** with the **DADA2 plugin**.

---
### Configuration
The script uses environment variables for configuration.
Main variables:
```bash
READS_DIR=/path/to/raw_reads
METADATA=/path/to/samples.tsv
PRIMER_DIR=/path/to/primers
OUTPUT_DIR=/path/to/output
THREADS=32
MODE=both
```
Taxonomy variables:
```bash
BLAST_DB=/path/to/blast_db/custom_silva_ncbi_metabarcoding
CUSTOM_DB_FASTA=/path/to/custom_database.fasta
TAXONKIT_DB=/path/to/taxdump
```
If `CUSTOM_DB_FASTA` is provided, the pipeline will build the BLAST database using `makeblastdb`.

---
### Running the pipeline
Run a single step
```bash
./edna_metabarcoding_pipeline.sh -s qc_raw
```
Run several steps
```bash
./edna_metabarcoding_pipeline.sh -s concat,qc_raw,fastp,qc_filtered,demux
```
Run the ASV route
```bash
MODE=asv ./edna_metabarcoding_pipeline.sh -s all
```
Run the OTU route
```bash
MODE=otu ./edna_metabarcoding_pipeline.sh -s all
```
Run both ASV and OTU routes
```bash
MODE=both ./edna_metabarcoding_pipeline.sh -s all
```
---
Example full run
```bash
READS_DIR=/staging/leuven/stg_00100/Rodgee/metabarcoding_data/raw_reads \
METADATA=/staging/leuven/stg_00100/Rodgee/metabarcoding_data/metadata/samples.tsv \
PRIMER_DIR=/staging/leuven/stg_00100/Rodgee/metabarcoding_data/primers \
OUTPUT_DIR=/staging/leuven/stg_00100/Rodgee/metabarcoding_data/output \
THREADS=32 \
MODE=both \
BLAST_DB=/staging/leuven/stg_00100/Rodgee/databases/metabarcoding/custom_silva_ncbi_metabarcoding \
TAXONKIT_DB=/staging/leuven/stg_00100/Rodgee/databases/taxdump \
./edna_metabarcoding_pipeline.sh -s all
```
---
### SLURM array usage
Several sample-level steps support SLURM array mode:
```text
concat
qc_raw
fastp
qc_filtered
demux
```
Example SLURM script:
```bash
#!/bin/bash -l
#SBATCH --job-name=edna_metabarcoding
#SBATCH --output=logs/edna_metabarcoding_%A_%a.out
#SBATCH --error=logs/edna_metabarcoding_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1-24

set -euo pipefail

READS_DIR=/path/to/raw_reads \
METADATA=/path/to/samples.tsv \
PRIMER_DIR=/path/to/primers \
OUTPUT_DIR=/path/to/output \
THREADS=${SLURM_CPUS_PER_TASK:-16} \
./edna_metabarcoding_pipeline.sh -s concat,fastp,demux
```
For global steps such as ASV inference, OTU clustering and taxonomy, it is usually better to run per marker or as separate non-array jobs.

---
### Output structure
```text
output/
├── concat/
├── qc/
│   ├── raw/
│   └── filtered/
├── fastp/
├── demux/
│   ├── COI/
│   ├── 12S/
│   ├── 16S/
│   └── 18S/
├── qiime/
│   ├── manifests/
│   ├── imported/
│   ├── dada2/
│   └── exported/
├── otu/
│   ├── COI/
│   ├── 12S/
│   ├── 16S/
│   └── 18S/
├── taxonomy/
├── multiqc/
└── logs/
```
---
### fastp settings
The quality-filtering step uses:
```bash
--detect_adapter_for_pe
--trim_poly_g
--cut_front
--cut_front_window_size 20
--cut_front_mean_quality 30
--average_qual 35
--length_required 50
--low_complexity_filter
--overrepresentation_analysis
```
These settings can be changed with:
```bash
AVERAGE_QUAL=35
LENGTH_REQUIRED=50
CUT_FRONT_WINDOW_SIZE=20
CUT_FRONT_MEAN_QUAL=30
```
---
### cutadapt settings
Marker demultiplexing uses marker-specific primer pairs and retains only reads matching both primers:
```bash
-e 0.10
-O 20
--discard-untrimmed
```
The default also includes:
```bash
--revcomp
```
If your cutadapt version does not support `--revcomp`, run with:
```bash
CUTADAPT_EXTRA="" ./edna_metabarcoding_pipeline.sh -s demux
```
---
### ASV workflow
The ASV workflow uses QIIME2/DADA2 separately for each marker.
DADA2 is run in paired-end mode without additional truncation:
```bash
--p-trunc-len-f 0
--p-trunc-len-r 0
```
This assumes that adapter trimming and quality filtering were already performed with fastp.
Main ASV outputs:
```text
output/qiime/dada2/<MARKER>_table.qza
output/qiime/dada2/<MARKER>_rep-seqs.qza
output/qiime/dada2/<MARKER>_stats.qza
output/qiime/exported/<MARKER>/dna-sequences.fasta
```
---
### OTU workflow
The OTU workflow uses VSEARCH.
For each marker, the pipeline:
- Merges paired-end reads.
- Converts merged reads to FASTA.
- Dereplicates full-length reads.
- Clusters OTUs at 97% identity.
- Generates an OTU table.

Main OTU outputs:
```text
output/otu/<MARKER>/<MARKER>_all_samples.fasta
output/otu/<MARKER>/<MARKER>_uniques.fasta
output/otu/<MARKER>/<MARKER>_otus97.fasta
output/otu/<MARKER>/<MARKER>_otu_table.tsv
```
---
### Taxonomic assignment
Taxonomic assignment is performed with BLAST.
Default thresholds:
```bash
-perc_identity 80
-evalue 0.001
-max_target_seqs 10
```
The best hit per ASV or OTU is selected by highest bitscore. If TaxonKit is available, full lineages are appended to the taxonomy table.
Main outputs:
```text
output/taxonomy/<MARKER>_ASV_blast.tsv
output/taxonomy/<MARKER>_ASV_best_hits.tsv
output/taxonomy/<MARKER>_ASV_taxonomy.tsv

output/taxonomy/<MARKER>_OTU97_blast.tsv
output/taxonomy/<MARKER>_OTU97_best_hits.tsv
output/taxonomy/<MARKER>_OTU97_taxonomy.tsv
```
---
### Reproducibility checklist
Record the following for each analysis:
```text
sample metadata file
primer sequences
software versions
pipeline version or Git commit
reference database names
reference database versions
database download dates
filtering thresholds
taxonomic assignment thresholds
```
Useful version commands:
```bash
fastqc --version
multiqc --version
fastp --version
cutadapt --version
vsearch --version
blastn -version
qiime --version
taxonkit version
```
---
## Citation and licence

This workflow was developed by Rodgee Mae Guden for eDNA metabarcoding analyses of African freshwater systems, in collaboration with KU Leuven and the Royal Museum for Central Africa.

The workflow and accompanying documentation are licensed under the **Creative Commons Attribution 4.0 International Licence (CC BY 4.0)**.

Under this licence, users may share and adapt the material for any purpose, provided that appropriate credit is given, a link to the licence is included, and any modifications are clearly indicated.

When using or adapting this workflow, please cite the associated Figshare repository, archived pipeline version, or related publication.

The external software tools and reference databases used by the pipeline remain subject to their respective licences and terms of use.

---

## Contact

**Rodgee Mae Guden**  
Postdoctoral researcher  
Laboratory of Ecology, Evolution and Biodiversity Conservation  
Department of Biology, KU Leuven  
Leuven, Belgium  

Email: `rodgee.guden@kuleuven.be`
