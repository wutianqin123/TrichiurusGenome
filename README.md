# TrichiurusGenome
The script commands and usage methods that appear in the article "Chromosome-level Genome Assembly Resolves Taxonomic Conflicts in Trichiurus and Informs Its Conservation Strategy in the Northwest Pacific" are stored in this document.
———————————————————————————————————
# Chromosome-Level Genome Assembly Using Docker

1. Basic Commands

$ docker image ls # List all available Docker images

$ docker ps # List all running containers

$ docker ps -a # List all containers (both running and stopped)

$ exit # Exit the current container shell

2. Creating a New Container

$ docker run --gpus all --name your_name_assemble -idt centos:centos7 /bin/bash # Create a new container with GPU support

--name Assign a name to the container for easier management

$ docker exec -it your_name_assemble /bin/bash # Open an interactive shell inside a running container

$ docker start your_name_assemble # Start a stopped container

$ docker stop your_name_assemble # Stop a running container

3. Transferring Data

$ docker cp your_name_assemble:/root/test.text /home/vagrant/test.txt # Copy a file from the Docker container to the host machine

$ docker cp /home/vagrant/test.txt your_name_assemble:/root/test.text # Copy a file from the host machine into the Docker container

Note: Within the Docker environment, it is recommended to store data in the /data directory.

4. Running NextDenovo for Assembly
4.0. Estimate Genome Size (Execute this on the host system, NOT inside the Docker container)

$ jellyfish count -m 21 -s 20G -t 30 -o 21mer_out -C <(zcat sample1_1.fq.gz) <(zcat sample1_2.fq.gz) # Count k-mers in the sequencing files

-m K-mer length

-t Number of threads to use

-o Output file prefix

-C Count canonical k-mers (consider both strands)

The syntax < (zcat *.fq.gz) streams decompressed FASTQ files directly into Jellyfish.

$ jellyfish histo -o 21mer_out.histo 21mer_out # Generate a k-mer frequency histogram

$ /home/software/genomescope2.0/genomescope.R -i 21mer_out.histo -o 21mer_out_plot -k 21 # Estimate genome size and heterozygosity, and generate plots

-o Output directory; the primary visualization is typically named linear_plot.png

4.1. Obtain the full path to the sequencing data

$ realpath seq.fq

4.2. Create an input file of sequence paths (FOFN)

$ touch input_seq.txt

$ vi input_seq.txt

Add the full path obtained in step 4.1 to this file and save.

4.3. Configure the run.cfg file (Copy a template from /home/software/NextDenovo/test_data/run.cfg to your working directory)

$ vi run.cfg

Modify the following key parameters:

input_fofn = # Full path to the input_seq.txt file

workdir = /data/1_assemble # Full path for the working directory

genome_size = # Estimated genome size (obtained from Jellyfish/GenomeScope)

parallel_jobs = # Number of parallel jobs to run (each uses 8 threads; typically set between 5-10)

4.4. Execute NextDenovo

$ /home/software/NextDenovo/nextDenovo run.cfg

4.5. Locate the assembly result

The primary assembly output file is: /data/01_assemble/03.ctg_graph/nd.asm.fasta

5. Running NextPolish for Error Correction
5.1. Prepare the Short-Read (NGS) data paths file

$ realpath seq_1.fq >> sgs.fofn && realpath seq_2.fq >> sgs.fofn

5.2. Configure the run.cfg file (Copy a template from /home/software/NextPolish/test_data/run.cfg)

Modify the following lines marked as crucial:

text
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6                           # Number of parallel jobs (5-10)
multithread_jobs = 5
genome = ./raw.genome.fasta                 # Input genome file (from step 4.5)
genome_size = auto
workdir = ./01_rundir                       # Output directory for this run
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./sgs.fofn                       # Path to the FOFN created in step 5.1
sgs_options = -max_depth 100 -bwa
5.3. Execute NextPolish

$ /home/software/NextPolish/nextPolish run.cfg

The polished genome will be located in the 01_rundir directory: genome.nextpolish.fasta

6. Purge Haplotigs for Haplotypic Purification
6.1. Step 1: Generate initial histograms

$ conda activate purgehap # Switch to the appropriate Conda environment

$ /home/software/purge_hap/purge_haplogs1.sh run_path refGenome.fa nanopore_reads.fq threads_num

run_path # Path to a new directory for this analysis

refGenome.fa # Path to the assembled genome (polished or unpolished)

nanopore_reads.fq # Path to the uncompressed Nanopore long reads

threads_num # Number of threads to use (e.g., 30)

6.2. Determine parameters from the histogram

A PNG histogram plot will be generated in the run directory.

Identify the two main peaks and determine the following three values from the plot (depth is on the X-axis):

LowCutOff # Depth value at the left boundary of the left peak

MidPoint # Depth value at the midpoint (valley) between the two peaks

HighCutOff # Depth value at the right boundary of the right peak

6.3. Step 2: Run the main purification and cleanup

$ /home/software/purge_hap/purge_haplogs2.sh run_path refGenome.fa nanopore_reads.fq threads_num LowCutOff MidPoint HighCutOff

The last three parameters are from step 6.2; others remain the same as in 6.1.
$ conda activate base # Switch back to the base Conda environment

7. HiC-Pro Data QC (Optional - Can be skipped if high-quality Hi-C data is available)
Requires an assembled genome and Hi-C sequencing data.
7.1. File Preparation

$ mkdir reference && cp genome.fa ./reference && cd reference # Create a reference directory and copy the genome

$ bowtie2-build -f genome.fa genome_idx # Build a Bowtie2 index for the genome

$ digest_genome.py genome.fa -r dpnii -o genome_DpnII.bed # Generate a BED file of restriction fragment sites (for DpnII)

$ samtools faidx genome.fa # Index the genome

$ awk '{print $1 "\t" $2}' genome.fa.fai > genome.chrom.sizes # Create a chromosome sizes file

Create the input directory structure:

text
InputPath/
├── sample1
│   ├── sample1_R1.fastq.gz
│   └── sample1_R2.fastq.gz
└── sample2
    ├── sample2_R1.fastq.gz
    └── sample2_R2.fastq.gz
Filenames must end with _R1.fastq.gz / _R2.fastq.gz.

7.2. Modify the config-hicpro.txt file

$ cp /home/software/HiC-Pro-master/config-hicpro.txt ./ # Copy the config template

Edit the following parameters:

N_CPU # Number of CPUs/threads to use

BOWTIE2_IDX_PATH # Path to the directory containing the Bowtie2 index (just the directory path)

REFERENCE_GENOME # Basename of the Bowtie2 index (e.g., genome_idx, without path or extension)

GENOME_SIZE # Path to the genome.chrom.sizes file

GENOME_FRAGMENT # Path to the restriction fragment BED file (genome_DpnII.bed)

LIGATION_SITE # Ligation junction sequence (for DpnII: GATCGATC)

JOB_MEM # Memory allocation per job (e.g., 100GB)

7.3. Run HiC-Pro

$ HiC-Pro -c config-hicpro.txt -i InputPath -o Output_Directory_HiCPro

7.4. Interpreting Results

Key result PDFs are found in: Output_Directory_HiCPro/hic_results/pictures/

*Fragment* statistics:

Valid_interaction_pairs should be > 90%. Higher is better.

*ContactRanges* statistics:

Valid Interactions should be > 80%.

Cis long-range contacts should be > 40%.

8. 3D-DNA for Chromosome Scaffolding
8.1. File Preparation (Note: Corrected section number from 5.1 to 8.1)

$ mkdir fastq && cp reads_R1.fastq.gz reads_R2.fastq.gz ./fastq # Create dir and copy Hi-C FASTQs

$ mkdir reference && cp genome.fa ./reference && cd reference # Create reference dir and copy genome

$ bwa index genome.fa # Build a BWA index for the genome

$ python /home/software/juicer/misc/generate_site_positions.py DpnII genome genome.fa # Generate DpnII site info file (genome_DpnII.txt)

$ awk 'BEGIN{OFS="\t"}{print $1, $NF}' genome_DpnII.txt > genome.chrom.sizes # Create chromosome sizes file from site info

8.2. Run Juicer

$ nohup /home/software/juicer/scripts/juicer.sh \

-g genome \ # Genome name (prefix for output files)

-s DpnII \ # Restriction enzyme used

-t 40 \ # Number of threads to use

-D /home/software/juicer \ # Path to the Juicer directory

-z reference/genome.fa \ # Path to the reference genome FASTA

-y reference/genome_DpnII.txt \ # Path to the restriction site file

-p reference/genome.chrom.sizes \ # Path to the chromosome sizes file

&> juicer.log & # Run in background, redirecting output to juicer.log

The file aligned/merged_nodups.txt is the main output needed for the next step.

8.3. Run 3D-DNA

$ nohup /home/software/3d-dna-master/run-asm-pipeline.sh -r 3 reference/genome.fa aligned/merged_nodups.txt &> 3d.log &

This will generate chromosome-scale scaffolds.
————————————————————————————————————————————————————————————————————
# Genome Annotation Workflow using MAKER (with Docker)
This document outlines the steps for de novo genome annotation using the MAKER pipeline, including preparatory steps for repeat masking and transcriptome assembly, all within Docker containers for reproducibility.

Prerequisites before running MAKER
1. Obtaining and Preparing a Protein Database from Related Species
Search databases like NCBI RefSeq or Ensembl using the taxonomic group of your target species.

Download available protein sequences and merge them into a single file.

Remove identical duplicate protein sequences using the script rm_dump_by_fa_seq.py.

Example command (adjust paths): python /path/to/rm_dump_by_fa_seq.py -i merged_proteins.fa -o protein.fasta

The final cleaned protein database should be named protein.fasta.

2. Building and Filtering a De Novo Repeat Database with RepeatModeler
This step is performed in a dedicated Docker container.

2.0 Container Setup and Data Transfer

bash

docker run --name your_name_rep -idt repeat:v1 /bin/bash

docker cp /path/to/your/genome.fa your_name_rep:/data/

docker exec -it your_name_rep /bin/bash
2.1 De Novo Repeat Library Construction

bash
cd /data

BuildDatabase -name your_genome_db -engine ncbi genome.fa

RepeatModeler -pa 4 -database your_genome_db -engine ncbi -LTRStruct

The primary output file is RM_.../consensi.fa.classified.

2.2 (Optional) Filtering the Repeat Library
This step removes repeats that might be real genes by blasting them against the protein database.

bash

makeblastdb -dbtype prot -in protein.fasta -out protein_db

blastx -query RM_*/consensi.fa.classified -db protein_db -outfmt 6 -evalue 1e-5 -num_threads 20 -out repeats_vs_protein.blast

python /path/to/Del_seq_in_blast.py -q RM_*/consensi.fa.classified -blast repeats_vs_protein.blast -l 30 -o filtered_repDB
The final repeat library is either the original consensi.fa.classified or the filtered filtered_repDB_accept.fa.

2.3 Extract Data and Clean Up Container

bash

realpath filtered_repDB_accept.fa  # or consensi.fa.classified

exit

docker cp your_name_rep:/data/filtered_repeats.fa /path/on/your/host/

docker stop your_name_rep
docker rm your_name_rep
MAKER Run Preparation
a. Create MAKER Container

bash
docker run --name your_name_maker -idt maker_image:tag /bin/bash
b. Copy All Necessary Data into the Container

bash
docker cp /host/path/to/genome.fa your_name_maker:/data/
docker cp /host/path/to/protein.fasta your_name_maker:/data/
docker cp /host/path/to/filtered_repeats.fa your_name_maker:/data/ # From step 2.3
docker cp /host/path/to/RNA-Seq/sample1_R1.fq your_name_maker:/data/
docker cp /host/path/to/RNA-Seq/sample1_R2.fq your_name_maker:/data/
c. Enter the MAKER Container

bash
docker exec -it your_name_maker /bin/bash
NOTE: All subsequent operations are performed INSIDE the your_name_maker container.

3. Transcriptome Assembly and Processing
3.1 Genome-Guided Transcriptome Assembly with Trinity
Perform this for each RNA-Seq sample/library, in separate directories if multiple exist.

bash
bwa index genome.fa

bwa mem -t 15 genome.fa sample1_R1.fq sample1_R2.fq | samtools sort -@ 10 -m 4G -O bam -o sample1.sorted.bam -

Trinity --genome_guided_bam sample1.sorted.bam --max_memory 50G --genome_guided_max_intron 10000 --CPU 6
The main assembly output is trinity_out_dir/Trinity-GG.fasta.

3.2 Deduplication of Assembled Transcripts

bash
cat trinity_out_dir1/Trinity-GG.fasta trinity_out_dir2/Trinity-GG.fasta > all_trinity_assembled.fasta

cd-hit-est -i all_trinity_assembled.fasta -o all_trinity_95.fasta -c 0.95 -n 10 -d 0 -M 32000 -T 20
The final EST evidence file for MAKER is all_trinity_95.fasta.

Running MAKER
The MAKER pipeline is run iteratively to improve annotation quality.

4. MAKER - Run 1 (Initial Annotation)
4.1 Configure MAKER

bash
python /home/tools/MAKER_ctl_create.py -run_type 1 -g genome.fa -pep protein.fasta -est all_trinity_95.fasta -rep filtered_repeats.fa -augustus_species zebrafish
-augustus_species: Critical parameter. Use a species from the Augustus library that is phylogenetically close to your target species.

zebrafish (fish), chicken (birds), fly (diptera), human (mammals) are common examples.

Optimization Note (maker_opts.ctl):

The pred_flank value controls how much genomic sequence flanking a gene prediction is considered. Adjust based on genome size to prevent fragmented or merged genes:

100MB - 500MB: pred_flank=80~100

500MB - 1.5GB: pred_flank=100~300

1.5GB: pred_flank>300

This parameter should be changed consistently in maker_opts.ctl for Run 2 and Run 3 if adjusted here.

4.2 Execute MAKER Run 1

bash
nohup python3 /home/tools/MAKER_multrun.py -g genome.fa -o /data/run1 -opt /data/ -p 30 1> /data/run1/maker_run1.log 2>&1 &
The main output GFF file is /data/run1/all.gff. This is used for training in the next step.

5. MAKER - Run 2 (Annotation with Trained Ab Initio Predictors)
5.1 Train Augustus

5.1.1 Extract Training Sequences
Uses the GFF from Run 1 (all.gff) to get genomic sequences around predicted genes.

bash
samtools faidx genome.fa # Ensure index exists
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' /data/run1/all.gff | awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | bedtools getfasta -fi genome.fa -bed - -fo augustus_training_data.fa
5.1.2 Train Augustus using BUSCO

bash
nohup /home/software/busco/bin/busco -i augustus_training_data.fa -o target_species_busco -l /home/software/busco/odb10/actinopterygii_odb10/ -m genome -c 20 --long --augustus_species zebrafish --augustus_parameters='--progress=true' 1> busco.log 2>&1 &

cp -r target_species_busco/run_/augustus_output/retraining_parameters/BUSCO_target_species_busco $AUGUSTUS_CONFIG_PATH/species/
-l: Choose the appropriate BUSCO database lineage (actinopterygii_odb10, diptera_odb10, insecta_odb10, vertebrata_odb10).

5.2 Train SNAP
Uses the GFF from Run 1 to train the SNAP gene predictor.

bash
mkdir snap && cd snap
maker2zff -c 0.8 -e 0.8 -o 0.8 -x 0.2 /data/run1/all.gff
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
mkdir params && cd params
forge ../export.ann ../export.dna
cd ../
hmm-assembler.pl genome params > genome.hmm
The final SNAP HMM model is snap/genome.hmm.

5.3 Configure MAKER for Run 2

bash
python /home/tools/MAKER_ctl_create.py -run_type 2 -g genome.fa -reanno_gff /data/run1/all.gff -augustus_new BUSCO_target_species_busco -snap_hmm /data/snap/genome.hmm
5.4 Execute MAKER Run 2

bash
nohup python3 /home/tools/MAKER_multrun.py -g genome.fa -o /data/run2 -opt /data/ -p 30 1> /data/run2/maker_run2.log 2>&1 &
The output GFF file is /data/run2/all.gff.

6. MAKER - Run 3 (Final Annotation)
6.1 Configure MAKER for Run 3

bash
python /home/tools/MAKER_ctl_create.py -run_type 3 -g genome.fa -reanno_gff /data/run2/all.gff -augustus_new BUSCO_target_species_busco -snap_hmm /data/snap/genome.hmm -pep protein.fasta -est all_trinity_95.fasta -rep filtered_repeats.fa
6.2 Execute MAKER Run 3

bash
nohup python3 /home/tools/MAKER_multrun.py -g genome.fa -o /data/run3 -opt /data/ -p 30 1> /data/run3/maker_run3.log 2>&1 &
The final annotation file is /data/run3/all.gff.

Post-MAKER Processing
6.3 Filtering the Final GFF

bash
python3 /home/tools/gff_filter_fromMAKER.py -i /data/run3/all.gff -o /data/run3/final_annotation_filtered.gff
6.4 Extract Sequences

bash
gffread /data/run3/final_annotation_filtered.gff -g genome.fa -x final_annotation.cds -y final_annotation.pep
6.5 Functional Annotation with InterProScan

bash
/home/software/interproscan/interproscan.sh -i final_annotation.pep -iprlookup -goterms -cpu 20 -t p -f tsv -o final_annotation_interpro.tsv
——————————————————————————————————————————————————————————
# Population Genetic Analysis Pipeline from Variant Calls
1. Purpose and Scope
This document outlines a standardized bioinformatics pipeline for performing population genetic analysis starting from a VCF file containing variant calls from multiple individuals. The pipeline includes quality control filtering, principal component analysis (PCA), and population structure inference using ADMIXTURE. Results are visualized using R and TBtools.

2. Prerequisites & Software

Input Data: A VCF file generated from a variant caller (e.g., GATK).

Software:

VCFtools (v0.1.16+)

PLINK (v1.90+)

ADMIXTURE (v1.3.0+)

R (v4.0+) with packages ggplot2 and readr

TBtools (for visualizing ADMIXTURE results)

3. Detailed Workflow
3.1. Quality Control Filtering with VCFtools
This step removes low-quality variants to reduce noise and potential errors in downstream analyses.

Command:

bash
vcftools --vcf input.vcf \
         --maf 0.01 \
         --hwe 1e-6 \
         --max-missing 0.8 \
         --recode \
         --recode-INFO-all \
         --out filtered_maf001_hwe1e6_miss08
Parameter Explanation:

--vcf input.vcf: Specifies the input VCF file.

--maf 0.01: Retains only variants with a Minor Allele Frequency (MAF) greater than or equal to 1%. This removes very rare variants.

--hwe 1e-6: Filters out variants that significantly deviate from Hardy-Weinberg Equilibrium (HWE) at a p-value threshold of 10⁻⁶. This helps remove genotyping errors.

--max-missing 0.8: Retains only variants that are called in at least 80% of individuals (--max-missing 0.8 means a maximum of 20% missing data is allowed).

--recode: Creates a new VCF file after applying the filters.

--recode-INFO-all: Preserves all INFO fields from the original VCF in the output file.

--out filtered_maf001_hwe1e6_miss08: Defines the prefix for the output files.

Output:

filtered_maf001_hwe1e6_miss08.recode.vcf: The filtered VCF file.

3.2. Data Conversion and PCA with PLINK
PLINK is used to convert the VCF into its binary format (for efficiency) and to perform Principal Component Analysis (PCA).

Step 2.1: Convert VCF to PLINK's binary format (BED/BIM/FAM)

bash
plink --vcf filtered_maf001_hwe1e6_miss08.recode.vcf \
      --make-bed \
      --out filtered_plink
Step 2.2: Perform PCA on the genotype data
This command generates the covariance matrix and calculates the eigenvectors.

bash
plink --bfile filtered_plink \
      --pca 15 header \
      --out filtered_pca
--bfile filtered_plink: Uses the binary PLINK files as input.

--pca 15 header: Calculates the first 15 principal components and includes a header in the output file. The number 15 is often chosen to capture most population structure, with the first 2-3 typically being the most important for visualization.

--out filtered_pca: Defines the output prefix.

Output:

filtered_pca.eigenval: The eigenvalues for each principal component, indicating their relative importance.

filtered_pca.eigenvec: The eigenvectors (the PCA coordinates) for each individual.

3.3. Visualizing PCA Results in R
The provided R script (PCA.R) creates a publication-quality PCA plot using the ggplot2 package.

R Script (PCA.R):

r

library(ggplot2)

 1. Read in the data
pca_data <- read.table("filtered_pca.eigenvec", header = TRUE)

eigenvals <- read.table("filtered_pca.eigenval", header = FALSE)

 2. (CRITICAL) Add a 'group' column to pca_data.

groups <- read.csv("sample_groups.csv", header = TRUE)
pca_data <- merge(pca_data, groups, by = c("FID", "IID")) # Merge group info

 3. Calculate percentage of variance explained by each PC
total_variance <- sum(eigenvals$V1)
variance_prop <- (eigenvals$V1 / total_variance) * 100

 4. Create axis labels with the variance percentage
x_label <- paste0("PC1 (", round(variance_prop[1], 2), "%)")
y_label <- paste0("PC2 (", round(variance_prop[2], 2), "%)")

 5. Generate the PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +                 # Scatter plot points
  theme_bw() +                           # Use a clean black-and-white theme
  theme(panel.grid = element_blank()) +  # Remove grid lines
  labs(x = x_label, y = y_label) +       # Add informative axis labels
  stat_ellipse(level = 0.95, linetype = 2) + # Add 95% confidence ellipses
  scale_color_manual(values = c("#FF9900", "#68cc99", "#EF566B", "#3E91D2")) # Custom colors
Execution:
Run the script in RStudio or from the command line:

bash
Rscript PCA.R
Output: A PDF or on-screen plot showing the first two principal components, with samples colored by their pre-defined group and surrounded by confidence ellipses.

3.4. Population Structure Analysis with ADMIXTURE
ADMIXTURE estimates ancestry proportions assuming K hypothetical ancestral populations.

Step 4.1: Convert PLINK files to a format suitable for ADMIXTURE
ADMIXTURE requires the .bed file but also needs the corresponding .fam and .bim files in the same directory.

Step 4.2: Run ADMIXTURE for different values of K
Run ADMIXTURE in cross-validation mode to estimate the quality of the fit for each K.

bash
for K in 3 4 5; do
    admixture --cv filtered_plink.bed $K | tee log${K}.out
done
--cv: Performs 5-fold cross-validation. The CV error is printed at the end of the log file and helps choose the best K (lower is generally better).

filtered_plink.bed: The input genotype file in PLINK binary format.

$K: The number of assumed ancestral populations.

tee log${K}.out: Saves the output to a log file for later inspection.

Output for each K:

filtered_plink.${K}.Q: The estimated ancestry proportions for each individual. This is the key file for visualization.

filtered_plink.${K}.P: The estimated allele frequencies in the ancestral populations.

log${K}.out: The log file containing the cross-validation error.

Step 4.3: Visualize ADMIXTURE results in TBtools

Prepare the .Q file (e.g., filtered_plink.3.Q, filtered_plink.4.Q) and a file listing the sample names in the same order as the .Q file (this can be extracted from the .fam file: cut -f1,2 -d' ' filtered_plink.fam > sample_list.txt).

Open TBtools.

Navigate to the visualization utilities.

Use the appropriate function (often called "ADMIXTURE Plot" or "Stacked Bar Chart Plot") to load the .Q file and the sample list.

Customize the plot (colors, sample order, labels) and export the figure.

4. Interpretation of Results

PCA: Clustering of samples along PC1 and PC2 indicates genetic similarity. Distinct clusters often represent different populations. Outliers can indicate admixed individuals or sample contamination.

ADMIXTURE: Each vertical bar represents an individual. The proportion of colors in each bar represents the estimated proportion of ancestry from K ancestral populations. The value of K with the lowest cross-validation error is often considered the most likely, but biological interpretation is also crucial.

Integration: Results from PCA and ADMIXTURE should be consistent and provide complementary views of the population structure within your samples.
———————————————————————————————————————————————————————————————————————————
# Targeted Gene Enrichment Sequencing Data Analysis Pipeline
Pipeline Overview:
(Raw sequencing data -> Decompression -> Adapter trimming) ->
Assembly -> {Merging different batches -> Removing low-quality/unneeded sequences} ->
Sequence Alignment -> Filtering -> {Removing reference/contaminant sequences -> Re-aligning} ->
Summary Statistics -> Locus Concatenation -> Gene Tree Construction (IQ-TREE/RAxML) -> Download and visualize tree
                                                                  -> Constructing multiple gene trees -> Coalescent-based Species Tree (ASTRAL) ->

Note: Steps in () are typically handled by the sequencing provider. Steps in {} are optional.

1. Assembly of Enriched Sequencing Data
Original Tool Repository: https://github.com/yhadevol/Assexon

1.1 Adapter Trimming
View usage instructions by running the following command in any directory:

bash
trim_adaptor.pl
Required Directory Structure (directrim):

1) A folder containing demultiplexed sequencing files (each file represents all sequences from one sample).

2) Text file: inlineindex.txt (contains fixed inline index sequences).

3) Text file: indexpair.txt (contains IS1 and IS2 barcode pairs for this dataset, self-created).

Format: {sample<TAB>IS1<TAB>IS2}

4) The Perl script: trim_adaptor.pl.

Execution:

bash
cd /path/to/directrim

perl trim_adaptor.pl --raw_reads demultiplexed --inline_index inlineindex.txt --index_pair indexpair.txt --trimmed trimmed

nohup perl trim_adaptor.pl --raw_reads demultiplexed --inline_index inlineindex.txt --index_pair indexpair.txt --trimmed trimmed > trim.log &

ps -ef | grep trim_adaptor.pl

kill <process_id>
Output (directrim directory):

1) trimmed/: Folder containing adapter-trimmed sequences (same format as input).

2) trimming_report/: Folder containing processing reports.

3) trimmed_reads_bases_count.txt: Text file summarizing total number of reads and total bases after trimming. View with less or a text editor.

4) trim.log: Log file.

1.2 Assembly (Using Modular Scripts)
View usage instructions:

bash
assemble.pl
Required Directory Structure (directassem):

1) The trimmed/ folder from step 1.1.

2) The modular assembly script assemble.sh.

Execution:

bash
cd /path/to/directassem

bash assemble.sh

ps -ef | grep assemble.sh
Output (directassem directory):

1) assemble_results/: Folder containing assembled sequences.

Subfolders: f (filtered?), nf (likely non-filtered/primary output, used downstream), p (potential?).

2) enriched_gene.txt: Text file listing the number of loci enriched per sample and the percentage of total reference loci captured. View with less.

3) assemb.log: Log file.

4) Other less critical files and folders.

1.3 Post-Assembly Processing
Merging results from multiple batches & deduplication:

bash
mkdir sample_merge
cd sample_merge


3. Edit 'merge.sh' (e.g., using a text editor like vim or nano). Modify the '--indir' paths to point to the 'nf' folders from different assembly batches. Example content for merge.sh:
    ------------------------------------------------------
!/bin/bash
    merge_loci.pl \
        --indir "/path/to/batch1/assemble_results/nf /path/to/batch2/assemble_results/nf" \
        --outdir merged_nf \
        --min_seq 3
    echo "Done!!!"
    ------------------------------------------------------

 4. Execute the merge script
bash merge.sh # This creates 'merged_nf/'

 5. Remove duplicate sequences resulting from the merge
pick_taxa.pl --indir merged_nf --outdir deduplicated --rm_dup_taxa # Creates 'deduplicated/'
Sequence Alignment:

Required Directory (direcalign):

deduplicated/ (or merged_nf/) folder from merging.

Aligning.sh module.

Edit Aligning.sh: Set --dna_unaligned to your input dir (e.g., deduplicated) and --dna_aligned to your desired output dir name (e.g., samples_nf_aligned).

Execution:

bash
cd /path/to/direcalign
nohup bash Aligning.sh > alignment.log &
# Check status: ps -ef | grep Aligning.sh or tail -f alignment.log
Output: Aligned sequences directory (e.g., samples_nf_aligned).

Filtering Poor Alignments:

Required Directory (direcfilter):

Aligned sequences directory (e.g., samples_nf_aligned).

filtering.sh module.

Execution:

bash
cd /path/to/direcfilter
nohup bash filtering.sh > filter.log &
 Check status: ps -ef | grep filtering.sh
Output: Filtered alignments directory (e.g., samples_nf_filtered).

Optional: Removing Reference/Unneeded Sequences:

Edit the Delete_unneeded_sequences.sh module to specify input, output, and sequence names/patterns to remove.

Realign the resulting dataset (repeat the Alignment step above).

Final output: Filtered and re-aligned sequences directory (e.g., samples_nf_aligned_final).

Summary Statistics for Enrichment Efficiency:

Required Directory (direcstatistic):

Filtered alignments directory (e.g., samples_nf_filtered or samples_nf_aligned_final).

statistics.sh module.

Execution:

bash
cd /path/to/direcstatistic
nohup bash statistics.sh > statistics.log &
 Check status: ps -ef | grep statistics.sh
Output:

sample_summary.txt: Number of loci enriched and average GC content per sample.

loci_summary.txt: Enrichment information per locus.

2. Phylogenetic Tree Construction
This section concatenates all aligned loci for each species into a supermatrix and builds a phylogenetic tree.

2.1 Locus Concatenation
Required Directory:

Directory containing the final aligned sequences (e.g., samples_nf_aligned_final).

Script: concat_loci.sh.

Edit concat_loci.sh: Set the input directory (--indir) and the desired output filename prefix (--outfile).

Execution:

bash
nohup bash concat_loci.sh > concat.log &
 Check concat.log for progress
Output: Three concatenated sequence files in different formats: .fas (FASTA), .nex (NEXUS), .phy (PHYLIP). The .phy file is often used for tree building.

2.2 Partitioning by Codon Position
Required Directory (on server):

The concatenated NEXUS file from step 2.1 (e.g., concat_seq.nex).

Perl script: extract_DNAblocks_iqtree.pl.

Edit extract_DNAblocks_iqtree.pl: Point my $file = "filename.nex"; to your NEXUS file.

Execution (on server):

bash
nohup perl extract_DNAblocks_iqtree.pl > dnabl.log &
Output: DNA_blocks.txt (defines codon positions).

Local Computer Processing:

Create a folder partitionfolder on your local machine.

Place the Python script transform_to_raxmlformat.py inside.

Create a new, empty text file DNA_partition.txt inside.

Download DNA_blocks.txt from the server to partitionfolder.

Run the Python script locally:

bash
 Navigate to the folder in your local terminal/command prompt
cd path/to/partitionfolder
python transform_to_raxmlformat.py
 Follow the prompts: Use 'DNA_blocks.txt' as input and 'DNA_partition.txt' as output.
Upload the newly created DNA_partition.txt file back to the server analysis directory.

2.3 Tree Inference (IQ-TREE)
Required Directory:

The concatenated PHYLIP file (concat_seq.phy).

The partition file (DNA_partition.txt).

Script: iqtree.sh.

Edit iqtree.sh: Ensure the parameters are correct:

-s concat_seq.phy (Input sequences)

-spp DNA_partition.txt (Partition file)

-prefix <desired_output_prefix> (e.g., MyGeneTree)

Execution:

bash
nohup bash iqtree.sh > iqtree.log &
Output: The main tree file is {prefix}.treefile (e.g., MyGeneTree.treefile). Other files include log, model information, and support values.
—————————————————————————————————————————————————————————————
# Mitochondrial Genome Assembly using GetOrganelle and Subsequent Annotation
1. Purpose and Scope
This document outlines a standard operating procedure for the de novo assembly of animal mitochondrial (mt) genomes from whole-genome sequencing (WGS) paired-end reads using the get_organelle toolkit, followed by annotation and extraction of specific genes.

2. Prerequisites & Software

Operating System: Unix/Linux (e.g., Ubuntu, CentOS) or macOS.

Software:

GetOrganelle (v1.7.5+ recommended). Install via pip install getorganelle or from GitHub.

Bowtie2 (bundled with GetOrganelle but must be in $PATH).

SPAdes (bundled with GetOrganelle).

Input Data: Illumina paired-end sequencing reads in FASTQ format (_R1.fq, _R2.fq).

Reference Sequence: A closely related mitochondrial genome sequence (in FASTA format) downloaded from NCBI GenBank.

3. Workflow Overview
The process involves three main stages:

Assembly: Using get_organelle_from_reads.py to assemble the mitochondrial genome for each sample.

Annotation: Uploading the assembled genome to the MITOS2 web server for structural annotation.

Gene Extraction: Using MEGA software to visually inspect the annotation and extract the nucleotide sequence of the target gene (e.g., COX1).

4. Detailed Methodology
4.1. Obtaining a Reference Sequence
Navigate to the NCBI Nucleotide database (https://www.ncbi.nlm.nih.gov/nucleotide).

Search for the mitochondrial genome of a species phylogenetically closest to your sample(s).

Download the complete sequence in FASTA format.

Note the file path for use in the assembly script (e.g., /path/to/reference_mt.fasta).

4.2. Automated Assembly Script
The provided Bash script automates the assembly process for multiple samples.

4.2.1. Script Explanation (assembly_script.sh)

bash
!/bin/bash

 Define directory variables
input_dir="/path/to/trimmed/reads/" # Directory containing input FASTQ files
output_dir="/path/to/assembly/output/" # Directory for assembly results
reference_path="/path/to/reference_mt.fasta" # Path to the NCBI reference FASTA

 Create the output directory if it doesn't exist
mkdir -p "$output_dir"

 Loop through all _R1.fq files in the input directory
for r1_file in "$input_dir"*_R1.fq; do
    # Extract the sample base name (assumes naming convention: 'sample_R1.fq')
    sample_name=$(basename "$r1_file" _R1.fq)
    
    # Construct the path to the corresponding reverse read file
    r2_file="$input_dir${sample_name}_R2.fq"
    
    # Check if the reverse read file exists before proceeding
    if [[ -f "$r2_file" ]]; then
        # Define a unique output directory for each sample
        output_sample_dir="$output_dir${sample_name}_out"
        
        # Execute the GetOrganelle command
        get_organelle_from_reads.py \
          -1 "$r1_file" \          # Forward reads
          -2 "$r2_file" \          # Reverse reads
          -R 10 \                  # Initial seed extension rounds
          -s "$reference_path" \   # Reference sequence for seed generation
          -k 21,45,65,85,105 \     # K-mer sizes for assembly graph
          -F animal_mt \           # Target genome type (animal mitochondrial)
          -o "$output_sample_dir" \ # Output directory
          -t 5 &                   # Number of threads; '&' runs jobs in background
    else
        echo "Error: Reverse reads $r2_file not found for $sample_name"
    fi
done

 Wait for all background processes to finish
wait
echo "All assembly jobs are complete."
4.2.2. Key Parameters:

-1, -2: Paths to forward and reverse read files.

-s: Path to the reference sequence. This is used to "bait" and recruit mitochondrial reads.

-F animal_mt: Specifies the target is an animal mitochondrial genome, optimizing the assembly strategy.

-k: A list of k-mer sizes. Using multiple k-mers improves the ability to resolve repeats and assemble complete circles.

-t: Number of CPU threads to use per sample.

&: The ampersand runs each sample's assembly as a background job, processing all samples in parallel. The final wait command ensures the script doesn't exit until all jobs finish.

4.2.3. Execution:

Save the script as run_assembly.sh.

Modify the input_dir, output_dir, and reference_path variables to match your system's paths.

Make the script executable: chmod +x run_assembly.sh.

Run the script: ./run_assembly.sh.

4.2.4. Output:
For each sample, a directory (${sample_name}_out) is created. The final assembly is typically found in the *.fasta file within this directory (e.g., sample1_out/extended_K105.complete.graph1.1.fasta). The *.fastg file contains the assembly graph for visualization in Bandage.

4.3. Annotation using MITOS2 Web Server
Access the MITOS2 Web Server: http://mitos2.bioinf.uni-leipzig.de/index.py.

Submit Assembly: Upload your assembled FASTA file.

Select Parameters:

Genetic Code: Select the appropriate code (e.g., 5: Invertebrate Mitochondrial or 2: Vertebrate Mitochondrial).

Reference: RefSeq 89 Metazoa (or latest available).

Keep other parameters at their default settings.

Run Annotation: Submit the job. You will receive results via email.

4.4. Gene Extraction using MEGA
Open MEGA software on your local computer.

Open Annotated Sequence: Load the FASTA file of your assembled and MITOS2-annotated genome.

View Annotation: Navigate to the Display menu to view the annotated features (genes, tRNAs, etc.).

Extract COX1:

Locate the COX1 (Cytochrome c oxidase subunit 1) gene in the feature table or graphical view.

Select the sequence corresponding to the COX1 gene.

Use MEGA's export or copy function to extract the selected nucleotide sequence into a new FASTA file. This file can now be used for downstream phylogenetic analysis or BLAST searches.

5. Notes and Troubleshooting

Reference Choice: The quality and phylogenetic proximity of the reference sequence are critical for successful seed-based assembly.

Read Depth: Low coverage of mitochondrial reads may result in incomplete assemblies.

Contamination: If the assembly is poor, the sample may have contamination from other species. Consider using --reduce-reads for large datasets or pre-filtering reads with bowtie2 against the reference.

Complexity: For genomes with long repeats, the assembly graph (*.fastg) should be visualized using Bandage to manually check for circularity and potential alternative assemblies.
—————————————————————————————————————————————————————————
# MSMC2 and MSMC-IM Analysis Pipeline Documentation

1. Data Preparation: Generating Mask and VCF Files from BAM Files

This section describes how to generate per-chromosome mask files and VCF files from high-depth sequencing BAM files for four samples.

1.1. Calculate Average Depth per Chromosome
The average depth is required for generating the sample-specific mask. While the example uses chr1, this should be done for each chromosome (chr1, chr2, ..., chrN) in a loop for each sample's BAM file.

bash
Example for chromosome 1 of a sample
samtools depth -r chr1 <sample.bam> | awk '{sum += $3} END {print sum / NR}'
Replace <sample.bam> with the path to your BAM file.The output is the mean coverage (depth) for chr1.
1.2. Generate Sample Mask and VCF Files
Using the mean coverage value calculated in the previous step, generate a mask file (in BED format) and a VCF file for each chromosome of each sample. The bamCaller.py script (part of the msmctools utilities) is used for this purpose.

bash
Example command for a single chromosome (e.g., chr1) of a single sample
bcftools mpileup -C 50 -u -r chr1 -f <reference_genome.fa> <sample.bam> --threads 16 | \
bcftools call -c -V indels --threads 16 | \
./bamCaller.py <mean_coverage> <sample_chr1_mask.bed.gz> | \
bgzip -c > <sample_chr1.vcf.gz>

Replace the placeholders:
<reference_genome.fa>: Path to your reference genome FASTA file.
<sample.bam>: Path to the input BAM file for the sample.
<mean_coverage>: The average depth calculated for this chromosome and sample.
<sample_chr1_mask.bed.gz>: Output path for the sample's mask file for chr1 (gzipped).
<sample_chr1.vcf.gz>: Output path for the sample's VCF file for chr1 (gzipped).
Note: This process must be repeated for each chromosome and for each of the four samples, resulting in files like sample1_chr1.vcf.gz, sample1_chr1_mask.bed.gz, sample2_chr1.vcf.gz, etc.

2. Generating the Genome Mappability Mask

A general mappability mask for the reference genome is also required to exclude low-complexity or unmappable regions.

2.1. Run run_snpable2.sh
This script (often part of the MSMC tools ecosystem) processes the reference genome to identify mappable regions.

bash
Assuming run_snpable2.sh is in your path or current directory
./run_snpable2.sh <reference_genome.fa>
This will generate output files used in the next step.
2.2. Process with makeMappabilityMask.py
This Python script converts the output from run_snpable2.sh into the mask format required by MSMC.

bash
Example command
python makeMappabilityMask.py <output_from_snpable> > <genome_mask.bed>
Replace <output_from_snpable> with the relevant file generated by run_snpable2.sh (often a .mask file).
Replace <genome_mask.bed> with the desired output file name for the genome-wide mask.This genome mask must also be split by chromosome for use in step 3.
3. Phasing VCFs with WhatsHap

Phase the per-chromosome VCF files for each sample using WhatsHap to improve haplotype resolution for MSMC2.

bash
Example command for a single chromosome of a single sample
whatshap phase --output <phased_sample_chr1.vcf.gz> --reference <reference_genome.fa> <sample_chr1.vcf.gz> <sample.bam>
Replace <phased_sample_chr1.vcf.gz> with the desired output name for the phased VCF. The original <sample_chr1.vcf.gz> and <sample.bam> are used as inputs.
Note: Repeat this for all chromosomes and all four samples. The phased VCFs (phased_*.vcf.gz) will be used as input for the next step.

4. Preparing MSMC2 Input Files

MSMC2 requires a specific input format (multihetsep). The generate_multihetsep.py script combines the sample masks, the genome mappability mask, and the phased VCFs to create this input for each chromosome.

bash
Example command for chromosome 1
generate_multihetsep.py --chr chr1 \
    --mask sample1_chr1_mask.bed.gz \
    --mask sample2_chr1_mask.bed.gz \
    --mask sample3_chr1_mask.bed.gz \
    --mask sample4_chr1_mask.bed.gz \
    --mask genome_chr1_mask.bed \ # The chromosome-specific genome mappability mask
    phased_sample1_chr1.vcf.gz \
    phased_sample2_chr1.vcf.gz \
    phased_sample3_chr1.vcf.gz \
    phased_sample4_chr1.vcf.gz \
    > chr1.multihetsep.txt
Note: This command must be run for each chromosome independently, generating a *.multihetsep.txt file for each (chr1.multihetsep.txt, chr2.multihetsep.txt, ...).

5. Running MSMC2

MSMC2 is run twice: once for analysis within a species/population ("within") and once for analysis between species/populations ("cross").

5.1. Within-Species Analysis
This command runs MSMC2 on haplotypes from the same species (e.g., two individuals from Species 1). The -I flag specifies which haplotypes to use (0,1 = first two haplotypes; 2,3 = next two, etc.).

bash
Example for Species 1 (using haplotypes 0 and 1 from the input)
msmc2 -I 0,1 -o species1_msmc2 chr1.multihetsep.txt chr2.multihetsep.txt chr3.multihetsep.txt ...
The input files are all the chromosome-specific multihetsep files.

Example for Species 2 (using haplotypes 2 and 3)
msmc2 -I 2,3 -o species2_msmc2 chr1.multihetsep.txt chr2.multihetsep.txt chr3.multihetsep.txt ...
5.2. Cross-Species Analysis
This command runs MSMC2 on all possible pairs of haplotypes between the two species. The -I flag defines the cross-haplotype pairs (0-2, 0-3, 1-2, 1-3).

bash
msmc2 -I 0-2,0-3,1-2,1-3 -o cross_species_msmc2 chr1.multihetsep.txt chr2.multihetsep.txt chr3.multihetsep.txt ...
6. Combining Results for MSMC-IM

The combineCrossCoal.py script combines the output from the within and cross analyses into a single file suitable for MSMC-IM.

bash
python combineCrossCoal.py cross_species_msmc2.final.txt species1_msmc2.final.txt species2_msmc2.final.txt > combined_crosspecies.msmc2.final.txt
Replace the .final.txt filenames with the actual outputs from your MSMC2 runs.
7. Running MSMC-IM

Finally, run MSMC-IM using the combined file to infer isolation and migration history.

bash
python MSMC_IM.py \
    -beta 1e-8,1e-6 \ # Defines the grid search range for the migration rate (M)
    -o result \        # Output prefix
    -mu 5.97e-9 \      # Assumed mutation rate per base per generation
    --printfittingdetails \
    --plotfittingdetails \
    --xlog \           # Plot time axis on a logarithmic scale
    combined_crosspecies.msmc2.final.txt
—————————————————————————————————————————————————————————————————————————————————————
# FastSimCoal (fsc) Automation Tool Documentation
Overview
This documentation describes two complementary scripts that automate running FastSimCoal (fsc) multiple times and selecting the best run based on maximum likelihood values.

Script 1: run_fsc_4D.py
Purpose
A Python script that automates running fsc 100 times with the same demographic model but different random seeds.

Usage
bash
python run_fsc_4D.py --fsc <fsc_directory> --dir <working_directory> --pre <file_prefix>
Parameters
--fsc or -F: Path to the directory containing the fsc executable (required)

--dir or -D: Working directory where input files are located (required)

--pre or -P: Prefix of the fsc input files (required)

Input Files Required
The script expects the following files in the working directory:

<prefix>.tpl - Template file

<prefix>.est - Estimation file

Multiple joint MAF files: <prefix>_jointMAFpopX_Y.obs

Functionality
Creates 100 subdirectories (run1, run2, ..., run100)

Copies all necessary input files to each subdirectory

Runs fsc in each directory with specific parameters:

-t - template file

-e - estimation file

-m -0 - use 0 migration matrix

-n 100000 - number of simulations

-L 30 - number of loops

-s 0 - random seed

-M - perform parameter estimation

-c 15 - number of cores to use

-B 15 - number of bootstrap replicates

After completing all runs, executes the selection script

Script 2: selectbestrun.sh
Purpose
A Bash script that identifies the best fsc run from multiple executions by comparing likelihood values.

Usage
bash
./selectbestrun.sh <file_prefix>
Parameters
<file_prefix>: Prefix of the fsc input files (without .tpl extension)

Functionality
Searches through all run* directories

Looks for .bestlhoods files in each run directory

Compares the likelihood values from all runs

Identifies the run with the highest likelihood value

Copies all files from the best run to a bestrun directory

Reports the number of runs analyzed and the best run found

Complete Workflow
Step 1: Prepare Input Files
Ensure all required fsc input files are in your working directory:

.tpl file

.est file

.obs files (multiple joint MAF files)

Step 2: Make Scripts Executable
bash
chmod +x selectbestrun.sh
Step 3: Run the Automation Script
bash
python run_fsc_4D.py --fsc /path/to/fsc/ --dir /path/to/working/directory/ --pre your_file_prefix
Step 4: Review Results
After completion:

The best run files will be copied to bestrun/ directory

Summary information will be displayed in the console

Detailed logs are available in bestrun.log

Notes
The scripts assume a specific naming convention for joint MAF files

The fsc version used is fsc28 (can be modified in the Python script if needed)

The number of cores (-c) and bootstrap replicates (-B) are set to 15

The number of runs is fixed at 100 but can be modified in the Python script

Ensure you have appropriate permissions to execute scripts and create directories
———————————————————————————————————————————————————————————————=
# Population Genetic Analysis Pipeline: Fst, Dxy, and Pi using pixy and CMplot
1. Overview
This pipeline calculates key population genetics statistics—Fst (genetic differentiation), Dxy (absolute genetic divergence), and Pi (nucleotide diversity)—from a VCF file using the pixy software. The raw output from pixy is then processed by custom Python scripts to clean, split, and format the data, preparing it for final visualization using the CMplot package in R to generate Manhattan-style plots.

2. Core Analysis: Calculating Statistics with pixy
The first step is to run pixy to compute Fst, Dxy, and Pi statistics in windows across the genome.

Command:

bash
pixy --stats pi fst dxy \
     --vcf your_data.vcf.gz \        # Input VCF file (must be bgzipped and indexed)
     --populations popfile.txt \     # Tab-delimited file defining populations: `sample_name<TAB>population_name`
     --window_size 50000 \           # Size of the sliding window (50 kb in this case)
     --n_cores 8 \                   # Number of CPU cores to use (recommended)
     --output_prefix pixy_results    # Base name for output files
Output Files:

pixy_results_pi.txt: Nucleotide diversity (π) per population per window.

pixy_results_fst.txt: Genetic differentiation (Fst) per pair of populations per window.

pixy_results_dxy.txt: Absolute genetic divergence (Dxy) per pair of populations per window.

Input File (popfile.txt) Format:

text
sample1    popA
sample2    popA
sample3    popB
sample4    popB
sample5    popC
...
3. Data Processing with Custom Python Scripts
The raw pixy output files require processing to handle missing values, split data by population pairs, and reformat them for visualization.

3.1. Script 1: Fst_Dxy_1st.py (for Fst and Dxy data)
Purpose: This script serves as the primary cleaner and splitter for the pixy_results_fst.txt and pixy_results_dxy.txt files.

Key Functions:

Data Cleaning: Reads the input file and targets the 6th column (which contains the Fst or Dxy value). It replaces NA, NaN, na, nan values and any negative values with 0.

Data Splitting: Splits the combined data into separate files for each pre-defined population pair (e.g., popA_popB_fst.txt). It ensures the population names in the first two columns are consistently ordered (alphabetically) for each pair.

Sorting: Sorts the data in each output file by chromosome and then by genomic position.

Usage:

bash
python Fst_Dxy_1st.py
Then enter the path to your file, e.g., `pixy_results_fst.txt`
Output: A directory (split_fst_files/) containing individual files for each population pair (e.g., popA_popB_fst.txt, popA_popC_fst.txt).

3.2. Script 2: Fst_Dxy_2nd.py (for Fst and Dxy data)
Purpose: This script acts as the formatter for the files generated by Script 1. It converts them into the standardized format required by CMplot.

Key Functions:

Format Conversion: Extracts the chromosome (column 3), position (column 4), and statistic value (column 6, the cleaned Fst or Dxy) from the split files.

Adds SNP Identifier: Creates a simple sequential SNP ID number for each data point, which is required by CMplot.

Outputs a clean, four-column table with the columns: SNP, Chromosome, Position, trait1 (which contains the Fst or Dxy value).

Usage:

bash
python Fst_Dxy_2nd.py popA_popB_fst.txt popA_popB_fst_CMplot.txt
Output: A formatted file (e.g., popA_popB_fst_CMplot.txt) ready for visualization.

3.3. Script 3: pi_1st.py (for Pi data)
Purpose: This script is the primary processor for the pixy_results_pi.txt file.

Key Functions:

Data Cleaning & Splitting: Reads the input Pi file. It replaces NA values in the 5th column (the Pi value) with 0 and filters the data to only include chromosomes 1 through 24.

Population Separation: Splits the single Pi file into four separate files, one for each population (popA, popB, popC, popD).

Sorting: Sorts the data in each population-specific file by chromosome and position.

Usage:

bash
python pi_1st.py pixy_results_pi.txt
Output: Four files: popA_output.txt, popB_output.txt, popC_output.txt, popD_output.txt.

3.4. Script 4: pi_2nd.py (for Pi data)
Purpose: This script is the formatter for the population-specific Pi files generated by Script 3.

Key Functions:

Format Conversion: Extracts the chromosome (column 2), position (column 3), and Pi value (column 5) from the split Pi files.

Adds SNP Identifier: Creates a sequential SNP ID number.

Outputs a clean, four-column table with the columns: SNP, Chromosome, Position, trait1 (which contains the Pi value).

Usage:

bash
python pi_2nd.py popA_output.txt popA_pi_CMplot.txt
Output: A formatted file (e.g., popA_pi_CMplot.txt) ready for visualization.

4. Visualization with CMplot in R
The final step is to create Manhattan plots to visualize the distribution of the statistics across the genome.

R Code Example (for Pi):

r
Load the CMplot library
library(CMplot)

Read the formatted data
pi_data <- read.table("popA_pi_CMplot.txt", header = TRUE)

Generate a Manhattan plot
CMplot(pi_data,
       type = "p",           # Type of plot: "p" for points
       plot.type = "m",      # Plot type: "m" for Manhattan
       threshold = NULL,     # No significance threshold line
       file = "jpg",         # Output format: "jpg", "pdf", etc.
       file.name = "popA_Pi", # Base name for output file
       dpi = 300,            # High resolution output
       file.output = TRUE,   # Save to file
       verbose = TRUE,       # Print progress
       width = 14,           # Plot width (inches)
       height = 3,           # Plot height (inches)
       chr.labels.angle = 45,# Angle of chromosome labels
       LOG10 = FALSE,        # Do NOT plot -log10(values)
       ylab = "pi value",    # Y-axis label
       ylim = c(0, 0.6),     # Y-axis limits
       cex = 0.2             # Size of the points
)
Explanation of Key CMplot Arguments:

plot.type="m": Specifies a Manhattan plot.

LOG10=FALSE: Plots the raw values (e.g., 0.05) instead of their negative log-transformed value. This is crucial for population genetics statistics like Fst, Dxy, and Pi.

ylab and ylim: Set the Y-axis label and limits, respectively. Adjust ylim based on the range of your data.

cex: Controls the size of the points on the plot.

To create plots for Fst or Dxy, simply replace the input data and adjust the ylab and ylim parameters accordingly (e.g., ylab="Fst", ylim=c(0, 1)).

5. Summary of the Full Pipeline
Run pixy to calculate genome-wide Fst, Dxy, and Pi.

Process Fst/Dxy:

Use Fst_Dxy_1st.py to clean and split the combined Fst/Dxy file by population pair.

Use Fst_Dxy_2nd.py on each split file to format it for CMplot.

Process Pi:

Use pi_1st.py to clean, filter, and split the Pi file by population.

Use pi_2nd.py on each population file to format it for CMplot.

Visualize in R using the CMplot function for each formatted file to generate publication-ready Manhattan plots.
