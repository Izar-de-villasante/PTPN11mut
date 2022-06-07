#!/bin/bash
####### description:
# Template for distributed task, where we want to perform a single run
# of our program with multiple CPUs and shared memory.

# Job description:
#SBATCH --job-name="mut_gatk_ptpn11"
#SBATCH --comment="Prueba de SBATCH"

# ******** Main parameters: *********
#SBATCH --mem=64gb          # Total amount of memory 
# ---- distributed ----> (doParallel,foreach,mclapply) "usually single task multiple workers"
#SBATCH --cpus-per-task=32 # number of cores per task <= 128 
#SBATCH --ntasks=1          # number of tasks, (1cpu x task is default)unles MPI comunication  
#SBATCH --nodes=3           # number of nodes ( <= ntasks since no communication between tasks)
#SBATCH --time=1-23:10:00   # walltime dd-hh:mm:ss. current max  = 07-00:00:00 (7 days)
#SBATCH --array=1-99
#SBATCH -o slurm.%N.%J.%u.%a.out # STDOUT
#SBATCH -e slurm.%N.%J.%u.%a.err # STDERR
#----- Info: 

echo "Starting at $(date)"
echo "Job submitted to the ${SLURM_JOB_PARTITION} partition, the default partition on ${SLURM_CLUSTER_NAME}"
echo "Job name: ${SLURM_JOB_NAME}, Job ID: ${SLURM_JOB_ID}"
echo "  I have ${SLURM_CPUS_ON_NODE} CPUs on compute node $(hostname)"
pwd
# ------ Environment Configuration ------
cd $SLURM_SUBMIT_DIR
n=$(echo "scale=2 ; ${SLURM_ARRAY_TASK_ID}/100" | bc) 
data_raw="/mnt/beegfs/idevillasante/Projects/PTPN11/90-612199113/00_fastq/"
fastq_dir="/mnt/beegfs/idevillasante/Projects/PTPN11/mut/fastq/"
bam_dir="/mnt/beegfs/idevillasante/Projects/PTPN11/mut/bam/"
name="PTPN11_sub$n"
input="${bam_dir}${name}.bam"
#cores=124
# ------ Load modules ------
module load spack
module load Java
module load bwa
#module load R
# ------ Run program ------

#bwa mem -t 64 -R "@RG\tID:2\tSM:pool1" refgene/PTPN11.fasta 90-612199113/00_fastq/out_R1.fastq.gz 90-612199113/00_fastq/out_R2.fastq.gz > trimmed_PTPN11.bam
#gatk --java-options "-Xmx62G" AnalyzeSaturationMutagenesis -I $1  -R refgene/gene.fasta --orf 1-1783 -O "./mut3/$1"

#bwa mem -M -t 64 -R "@RG\tID:ptpn11amp\tSM:ptpn11amp\tPL:ILLUMINA" ./refgene/gene.fasta 90-612199113/00_fastq/PTPN11_sub_R1.fastq.gz 90-612199113/00_fastq/PTPN11_sub_R2.fastq.gz > sub30M_trimmed.bam
echo $n
seqtk sample -s12345 ${data_raw}PTPN11-SeSaM-Library_R1_001.fastq.gz $n > ${fastq_dir}${name}_R1.fastq
seqtk sample -s12345 ${data_raw}PTPN11-SeSaM-Library_R2_001.fastq.gz $n > ${fastq_dir}${name}_R2.fastq
bwa mem -M -t 32 -R "@RG\tID:ptpn11amp\tSM:ptpn11amp\tPL:ILLUMINA" ./refgene/gene.fasta ${fastq_dir}${name}_R1.fastq ${fastq_dir}${name}_R2.fastq > "${bam_dir}${name}.bam"
gatk --java-options "-Xmx64G" AnalyzeSaturationMutagenesis -I $input  -R /mnt/beegfs/idevillasante/Projects/PTPN11/refgene/gene.fasta --orf 1-1782 -O "/mnt/beegfs/idevillasante/Projects/PTPN11/mut/out/${name}"
Rscript -e "rmarkdown::render('./R/codons2.Rmd',params=list(prefix = '../data-raw/mut4/PTPN11_trimmed',gene_file='../data-raw/gene.fasta'))"
# ------ Job stats:
sstat  -j   $SLURM_JOB_ID.batch   --format=JobID,MaxVMSize
seff $SLURM_JOB_ID.batch
