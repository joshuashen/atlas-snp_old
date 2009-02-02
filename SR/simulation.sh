# input :  a sample fastq file, reference fasta


sample=$1
ref=$2
dist=$3
var=$4
num=$5

maq simutrain simu_train.dat $sample
maq simulate -d $dist -s $var -N $num -h simu_1.fastq simu_2.fastq $ref simu_train.dat

mkdir batches

split -d -l 2000000 simu_1.fastq batches/F_reads_
split -d -l 2000000 simu_2.fastq batches/R_reads_

