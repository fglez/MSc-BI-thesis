#!/bin/bash

prefix=$1 #An assigned name to identify your sample
path=$(pwd)
sign='$'

cat > 1clean-${prefix}.sh << EOF
#PBS -N clean${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=16g,vmem=16g,walltime=9:00:00
#PBS -V
#PBS -e 1clean${prefix}.err
#PBS -o 1clean${prefix}.out

cd $path

module load fastp/0.20 FastQC/0.11.2

mkdir RAWS

fastp -i FC-${prefix}-Cult_1.fastq.gz -I FC-${prefix}-Cult_2.fastq.gz -o Clean-${prefix}-Cult_1.fastq.gz -O Clean-${prefix}-Cult_2.fastq.gz -q 20 -l 150

fastqc Clean-${prefix}-Cult_1.fastq.gz Clean-${prefix}-Cult_2.fastq.gz

mv FC-${prefix}-Cult_*.fastq.gz RAWS/

EOF

cat > 2kraken-${prefix}.sh << EOF
#PBS -N kraken${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=12,mem=48g,vmem=48g,walltime=99:99:99
#PBS -V
#PBS -e 2kraken${prefix}.err
#PBS -o 2kraken${prefix}.out

module load kraken/2.0.7

cd $path

kraken2 --db /LUSTRE/usuario/mcontreras/kraken --threads 12 --report ${prefix}.kreport --paired Clean-${prefix}-Cult_1.fastq.gz Clean-${prefix}-Cult_2.fastq.gz > ${prefix}.kraken

mkdir kraken/

mv ${prefix}.kreport ${prefix}.kraken kraken/

EOF

cat > 3metaidba-${prefix}.sh << EOF
#PBS -N idba${prefix}
#PBS -q ensam
#PBS -l nodes=1:ppn=16,mem=128g,vmem=256g,walltime=99:00:00
#PBS -e 3idba${prefix}.err
#PBS -o 3idba${prefix}.out
#PBS -V

cd $path 

module load idba/2.0

gunzip Clean-${prefix}-Cult_*.fastq.gz

fq2fa --merge Clean-${prefix}-Cult_1.fastq Clean-${prefix}-Cult_2.fastq CLEAN_${prefix}_merged.fa

idba_ud -r CLEAN_${prefix}_merged.fa -o idba-ud

EOF

cat > 4bowsam-${prefix}.sh << EOF
#PBS -N bowtie${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=24g,vmem=24g,walltime=99:00:00
#PBS -e 4bowsam${prefix}.err
#PBS -o 4bowsam${prefix}.out
#PBS -V

cd $path

module load bowtie2/2.3.5.1 samtools/1.9

mkdir bowtie/

bowtie2-build idba-ud/contig.fa bowtie/${prefix}.ind

bowtie2 -p 8 -x bowtie/${prefix}.ind -1 Clean-${prefix}-Cult_1.fastq -2 Clean-${prefix}-Cult_2.fastq -S bowtie/${prefix}.map.sam

samtools sort -o bowtie/${prefix}.map.sorted.bam -O bam bowtie/${prefix}.map.sam

EOF

cat > 5metabat-${prefix}.sh << EOF
#PBS -N MetaBat${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=32g,vmem=32g,walltime=99:00:00
#PBS -e 5metabat${prefix}.err
#PBS -o 5metabat${prefix}.out
#PBS -V

cd $path

module load metabat/2.13

mkdir metabat/

jgi_summarize_bam_contig_depths --outputDepth metabat/${prefix}depth.txt bowtie/${prefix}.map.sorted.bam

metabat2 --saveCls -i idba-ud/contig.fa -a metabat/${prefix}depth.txt -o metabat/${prefix}bin -m 2500 --maxEdges 600 --noAdd --seed 4

EOF

cat > 6checkm-${prefix}.sh << EOF
#PBS -q default
#PBS -l nodes=1:ppn=8,walltime=9:59:59,vmem=32gb,mem=32gb
#PBS -N check${prefix}
#PBS -V
#PBS -o 6check${prefix}.out
#PBS -e 6check${prefix}.err

cd $path

module load CheckM/1.1.3 Prodigal/2.6.2 hmmer/3.1b2

mkdir checkm/

checkm taxonomy_wf domain Bacteria -t 8 -x fa metabat/ checkm/

EOF

cat > 7mtphln-${prefix}.sh << EOF
#PBS -q default
#PBS -l nodes=1:ppn=16,walltime=9:59:59,vmem=32gb,mem=32gb
#PBS -N metaphlan${prefix}
#PBS -V
#PBS -o 7mtphln${prefix}.out
#PBS -e 7mtphln${prefix}.err

cd $path

module load metaphlan2/2.7.7 bowtie2/2.3.5.1
 
mkdir metaphlan/

metaphlan2.py  CLEAN_${prefix}_merged.fa --input_type fasta --nproc 16 --ignore_viruses --ignore_eukaryotes > metaphlan/${prefix}_mtphln_prfl.txt

EOF

