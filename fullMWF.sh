#!/bin/bash

prefix=$1 #An assigned name to identify your sample
path=$(pwd)
sign='$'

cat > 0clean-${prefix}.sh << EOF

#PBS -N fastp${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=4,mem=16g,vmem=16g
#PBS -V
#PBS -o fastp${prefix}.out
#PBS -e fastp${prefix}.err

module load fastp/0.20 FastQC/0.11.2

cd $path

mkdir fastp/

fastp -i RAWS/AC1ME2S${prefix}_${prefix}_L002_R1_001.fastq.gz -I RAWS/AC1ME2S${prefix}_${prefix}_L002_R2_001.fastq.gz -o CLEAN_${prefix}_L002_R1_001.fastq.gz -O CLEAN_${prefix}_L002_R2_001.fastq.gz -q 20 -l 98 -y --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

fastqc CLEAN_${prefix}_L002_R1_001.fastq.gz CLEAN_${prefix}_L002_R2_001.fastq.gz

mv CLEAN_${prefix}_L002_R1_001.fastq.gz CLEAN_${prefix}_L002_R2_001.fastq.gz fastp/

EOF

cat > 1rmvhst-${prefix}.sh << EOF

#PBS -N bowtiehost${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=24g,vmem=24g,walltime=99:00:00
#PBS -e bthost${prefix}.err
#PBS -o bthost${prefix}.out
#PBS -V

cd $path

module load bowtie2/2.3.5.1

mkdir bthost nohost host

bowtie2-build /LUSTRE/usuario/mcontreras/thesisseq/eatala20/genomahost/GCA_017140195.1_ASM1714019v1_genomic.fna bthost/${prefix}eata.ind

bowtie2 -p 8 -x bthost/${prefix}eata.ind -1 fastp/CLEAN_${prefix}_L002_R1_001.fastq.gz -2 fastp/CLEAN_${prefix}_L002_R2_001.fastq.gz --un-conc nohost/nohost${prefix}.fastq --al-conc host/host${prefix}.fastq

EOF

cat > 2kraken-${prefix}.sh << EOF

#PBS -N kraken${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=12,mem=48g,vmem=48g,walltime=99:99:99
#PBS -e kraken${prefix}.err
#PBS -o kraken${prefix}.out
#PBS -V

cd $path

module load kraken/2.0.7

kraken2 --db /LUSTRE/usuario/mcontreras/kraken --threads 12 --report ${prefix}.kreport --paired nohost/nohost${prefix}.1.fastq nohost/nohost${prefix}.2.fastq > ${prefix}.kraken

mkdir kraken/

mv ${prefix}.kreport ${prefix}.kraken kraken/

EOF

cat > 3metaidba-${prefix}.sh << EOF

#PBS -N idba${prefix}
#PBS -q ensam
#PBS -l nodes=1:ppn=16,mem=128g,vmem=256g,walltime=99:00:00
#PBS -e idba${prefix}.err
#PBS -o idba${prefix}.out
#PBS -V

cd $path 

module load idba/2.0

fq2fa --merge nohost/nohost${prefix}.1.fastq nohost/nohost${prefix}.2.fastq nohost_${prefix}_merged.fa

idba_ud -r nohost_${prefix}_merged.fa -o idba-ud

EOF

cat > 4bowsam-${prefix}.sh << EOF
#PBS -N bowsam${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=24g,vmem=24g,walltime=99:00:00
#PBS -e bowsam${prefix}.err
#PBS -o bowsam${prefix}.out
#PBS -V

cd $path

module load bowtie2/2.3.5.1 samtools/1.9

mkdir bowsam/

bowtie2-build idba-ud/contig.fa bowsam/${prefix}.ind

bowtie2 -p 8 -x bowsam/${prefix}.ind -1 nohost/nohost${prefix}.1.fastq -2 nohost/nohost${prefix}.2.fastq -S bowsam/${prefix}.map.sam

samtools sort -o bowsam/${prefix}.map.sorted.bam -O bam bowsam/${prefix}.map.sam

EOF

cat > 5metabat-${prefix}.sh << EOF
#PBS -N MetaBat${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=32g,vmem=32g,walltime=99:00:00
#PBS -e 2metabat${prefix}.err
#PBS -o 2metabat${prefix}.out
#PBS -V

cd $path

module load metabat/2.13

mkdir metabat/

jgi_summarize_bam_contig_depths --outputDepth metabat/${prefix}depth.txt bowsam/${prefix}.map.sorted.bam

metabat2 --saveCls -i idba-ud/contig.fa -a metabat/${prefix}depth.txt -o metabat/${prefix}bin -m 2500 --maxEdges 600 --noAdd --seed 4

EOF

cat > 6checkm-${prefix}.sh << EOF
#PBS -q default
#PBS -l nodes=1:ppn=8,walltime=9:59:59,vmem=32gb,mem=32gb
#PBS -N check${prefix}
#PBS -o check${prefix}.out
#PBS -e check${prefix}.err
#PBS -V

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

metaphlan2.py  nohost_${prefix}_merged.fa --input_type fasta --nproc 16 --ignore_viruses --ignore_eukaryotes > metaphlan/${prefix}_prfl.txt

EOF