#! /bin/sh

#### indexing the human genome using ENCODE pri fasta and gtf 
# gencode.v39.primary_assembly.annotation -> gtf
# GRCh38.primary_assembly.genome -> fasta
# ENCODE version 5 (14 Feb 2022)

# cd /home/parnika/Documents/Projects/Monocyte_bulRNAseq/

zcat Data/GRCh38.primary_assembly.genome.fa.gz > human.fa
zcat Data/gencode.v39.primary_assembly.annotation.gtf.gz > human.gtf

STAR \
-- runThreadN 12 \
-- runMode genomeGenerate \
-- genomeDir Data/human_STARindex/ \
-- genomeFastaFiles human.fa \
-- sjdbGTFfile human.gtf \
-- sjdbOverhang 99

#rm human.*
#### end

STAR \
-- runThreadN 12 \
-- runMode genomeGenerate \
-- genomeDir human_STARindex/ \
-- genomeFastaFiles human.fa \
-- sjdbGTFfile human.gtf \
-- sjdbOverhang 99
