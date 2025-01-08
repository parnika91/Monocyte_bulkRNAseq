#! /bin/sh -x

#### take fq name from call_samples_to_STAR.sh

fastq=$1
echo $fastq

ulimit -n 5000

## Local alignment with STAR
	STAR --runThreadN 20 \
	--runMode alignReads \
	--genomeDir Data/human_STARindex \
	--readFilesIn Data/fastq/$fastq*\_R1_001.fastq.gz Data/fastq/$fastq*\_R2_001.fastq.gz \
	--readFilesCommand zcat \
	--sjdbGTFfile human.gtf \
	--outFileNamePrefix Results/bam/$fastq \
	--outSAMtype BAM SortedByCoordinate
	#--readMapNumber 100000 \
	#--limitGenomeGenerateRAM 210000000000 \
	#--limitBAMsortRAM 210000000000

#Input reads can be name-sorted or location-sorted. 
#Users do not need to resort the reads before feeding them to featureCounts.
#from featureCounts manual, recent update

#### end

	
