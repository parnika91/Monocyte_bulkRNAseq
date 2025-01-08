#! /bin/sh/ -x

### script to run the bulk RNAseq pipeline
file="Data/samplelist.csv"
samples="$(cat $file | cut -f1)"

for sample in $samples; 
do 
  #parallel --eta -j 2 --link bash --verbose echo ::: $(cat $line);
  sh Scripts/bulkRNA_noQC.sh $sample
  #echo $sample
  done ##< $samples ## check

#Scripts/bulkRNA_noQC.sh