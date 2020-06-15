while read i
do

echo "Processing sample $i ..." 

java -jar /ohta/haoran.xue/programs/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 10 -phred33 -trimlog trimlog $i.fastq.gz $i.trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done < /ohta/haoran.xue/PVE_Eichhornia_project/samplelist

