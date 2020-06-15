while read i
do

echo "Processing sample $i ..." 

nice -n 20 java -jar /ohta/haoran.xue/programs/picard.jar FixMateInformation \
 VALIDATION_STRINGENCY=LENIENT \
       I=/ohta/haoran.xue/PVE_Eichhornia_project/align_PGAbf_pilon/$i.PGAbfr_pilon.align.bam

       
nice -n 20 java -jar /ohta/haoran.xue/programs/picard.jar CleanSam \
 VALIDATION_STRINGENCY=LENIENT \
INPUT=/ohta/haoran.xue/PVE_Eichhornia_project/align_PGAbf_pilon/$i.PGAbfr_pilon.align.bam \
OUTPUT=/ohta/haoran.xue/PVE_Eichhornia_project/align_PGAbf_pilon/$i.tmp_clean.PGAbfr_pilon.bam

nice -n 20 java -jar /ohta/haoran.xue/programs/picard.jar SortSam \
 INPUT=/ohta/haoran.xue/PVE_Eichhornia_project/align_PGAbf_pilon/$i.tmp_clean.PGAbfr_pilon.bam \
 OUTPUT=/ohta/haoran.xue/PVE_Eichhornia_project/align_PGAbf_pilon/$i.tmp_clean_sorted.PGAbfr_pilon.bam \
 SORT_ORDER=coordinate CREATE_INDEX=true \
 VALIDATION_STRINGENCY=LENIENT \
 TMP_DIR=`pwd`/tmp 

 
nice -n 20 java -jar /ohta/haoran.xue/programs/picard.jar AddOrReplaceReadGroups \
      VALIDATION_STRINGENCY=LENIENT \
      I=/ohta/haoran.xue/PVE_Eichhornia_project/align_PGAbf_pilon/$i.tmp_clean_sorted.PGAbfr_pilon.bam  \
      O=/ohta/haoran.xue/PVE_Eichhornia_project/align_PGAbf_pilon/$i.PGAbfr_pilon.sorted.bam  \
      RGID=$i \
      RGLB=lib \
      RGPL=illumina \
      RGPU=$i \
      RGSM=sample$i
	  
	    
done < /ohta/haoran.xue/PVE_Eichhornia_project/samplelist

