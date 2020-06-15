while read i
do

echo "Processing sample $i ..."

ngm -t 20 -r /ohta/haoran.xue/epaniculata_pacbio/PGA_pilon/PGAbf_pilon.fasta -q /ohta/haoran.xue/PVE_Eichhornia_project/eichh/$i\.trimmed.fastq -o /ohta/haoran.xue/PVE_Eichhornia_project/align_PGAbf_pilon/$i\.PGAbfr_pilon.align.bam

done < /ohta/haoran.xue/PVE_Eichhornia_project/samplelist
