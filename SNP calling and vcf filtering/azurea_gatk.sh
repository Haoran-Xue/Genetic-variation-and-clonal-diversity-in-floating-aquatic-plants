#!/bin/bash

#install GATK: https://software.broadinstitute.org/gatk/guide/quickstart
#install Picard: https://broadinstitute.github.io/picard/ , https://github.com/broadinstitute/picard/releases/tag/2.8.1

project=PVE_azurea
prefix=/ohta/haoran.xue/PVE_Eichhornia_project/align_PGAbf_pilon/azurea
picard=/ohta/haoran.xue/programs/picard.jar
gatk=/ohta/haoran.xue/programs/gatk-4.1.2.0/gatk
names=/ohta/haoran.xue/PVE_Eichhornia_project/align_PGAbf_pilon/azurea/azurea_samples
ref=/ohta/haoran.xue/epaniculata_pacbio/PGA_pilon/PGAbf_pilon.fasta
tmp=/ohta/haoran.xue/tmp

mkdir -p ${tmp}

#TK: set hard paths as variables so you don't need to touch the script below this section


#TK: loop over names in a file...more reproducible
while read sample
do

###################
#VARIANT DISCOVERY
###################

#Run haplotype caller with GVCF function
$gatk HaplotypeCaller -R ${ref} -I ${prefix}/${sample}.PGAbfr_pilon.sorted.bam --dont-use-soft-clipped-bases true --stand-call-conf 20 -o ${prefix}/${sample}.g.vcf -ERC GVCF >${prefix}/${sample}_hapl.out 2>${prefix}/${sample}_hapl.err


done < $names

# Nicolay: I will work better on this loop to combine the vcf files.
#up until now running commands over a loop for each file, the next step joins all files into one.
#CombineGVCFs



all_variants=""
for variant in $(ls ${prefix}/*.g.vcf);do
    variant=$(basename $variant)
    all_variants="$all_variants --variant $variant"
done


#GenotypeGVCF
java -Djava.io.tmpdir=${tmp} -jar $gatk -T GenotypeGVCFs -R ${ref} --variant ${project}.g.vcf --includeNonVariantSites -o ${project}_allsites.vcf
