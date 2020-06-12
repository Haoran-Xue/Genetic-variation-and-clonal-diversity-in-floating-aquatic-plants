#Julia Kreiner, June 2018
#STEPS FOR FILTERING WHOLE GENOME SEQUENCE DATA, SNPCALLED BY FREEBAYES (modified from http://ddocent.com/filtering/ for more details)

######################################################################
#only keep variant sites that have been successfully genotyped in 50% of individuals, a minimum quality score of 30, and a minor allele frequency of less than 0.45 (to remove possible tetraploid sites)
#only keep genotypes with 3 or more reads
######################################################################

vcftools --vcf PVE_crassipes_allsites.vcf --recode --recode-INFO-all --max-missing 0.5 --minQ 30 --max-maf 0.45 --minDP 3 --out PVE_crassipes_PGAbfp_HC_50miss

#If --max-missing 0.25:
#After filtering, kept 215327 out of a possible 816445568 Sites

#If --max-missing 0.5:
#After filtering, kept 200279 out of a possible 816445568 Sites

#If --max-missing 0.75:
#After filtering, kept 188372 out of a possible 816445568 Sites

######################################################################
#estimate missing data for each individual
######################################################################

vcftools --vcf PVE_crassipes_PGAbfp_HC_50miss.recode.vcf --missing-indv

######################################################################
#use cutoff from missing data distribution to drop individuals with more than 50% missing data
######################################################################

/ohta/julia.kreiner/software/bin/bin/mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv

#After filtering, kept 150 out of 154 Individuals

vcftools --vcf PVE_crassipes_PGAbfp_HC_50miss.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out PVE_crassipes_PGAbfp_HC_50miss_indmiss

######################################################################
#mask low complexity regions
#maybe we should not do it
#if we do it, maybe repeatmasker is better than dustmasker because repeatmasker masks 64% of the genome while dustmasker only mask 5% (too low, most of the TEs are not included).
#also, bedtools does not work for this vcf because some sites have format issue
######################################################################

#commands for dustmasker only:
#/ohta/haoran.xue/programs/ncbi-blast-2.10.0+/bin/dustmasker -in /ohta/haoran.xue/epaniculata_pacbio/PGA_pilon/PGAbf_pilon.fasta -out /ohta/haoran.xue/epaniculata_pacbio/PGA_pilon/PGAbf_pilon.fasta.dm.list -outfmt acclist
#awk 'split($1,a,">") {print a[2]"\t"$2-1"\t"$3}' /ohta/haoran.xue/epaniculata_pacbio/PGA_pilon/PGAbf_pilon.fasta.dm.list > PGAbfp.dustmasked.bed
#also need to remove -1 values
#sed 's/-1/0/g' PGAbfp.dustmasked.bed

#commands for both maskers:
#bedtools subtract -a PVE_crassipes_PGAbfp_HC_50miss_indmiss.recode.vcf -b /ohta/haoran.xue/epaniculata_pacbio/PGA_pilon_rm/softMasked/PGAbf_pilon.rm.bed > PVE_crassipes_PGAbfp_HC_50miss_indmiss_rm.vcf
#grep "#" PVE_crassipes_PGAbfp_HC_50miss_indmiss.recode.vcf > header
#cat header PVE_crassipes_PGAbfp_HC_50miss_indmiss_rm.vcf > PVE_crassipes_PGAbfp_HC_50miss_indmiss_rm_rehdrs.vcf

######################################################################
#filter by QUAL/DEPTH (error prone and high copy number loci)
#i#####################################################################

vcffilter -f "QUAL / DP > 0.25" PVE_crassipes_PGAbfp_HC_50miss_indmiss.recode.vcf > PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp.vcf

cut -f8 PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp.vcf | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp.DEPTH

grep -v "#" PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp.vcf | cut -f1,2,6 > PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp.vcf.loci.qual

# Run this command in Terminal #
/ohta/julia.kreiner/software/bin/bin/mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp.DEPTH
#This command calculates the mean depth and print it on the screen

# Run this command in Terminal #
python -c "print int(4640.22+(3*(4640.22**0.5)))" = 4844

paste PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp.vcf.loci.qual PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp.DEPTH | /ohta/julia.kreiner/software/bin/bin/mawk -v x=4844 '$4 > x' | /ohta/julia.kreiner/software/bin/bin/mawk '$3 < 2 * $4' > PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp.lowQDloci

vcftools --vcf PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp.vcf --recode-INFO-all --exclude-positions PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp.lowQDloci --out PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth --recode

#After filtering, kept 20275 out of a possible 21985 Site

######################################################################
#filter by max mean DEPTH (not sure whether we should do it)
######################################################################

# Run this command in Terminal #
#cut -f3 PVE_crassipes_PGAbfp_HC_50miss_indmiss_filtqualbydepth.ldepth > PVE_crassipes_PGAbfp_HC_50miss_indmiss_filtqualbydepth.justdepth

#/ohta/julia.kreiner/software/bin/bin/mawk '!/D/' PVE_crassipes_PGAbfp_HC_50miss_indmiss_filtqualbydepth.justdepth  | /ohta/julia.kreiner/software/bin/bin/mawk -v x=121 '{print $1/x}' > PVE_crassipes_PGAbfp_HC_50miss_indmiss_filtqualbydepth.meandepthpersite

#vcftools --vcf PVE_crassipes_PGAbfp_HC_50miss_indmiss_filtqualbydepth.recode.vcf --recode-INFO-all --out PVE_crassipes_PGAbfp_HC_50miss_indmiss_filtqualbydepth_maxmeanDP --max-meanDP --recode

######################################################################
#filter on allelic balance for heterozygous calls
#filter sites that have over 100 times more forward alternate reads than reverse alternate reads and 100 times more forward reference reads than reverse reference reads along with the reciprocal
#filter on whether or not their is a discrepancy in the properly paired status of for reads supporting reference or alternate alleles
#since freebayes, drop sites where QUAL < 30 (1/1000 error rate)
#not going to use them as vcf files called with GATK HaplotypeCaller do not have these tags
######################################################################
#vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth.recode.vcf > PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_AB.vcf
#vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth.recode.vcf > PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_AB_RB.vcf
#vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth.recode.vcf > PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_PS.vcf
#vcffilter -f "QUAL > 30" PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth.vcf > PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_QUAL.vcf

######################################################################
#vcfallelicprimitives
#If multiple alleleic primitives (gaps or mismatches) are specified in a single VCF record, split the record into multiple lines, but drop all INFO fields
######################################################################

vcfallelicprimitives PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth.recode.vcf --keep-info --keep-geno > PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_prim.vcf

######################################################################
#HW equilibrium
#not doing this since samples are not in a panmictic population
######################################################################
./filter_hwe_by_pop.pl -v PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_AB_PS_QUAL_noindels.recode.vcf -p popmap -o PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_AB_PS_QUAL_HWE

######################################################################
#recode
######################################################################

vcftools --vcf PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_prim.vcf --recode --recode-INFO-all --out PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_prim

######################################################################
#remove indels
######################################################################

vcftools --vcf PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_prim.vcf --remove-indels --recode --recode-INFO-all --out PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_prim_noindels

#After filtering, kept 17796 out of a possible 20275 Sites

######################################################################
#only keep variant sites
######################################################################

vcftools --vcf PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_prim.recode.vcf --recode --recode-INFO-all --mac 1 --out PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_filtqualbydepth_prim_variantsitesonly

#After filtering, kept 6083 out of a possible 20275 Sites