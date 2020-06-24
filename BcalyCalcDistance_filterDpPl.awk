#Read vcf file, find average distance between samples; only compare sites that both individuals possess a Gt that PASSES pl filter
BEGIN {print "-----Begin BcalyCalcDistance_filterPl.awk-----";
	print "OutputFile: " outputFileName;
	print "thresholdDp: " thresholdDp
	print "thresholdPl: " thresholdPl;
	tableFileName1 = outputFileName ".distance1Table.Dp" thresholdDp "Pl" thresholdPl ".tbl";
	tableFileName2 = outputFileName ".distance2Table.Dp" thresholdDp "Pl" thresholdPl ".tbl";
	outFileName = outputFileName ".distance.Dp" thresholdDp "Pl" thresholdPl ".txt";
	
	numSite=0; numSample = 0;
	countLine=0;switch1=0;}
{	
	#Print to log file only when reading first record
	if(FNR==1){
		print "thresholdPl\tSample1\tSample2\tnumSiteCompared\tdistance1\tdistance2\tnumRefHet\tnumAltHet\tnumHetHet\tnumRefAlt\tnumRefRef\tnumAltAlt" > outFileName ;
	}
	
	#to display while processing file
	countLine++;
	if(countLine%1000==0){
		print "Read " countLine " lines";
	}
	
	
	#if header line
	if($0 ~ /\#/){

		#Last header row, contains column names of site info, used to find field that contains sample being removed
		if($0~/\#CHROM/){
			#look from field 10 and onward for sample names
			for(i = 10; i <= NF; i++){
				numSample++;
				sampleName[i-9]=$i;
			}
		}
	}
	#else all other lines (i.e. sites)
	else{
	
		#Initialize stuff after going through all header lines
		if(switch1==0){
			#Initialize site compared count
			for(i = 1; i <=numSample; i++){
				for(j = i; j <= numSample; j++){
					numSiteCompared[i,j]=0;
					distance1[i,j]=0;
					distance2[i,j]=0;
					numRefHet[i,j]=0;
					numAltHet[i,j]=0;
					numHetHet[i,j]=0;
					numRefAlt[i,j]=0;
					numRefRef[i,j]=0;
					numAltAlt[i,j]=0;
					
					cRefRef[i,j]=0;
					cHetHet[i,j]=0;
					cAltAlt[i,j]=0;
					cRefHet[i,j]=0;
					cAltHet[i,j]=0;
					cRefAlt[i,j]=0;
				}
			}
			#change switch1 to 1 so will not run initialization again
			switch1=1;
		}
	
		
		#increase count for number of sites
		numSite++;
		
		#####Find position PL and DP in FORMAT (field 9)
		#split FORMAT field by ":"
		nFormat=split($9,format,":");
		positionDp=0;
		positionPl=0;
		#find element in format that is "DP"
		for(f=1;f<=nFormat;f++){
			if(format[f]~/DP/){positionDp=f}
			if(format[f]~/PL/){positionPl=f}
		}
		#if no Dp field, then skip to next line; else continue forward to determine if site passes filter
		if(positionDp==0){nSiteNoDp++;print$0;next;}else{nSiteYesDp++;}
		#####Find position PL and DP in FORMAT (field 9)
		
		#loop through all pair-wise comparisons for individuals in a given site
		#Start at field 10 since field 1-9 are not genotype calls of samples (it contains summary stats of site)
		#i is focal individual
		for(i = 10; i <= NF; i++){
			
			#j is comparison to individual i
			for(j = i; j <=NF; j++){
			
				#indicate of ind i and j passes filter
				pass_i = 0;
				pass_j = 0;
				#indicate genotype on ind, if called; 1 = homo ref, 2 = getero, 3 = homo alt
				gt1 = 0;
				gt2 = 0;
				
				#examine if genotype is called for ind i
				if($i~/0\/0/ || $i~/0\/1/ || $i~/1\/0/ || $i~/1\/1/){
					
					##########Determine if ind i pass pl filter
					#####Split field indicated by samplePosition with ":";
					nFields=split($i,info,":");
					
					#####record dp of Gt
					gtDp=info[positionDp];
					
					#####Find second largest phred likelihood of genotype
					#split pl field by ","; record number of fields in nPl; record elements into pl
					#pl is in element positionPl of info
					nPl=split(info[positionPl],pl,",");
	 
					#find the two likelihoods larger than 0 (format always has 0 and two non-zeros)
					num1=0;num2=0;
					for(k=1;k<=length(pl);k++){
					 if(pl[k] > 0 && num1 == 0){
					  num1=pl[k];
					 }
					 else if(pl[k] > 0 && num1 > 0){
					  num2 = pl[k];
					 }
					}

					#determine second largest pl (the smaller of num1,num2)
					midPl=0;
					if(num1>num2){midPl=num2;}else{midPl=num1;}
					#####Find second largest phred likelihood of genotype
					
					#####Pass threshold: indicate and record which Gt
					if(midPl>=thresholdPl && gtDp>=thresholdDp){
						
						pass_i = 1;
						
						if($i~/0\/0/){
							gt1 = 1;
						}
						else if($i~/1\/1/){
							gt1 = 3;
						}
						else{
							gt1 = 2;
						}
					 }
					
				}
				
				#examine if genotype is called for ind j
				if($j~/0\/0/ || $j~/0\/1/ || $j~/1\/0/ || $j~/1\/1/){
					
					##########Determine if ind j pass Pl filter
					#####Split field indicated by samplePosition with ":";
					nFields=split($j,info,":");
					
					#####record dp of Gt
					gtDp=info[positionDp];
					
					#####Find second largest phred likelihood of genotype
					#split pl field by ","; record number of fields in nPl; record elements into pl
					#pl is in element positionPl of info
					nPl=split(info[positionPl],pl,",");
	 
					#find the two likelihoods larger than 0 (format always has 0 and two non-zeros)
					num1=0;num2=0;
					for(k=1;k<=length(pl);k++){
					 if(pl[k] > 0 && num1 == 0){
					  num1=pl[k];
					 }
					 else if(pl[k] > 0 && num1 > 0){
					  num2 = pl[k];
					 }
					}

					#determine second largest pl (the smaller of num1,num2)
					midPl=0;
					if(num1>num2){midPl=num2;}else{midPl=num1;}
					#####Find second largest phred likelihood of genotype
					
					#####Pass threshold: indicate and record which Gt
					if(midPl>=thresholdPl && gtDp>=thresholdDp){
						
						pass_j = 1;
						
						if($j~/0\/0/){
							gt2 = 1;
						}
						else if($j~/1\/1/){
							gt2 = 3;
						}
						else{
							gt2 = 2;
						}
					 }
				}
				
				#Site in BOTH ind are called if pass_i+pass_j==2 (both sites pass midPl threshold)
				if(pass_i+pass_j==2){
				
					numSiteCompared[i-9,j-9]++;
					
					#####Calc genetic distance
					#if genotypes the same 
					if(gt1==gt2){
						distance1[i-9,j-9] += 0;
						
						#distance2 different from distance1 when both genotype heterozygous
						if(gt1==2 && gt2==2){
							distance2[i-9,j-9] += 0.5
						}
						else{
							distance2[i-9,j-9] += 0
						}
					}
					#else if homo ref (gt=1) or homo alt (gt=3) vs hetero (gt=2); 1*2==2 and 3*2==6
					else if(gt1*gt2==2 || gt1*gt2==6){
						distance1[i-9,j-9] += 0.5;
						distance2[i-9,j-9] += 0.5;
					}
					#else if homo ref (gt=1) vs homo alt (gt=3); 1*2==2
					else if(gt1*gt2==3){
						distance1[i-9,j-9] += 1;
						distance2[i-9,j-9] += 1;
					}
					#####Calc genetic distance
					
					#print i-9 "," j-9 ": " $i " , " $j ": " distance1[i-9,j-9] "," distance2[i-9,j-9];
					
					#####Count genotype comparisons
					if(gt1*gt2==1){
						numRefRef[i-9,j-9]++;
					}
					else if(gt1*gt2==4){
						numHetHet[i-9,j-9]++;
					}
					else if(gt1*gt2==9){
						numAltAlt[i-9,j-9]++;
					}
					else if(gt1*gt2==2){
						numRefHet[i-9,j-9]++;
					}
					else if(gt1*gt2==6){
						numAltHet[i-9,j-9]++;
					}
					else if(gt1*gt2==3){
						numRefAlt[i-9,j-9]++;
					}
					#####Count genotype comparisons
				}
			}
		}	
	}
}
END {
	print "ENDING";
	print numSample;
	#Calc average distance
	for(i = 1; i <=numSample; i++){
		for(j = i; j <= numSample; j++){
			#calc average distance if more than 0 sites compared
			if(numSiteCompared[i,j] > 0){
				distance1[i,j] /= numSiteCompared[i,j];
				distance2[i,j] /= numSiteCompared[i,j];
			}
			#if no sites compared, give -1
			else{
				distance1[i,j] = -1;
				distance2[i,j] = -1
			}			
		}
	}
	
	#Print list of distances
	for(i = 1; i <=numSample; i++){
		for(j = i; j <= numSample; j++){
			print thresholdPl "\t" sampleName[i] "\t" sampleName[j] "\t" numSiteCompared[i,j] "\t" distance1[i,j] "\t" distance2[i,j] "\t" numRefHet[i,j] "\t" numAltHet[i,j] "\t" numHetHet[i,j] "\t" numRefAlt[i,j] "\t" numRefRef[i,j] "\t" numAltAlt[i,j] > outFileName;
			
			distance1Table[i,j]=distance1[i,j]
			distance1Table[j,i]=distance1[i,j]
			distance2Table[i,j]=distance2[i,j]
			distance2Table[j,i]=distance2[i,j]
		}
	}
	
	#Print table of distances
	printf "Distance1\t" > tableFileName1
	printf "Distance2\t" > tableFileName2
	for(i = 1; i <=numSample; i++){
		printf sampleName[i] > tableFileName1
		printf sampleName[i] > tableFileName2
		if(i<numSample){
			printf "\t" > tableFileName1
			printf "\t" > tableFileName2
		}
		else{
			printf "\n" > tableFileName1
			printf "\n" > tableFileName2
		}
	}
		
	for(i = 1; i <=numSample; i++){
		
		printf sampleName[i] "\t" > tableFileName1
		printf sampleName[i] "\t" > tableFileName2
		
		for(j = 1; j <= numSample; j++){
			
			printf distance1Table[i,j] > tableFileName1
			printf distance2Table[i,j] > tableFileName2
			
			if(j<numSample){
			printf "\t" > tableFileName1
			printf "\t" > tableFileName2
			}
			else{
				printf "\n" > tableFileName1
				printf "\n" > tableFileName2
			}
		}
	}
}