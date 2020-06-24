BEGIN{print "-----Begin BcalySiteDepth.awk-----"
	print "Output file Name: " outputFileName
	logFileName = outputFileName ".siteDepth.txt"
	logFileName2 = outputFileName ".siteDepthPerInd.txt"
	nLines=0; nSite=0}
{	#count num of lines read
	nLines++
	if(nLines%100000==0){print "Lines read: " nLines;}
	#header info
	if($0~/\#/){
		#Last header row, print output file header (need sample names)
		if($0~/\#CHROM/){
			
			#print column names
			print "CHROM\tPOS\tnCalledGt\ttotalDp\tmeanDp">logFileName
			printf "CHROM\tPOS\t" > logFileName2
		
			#look from field 10 and onward for sample names
			for(i = 10; i <= NF; i++){
				numSample++;
				sampleName[i-9]=$i;
				
				#print column names
				printf "%s", $i > logFileName2
				if(i < NF){printf "\t" > logFileName2}
				if(i == NF){printf "\n" > logFileName2}
			}
		}
	}
	#line is a site if not a header
	else{
		#count num of sites
		nSite++
		
		#initialize variables
		totalDp=0
		nCalledGt=0

		#####Find position PL and DP in FORMAT (field 9)
		#split FORMAT field by ":"
		nFormat=split($9,format,":")
		positionDp=0;
		positionPl=0;
		positionAd=0;
		#find element in format that is "DP"
		for(f=1;f<=nFormat;f++){
			if(format[f]~/DP/){positionDp=f}
			if(format[f]~/PL/){positionPl=f}
			if(format[f]~/AD/){positionAd=f}
		}
		#if no Dp field, then skip to next line; else continue forward to determine if site passes filter
		if(positionDp==0){nSiteNoDp++;print$0;next;}else{nSiteYesDp++;}
		#####Find position PL and DP in FORMAT (field 9)

		#####Record depth about each site
		for(i = 10; i <= NF; i++){
			#examine if genotype is called for ind i
			if($i~/0\/0/ || $i~/0\/1/ || $i~/1\/0/ || $i~/1\/1/){
				
				#count num of gt regardless of quality
				nCalledGt++
			
				#Split field by ":"; record number of fields in nFields; record elements into info
				nFields=split($i,info,":")
				
				#Record total Dp of site
				totalDp += info[positionDp]
				
				#Record dp per ind
				dp[i-9]=info[positionDp]
			}
			else{
				dp[i-9]=-1
			}
		}
		
		#calc mean dp
		if(nCalledGt>0){
			meanDp = totalDp / nCalledGt
		}
		else{
			meanDp
		}
		
		#####Output data
		print $1"\t"$2"\t"nCalledGt"\t"totalDp"\t"meanDp>logFileName
		
		#####Print info and gt for each sample
		printf $1 "\t" $2 "\t" > logFileName2
		
		for(i = 1; i <= numSample; i++){
			printf "%s", dp[i] > logFileName2
			if(i < numSample){printf "\t" > logFileName2}
			if(i == numSample){printf "\n" > logFileName2}
		}	
	}
}
END{
	print "nSites = " nSite
}