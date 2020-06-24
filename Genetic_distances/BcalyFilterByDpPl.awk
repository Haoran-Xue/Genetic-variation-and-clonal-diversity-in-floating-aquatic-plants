#Filters vcf
#Keep sites with >= thresholdInd idividuals with >= thresholdPl midPl
#Only looks at BIALLELIC sites
BEGIN{print "-----Begin BcalyFilterByDpPl.awk-----"
	print "Output file Name: " outputFileName
	print "thresholdInd: " thresholdInd
	print "minDp: " minDp
	print "maxDp: " maxDp
	print "thresholdPl: " thresholdPl
	var_vcf = outputFileName ".minDp" minDp "maxDp" maxDp "Pl" thresholdPl "Ind" thresholdInd "_var.vcf"
	invar_vcf = outputFileName ".minDp" minDp "maxDp" maxDp "Pl" thresholdPl "Ind" thresholdInd "_invar.vcf"
	nLines=0; nSite=0; nSite_bi=0; nSite_mult=0; nSitePass_var=0; nSitePass_invar=0; nSiteFail=0;
	nPassGt_var=0; nPassGt_invar=0; meanDp_var=0; meanDp_invar=0;
}
{	#count num of lines read
	nLines++;
	if(nLines%10000==0){print "Lines read: " nLines};
 
	#header info
	if($0~/\#/){
		print $0 > var_vcf
		print $0 > invar_vcf
	}
	#line is a site if not a header
	else{
		#count num of sites
		nSite++
		
		#Check if site biallelic or multiallelc
		if($5~/\,/){
			nSite_mult++
		}
		else{
			nSite_bi++
			
			#initialize Gt count
			nHomoRef=0
			nHomoAlt=0
			nHetero=0
		
			nCalledGt=0
			nPassingGt=0
			totalDp=0

			#####Find position PL and DP in FORMAT (field 9)
			#split FORMAT field by ":"
			nFormat=split($9,format,":")
			positionDp=0
			positionPl=0
			positionAd=0
			#find element in format that is "DP"
			for(f=1;f<=nFormat;f++){
				if(format[f]~/DP/){positionDp=f}
				if(format[f]~/PL/){positionPl=f}
				if(format[f]~/AD/){positionAd=f}
			}
			#if no Dp field, then skip to next line; else continue forward to determine if site passes filter
			if(positionDp==0){nSiteNoDp++;print$0;next;}else{nSiteYesDp++;}
			#####Find position PL and DP in FORMAT (field 9)

			#####Record info about each site
			for(i = 10; i <= NF; i++){
				#examine if genotype is called for ind i
				if($i~/0\/0/ || $i~/0\/1/ || $i~/1\/0/ || $i~/1\/1/){
					#count num of sites regardless of quality
					nCalledGt++
			
					##########Determine if ind pass pl filter
					#####Split field by ":"; record number of fields in nFields; record elements into info
					nFields=split($i,info,":")
				
					#####record dp of Gt
					gtDp=info[positionDp]
	
					#####Find second largest phred likelihood of genotype
					#split pl field by ","; record number of fields in nPl; record elements into pl
					#pl is in element positionPl of info
					nPl=split(info[positionPl],pl,",")
 
					#find the two likelihoods larger than 0 (format always has 0 and two non-zeros)
					num1=0;num2=0;
					for(k=1;k<=length(pl);k++){
						if(pl[k] > 0 && num1 == 0){
							num1=pl[k]
						}
						else if(pl[k] > 0 && num1 > 0){
							num2 = pl[k]
						}
					}

					#determine second largest pl (the smaller of num1,num2)
					midPl=0;
					if(num1>num2){midPl=num2}else{midPl=num1}
					
					####Find allelic depth and frequency of Ref allele
					split(info[positionAd],ad,",")
					nRef=ad[1]
					nAlt=ad[2]
					if((nRef+nAlt)>0){pRef = nRef/(nRef+nAlt)}
				
					#####Pass threshold: indicate and record which Gt
					if(midPl>=thresholdPl && gtDp>=minDp && gtDp <= maxDp){
						#count num of sites passing thresholdPl
						nPassingGt++
						#increase totalDp (only for Gt that pass filter)
						totalDp+=gtDp
						#record Gt of passing site
						if($i~/0\/0/){
							nHomoRef++
						}
						else if($i~/1\/1/){
							nHomoAlt++
						}
						else if($i~/0\/1/ || $i~/1\/0/){
							nHetero++
						}
					}
					##########Determine if ind pass pl filter
				}#end Gt called check
			}#end for loop
			#####Record info about each site
			
			
			
			#####Calc allele frequency
			totalAllele=2*(nHomoRef+nHomoAlt+nHetero)
			if(totalAllele>0){
				freqRef=(2*nHomoRef+nHetero)/totalAllele
				freqAlt=(2*nHomoAlt+nHetero)/totalAllele
			}
			#####Calc allele frequency
			
		
			#####Determine if site pass threshold and if site variable
			passQualityFilter=0
			passVariantFilter=0
  
			if(nPassingGt>=thresholdInd){
				passQualityFilter=1
			}
		
			#site variable both ref and alt alleles exist
			if(freqRef>0 && freqAlt>0){
				passVariantFilter=1
			}
			#####Determine if site pass threshold and if site variable
		
			#site pass quality filter and is variant
			if(passQualityFilter==1 && passVariantFilter==1){
				#count
				nSitePass_var++

				#data for calc avg Dp
				nPassGt_var+=nPassingGt
				meanDp_var+=totalDp; 

				print $0 > var_vcf
			}
			#site pass quality filter and but is invariant
			else if(passQualityFilter==1 && passVariantFilter==0){
				#count
				nSitePass_invar++

				#data for calc avg Dp
				nPassGt_invar+=nPassingGt
				meanDp_invar+=totalDp

				print $0 > invar_vcf
			}
			#site fails quality filter
			else{
				nSiteFail++
			}
			
		}#end biallellic/multiallelic site check
	}#end line is non-header check
}#end body
END{
	print "nSite: " nSite; print "nSite_bi: " nSite_bi; print "nSite_mult: " nSite_mult
	print "nSitePass_var: " nSitePass_var; print "nSitePass_invar: " nSitePass_invar; print "nSiteFail: " nSiteFail
	if(nPassGt_var>0){
		meanDp_var/=nPassGt_var
	}
	if(nPassGt_invar>0){
		meanDp_invar/=nPassGt_invar
	}
	print "meanDp_var: " meanDp_var; print "meanDp_invar: " meanDp_invar
}
