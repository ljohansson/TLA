# TLA
This repository contains scripts associated with the analysis of data created with the Cergentis TLA technology

##Dependencies and used scripts
The script splitfqbyCATG.pl can be found on https://gist.github.com/mmterpstra/2417200a96f841862c82220893490202

The script in remove_duplicates_from_fastq.py was created by Pengfei Yu and can also be copied from https://www.biostars.org/p/2733/

The Regions_of_interest.txt file should have the following format. Tab delimited. First column gene name; second column chr:start-stop. For instance: 

> RUNX1   21:36150098-36421595

The scripts 
> windowed_coverage.pl
> in_silico_enrich.pl
> human-mouse_singleGenomePlot_normalized.R
have been created by Cergentis B.V.

> Cergentis B.V.
> Yalelaan 62
> 3584 CM Utrecht
> The Netherlands
> +31 (0)30 760 16 36
> info@cergentis.com

The file human-mouse_singleGenomePlot_normalized_fix_20160920_abscut100_Inhousefiltered.R is an adaptation of the human-mouse_singleGenomePlot_normalized.R script. The following changes are made:

Function removeSmall: 
> n <- names(cnt[cnt >= 0])        		is changed to          	n <- names(cnt[cnt > 100])

Function translocationPlot: 
> abs.cut <-3       						is changed to  		   	abs.cut <-100        

Added lines to section "#remove 'chr' string"
> \# remove chr names that contain "GL"                                                            
> td <- td[grep("GL",td[,1], invert=T),]                                                         
> \# remove chr names that contain "25"                                                         
> td <- td[grep("25",td[,1], invert=T),]                                                       

Changed in section "	# check if chrY exists"
> if(sum(chroms == 'chrY')){          	is changed to     		if(sum(chroms == 'chr24')){     
> td[td[,1]=='X',1] <- length(chroms)-1	is changed to			td[td[,1]=='chr23',1] <- length(chroms)-1
> td[td[,1]=='Y',1] <- length(chroms)		is changed to			td[td[,1]=='chr24',1] <- length(chroms)
> td[td[,1]=='X',1] <- length(chroms)		is changed to			td[td[,1]=='chr23',1] <- length(chroms)


###Further dependencies: 
cutadapt (v1.8.1): https://pypi.python.org/pypi/cutadapt

samtools (v1.2 and v1.3): http://www.htslib.org/

bwa (v0.7.12): http://bio-bwa.sourceforge.net/

Perl (v5.22.0): https://www.perl.org/

R (V3.2.2 or higher): https://www.r-project.org/
