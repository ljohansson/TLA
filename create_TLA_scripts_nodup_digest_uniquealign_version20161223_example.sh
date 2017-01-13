#!/bin/bash

#Runname
NUMBEROFSAMPLES="8"
NUMBEROFLANES="4" #NOTE THE SCRIPT 3_MERGE_PAIRS IS HARDCODED TO FOUR LANES
CNTMAX=$(($NUMBEROFSAMPLES*$NUMBEROFLANES+1))
CNTSAMPMAX=$(($NUMBEROFSAMPLES+1))

RUNNAME="STATE_RUNNAME"						#Change per run
BARCODESAMPLE1="ATTACTCG-AGGCTATA"			#Change per run
BARCODESAMPLE2="ATTACTCG-GCCTCTAT"			#Change per run
BARCODESAMPLE3="ATTACTCG-AGGATAGG"			#Change per run
BARCODESAMPLE4="ATTACTCG-TCAGAGCC"			#Change per run
BARCODESAMPLE5="ATTACTCG-CTTCGCCT"			#Change per run
BARCODESAMPLE6="ATTACTCG-TAAGATTA"			#Change per run
BARCODESAMPLE7="ATTACTCG-ACGTCCTG"			#Change per run
BARCODESAMPLE8="ATTACTCG-GTCAGTAC"			#Change per run

#IN OUR SAMPLES WE HAVE TWO MULTIPLEXES THAT ARE MERGED
PREFIXSAMPLE1="1_MP1"  						#Change per run			#NOTE THAT STEP FOUR ONLY WORKS IF PREXIXSAMPLE1 and 2; 3 and 4; 5 and 6; 7 and 8 are each paired MP1 and 2 of the same sample
PREFIXSAMPLE2="1_MP2"  						#Change per run			#STEP 5 AND HIGHER ARE DEPENDENT ON STEP 4
PREFIXSAMPLE3="2_MP1"						#Change per run
PREFIXSAMPLE4="2_MP2"						#Change per run
PREFIXSAMPLE5="3_MP1"						#Change per run
PREFIXSAMPLE6="3_MP2"						#Change per run
PREFIXSAMPLE7="4_MP1"						#Change per run
PREFIXSAMPLE8="4_MP2"						#Change per run



#ANALYSIS FOLDERS
MAINPATH="/PATH/TO/TLA/DIR/"				#Change per run
RUNPATH="$MAINPATH/RUNFOLDER/"				#Change per run
OUTPUTPATH="$RUNPATH/analysis/"
DATAPATH="$RUNPATH/data"  #should containt original fastq.gz files
SCRIPTSFOLDER="$MAINPATH/TLA_scripts/"
LOCALOUTPUTPATH="/cygdrive/c/TLA/analysis/"						#For local processing if needed
LOCALSCRIPTSFOLDER="/cygdrive/c/TLA/scripts/"						#For local processing if needed
LOCALOUTPUTPATHIGV="C:\\\\TLA/analysis/"						#For local processing if needed
LOCALREFERENCEFOLDER="C:\\\\TLA/reference_genome/human_g1k_v37.fasta"			#For local processing if needed
LOCALIGVTOOLS="/cygdrive/c/Program\ Files/igvtools_2.3.63/IGVTools/igvtools.bat"	#For local processing if needed
REGIONSOFINTEREST="$RUNPATH/scripts/Regions_AML_TLA.txt"				#Change per Multiplex design. Make sure this file is in the correct location.

#Check if regionsfile exists or die
if [ ! -f $REGIONSOFINTEREST ]; then
    	echo "ERROR: No ROI file found"
	exit 1
fi


#Create scripts
#without extensions. Expected extensions _1.fz.gz for read 1 and _2.fq.gz for read 2. Expected 8 samples in 4 runs
FILE1="$PREFIXSAMPLE1""_""$RUNNAME""_L1_""$BARCODESAMPLE1"
FILE9="$PREFIXSAMPLE1""_""$RUNNAME""_L2_""$BARCODESAMPLE1"
FILE17="$PREFIXSAMPLE1""_""$RUNNAME""_L3_""$BARCODESAMPLE1"
FILE25="$PREFIXSAMPLE1""_""$RUNNAME""_L4_""$BARCODESAMPLE1"
FILE2="$PREFIXSAMPLE2""_""$RUNNAME""_L1_""$BARCODESAMPLE2"
FILE10="$PREFIXSAMPLE2""_""$RUNNAME""_L2_""$BARCODESAMPLE2"
FILE18="$PREFIXSAMPLE2""_""$RUNNAME""_L3_""$BARCODESAMPLE2"
FILE26="$PREFIXSAMPLE2""_""$RUNNAME""_L4_""$BARCODESAMPLE2"
FILE3="$PREFIXSAMPLE3""_""$RUNNAME""_L1_""$BARCODESAMPLE3"
FILE11="$PREFIXSAMPLE3""_""$RUNNAME""_L2_""$BARCODESAMPLE3"
FILE19="$PREFIXSAMPLE3""_""$RUNNAME""_L3_""$BARCODESAMPLE3"
FILE27="$PREFIXSAMPLE3""_""$RUNNAME""_L4_""$BARCODESAMPLE3"
FILE4="$PREFIXSAMPLE4""_""$RUNNAME""_L1_""$BARCODESAMPLE4"
FILE12="$PREFIXSAMPLE4""_""$RUNNAME""_L2_""$BARCODESAMPLE4"
FILE20="$PREFIXSAMPLE4""_""$RUNNAME""_L3_""$BARCODESAMPLE4"
FILE28="$PREFIXSAMPLE4""_""$RUNNAME""_L4_""$BARCODESAMPLE4"
FILE5="$PREFIXSAMPLE5""_""$RUNNAME""_L1_""$BARCODESAMPLE5"
FILE13="$PREFIXSAMPLE5""_""$RUNNAME""_L2_""$BARCODESAMPLE5"
FILE21="$PREFIXSAMPLE5""_""$RUNNAME""_L3_""$BARCODESAMPLE5"
FILE29="$PREFIXSAMPLE5""_""$RUNNAME""_L4_""$BARCODESAMPLE5"
FILE6="$PREFIXSAMPLE6""_""$RUNNAME""_L1_""$BARCODESAMPLE6"
FILE14="$PREFIXSAMPLE6""_""$RUNNAME""_L2_""$BARCODESAMPLE6"
FILE22="$PREFIXSAMPLE6""_""$RUNNAME""_L3_""$BARCODESAMPLE6"
FILE30="$PREFIXSAMPLE6""_""$RUNNAME""_L4_""$BARCODESAMPLE6"
FILE7="$PREFIXSAMPLE7""_""$RUNNAME""_L1_""$BARCODESAMPLE7"
FILE15="$PREFIXSAMPLE7""_""$RUNNAME""_L2_""$BARCODESAMPLE7"
FILE23="$PREFIXSAMPLE7""_""$RUNNAME""_L3_""$BARCODESAMPLE7"
FILE31="$PREFIXSAMPLE7""_""$RUNNAME""_L4_""$BARCODESAMPLE7"
FILE8="$PREFIXSAMPLE8""_""$RUNNAME""_L1_""$BARCODESAMPLE8"
FILE16="$PREFIXSAMPLE8""_""$RUNNAME""_L2_""$BARCODESAMPLE8"
FILE24="$PREFIXSAMPLE8""_""$RUNNAME""_L3_""$BARCODESAMPLE8"
FILE32="$PREFIXSAMPLE8""_""$RUNNAME""_L4_""$BARCODESAMPLE8"


#remove_duplicates_from_fastq.py

#Script 0a_Remove_duplicates_from_fastq
COUNTER=1
while [ $COUNTER -lt $CNTMAX ]
do

FILENAME="FILE"$COUNTER

TOTAL_SAMPLENAME_R1=$DATAPATH"/"${!FILENAME}"_1.fq.gz"
TOTAL_SAMPLENAME_R2=$DATAPATH"/"${!FILENAME}"_2.fq.gz"
TOTAL_SAMPLENAME_R1B=$OUTPUTPATH"/"${!FILENAME}"_1.fq"
TOTAL_SAMPLENAME_R2B=$OUTPUTPATH"/"${!FILENAME}"_2.fq"
SAMPLENAME_R1=$(basename "$TOTAL_SAMPLENAME_R1B")
SAMPLENAME_R2=$(basename "$TOTAL_SAMPLENAME_R2B")
OUTPUTNAME="0a_Remove_dup_from_fq_"${!FILENAME}".sh"

echo "Creating Script 0a_Remove_duplicates_from_fq for sample" ${!FILENAME}

echo "#!/bin/bash
#SBATCH --job-name=0a_Remove_dup_from_fq_${!FILENAME}
#SBATCH --output=0a_Remove_dup_from_fq_${!FILENAME}.out
#SBATCH --error=0a_Remove_dup_from_fq_${!FILENAME}.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Python/2.7.11-foss-2015b

gunzip -c $TOTAL_SAMPLENAME_R1 > $TOTAL_SAMPLENAME_R1B
gunzip -c $TOTAL_SAMPLENAME_R2 > $TOTAL_SAMPLENAME_R2B

cd $OUTPUTPATH
python $SCRIPTSFOLDER/remove_duplicates_from_fastq.py \\
$SAMPLENAME_R1 \\
$SAMPLENAME_R2

rm $TOTAL_SAMPLENAME_R1B
rm $TOTAL_SAMPLENAME_R2B
mv $OUTPUTPATH/Rm_dupPE_$SAMPLENAME_R1 $OUTPUTPATH/${SAMPLENAME_R1%.*}_nodup.fq
mv $OUTPUTPATH/Rm_dupPE_$SAMPLENAME_R2 $OUTPUTPATH/${SAMPLENAME_R2%.*}_nodup.fq
gzip ${TOTAL_SAMPLENAME_R1B%.*}_nodup.fq
gzip ${TOTAL_SAMPLENAME_R2B%.*}_nodup.fq

" > $OUTPUTNAME
let COUNTER=COUNTER+1
done



#Script 0_Trim_adapters
COUNTER=1
while [ $COUNTER -lt $CNTMAX ]
do

FILENAME="FILE"$COUNTER

TOTAL_SAMPLENAME_R1=$OUTPUTPATH"/"${!FILENAME}"_1_nodup.fq.gz"
TOTAL_SAMPLENAME_R2=$OUTPUTPATH"/"${!FILENAME}"_2_nodup.fq.gz"
SAMPLENAME_R1=$(basename "$TOTAL_SAMPLENAME_R1")
SAMPLENAME_R2=$(basename "$TOTAL_SAMPLENAME_R2")
OUTPUTNAME="0_Trim_adapters_"${!FILENAME}".sh"

echo "Creating Script 0_Trim_adapters for sample" ${!FILENAME}

echo "#!/bin/bash
#SBATCH --job-name=0_Trim_adapters_${!FILENAME}
#SBATCH --output=0_Trim_adapters_${!FILENAME}.out
#SBATCH --error=0_Trim_adapters_${!FILENAME}.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


module load cutadapt/1.8.1-goolf-1.7.20-Python-2.7.9

cutadapt --format=fastq \
-a CTGTCTCTTATACACATCTGACGCTGCCGACGA \
-a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
-g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
-g GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
-o $OUTPUTPATH/${SAMPLENAME_R1%%.*}_trimmed.fq.gz $TOTAL_SAMPLENAME_R1

cutadapt --format=fastq \
-a CTGTCTCTTATACACATCTGACGCTGCCGACGA \
-a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
-g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
-g GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
-o $OUTPUTPATH/${SAMPLENAME_R2%%.*}_trimmed.fq.gz $TOTAL_SAMPLENAME_R2

" > $OUTPUTNAME
let COUNTER=COUNTER+1
done


#Script 1_In_Silico_digest
COUNTER=1
while [ $COUNTER -lt $CNTMAX ]
do

FILENAME="FILE"$COUNTER

TOTAL_SAMPLENAME_R1=$OUTPUTPATH"/"${!FILENAME}"_1_nodup_trimmed.fq.gz"
TOTAL_SAMPLENAME_R2=$OUTPUTPATH"/"${!FILENAME}"_2_nodup_trimmed.fq.gz"
SAMPLENAME_R1=$(basename "$TOTAL_SAMPLENAME_R1")
SAMPLENAME_R2=$(basename "$TOTAL_SAMPLENAME_R2")
OUTPUTNAME="1_In_silico_digest_"${!FILENAME}".sh"

echo "Creating Script 1_In_silico_digest for sample" ${!FILENAME}

echo "#!/bin/bash
#SBATCH --job-name=1_In_silico_digest_${!FILENAME}
#SBATCH --output=1_In_silico_digest_${!FILENAME}.out
#SBATCH --error=1_In_silico_digest_${!FILENAME}.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Perl/5.22.0-foss-2015b-bare

perl $SCRIPTSFOLDER/splitfqbyCATG.pl \\
$TOTAL_SAMPLENAME_R1 \\
${TOTAL_SAMPLENAME_R1%%.*}_digested.fq.gz

perl $SCRIPTSFOLDER/splitfqbyCATG.pl \\
$TOTAL_SAMPLENAME_R2 \\
${TOTAL_SAMPLENAME_R2%%.*}_digested.fq.gz

" > $OUTPUTNAME
let COUNTER=COUNTER+1
done


#Script 2_bwa_sw_alignment
COUNTER=1
while [ $COUNTER -lt $CNTMAX ]
do

FILENAME="FILE"$COUNTER

TOTAL_SAMPLENAME_R1=$OUTPUTPATH"/"${!FILENAME}"_1_nodup_trimmed_digested.fq.gz"
TOTAL_SAMPLENAME_R2=$OUTPUTPATH"/"${!FILENAME}"_2_nodup_trimmed_digested.fq.gz"
OUTPUTNAME="2_bwa_sw_alignment_"${!FILENAME}".sh"

echo "Creating Script 2_bwa_sw_alignment for sample" ${!FILENAME}

echo "#!/bin/bash
#SBATCH --job-name=2_bwa_sw_alignment_${!FILENAME}
#SBATCH --output=2_bwa_sw_alignment_${!FILENAME}.out
#SBATCH --error=2_bwa_sw_alignment_${!FILENAME}.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 12gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load BWA/0.7.12-goolf-1.7.20
module load SAMtools/1.2-foss-2015b #DOES NOT WORK WITH VERSION 1.3

bwa bwasw \
-t 20 -b 7 \
/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta $TOTAL_SAMPLENAME_R1 | \
samtools view -bS - | \
samtools sort - ${TOTAL_SAMPLENAME_R1%%.*}

bwa bwasw \
-t 20 -b 7 \
/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta $TOTAL_SAMPLENAME_R2 | \
samtools view -bS - | \
samtools sort - ${TOTAL_SAMPLENAME_R2%%.*}

" > $OUTPUTNAME
let COUNTER=COUNTER+1
done




#Script 3_Merge_pairs
COUNTER=1
while [ $COUNTER -lt $CNTSAMPMAX ]
do

FILENAME_int="FILE"$COUNTER
FILENAME="${!FILENAME_int}"
FILENAME2="${FILENAME/_L1_/_L1234_}"


TOTAL_SAMPLENAME_L1_R1="$OUTPUTPATH"/"$FILENAME""_1_nodup_trimmed_digested.bam"
TOTAL_SAMPLENAME_L1_R2="$OUTPUTPATH"/"$FILENAME""_2_nodup_trimmed_digested.bam"
TOTAL_SAMPLENAME_L2_R1="${TOTAL_SAMPLENAME_L1_R1/_L1_/_L2_}"
TOTAL_SAMPLENAME_L2_R2="${TOTAL_SAMPLENAME_L1_R2/_L1_/_L2_}"
TOTAL_SAMPLENAME_L3_R1="${TOTAL_SAMPLENAME_L1_R1/_L1_/_L3_}"
TOTAL_SAMPLENAME_L3_R2="${TOTAL_SAMPLENAME_L1_R2/_L1_/_L3_}"
TOTAL_SAMPLENAME_L4_R1="${TOTAL_SAMPLENAME_L1_R1/_L1_/_L4_}"
TOTAL_SAMPLENAME_L4_R2="${TOTAL_SAMPLENAME_L1_R2/_L1_/_L4_}"
OUTPUTNAME="3_Merge_Pairs_"$FILENAME2".sh"
SAMPLENAME_INT="${TOTAL_SAMPLENAME_L1_R1/_1_nodup_trimmed_digested.bam/_nodup_trimmed_digested.bam}"
SAMPLENAME="${SAMPLENAME_INT/_L1_/_L1234_}"


echo "Creating Script 3_Merge_pairs for sample" $FILENAME2

echo "#!/bin/bash
#SBATCH --job-name=3_Merge_pairs_$FILENAME2
#SBATCH --output=3_Merge_pairs_$FILENAME2.out
#SBATCH --error=3_Merge_pairs_$FILENAME2.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=2:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load SAMtools/1.3-foss-2015b

echo \"Merging bam files\"
samtools merge -f $SAMPLENAME $TOTAL_SAMPLENAME_L1_R1 $TOTAL_SAMPLENAME_L1_R2 $TOTAL_SAMPLENAME_L2_R1 $TOTAL_SAMPLENAME_L2_R2 $TOTAL_SAMPLENAME_L3_R1 $TOTAL_SAMPLENAME_L3_R2 $TOTAL_SAMPLENAME_L4_R1 $TOTAL_SAMPLENAME_L4_R2

echo \"Creating bam index\"
samtools index $SAMPLENAME > $SAMPLENAME.bai

" > $OUTPUTNAME
let COUNTER=COUNTER+1
done



#Script 4a_Merge_pairs_MP1and2 #NOTE THAT THIS STEP SHOULD ONLY BE RUN WHEN SAMPLES 1 and 2; 3 and 4; 5 and 6; 7 and 8 are paired MP1 and MP2 sets
COUNTER=1
while [ $COUNTER -lt $CNTSAMPMAX ]
do

#MP1 file set names
FILENAME_int="FILE"$COUNTER
FILENAME="${!FILENAME_int}"
FILENAME2="${FILENAME/_L1_/_L1234_}"

TOTAL_SAMPLENAME_L1_R1_MP1="$OUTPUTPATH"/"$FILENAME""_1_nodup_trimmed_digested.bam"
OUTPUTNAME="4a_Merge_MP1and2_"$FILENAME2".sh"
SAMPLENAME_INT_MP1="${TOTAL_SAMPLENAME_L1_R1_MP1/_1_nodup_trimmed_digested.bam/_nodup_trimmed_digested.bam}"
SAMPLENAME_MP1="${SAMPLENAME_INT_MP1/_L1_/_L1234_}"

#MP2 file set names
let COUNTER2=COUNTER+1
FILENAME_int="FILE"$COUNTER2
FILENAME3="${!FILENAME_int}"
FILENAME4="${FILENAME/_L1_/_L1234_}"

TOTAL_SAMPLENAME_L1_R1_MP2="$OUTPUTPATH"/"$FILENAME3""_1_nodup_trimmed_digested.bam"
SAMPLENAME_INT_MP2="${TOTAL_SAMPLENAME_L1_R1_MP2/_1_nodup_trimmed_digested.bam/_nodup_trimmed_digested.bam}"
SAMPLENAME_MP2="${SAMPLENAME_INT_MP2/_L1_/_L1234_}"



echo "Creating Script 4a_Merge_pairs MP1and2 for sample" $FILENAME2

echo "#!/bin/bash
#SBATCH --job-name=4a_Merge_pairs_MP1and2_$FILENAME2
#SBATCH --output=4a_Merge_pairs_MP1and2_$FILENAME2.out
#SBATCH --error=4a_Merge_pairs_MP1and2_$FILENAME2.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=2:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load SAMtools/1.3-foss-2015b

echo \"Merging bam files\"
samtools merge -f ${SAMPLENAME_MP1%.*}_MP1and2.bam $SAMPLENAME_MP1 $SAMPLENAME_MP2

echo \"Creating bam index\"
samtools index ${SAMPLENAME_MP1%.*}_MP1and2.bam > ${SAMPLENAME_MP1%.*}_MP1and2.bam.bai

" > $OUTPUTNAME
let COUNTER=COUNTER+2
done


#Script 4c_Filter_unique_MP1and2
COUNTER=1
while [ $COUNTER -lt $CNTSAMPMAX ]
do

#MP1 file set names
FILENAME_int="FILE"$COUNTER
FILENAME="${!FILENAME_int}"
FILENAME2="${FILENAME/_L1_/_L1234_}"

TOTAL_SAMPLENAME_L1_R1_MP1="$OUTPUTPATH"/"$FILENAME""_1_nodup_trimmed_digested_MP1and2.bam"
OUTPUTNAME="4c_Filter_unique_MP1and2_"$FILENAME2".sh"
SAMPLENAME_INT_MP1="${TOTAL_SAMPLENAME_L1_R1_MP1/_1_nodup_trimmed_digested_MP1and2.bam/_nodup_trimmed_digested_MP1and2.bam}"
SAMPLENAME_MP1="${SAMPLENAME_INT_MP1/_L1_/_L1234_}"


echo "Creating Script 4c_Filter_unique_MP1and2 for sample" $FILENAME2

echo "#!/bin/bash
#SBATCH --job-name=4c_Filter_unique_MP1and2_$FILENAME2
#SBATCH --output=4c_Filter_unique_MP1and2_$FILENAME2.out
#SBATCH --error=4c_Filter_unique_MP1and2_$FILENAME2.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=2:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load SAMtools/1.3-foss-2015b

echo \"Filtering unique reads from bam files\"
samtools view -bSq 3 $SAMPLENAME_MP1 > ${SAMPLENAME_MP1%.*}_uniquely_aligned.bam

echo \"Creating bam index\"
samtools index ${SAMPLENAME_MP1%.*}_uniquely_aligned.bam > ${SAMPLENAME_MP1%.*}_uniquely_aligned.bam.bai

" > $OUTPUTNAME
let COUNTER=COUNTER+2
done


#Script 5a_generate_genome_coverage_plots

COUNTER=1
while [ $COUNTER -lt $CNTSAMPMAX ]
do

FILENAME_int="FILE"$COUNTER
FILENAME="${!FILENAME_int}"
FILENAME2="${FILENAME/_L1_/_L1234_}"
OUTPUTNAME="5a_generate_genome_coverage_plots_"$FILENAME2".sh"
TOTAL_SAMPLENAME="$OUTPUTPATH"/"$FILENAME2""_nodup_trimmed_digested_MP1and2_uniquely_aligned.bam"

echo "Creating Script 5a_generate_genome_coverage_plots" $FILENAME2

echo "#!/bin/bash
#SBATCH --job-name=5a_generate_genome_coverage_plots_$FILENAME2
#SBATCH --output=5a_generate_genome_coverage_plots_$FILENAME2.out
#SBATCH --error=5a_generate_genome_coverage_plots_$FILENAME2.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

#The cluster installation does not have X11 available for R to create png images. Run script 5b locally to create the png

module load SAMtools/1.3-foss-2015b
#module load R/3.0.2

samtools depth -aa -Q 30 $TOTAL_SAMPLENAME| perl $SCRIPTSFOLDER/windowed_coverage.pl > ${TOTAL_SAMPLENAME%.*}_coverage.txt      ##ADDED -aa flag to fill in all bins, even if the coverage is 0. DEPENDENT ON SAMtools v 1.3

#The png image cannot be created on the cluster (X11 is not available). Use a local R installation and run step 5b locally to create the image
#Rscript $SCRIPTSFOLDER/singleGenomePlot.R ${TOTAL_SAMPLENAME%.*}_coverage.txt

#USE THIS SCRIPT ON A LOCAL COMPUTER WITH IGVTOOLS INSTALLED

" > $OUTPUTNAME
let COUNTER=COUNTER+2
done



#scripts 7a in silico enrichment all reads
COUNTER=1
while [ $COUNTER -lt $CNTSAMPMAX ]
do

	FILENAME_int="FILE"$COUNTER
	FILENAME="${!FILENAME_int}"
	FILENAME2="${FILENAME/_L1_/_L1234_}"
	TOTAL_SAMPLENAME="$OUTPUTPATH"/"$FILENAME2""_nodup_trimmed_digested_MP1and2_uniquely_aligned.bam"

	NUMBEROFPROBES=$(wc -l $REGIONSOFINTEREST | awk '{print $1}')



	COUNTERB=1
	while [ $COUNTERB -le $NUMBEROFPROBES ]
	do

		FILENAME3=$FILENAME2"_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_"$COUNTERB
		OUTPUTNAME="7a_in_silico_enrich_"$FILENAME3".sh"
		ROI=$(head -$COUNTERB $REGIONSOFINTEREST | tail -1 | awk '{print $2}' )
		TOTAL_SAMPLENAME_ROI="$OUTPUTPATH"/"$FILENAME3"

		echo "Creating Script 7a_in_silico_enrich_"$FILENAME3

		echo "#!/bin/bash
#SBATCH --job-name=7a_in_silico_enrich_$FILENAME3
#SBATCH --output=7a_in_silico_enrich_$FILENAME3.out
#SBATCH --error=7a_in_silico_enrich_$FILENAME3.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=2:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load SAMtools/1.3-foss-2015b

perl $SCRIPTSFOLDER/in_silico_enrich.pl \\
$TOTAL_SAMPLENAME \\
$ROI \\
$TOTAL_SAMPLENAME_ROI

samtools index $TOTAL_SAMPLENAME_ROI.bam $TOTAL_SAMPLENAME_ROI.bam.bai
" > $OUTPUTNAME
		let COUNTERB=COUNTERB+1
	done


	let COUNTER=COUNTER+2
done




#Script 8a_generate_genome_coverage_plots_ROI

COUNTER=1
while [ $COUNTER -lt $CNTSAMPMAX ]
do

FILENAME_int="FILE"$COUNTER
FILENAME="${!FILENAME_int}"
FILENAME2="${FILENAME/_L1_/_L1234_}"
TOTAL_SAMPLENAME="$OUTPUTPATH"/"$FILENAME2""_nodup_trimmed_digested_MP1and2_uniquely_aligned.bam"


        COUNTERB=1
        while [ $COUNTERB -le $NUMBEROFPROBES ]
        do

                FILENAME3=$FILENAME2"_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_"$COUNTERB
                OUTPUTNAME="8a_generate_genome_coverage_plots_"$FILENAME3".sh"
                ROI=$(head -$COUNTER $REGIONSOFINTEREST | tail -1 | awk '{print $2}' )
                TOTAL_SAMPLENAME_ROI="$OUTPUTPATH"/"$FILENAME3"

                echo "Creating Script 8a_generate_genome_coverage_plots_"$FILENAME3

                echo "#!/bin/bash
#SBATCH --job-name=8a_generate_genome_coverage_plots_$FILENAME3
#SBATCH --output=8a_generate_genome_coverage_plots_$FILENAME3.out
#SBATCH --error=8a_generate_genome_coverage_plots_$FILENAME3.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

#The cluster installation does not have X11 available for R to create png images. Run script 5b locally to create the png

module load SAMtools/1.3-foss-2015b
#module load R/3.0.2

samtools depth -aa -Q 30 $TOTAL_SAMPLENAME_ROI.bam | perl $SCRIPTSFOLDER/windowed_coverage.pl > $TOTAL_SAMPLENAME_ROI"_coverage.txt"

#The png image cannot be created on the cluster (X11 is not available). Use a local R installation and run step 5b locally to create the image
#Rscript $SCRIPTSFOLDER/singleGenomePlot.R ${TOTAL_SAMPLENAME%.*}_coverage.txt

#USE THIS SCRIPT ON A LOCAL COMPUTER WITH IGVTOOLS INSTALLED

" > $OUTPUTNAME
                let COUNTERB=COUNTERB+1
        done

let COUNTER=COUNTER+2
done



#Script 8c_generate_genome_coverage_plots_ROI

COUNTER=1
while [ $COUNTER -lt $CNTSAMPMAX ]
do

FILENAME_int="FILE"$COUNTER
FILENAME="${!FILENAME_int}"
FILENAME2="${FILENAME/_L1_/_L1234_}"
LOCALOUTPUTPATH="$LOCALOUTPUTPATH"
LOCALSCRIPTSFOLDER="$LOCALSCRIPTSFOLDER"
TOTAL_SAMPLENAME="$LOCALOUTPUTPATH"/"$FILENAME2""_nodup_trimmed_digested_MP1and2_uniquely_aligned.bam"

        COUNTERB=1
        while [ $COUNTERB -le $NUMBEROFPROBES ]
        do

                FILENAME3=$FILENAME2"_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_"$COUNTERB
                OUTPUTNAME="8c_generate_genome_coverage_plots_LOCAL_"$FILENAME3".sh"
                ROI=$(head -$COUNTER $REGIONSOFINTEREST | tail -1 | awk '{print $2}' )
                TOTAL_SAMPLENAME_ROI="$LOCALOUTPUTPATH"/"$FILENAME3"".bam"

                echo "Creating Script 8c_generate_genome_coverage_plots_"$FILENAME3

                echo "#!/bin/bash


#The cluster installation does not have X11 available for R to create png images. Run script 5b locally to create the png

#The png image can be created using a local R installation
Rscript $LOCALSCRIPTSFOLDER/human-mouse_singleGenomePlot_normalized_fix_20160920_abscut100_Inhousefiltered.R $TOTAL_SAMPLENAME_ROI"_coverage.txt"

" > $OUTPUTNAME
                let COUNTERB=COUNTERB+1
        done

let COUNTER=COUNTER+2
done



#Script 9a_generate_tdf_files_ROI

COUNTER=1
while [ $COUNTER -lt $CNTSAMPMAX ]
do

FILENAME_int="FILE"$COUNTER
FILENAME="${!FILENAME_int}"
FILENAME2="${FILENAME/_L1_/_L1234_}"
TOTAL_SAMPLENAME="$LOCALOUTPUTPATHIGV"/"$FILENAME2""_nodup_trimmed_digested_MP1and2_uniquely_aligned.bam"


        COUNTERB=1
        while [ $COUNTERB -le $NUMBEROFPROBES ]
        do

                FILENAME3=$FILENAME2"_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_"$COUNTERB
                OUTPUTNAME="9a_generate_tdf_files_ROI_LOCAL_"$FILENAME3".sh"
                ROI=$(head -$COUNTER $REGIONSOFINTEREST | tail -1 | awk '{print $2}' )
                TOTAL_SAMPLENAME_ROI="$LOCALOUTPUTPATHIGV"/"$FILENAME3"".bam"

                echo "Creating Script 9a_generate_tdf_files_"$FILENAME3

                echo "#!/bin/bash

$LOCALIGVTOOLS count -w 25 -f mean,median $TOTAL_SAMPLENAME_ROI $TOTAL_SAMPLENAME_ROI.tdf $LOCALREFERENCEFOLDER

" > $OUTPUTNAME
                let COUNTERB=COUNTERB+1
        done
let COUNTER=COUNTER+2
done


#Script 9b TLA_translocation_Inhouse_peak_detection_ROI
NUMBEROFPROBES=$(wc -l $REGIONSOFINTEREST | awk '{print $1}')

COUNTER=1
while [ $COUNTER -lt $CNTSAMPMAX ]
do

FILENAME_int="FILE"$COUNTER
FILENAME="${!FILENAME_int}"
FILENAME2="${FILENAME/_L1_/_L1234_}"
TOTAL_SAMPLENAME="$OUTPUTPATH"/"$FILENAME2""_nodup_trimmed_digested_MP1and2_uniquely_aligned.bam"


        COUNTERB=1
        while [ $COUNTERB -le $NUMBEROFPROBES ]
        do

                FILENAME3=$FILENAME2"_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_"$COUNTERB
                OUTPUTNAME="9b_peak_filter_"$FILENAME3".sh"
                ROI=$(head -$COUNTER $REGIONSOFINTEREST | tail -1 | awk '{print $2}' )
                TOTAL_SAMPLENAME_ROI="$OUTPUTPATH"/"$FILENAME3"


                echo "Creating Script 9b_peak_filter_"$FILENAME3

                echo "#!/bin/bash
#SBATCH --job-name=9b_peak_filter_$FILENAME3
#SBATCH --output=9b_peak_filter_$FILENAME3.out
#SBATCH --error=9b_peak_filter_$FILENAME3.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 1gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


module load R/3.2.1-foss-2015b

#Filtering noise peaks from the coverage files
Rscript /groups/umcg-gdio/tmp04/umcg-edeboer/TLA/TLA_scripts/TLA_translocation_peak_detection_general_20161116_part1.R \
${TOTAL_SAMPLENAME_ROI%.*}_coverage.txt \
$OUTPUTPATH
" > $OUTPUTNAME
        let COUNTERB=COUNTERB+1
        done

let COUNTER=COUNTER+2
done



#Script 9c TLA_translocation_Inhouse_peak_detection_ROI_part2
NUMBEROFPROBES=$(wc -l $REGIONSOFINTEREST | awk '{print $1}')

COUNTER=1
while [ $COUNTER -lt $CNTSAMPMAX ]
do

FILENAME_int="FILE"$COUNTER
FILENAME="${!FILENAME_int}"
FILENAME2="${FILENAME/_L1_/_L1234_}"
TOTAL_SAMPLENAME="$OUTPUTPATH"/"$FILENAME2""_nodup_trimmed_digested_MP1and2_uniquely_aligned.bam"


        COUNTERB=1
        while [ $COUNTERB -le $NUMBEROFPROBES ]
        do

                FILENAME3=$FILENAME2"_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_"$COUNTERB
                OUTPUTNAME="9c_peak_filter_part2_"$FILENAME3".sh"
                ROI=$(head -$COUNTER $REGIONSOFINTEREST | tail -1 | awk '{print $2}' )
                TOTAL_SAMPLENAME_ROI="$OUTPUTPATH"/"tmp"/"$FILENAME3"


                echo "Creating Script 9c_peak_filter_part2_"$FILENAME3

                echo "#!/bin/bash
#SBATCH --job-name=9c_peak_filter_part2_$FILENAME3
#SBATCH --output=9c_peak_filter_part2_$FILENAME3.out
#SBATCH --error=9c_peak_filter_part2_$FILENAME3.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 1gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


module load R/3.2.1-foss-2015b

#Filtering noise peaks from the coverage files
Rscript /groups/umcg-gdio/tmp04/umcg-edeboer/TLA/TLA_scripts/TLA_translocation_peak_detection_general_20161116_part2.R \
${TOTAL_SAMPLENAME_ROI%.*}_coverage.txt \
$OUTPUTPATH
" > $OUTPUTNAME
        let COUNTERB=COUNTERB+1
        done

let COUNTER=COUNTER+2
done


#Script 9d TLA_translocation_Inhouse_peak_detection_ROI_part3
NUMBEROFPROBES=$(wc -l $REGIONSOFINTEREST | awk '{print $1}')

COUNTER=1
while [ $COUNTER -lt $CNTSAMPMAX ]
do

FILENAME_int="FILE"$COUNTER
FILENAME="${!FILENAME_int}"
FILENAME2="${FILENAME/_L1_/_L1234_}"
TOTAL_SAMPLENAME="$OUTPUTPATH"/"$FILENAME2""_nodup_trimmed_digested_MP1and2_uniquely_aligned.bam"


        COUNTERB=1
        while [ $COUNTERB -le $NUMBEROFPROBES ]
        do

                FILENAME3=$FILENAME2"_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_"$COUNTERB
                OUTPUTNAME="9d_peak_filter_part3_"$FILENAME3".sh"
                ROI=$(head -$COUNTER $REGIONSOFINTEREST | tail -1 | awk '{print $2}' )
                TOTAL_SAMPLENAME_ROI="$OUTPUTPATH"/"$FILENAME3"


                echo "Creating Script 9d_peak_filter_part3_"$FILENAME3

                echo "#!/bin/bash
#SBATCH --job-name=9d_peak_filter_part3_$FILENAME3
#SBATCH --output=9d_peak_filter_part3_$FILENAME3.out
#SBATCH --error=9d_peak_filter_part3_$FILENAME3.err
#SBATCH --qos=regular
#SBATCH --constraint=tmp04
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 1gb
#SBATCH --nodes 1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


module load R/3.2.1-foss-2015b

#Filtering noise peaks from the coverage files
Rscript /groups/umcg-gdio/tmp04/umcg-edeboer/TLA/TLA_scripts/TLA_translocation_peak_detection_general_20161116_part3.R \
${TOTAL_SAMPLENAME_ROI%.*}_coverage.txt \
$OUTPUTPATH
" > $OUTPUTNAME
        let COUNTERB=COUNTERB+1
        done

let COUNTER=COUNTER+2
done


