#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#

#THIS CODE TAKES A SAM FILE (EXAMPLE SAM, sample.sam WITH 4nt UMI WRITTEN IN THE READS' NAME IS PROVIDED  ), REFORMAT IT AND RUN THE sscs_generator.py AND duplex_generator.py   
#USAGE INSTRACTIONS
# qsub Consensus_generator.sh -d "main_folder"/"sample_name" -n "sample_name" -c "path to the sscs_generator.py duplex_generator.py files"
#
#THE DIRECTORY STRUCTURE NEED TO BE AS FOLLOWS:
#"main folder"--- 
#               |
#		----"sample1"
#		----"sample2"
#		----"sample3"---
#				|
#				----sample3.sam
#

while getopts ":n:d:c:" opt; do
        case $opt in

                n)
                        OutPrefix=$OPTARG;
                        ((count+=1))
                        ;;
                d)
                        WorkingDir=$OPTARG;
                        ((count+=1))
                        ;;
		c)
                        c=$OPTARG;
                        ((count+=1))
                        ;;
                :)
                        echo
                        exit 1;
                        ;;
                \?)
                        exit;
                        ;;
	esac
done


STARTTIME=$(date +"%Y%m%d-%H%M%S")
echo "Started: " $STARTTIME
echo "Working on " ${OutPrefix}

## load modules
module load samtools/0.1.18
module load python/2.7.11 

cd ${WorkingDir}

#############################################################
###### SPLIT THE SAM TO SINGLETONS, 2 OR MORE (2OM), 3OM, 4OM
#############################################################
awk -F"[:|\t]" '{print$10":"$11":"$15":"$9":"$8}' ${OutPrefix}.sam > col1_${OutPrefix}
paste col1_${OutPrefix} ${OutPrefix}.sam | sort -T ${WorkingDir} -k1,1n > ${OutPrefix}_sorted_forSplit.sam
awk '{print $1}' ${OutPrefix}_sorted_forSplit.sam | uniq -c > ${OutPrefix}_sorted_col1_uniq
awk '{print $1}' ${OutPrefix}_sorted_col1_uniq > ${OutPrefix}_sorted_col1_uniqCounts 
awk '{ for (i = 1; i <= $1; i++) { print $0 } }' ${OutPrefix}_sorted_col1_uniqCounts > ${OutPrefix}_sorted_col1_uniqCounts_rep
paste ${OutPrefix}_sorted_col1_uniqCounts_rep ${OutPrefix}_sorted_forSplit.sam > ${OutPrefix}_sorted_Part_of_family_of_SizeX.sam

rm ${OutPrefix}_sorted_forSplit.sam ${OutPrefix}_sorted_col1_uniq ${OutPrefix}_sorted_col1_uniqCounts ${OutPrefix}_sorted_col1_uniqCounts_rep col1_${OutPrefix}

####################################
###### RUN CONSENSUS CALLER (python)
####################################

python ${c}/sscs_generator.py -i ${OutPrefix}_sorted_Part_of_family_of_SizeX.sam -outSin ${OutPrefix}_singletons.sam -out2OM ${OutPrefix}_2OM.sam
rm ${OutPrefix}_sorted_Part_of_family_of_SizeX.sam

#########################################
###### PREPARE THE FILES FOR DUPLEX
########################################

awk '{print $3":"$4":"$8":"$6}' ${OutPrefix}_2OM.sam > ${OutPrefix}_Col1_forDuplex_Chr_P_PNEXT_CIGAR
paste ${OutPrefix}_Col1_forDuplex_Chr_P_PNEXT_CIGAR ${OutPrefix}_2OM.sam | sort -k1,1 > ${OutPrefix}_PFDM_sorted.sam
cut -f1 ${OutPrefix}_PFDM_sorted.sam | uniq -c | awk '{print $1}' > ${OutPrefix}_PFDM_sorted_col1_uniqCounts
awk '{ for (i = 1; i <= $1; i++) { print $0 } }' ${OutPrefix}_PFDM_sorted_col1_uniqCounts > ${OutPrefix}_PFDM_sorted_col1_uniqCounts_rep
paste ${OutPrefix}_PFDM_sorted_col1_uniqCounts_rep ${OutPrefix}_PFDM_sorted.sam | awk '$1 > 1 {print}'  > ${OutPrefix}_PFDM_sorted_preDuplex.sam

rm ${OutPrefix}_Col1_forDuplex_Chr_P_PNEXT_CIGAR ${OutPrefix}_PFDM_sorted_col1_uniqCounts ${OutPrefix}_PFDM_sorted_col1_uniqCounts_rep ${OutPrefix}_PFDM_sorted.sam

####################################
###### RUN DUPLEX CALLER (python) to generate duplex from SSCS-SSCS
####################################

python ${c}/duplex_generator.py -i ${OutPrefix}_PFDM_sorted_preDuplex.sam -outDuplex ${OutPrefix}_Duplex.sam -n ${OutPrefix}
rm ${OutPrefix}*PFDM* ${OutPrefix}*usedForDuplex*

ENDTIME=$(date +"%Y%m%d-%H%M%S")
echo "DONE: " $ENDTIME
