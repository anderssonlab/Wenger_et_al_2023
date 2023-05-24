#!/bin/bash
## Script to process SCAR bam files and obtain
## Partition scores and CPMs in windows
## Commands by Maria Dalby (0310-2018)
## Modified by Nicolas Alcaraz (3003-2019)

## QC and mapping ##
## Install Kundaje pipeline and genome data from https://github.com/kundajelab/chipseq_pipeline#pipeline
## and use dedup bams as inputs
## requires scripts partition_smooth.pl and bdg2bw to be installed
## in addtion to samtools, bedtools, bigWigAverageOverBed and bedGraphToBigWig

BAM_FILE=""
PAIRED_END="false"
INPUT_PAIRED_END="false"
INPUT_BAM_FILE=""
OUT_BASE_DIR="MAPPING_mm10"
CHROM_SIZES="mm10.chrom.sizes"
WIN=1000
RADIUS=30
DRADIUS=30
ZRADIUS=1


usage() {
        echo
        echo ${0##/*}" Usage: SCAR_seq_process.sh [OPTIONS] <bam_file>"
        echo "Script that 1) Splits SCAR-seq <bam_file> by strand."
        echo "            2) Smooths the signal."
        echo "            3) Computes the breakpoints and the RFD scores."
        echo "requires partition_smooth.pl script to be in the path"
        echo
        echo "OPTIONS:"
        echo "[-b STRING]   BAM file of sample to analyze (required)"
        echo "[-p STRING]   if sample BAM is paired-end (default false)"
        echo "[-e STRING]   BAM file of matched input to sample (optional)"
	echo "[-i STRING]   if input BAM is paired-end (default false)"
        echo "[-o STRING]   Directory where all output will be written to"
        echo "[-g STRING]   Chromosome sizes file (required)"
        echo "[-w INTEGER]  Window size for smoothing (default 1000)"
        echo "[-r INTEGER]  Radius for smoothing (default 30)"
        echo "[-d INTEGER]  D-radius for smoothing (default 30)"
        echo "[-z INTEGER]  Z-raduus for smoothing (default 1)"
        exit
}

## show usage if '-h' or  '--help' is the first argument or no argument is given
case $1 in
        ""|"-h"|"--help") usage ;;
esac

echo "Running \"`basename $0`\" with parameters \"$*\""

## read the paramters
while getopts b:p:e:o:i:g:w:r:d:z: opt; do
        case ${opt} in
                b) BAM_FILE=${OPTARG};;
                p) PAIRED_END=${OPTARG};;
                e) INPUT_BAM_FILE=${OPTARG};;
                i) INPUT_PAIRED_END=${OPTARG};;
                o) OUT_BASE_DIR=${OPTARG};;
                g) CHROM_SIZES=${OPTARG};;
                w) WIN=${OPTARG};;
                r) RADIUS=${OPTARG};;
                d) DRADIUS=${OPTARG};;
                z) ZRADIUS=${OPTARG};;
                *) echo "unrecognized option ${opt}"; usage;;
        esac
done


if [ ! -d ${OUT_BASE_DIR} ]; then
  mkdir ${OUT_BASE_DIR}
  chmod -R g+ws ${OUT_BASE_DIR}
fi

OUT_BAM_DIR="${OUT_BASE_DIR}/bam_files"
OUT_BW_DIR="${OUT_BASE_DIR}/bw_files"
OUT_RFD_DIR="${OUT_BASE_DIR}/rfd_files"

if [ ! -d ${OUT_BAM_DIR} ]; then mkdir ${OUT_BAM_DIR}; chmod -R g+ws ${OUT_BAM_DIR}; fi
if [ ! -d ${OUT_BW_DIR} ]; then mkdir ${OUT_BW_DIR}; chmod -R g+ws ${OUT_BW_DIR}; fi
if [ ! -d ${OUT_RFD_DIR} ]; then mkdir ${OUT_RFD_DIR}; chmod -R g+ws ${OUT_RFD_DIR}; fi


PREFIX=$(basename ${BAM_FILE} .bam)
INPUT_PREFIX=$(basename ${INPUT_BAM_FILE} .bam)

TMP_DIR=${OUT_BASE_DIR}/${PREFIX}_tmp

if [ ! -d ${TMP_DIR} ]; then mkdir ${TMP_DIR}; chmod -R g+ws ${TMP_DIR}; fi


OUT_F_BAM=${OUT_BAM_DIR}/${PREFIX}_F.bam
OUT_R_BAM=${OUT_BAM_DIR}/${PREFIX}_R.bam

INPUT_F_BAM=${OUT_BAM_DIR}/${INPUT_PREFIX}_F.bam
INPUT_R_BAM=${OUT_BAM_DIR}/${INPUT_PREFIX}_R.bam



PREFIX=$(basename ${BAM_FILE} .bam)_SE
if [ ${PAIRED_END} == 'true' ]; then
        PREFIX=$(basename ${BAM_FILE} .bam)_PE
fi

INPUT_PREFIX=$(basename ${INPUT_BAM_FILE} .bam)_SE
if [ ${INPUT_PAIRED_END} == 'true' ]; then
        INPUT_PREFIX=$(basename ${INPUT_BAM_FILE} .bam)_PE
fi


OUT_F_BG=${OUT_BW_DIR}/${PREFIX}_F.bdg
OUT_R_BG=${OUT_BW_DIR}/${PREFIX}_R.bdg
OUT_F_BW=${OUT_BW_DIR}/${PREFIX}_F.bw
OUT_R_BW=${OUT_BW_DIR}/${PREFIX}_R.bw

INPUT_F_BG=${OUT_BW_DIR}/${INPUT_PREFIX}_F.bdg
INPUT_R_BG=${OUT_BW_DIR}/${INPUT_PREFIX}_R.bdg
INPUT_F_BW=${OUT_BW_DIR}/${INPUT_PREFIX}_F.bw
INPUT_R_BW=${OUT_BW_DIR}/${INPUT_PREFIX}_R.bw



if [[ ! -f ${OUT_F_BAM} ]] || [[ ! -f ${OUT_R_BAM} ]]; then
  echo "Splitting bam file in forward and reverse strands..."
  samtools view -F 20 -h ${BAM_FILE} | samtools view -Sb -h > ${OUT_F_BAM}
  samtools view -f 16 -h ${BAM_FILE} | samtools view -Sb -h > ${OUT_R_BAM}
  samtools index ${OUT_F_BAM}
  samtools index ${OUT_R_BAM}
fi

if [[ ! -f ${INPUT_F_BAM} ]] || [[ ! -f ${INPUT_R_BAM} ]]; then
  echo "Splitting bam file in forward and reverse strands..."
  samtools view -F 20 -h ${INPUT_BAM_FILE} | samtools view -Sb -h > ${INPUT_F_BAM}
  samtools view -f 16 -h ${INPUT_BAM_FILE} | samtools view -Sb -h > ${INPUT_R_BAM}
  samtools index ${INPUT_F_BAM}
  samtools index ${INPUT_R_BAM}
fi


if [[ ! -f ${OUT_F_BW} ]] || [[ ! -f ${OUT_R_BW} ]]; then
  echo "Creating bedgraph files...."
  if [ ${PAIRED_END} == 'false' ]; then
  	genomeCoverageBed -ibam ${OUT_F_BAM} -bg -fs 150 -g ${CHROM_SIZES} > ${OUT_F_BG}
  	genomeCoverageBed -ibam ${OUT_R_BAM} -bg -fs 150 -g ${CHROM_SIZES} > ${OUT_R_BG}
  else
	genomeCoverageBed -ibam ${OUT_F_BAM} -bg -pc -g ${CHROM_SIZES} > ${OUT_F_BG}
        genomeCoverageBed -ibam ${OUT_R_BAM} -bg -pc -g ${CHROM_SIZES} > ${OUT_R_BG}
  fi
  echo "Creating bigwig files.... ${OUT_F_BW}"
  bdg2bw ${OUT_F_BG} ${CHROM_SIZES} ${OUT_F_BW}
  bdg2bw ${OUT_R_BG} ${CHROM_SIZES} ${OUT_R_BW}
fi


if [[ ! -f ${INPUT_F_BW} ]] || [[ ! -f ${INPUT_R_BW} ]]; then
  echo "Creating bedgraph files...."
  if [ ${INPUT_PAIRED_END} == 'false' ]; then
        genomeCoverageBed -ibam ${INPUT_F_BAM} -bg -fs 150 -g ${CHROM_SIZES} > ${INPUT_F_BG}
        genomeCoverageBed -ibam ${INPUT_R_BAM} -bg -fs 150 -g ${CHROM_SIZES} > ${INPUT_R_BG}
  else
        genomeCoverageBed -ibam ${INPUT_F_BAM} -bg -pc -g ${CHROM_SIZES} > ${INPUT_F_BG}
        genomeCoverageBed -ibam ${INPUT_R_BAM} -bg -pc -g ${CHROM_SIZES} > ${INPUT_R_BG}
  fi
  echo "Creating bigwig files.... ${INPUT_F_BW}"
  bdg2bw ${INPUT_F_BG} ${CHROM_SIZES} ${INPUT_F_BW}
  bdg2bw ${INPUT_R_BG} ${CHROM_SIZES} ${INPUT_R_BW}
fi



### Smooth and partion
echo "Smoothing...."
bedtools makewindows -i srcwinnum -g ${CHROM_SIZES} -w ${WIN} -s ${WIN} > ${TMP_DIR}/windows.bed

echo "Splitting windows...."
## Split windows on chromosome
for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
	awk -v c=${CHR} '$1==c' ${TMP_DIR}/windows.bed > ${TMP_DIR}/${CHR}_windows.bed &
done
wait

echo "Aggregating counts..."
## Aggregate counts in windows
for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
	bigWigAverageOverBed ${OUT_F_BW} ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_F.tab &
	bigWigAverageOverBed ${OUT_R_BW} ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_R.tab &
	bigWigAverageOverBed ${INPUT_F_BW} ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_input_F.tab &
        bigWigAverageOverBed ${INPUT_R_BW} ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_input_R.tab &
done
wait




## CPM normalise strands individually (with pseudocount 1)
cat ${TMP_DIR}/*_windows_F.tab ${TMP_DIR}/*_windows_R.tab > ${TMP_DIR}/tmp_windows.bed
CPM=$(awk '{SUM += ($4+1)} END {print SUM/10^6}' ${TMP_DIR}/tmp_windows.bed)

cat ${TMP_DIR}/*_windows_input_F.tab ${TMP_DIR}/*_windows_input_R.tab > ${TMP_DIR}/tmp_windows_input.bed
CPM2=$(awk '{SUM += ($4+1)} END {print SUM/10^6}' ${TMP_DIR}/tmp_windows_input.bed)



echo "$CPM"
echo "CPM factors: $CPM"
echo "CPM input factors: $CPM2"

for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
	echo "Normalizing for chromosome ${CHR}..."
	awk -v c=${CPM} 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${TMP_DIR}/${CHR}_windows_F.tab > \
		${TMP_DIR}/${CHR}_windows_F_CPM.tab
	awk -v c=${CPM} 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${TMP_DIR}/${CHR}_windows_R.tab > \
		${TMP_DIR}/${CHR}_windows_R_CPM.tab
	awk -v c=${CPM2} 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${TMP_DIR}/${CHR}_windows_input_F.tab > \
		${TMP_DIR}/${CHR}_windows_input_F_CPM.tab
        awk -v c=${CPM2} 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${TMP_DIR}/${CHR}_windows_input_R.tab > \
		${TMP_DIR}/${CHR}_windows_input_R_CPM.tab
        paste ${TMP_DIR}/${CHR}_windows_F_CPM.tab ${TMP_DIR}/${CHR}_windows_input_F_CPM.tab | \
        	awk -v c=${CPM} 'BEGIN{OFS="\t"}{if (($4-$10) > 1/c) print $1, $2, $3, $4-$10, $5, $6; else print $1, $2, $3, 1/c, $5, $6}' - > \
		${TMP_DIR}/${CHR}_windows_F_CPM_minusinput.tab
	paste ${TMP_DIR}/${CHR}_windows_R_CPM.tab ${TMP_DIR}/${CHR}_windows_input_R_CPM.tab | \
        	awk -v c=${CPM} 'BEGIN{OFS="\t"}{if (($4-$10) > 1/c) print $1, $2, $3, $4-$10, $5, $6; else print $1, $2, $3, 1/c, $5, $6}' - > \
		${TMP_DIR}/${CHR}_windows_R_CPM_minusinput.tab
done
wait




echo "Calculating Partition scores"
## Calculate Partition, derivative and boundary scores
for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
        echo "doing for ${CHR}"
        partition_smooth.pl ${TMP_DIR}/${CHR}_windows_F_CPM_minusinput.tab \
          ${TMP_DIR}/${CHR}_windows_R_CPM_minusinput.tab \
          ${RADIUS} ${DRADIUS} ${ZRADIUS} > ${TMP_DIR}/${CHR}_windows_RFD_minusinput.txt &
done
wait

for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
        echo "doing for ${CHR}"
        partition_smooth.pl ${TMP_DIR}/${CHR}_windows_F_CPM.tab \
          ${TMP_DIR}/${CHR}_windows_R_CPM.tab \
          ${RADIUS} ${DRADIUS} ${ZRADIUS} > ${TMP_DIR}/${CHR}_windows_RFD.txt &
done
wait




## Gather output
OUT_INPUT=${OUT_RFD_DIR}/${PREFIX}_tmp1_smooth_results_w${WIN}_s${RADIUS}_d${DRADIUS}_z${ZRADIUS}.txt
if [ -f ${OUT_INPUT} ]; then
  rm ${OUT_INPUT}
fi

touch ${OUT_INPUT}
for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
        paste ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_F.tab \
          ${TMP_DIR}/${CHR}_windows_R.tab ${TMP_DIR}/${CHR}_windows_F_CPM_minusinput.tab \
          ${TMP_DIR}/${CHR}_windows_R_CPM_minusinput.tab ${TMP_DIR}/${CHR}_windows_RFD_minusinput.txt | \
          awk '$8>0 || $14>0' | cut -f -3,8,14,20,26,30- >> ${OUT_INPUT}
done

OUT_NORMAL=${OUT_RFD_DIR}/${PREFIX}_tmp2_smooth_results_w${WIN}_s${RADIUS}_d${DRADIUS}_z${ZRADIUS}.txt
if [ -f ${OUT_NORMAL} ]; then
  rm ${OUT_NORMAL}
fi

touch ${OUT_NORMAL}
for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
        paste ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_F.tab \
          ${TMP_DIR}/${CHR}_windows_R.tab ${TMP_DIR}/${CHR}_windows_F_CPM.tab \
          ${TMP_DIR}/${CHR}_windows_R_CPM.tab ${TMP_DIR}/${CHR}_windows_RFD.txt | \
          awk '$8>0 || $14>0' | cut -f -3,8,14,20,26,30- >> ${OUT_NORMAL}
done

OUT_FINAL=${OUT_RFD_DIR}/${PREFIX}_minusinput_smooth_results_w${WIN}_s${RADIUS}_d${DRADIUS}_z${ZRADIUS}.txt
if [ -f ${OUT_FINAL} ]; then
  rm ${OUT_FINAL}
fi

paste ${OUT_NORMAL} ${OUT_INPUT} | cut -f 1-12,18,19,21- > ${OUT_FINAL}

rm ${OUT_NORMAL} ${OUT_INPUT}

OUT_BASE=$(basename ${OUT_FINAL} .txt)
OUT_RFD_BW=${OUT_BW_DIR}/${OUT_BASE}_RFD.bw
OUT_CPM_F_BW=${OUT_BW_DIR}/${OUT_BASE}_CPM_F.bw
OUT_CPM_R_BW=${OUT_BW_DIR}/${OUT_BASE}_CPM_R.bw
OUT_RFD_BDG=${TMP_DIR}/${OUT_BASE}_RFD.bdg
OUT_CPM_F_BDG=${TMP_DIR}/${OUT_BASE}_CPM_F.bdg
OUT_CPM_R_BDG=${TMP_DIR}/${OUT_BASE}_CPM_R.bdg



awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$9}' ${OUT_FINAL} | \
  sort -k1,1 -k2,2n - > ${OUT_RFD_BDG}

bedGraphToBigWig ${OUT_RFD_BDG} ${CHROM_SIZES} ${OUT_RFD_BW}

awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$6}' ${OUT_FINAL} | \
  sort -k1,1 -k2,2n - > ${OUT_CPM_F_BDG}

bedGraphToBigWig ${OUT_CPM_F_BDG} ${CHROM_SIZES} ${OUT_CPM_F_BW}

awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$7}' ${OUT_FINAL} | \
  sort -k1,1 -k2,2n - > ${OUT_CPM_R_BDG}

bedGraphToBigWig ${OUT_CPM_R_BDG} ${CHROM_SIZES} ${OUT_CPM_R_BW}

rm ${OUT_F_BG} ${OUT_R_BG} ${INPUT_F_BG} ${INPUT_R_BG}
rm -rf ${TMP_DIR}









