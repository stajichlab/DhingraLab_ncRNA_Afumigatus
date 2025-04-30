#!/usr/bin/bash -l
#SBATCH -p epyc -c 24 -N 1 -n 1 --mem 64gb --out logs/make_bam.log

CPU=8
module load samtools
FOLDER=results
for STRAIN in A1163 Af293
do
	#parallel -j 6 samtools view --threads $CPU -O BAM -o {.}.bam {} ::: $(find $FOLDER/STAR_$STRAIN -name "*.sam")
	echo "$FOLDER/STAR_$STRAIN"
	mkdir -p $FOLDER/STAR_${STRAIN}_sort
	for file in $(find $FOLDER/STAR_$STRAIN -name "*.bam")
	do
		OUT=$(basename $file .Aligned.out.bam)
		if [ ! -f $FOLDER/STAR_${STRAIN}_sort/$OUT.bam.csi ]; then
			samtools sort -T $SCRATCH/$OUT --threads $CPU --write-index -o $FOLDER/STAR_${STRAIN}_sort/$OUT.bam $file
		fi
	done
done

