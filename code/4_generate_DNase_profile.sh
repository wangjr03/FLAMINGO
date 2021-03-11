#!/bin/bash
cd $1
module load bedtools
for i in *;do
	echo $i
	bedtools map -a $i -b $2 -c 4 -o mean > DNase_$i
done

