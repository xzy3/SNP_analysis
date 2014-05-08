#!/bin/sh

#  stats_bamtools.sh
#  This script will get bamtools stats.

for i in *.bam; do

n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`

echo "***Getting stats for $n"
bamtools stats -in $i >> $n.stats2.txt

echo "fastq.gz file sizes:" >> $n.stats2.txt
ls -lh ../Zips/ | awk '{print $5}' | egrep -v '^$' >> $n.stats2.txt

echo "Unmapped fastq file sizes:" >> $n.stats2.txt
ls -lh ../unmappedReads/*.fastq | awk '{print $5}' | egrep -v '^$' >> $n.stats2.txt

echo "Unmapped contig count:" >> $n.stats2.txt
grep -c ">" ../unmappedReads/${n}_abyss-3.fa >> $n.stats2.txt

#Space
echo "" >> $n.stats2.txt

bamtools coverage -in $i | awk '{sum+=$3} END { print "Average coverage: ",sum/NR"X"}' >> $n.stats2.txt

bamtools coverage -in $i | awk '{if ($3 < 1) ++b } END {print "Reference with coverage:  "((FNR-b)/FNR)*100"%"}' >> $n.stats2.txt

cat $n.stats2.txt | grep -v "Failed" | grep -v "Duplicates" | grep -v "Proper-pairs" >> $n.stats.txt

rm $n.stats2.txt

#bamtools coverage -in $i | awk '{if ($3 <= 10) ++b } END {print "Reference has coverage at or above 10X:  "((FNR-b)/FNR)*100"%"}' >> $n.stats.txt

#bamtools coverage -in $i | awk '{if ($3 <= 50) ++b } END {print "Reference has coverage at or above 50X:  "((FNR-b)/FNR)*100"%"}' >> $n.stats.txt

#bamtools coverage -in $i | awk '{if ($3 <= 100) ++b } END {print "Reference has coverage at or above 100X:  "((FNR-b)/FNR)*100"%"}' >> $n.stats.txt

done

#
#  Created by Stuber, Tod P - APHIS on 11/07/12.
#

