#!/bin/sh

#  bovis_processZips.sh
#  Working directory should contain paired-end fastq reads
#  Sample must be labeled with TB/Bruc number.
#  TB numbers must be labeled ##-####
#  Bruc numbers must be labeled B##-####
#  Reads must be included as _R1 and _R2
#  See loopfiles.sh and email_loopfiles for multiple samples.

#################################################################################
#  Dependencies ---
#   Tested and ran on OS X V10.9
#   -Xcode developer tools
#   bwa, http://bio-bwa.sourceforge.net/bwa.shtml
#   samtools, http://samtools.sourceforge.net/samtools.shtml
#   picard, http://picard.sourceforge.net/command-line-overview.shtml
#   gatk, http://www.broadinstitute.org/gatk/
#   bamtools
#   File containing high quality SNPs, Volumes/Mycobacterium/Go_To_File/HighestQualitySNPs.vcf
#   Reference in fasta format, /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta
#################################################################################
date
echo "**************************************************************************"
echo "**************************** START ${PWD##*/} ****************************"
echo "**************************************************************************"

echo "current directory"
pwd

# Move zip files to their own directory
mkdir ./Zips
mv *.fastq* ./Zips
mkdir ./BWAmem-GATK
cd BWAmem-GATK/


# Make alias links in BWA-GATK directory to zip files
ls ../Zips/*.fastq* | while read file; do ln -s $file; done

if [ $1 == ab1 ]; then
    cp /Volumes/Data_HD/Brucella/processZips_dependencies/Abortus1/script_dependents/NC_00693c.fasta ./
    hqs="/Volumes/Data_HD/Brucella/processZips_dependencies/Abortus1/script_dependents/NC_00693cHighestQualitySNPs.vcf"
    bioinfo="/Volumes/TStuber/Results/_Brucella/Abortus1"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/Abortus1/coverageFiles"
    copyto="/Volumes/Data_HD/Brucella/processZips_dependencies/Abortus1/newFiles"

elif [ $1 == mel ]; then
    cp /Volumes/Data_HD/Brucella/processZips_dependencies/Melitensis/script_dependents/BmelitensisM5-90.fasta ./
    hqs="/Volumes/Data_HD/Brucella/processZips_dependencies/Melitensis/script_dependents/melHighestQualitySNPs.vcf"
    bioinfo="/Volumes/TStuber/Results/_Brucella/Melitensis"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    copyto="/Volumes/Data_HD/Brucella/processZips_dependencies/Melitensis/newFiles"

elif [ $1 == suis1 ]; then
    cp /Volumes/Data_HD/Brucella/processZips_dependencies/Suis1/script_dependents/NC_01725c.fasta ./
    hqs="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis1/script_dependents/NC_01725cHighestQualitySNPs.vcf"
    bioinfo="/Volumes/TStuber/Results/_Brucella/Suis1"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    copyto="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis1/newFiles"

elif [ $1 == suis2 ]; then
    cp /Volumes/Data_HD/Brucella/processZips_dependencies/Suis2/script_dependents/Bsuisbv2-94-11.fasta ./
    hqs="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis2/script_dependents/suis2HighestQualitySNPs.vcf"
    bioinfo="/Volumes/TStuber/Results/_Brucella/Suis2"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    copyto="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis2/newFiles"

elif [ $1 == suis3 ]; then
    cp /Volumes/Data_HD/Brucella/processZips_dependencies/Suis3/script_dependents/B-REF-BS3-686.fasta ./
    hqs="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis3/script_dependents/suis3HighestQualitySNPs.vcf"
    bioinfo="/Volumes/TStuber/Results/_Brucella/Suis3"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    copyto="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis3/newFiles"

elif [ $1 == suis4 ]; then
    cp /Volumes/Data_HD/Brucella/processZips_dependencies/Suis4/script_dependents/B-REF-BS4-40.fasta ./
    hqs="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis4/script_dependents/suis4HighestQualitySNPs.vcf"
    bioinfo="/Volumes/TStuber/Results/_Brucella/Suis4"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    copyto="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis4/newFiles"

elif [ $1 == canis ]; then
    cp /Volumes/Data_HD/Brucella/processZips_dependencies/Canis/script_dependents/BcanisATCC23365.fasta ./
    hqs="/Volumes/Data_HD/Brucella/processZips_dependencies/Canis/script_dependents/canisHighestQualitySNPs.vcf"
    bioinfo="/Volumes/TStuber/Results/_Brucella/Canis"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    copyto="/Volumes/Data_HD/Brucella/processZips_dependencies/Canis/newFiles"

elif [ $1 == ceti1 ]; then
    cp /Volumes/Data_HD/Brucella/processZips_dependencies/Ceti1/script_dependents/Bceti1Cudo.fasta ./
    hqs="/Volumes/Data_HD/Brucella/processZips_dependencies/Ceti1/script_dependents/ceti1HighestQualitySNPs.vcf"
    bioinfo="/Volumes/TStuber/Results/_Brucella/Ceti1"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    copyto="/Volumes/Data_HD/Brucella/processZips_dependencies/Ceti1/newFiles"

elif [ $1 == ceti2 ]; then
    cp /Volumes/Data_HD/Brucella/processZips_dependencies/Ceti2/script_dependents/Bceti2-TE10759.fasta ./
    hqs="/Volumes/Data_HD/Brucella/processZips_dependencies/Ceti2/script_dependents/ceti2HighestQualitySNPs.vcf"
    bioinfo="/Volumes/TStuber/Results/_Brucella/Ceti2"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    copyto="/Volumes/Data_HD/Brucella/processZips_dependencies/Ceti2/newFiles"

elif [ $1 == ovis ]; then
    cp /Volumes/Data_HD/Brucella/processZips_dependencies/Ovis/script_dependents/BovisATCC25840.fasta ./
    hqs="/Volumes/Data_HD/Brucella/processZips_dependencies/Ovis/script_dependents/BovisATCC25840HighestQualitySNPs.vcf"
    bioinfo="/Volumes/TStuber/Results/_Brucella/Ovis"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    copyto="/Volumes/Data_HD/Brucella/processZips_dependencies/Ovis/newFiles"

elif [ $1 == bovis ]; then
    cp /Volumes/Data_HD/Mycobacterium/script_dependents/NC_002945.fasta ./
    hqs="/Volumes/Data_HD/Mycobacterium/script_dependents/HighestQualitySNPs.vcf"
    bioinfo="/Volumes/TStuber/Results/_Mycobacterium"
    coverageFiles="/Volumes/Data_HD/Mycobacterium/script_dependents/coverageFiles-chrom"
    copyto="/Volumes/Data_HD/Mycobacterium/Analysis_new/2014-05-07/copytoTest"

    ###################################################################

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    spoligoSpacerFinder.sh &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################

else
    echo "Incorrect argument!  Must use one of the following arguments: ab1, mel, suis1, suis2, suis3, suis4, canis, ceti1, ceti2, ovis, bovis"
    exit 1
fi

# Grab reads and reference and place them in variables
ref=`ls | grep .fasta`
echo "Reference Input:  $ref"

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"

revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

#   retrieves reference name and name from sorted BAM file name
r=`echo $ref | sed 's/\..*//'`
n=`echo $revReads | sed 's/_.*//' | sed 's/\..*//'`

echo "***Reference naming convention:  $r"
echo "***Isolate naming convention:  $n"

samtools faidx $ref
java -Xmx4g -jar /Users/Shared/_programs/picard-tools-1.100/CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${r}.dict

if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
    echo "Index and dict are present, continue script"
    else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    samtools faidx $ref
    java -Xmx4g -jar /Users/Shared/_programs/picard-tools-1.100/CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${r}.dict
        if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
        read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
        fi
fi

# See echo comments
echo "***bwa index $r"
bwa index $ref
#echo "***bwa aln forward $n"
# -t sets the number of threads/cores
#please use option `-M' to flag extra hits as secondary.
#/Users/Shared/_programs/bwa-0.7.5a/bwa mem -M -t2 $ref $forReads > $n.sr1.sai
#echo "***bwa aln reverse $n"
#/Users/Shared/_programs/bwa-0.7.5a/bwa mem -M -t2 $ref $revReads > $n.sr2.sai
# -r STR	 Specify the read group in a format like ‘@RG\tID:foo\tSM:bar’ Needed for GATK
#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow more mismatch per read.
echo "***Making Sam file"
bwa mem -M -t 4 -R @RG"\t"ID:"$n""\t"PL:ILLUMINA"\t"PU:"$n"_RG1_UNIT1"\t"LB:"$n"_LIB1"\t"SM:"$n" $ref $forReads $revReads > $n.sam


# -b	 Output in the BAM format.
# -h	 Include the header in the output.
#-F INT	 Skip alignments with bits present in INT [0]
echo "***Making Bam file"
samtools view -bh -F4 -T $ref $n.sam > $n.raw.bam

####### unmapped reads #######
#Bam with mapped and unmapped reads
samtools view -bh -T $ref $n.sam > $n.all.bam
#Strip off the unmapped reads
samtools view -h -f4 $n.all.bam > $n.unmappedReads.sam
#Create fastqs of unmapped reads to assemble
java -Xmx2g -jar  /Users/Shared/_programs/picard-tools-1.100/SamToFastq.jar INPUT=$n.unmappedReads.sam FASTQ=${n}-unmapped_R1.fastq SECOND_END_FASTQ=${n}-unmapped_R2.fastq
rm $n.all.bam
rm $n.unmappedReads.sam
abyss-pe name=${n}_abyss k=64 in="${n}-unmapped_R1.fastq ${n}-unmapped_R1.fastq"

mkdir ../unmappedReads
mv ${n}-unmapped_R1.fastq ../unmappedReads
mv ${n}-unmapped_R2.fastq ../unmappedReads
mv ${n}_abyss-3.fa ../unmappedReads
mv ${n}_abyss-7.fa ../unmappedReads
mv ${n}_abyss-stats ../unmappedReads
mv *coverage* ../unmappedReads
rm *abyss*
######################

echo "***Sorting Bam"
samtools sort $n.raw.bam $n.sorted
echo "***Indexing Bam"
samtools index $n.sorted.bam
# Remove duplicate molecules

echo "***Marking Duplicates"
java -Xmx2g -jar  /Users/Shared/_programs/picard-tools-1.100/MarkDuplicates.jar INPUT=$n.sorted.bam OUTPUT=$n.dup.bam METRICS_FILE=$n.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index $n.dup.bam"
samtools index $n.dup.bam

# Creates file that is used in the next step
# locally realign reads such that the number of mismatching bases is minimized across all the reads
# http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_indels_RealignerTargetCreator.html
echo "***Realigner Target Creator"
java -Xmx2g -jar /Users/Shared/_programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T RealignerTargetCreator -I $n.dup.bam -R $ref -o $n.forIndelRealigner.intervals

# Uses the RealignerTargetCreator output file to improve BAM alignment
# http://www.broadinstitute.org/gatk/guide/tagged?tag=indelrealigner
echo "***Target Intervals"
java -Xmx2g -jar /Users/Shared/_programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T IndelRealigner -I $n.dup.bam -R $ref -targetIntervals $n.forIndelRealigner.intervals -o $n.realignedBam.bam

# Uses a .vcf file which contains SNP calls of known high value to recalibrates base quality scores
# http://www.broadinstitute.org/gatk/guide/tagged?tag=baserecalibrator
echo "***Base Recalibrator"
java -Xmx2g -jar /Users/Shared/_programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T BaseRecalibrator -I $n.realignedBam.bam -R $ref -knownSites ${hqs} -o $n.recal_data.grp

# Make the finished "ready" .bam file
echo "***Print Reads"
java -Xmx2g -jar /Users/Shared/_programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T PrintReads -R $ref -I $n.realignedBam.bam -BQSR $n.recal_data.grp -o $n.ready-mem.bam

# SNP calling and .vcf making
# Threads used can be changed
# http://www.broadinstitute.org/gatk/guide/tagged?tag=unifiedgenotyper
#echo "***UnifiedGenotyper, aka calling SNPs"
#java -Xmx2g -jar /Users/Shared/_programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper -I $n.ready-mem.bam -o $n.ready-mem.vcf -nt 8

# Add zero positions to vcf
java -Xmx2g -jar /Users/Shared/_programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper -out_mode EMIT_ALL_SITES -I ${n}.ready-mem.bam -o ${n}.allsites.vcf -nt 8
awk ' $0 ~ /#/ || $8 !~ /^AN=2;/ {print $0}' ${n}.allsites.vcf > $n.ready-mem.vcf

# SNP calling and .vcf making
# Threads used can be changed
# http://www.broadinstitute.org/gatk/guide/tagged?tag=unifiedgenotyper
echo "***HaplotypeCaller, aka calling SNPs"
java -Xmx2g -jar /Users/Shared/_programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller -I $n.ready-mem.bam -o $n.hapreadyAll.vcf -nct 8

echo "******Awk VCF leaving just SNPs******"
awk '/#/ || $4 ~ /^[ATGC]$/ && $5 ~ /^[ATGC]$/ {print $0}' $n.hapreadyAll.vcf > $n.hapreadyOnlySNPs.vcf

java -Xmx2g -jar /Users/Shared/_programs/IGVTools/igvtools.jar index $n.hapreadyOnlySNPs.vcf

echo "***Deleting Files"
rm $n.sr1.sai
rm $n.sr2.sai
rm $n.sam
rm $n.raw.bam
rm $n.dup.bam
rm $n.dup.bam.bai
rm $n.sorted.bam
rm $n.sorted.bam.bai
rm $n.realignedBam.bam
rm $n.realignedBam.bai
rm $forReads
rm $revReads
rm igv.log

###################################
# The next 6 steps collect metrics
###################################

#Collect Depth of coverage info
echo "***Collect Depth of Coverage"
java -jar /Users/Shared/_programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T DepthOfCoverage -R $ref -I $n.ready-mem.bam --omitDepthOutputAtEachBase > $n.DepthofCoverage.xls

#Quality Score Distribution
echo "***Quality Score Distribution"
java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.100/QualityScoreDistribution.jar REFERENCE_SEQUENCE=$ref INPUT=$n.ready-mem.bam CHART_OUTPUT=$n.QualityScorceDistribution.pdf OUTPUT=$n.QualityScoreDistribution ASSUME_SORTED=true

#Mean Quality by Cycle
echo "***Mean Quality by Cycle"
java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.100/CollectMultipleMetrics.jar REFERENCE_SEQUENCE=$ref INPUT=$n.ready-mem.bam OUTPUT=$n.Quality_by_cycle PROGRAM=MeanQualityByCycle ASSUME_SORTED=true

#Collect Alignment Summary Metrics
echo "***Collect Alignment Summary Metrics"
java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.100/CollectAlignmentSummaryMetrics.jar REFERENCE_SEQUENCE=$ref INPUT=$n.ready-mem.bam OUTPUT=$n.AlignmentMetrics ASSUME_SORTED=true

#Collect GC Bias Error
echo "***Collect GC Bias Error"
java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.100/CollectGcBiasMetrics.jar REFERENCE_SEQUENCE=$ref INPUT=$n.ready-mem.bam OUTPUT=$n.CollectGcBiasMetrics CHART_OUTPUT=$n.GC.PDF ASSUME_SORTED=true

#Collect Insert Size Metrics
echo "***Collect Insert Size Metrics"
java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.100/CollectInsertSizeMetrics.jar REFERENCE_SEQUENCE=$ref INPUT=$n.ready-mem.bam HISTOGRAM_FILE=$n.InsertSize.pdf OUTPUT=$n.CollectInsertSizeMetrics ASSUME_SORTED=true

# Make coverage file
#bamtools coverage -in *.ready-mem.bam > ${n}-coverage
#cat ${n}-coverage > ${n}-coveragetemp-orginal
#chroms=`awk '{a[$1]++} END{for (var in a) print var}' ${n}-coverage | sort -n`
#num=1
#echo "chroms are: $chroms"
#for i in $chroms; do
#    sed "s/$i/chrom${num}/g" ${n}-coverage > temp.vcf
#    mv temp.vcf ${n}-coverage
#    echo "$i was marked as chrom${num}"
#    num=$(( $num + 1 ))
#done

#mv ${n}-coverage ${coverageFiles}/${n}-coverage

###

cat $n.DepthofCoverage.xls >> $n.Metrics_summary.xls
cat $n.AlignmentMetrics >> $n.Metrics_summary.xls
cat $n.CollectInsertSizeMetrics >> $n.Metrics_summary.xls

echo "***Organizing files"

#Move to QualityValues subfolder
mkdir QualityValues
mv $n.GC.PDF QualityValues/
mv $n.QualityScorceDistribution.pdf QualityValues/
mv $n.InsertSize.pdf QualityValues/
mv $n.Quality_by_cycle.quality_by_cycle.pdf QualityValues/

#Remove files
rm $n.DepthofCoverage.xls
rm $n.CollectInsertSizeMetrics
rm $n.forIndelRealigner.intervals
rm $n.recal_data.grp
rm $n.FilteredReads.xls
rm $n.Quality_by_cycle.quality_distribution_metrics
rm $n.Quality_by_cycle.quality_by_cycle_metrics
rm $n.Quality_by_cycle.alignment_summary_metrics
rm $n.Quality_by_cycle.insert_size_histogram.pdf
rm $n.Quality_by_cycle.quality_distribution.pdf
rm $n.CollectGcBiasMetrics
rm $n.QualityScoreDistribution

###########################
echo "***Getting stats for $n"

bamtools stats -in $n.ready-mem.bam > $n.stats2.txt

echo "fastq.gz file sizes:" >> $n.stats2.txt
ls -lh ../Zips/ | awk '{print $5}' | egrep -v '^$' >> $n.stats2.txt

echo "Unmapped fastq file sizes:" >> $n.stats2.txt
ls -lh ../unmappedReads/*.fastq | awk '{print $5}' | egrep -v '^$' >> $n.stats2.txt

echo "Unmapped contig count:" >> $n.stats2.txt
grep -c ">" ../unmappedReads/${n}_abyss-3.fa >> $n.stats2.txt
echo "" >> $n.stats2.txt

# Split output to multiple processes using tee.
#bamtools coverage -in $n.ready-mem.bam | tee >((awk '{sum+=$3} END { print "Average coverage: ",sum/NR"X"}' >> $n.stats2.txt)) | awk '{if ($3 < 1) ++b } END {print "Reference with coverage:  "((FNR-b)/FNR)*100 "%"}' >> $n.stats2.txt

bamtools coverage -in $n.ready-mem.bam | awk '{sum+=$3} END { print "Average coverage: ",sum/NR"X"}' >> $n.stats2.txt
bamtools coverage -in $n.ready-mem.bam | awk '{if ($3 < 1) ++b } END {print "Reference with coverage:  "((FNR-b)/FNR)*100 "%"}' >> $n.stats2.txt

cat $n.stats2.txt | grep -v "Failed" | grep -v "Duplicates" | grep -v "Proper-pairs" >> $n.stats.txt
rm $n.stats2.txt
echo "" >> $n.stats.txt
###########################

#  Add Insert_Size and Read_Length to stats.txt file
echo 'Mean_Insert_Size  Standard_Deviation:' >> $n.stats.txt
awk 'BEGIN {OFS="\t"} { print $5,$6 }' $n.Quality_by_cycle.insert_size_metrics | awk 'FNR == 8 {print $0}' >> $n.stats.txt

echo 'Mean_Read_Length:' >> $n.stats.txt
awk 'BEGIN {OFS="\t"} { print $16 }' $n.AlignmentMetrics | awk 'FNR == 10 {print $0}' >> $n.stats.txt

echo "" >> $n.stats.txt

#  Add SNP call numbers to stats.txt file
echo "All SNPs called in ready-mem.vcf:" >> $n.stats.txt
egrep -v "#" $n.ready-mem.vcf | grep -c ".*" >> $n.stats.txt

echo "SNPs of AC2 and QUAL > 150:" >> $n.stats.txt
egrep -v "#" $n.ready-mem.vcf | egrep "AC=2" | awk '$6 > 150' | grep -c ".*" >> $n.stats.txt

#  Show Mean Coverage at Terminal and coverageReport
echo "Mean Coverage"
awk -v number="$n" 'BEGIN {OFS="\t"} $0 ~ number { print $1,$2,$3,$7 }' $n.Metrics_summary.xls | awk 'FNR == 2 {print $0}'

awk -v number="$n" 'BEGIN {OFS="\t"} $0 ~ number { print $1,$2,$3,$7 }' $n.Metrics_summary.xls | awk 'FNR == 2 {print $0}' >> /Users/Shared/_WGS/coverageReport.txt

echo "Sample identified and ran as:  $1" >> /Users/Shared/_WGS/dailyReport.txt

awk -v number="$n" 'BEGIN {OFS="\t"} $0 ~ number { print $1,$2,$3,$7 }' $n.Metrics_summary.xls | awk 'FNR == 2 {print $0}' >> /Users/Shared/_WGS/dailyReport.txt

mv $n.Metrics_summary.xls QualityValues/
mv $n.stats.txt QualityValues/
rm $n.Quality_by_cycle.insert_size_metrics
rm $n.AlignmentMetrics
cat ../*out1* ../*out2* > ../${n}-identification.txt
rm ../*identifier_out1*
rm ../*identifier_out2*
rm -r temp

#  Send files to the Network
echo "***Sending files to the Network"
mkdir ${bioinfo}/newFiles/$n/
mkdir ${bioinfo}/newFiles/$n/BWAmem-GATK
mkdir ${bioinfo}/newFiles/$n/Zips
mkdir ${bioinfo}/newFiles/$n/BWAmem-GATK/QualityValues
mkdir ${bioinfo}/newFiles/$n/unmappedReads

cp ../${n}-identification.txt ${bioinfo}/newFiles/$n
cp $n.* ${bioinfo}/newFiles/$n/BWAmem-GATK
cp $ref ${bioinfo}/newFiles/$n/BWAmem-GATK
cp $ref.fai ${bioinfo}/newFiles/$n/BWAmem-GATK
cp *-coverage ${bioinfo}/newFiles/$n/BWAmem-GATK
cp -r ../Zips/*.gz ${bioinfo}/newFiles/$n/Zips
cp -r ./QualityValues/* ${bioinfo}/newFiles/$n/BWAmem-GATK/QualityValues
cp -r ../unmappedReads/* ${bioinfo}/newFiles/$n/unmappedReads
cp -r ../../${n} $copyto
pwd
#Make dailyStats.txt for each stats.txt made for each isolate.
echo "" >> /Users/Shared/_WGS/dailyStats.txt
echo "" >> /Users/Shared/_WGS/dailyStats.txt
echo "          <--------- $n --------->" >> /Users/Shared/_WGS/dailyStats.txt
echo "*** IDENTIFIED AS $1 ***" >> /Users/Shared/_WGS/dailyStats.txt
cat QualityValues/$n.stats.txt >> /Users/Shared/_WGS/dailyStats.txt
cp QualityValues/$n.stats.txt /Volumes/Data_HD/nvsl/DBL/MB/Sequencing/WGS/_stats
cp QualityValues/$n.stats.txt /Users/Shared/_WGS/stats

#open -a textwrangler QualityValues/$n.stats.txt
echo "**************************** END $n ****************************"
date

#
#  Created by Stuber, Tod P - APHIS on 11/08/12.
#
