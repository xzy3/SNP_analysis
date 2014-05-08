#!/bin/sh

: <<'END'
This script is the second script in a two script workflow.  Script 2 genotypes Mycobacterium tuberculosis complex and Brucella species from SNP data contained in VCFs.  It operates on VCFs generated with the same reference output from script 1.  VCFs are collected into a single working directory.  Comparisons are output as SNP tables and alignment files to view as trees in your program of choice.

Script 2 will run and output tables and alignments from just the data generated from script 1, however the data will be more informative if additional files are provide.  Those files are:
1) A file that contains positions to cluster individual isolates/VCFs into groups, subgroups and clades.
2) Files that contain positions to remove from the analysis.

Paradigm
1) Once a SNP occurs and establishes in a population it does not revert back
2) Observations of homoplasy are rare
3) Group, subgroup and clade clusters only show parsimony informative SNPs for the isolates within that cluster
4) SNPs observed in a single isolate are less informative than SNPs seen in multiple isolates and therefore established in a population

END

#################################################################################
#  Dependencies ---
#   Tested and ran on OS X V10.9
#   -Xcode developer tools
#   -ClustalW
#################################################################################

#################################################################################
# Argument controls:

if [[ $1 == ab1 ]]; then

    genotypingcodes="/Volumes/Data_HD/Brucella/processZips_dependencies/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/Volumes/Data_HD/Brucella/processZips_dependencies/Abortus1/script_dependents/Abortus1_Defining_SNPs.txt"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/Abortus1/coverageFiles"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/Volumes/Data_HD/Brucella/processZips_dependencies/Abortus1/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/Volumes/TStuber/Results/_Brucella/Abortus1/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/Volumes/TStuber/Results/_Brucella/Abortus1/_VCFs"
    echo "vcftofasta.sh ran as Brucella abortus bv 1, 2 or 4"
    echo "Script vcftofasta.sh ran using Brucella abortus bv 1, 2 or 4 variables" > log

elif [[ $1 == mel ]]; then

    genotypingcodes="/Volumes/Data_HD/Brucella/processZips_dependencies/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/Volumes/Data_HD/Brucella/processZips_dependencies/Melitensis/script_dependents/Mel_Defining_SNPs.txt"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/Volumes/Data_HD/Brucella/processZips_dependencies/Melitensis/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/Volumes/TStuber/Results/_Brucella/Melitensis/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/Volumes/TStuber/Results/_Brucella/Melitensis/_VCFs"
    echo "vcftofasta.sh ran as B. melitensis"
    echo "Script vcftofasta.sh ran using B. melitensis variables" > log

elif [[ $1 == suis1 ]]; then

    genotypingcodes="/Volumes/Data_HD/Brucella/processZips_dependencies/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis1/script_dependents/Suis1_Defining_SNPs.txt"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis1/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/Volumes/TStuber/Results/_Brucella/Suis1/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/Volumes/TStuber/Results/_Brucella/Suis1/_VCFs"
    echo "vcftofasta.sh ran as B. suis bv1"
    echo "Script vcftofasta.sh ran using B. suis bv1 variables" > log

elif [[ $1 == suis2 ]]; then

    genotypingcodes="/Volumes/Data_HD/Brucella/processZips_dependencies/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis2/script_dependents/suis2_Defining_SNPs.txt"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis2/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/Volumes/TStuber/Results/_Brucella/Suis2/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF=""
    echo "vcftofasta.sh ran as B. suis bv2"
    echo "Script vcftofasta.sh ran using B. suis bv2 variables" > log

elif [[ $1 == suis3 ]]; then

    genotypingcodes="/Volumes/Data_HD/Brucella/processZips_dependencies/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis3/script_dependents/Suis3_Defining_SNPs.txt"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis3/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/Volumes/TStuber/Results/_Brucella/Suis3/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/Volumes/TStuber/Results/_Brucella/Suis3/_VCFs"
    echo "vcftofasta.sh ran as B. suis bv3"
    echo "Script vcftofasta.sh ran using B. suis bv3 variables" > log

elif [[ $1 == suis4 ]]; then

    genotypingcodes="/Volumes/Data_HD/Brucella/processZips_dependencies/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis4/script_dependents/Suis4_Defining_SNPs.txt"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/Volumes/Data_HD/Brucella/processZips_dependencies/Suis4/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/Volumes/TStuber/Results/_Brucella/Suis4/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/Volumes/TStuber/Results/_Brucella/Suis4/_VCF"
    echo "vcftofasta.sh ran as B. suis bv4"
    echo "Script vcftofasta.sh ran using B. suis bv4 variables" > log

elif [[ $1 == canis ]]; then

    genotypingcodes="/Volumes/Data_HD/Brucella/processZips_dependencies/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/Volumes/Data_HD/Brucella/processZips_dependencies/Canis/script_dependents/Canis_Defining_SNPs.txt"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/Volumes/Data_HD/Brucella/processZips_dependencies/Canis/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/Volumes/TStuber/Results/_Brucella/Canis/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/Volumes/TStuber/Results/_Brucella/Canis/_VCF"
    echo "vcftofasta.sh ran as B. canis"
    echo "Script vcftofasta.sh ran using B. canis variables" > log


elif [[ $1 == ceti1 ]]; then

    genotypingcodes="/Volumes/Data_HD/Brucella/processZips_dependencies/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/Volumes/Data_HD/Brucella/processZips_dependencies/Ceti1/script_dependents/Ceti1_Defining_SNPs.txt"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/Volumes/Data_HD/Brucella/processZips_dependencies/Ceti1/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/Volumes/TStuber/Results/_Brucella/Ceti1/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/Volumes/TStuber/Results/_Brucella/Ceti1/_VCF"
    echo "vcftofasta.sh ran as B ceti group 1"
    echo "Script vcftofasta.sh ran using B ceti group 1 variables" > log


elif [[ $1 == ceti2 ]]; then

    genotypingcodes="/Volumes/Data_HD/Brucella/processZips_dependencies/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/Volumes/Data_HD/Brucella/processZips_dependencies/Ceti2/script_dependents/Ceti2_Defining_SNPs.txt"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/Volumes/Data_HD/Brucella/processZips_dependencies/Ceti2/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/Volumes/TStuber/Results/_Brucella/Ceti2/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/Volumes/TStuber/Results/_Brucella/Ceti2/_VCF"
    echo "vcftofasta.sh ran as B ceti group 2"
    echo "Script vcftofasta.sh ran using B ceti group 2 variables" > log


elif [[ $1 == ovis ]]; then

    genotypingcodes="/Volumes/Data_HD/Brucella/processZips_dependencies/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/Volumes/Data_HD/Brucella/processZips_dependencies/Ovis/script_dependents/Ovis_Defining_SNPs.txt"
    coverageFiles="/Volumes/Data_HD/Brucella/processZips_dependencies/coverageFiles"
    FilterAllVCFs=no #(yes or no), Do you want to filter all VCFs?
    FilterGroups=no #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/Volumes/Data_HD/Brucella/processZips_dependencies/Ovis/script_dependents/FilterFiles" #Files containing positions to filter
    RemoveFromAnalysis="/Volumes/TStuber/Results/_Brucella/Ovis/script_dependents/RemoveFromAnalysis.txt"
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N
    bioinfoVCF="/Volumes/TStuber/Results/_Brucella/Ovis/_VCFs"
    echo "vcftofasta.sh ran as B. ovis"
    echo "Script vcftofasta.sh ran using B. ovis variables" > log

elif [[ $1 == bovis ]]; then

    genotypingcodes="/Users/tstuberadmin/Desktop/Untitled.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs="/Volumes/Data_HD/Mycobacterium/DefiningSNPsGroupDesignations.txt"
    coverageFiles="/Volumes/Data_HD/Mycobacterium/script_dependents/coverageFiles-chrom"
    FilterAllVCFs=yes #(yes or no), Do you want to filter all VCFs?
    FilterGroups=yes #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory="/Volumes/Data_HD/Mycobacterium/script_dependents/bovisGroups" #Files containing positions to filter
    RemoveFromAnalysis="/Volumes/TStuber/Results/_Mycobacterium/_TB-VCF/RemoveFromAnalysis.txt"
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N
    bioinfoVCF="/Volumes/TStuber/Results/_Mycobacterium/_TB-VCF"
    echo "vcftofasta.sh ran as M. bovis"
    echo "Script vcftofasta.sh ran using M. bovis variables" > log

else

    echo ""
    echo "Incorrect argument!  Must use one of the following arguments: ab1, mel, suis1, suis2, suis3, suis4, canis, ceti1, ceti2, ovis, bovis"
    echo "For example, type ~$ vcftofasta.sh bovis"
    echo ""
    exit 1

fi

#################################################################################
# Set variables:

# Sed searches put into variables
tbNumberV='s/_.*//' #Remove all charaters at and beyond "_"
tbNumberW='s/\..*//' #Remove all charaters at and beyond "."
tbNumberOnly='s/.*\([0-9]\{2\}-[0-9,FM]\{4,6\}\).*/\1/' #Only tb Number, *laboratory specific*
dropEXT='s/\(.*\)\..*/\1/' #Just drop the extention from the file

NR_CPUS=18 # Computer cores to use when analyzing

Ncov=1 # Coverage below this value will be changed to -

fulDir=$PWD # Current working directory, do not change.

#################################################################################

# Count the number of chromosomes used in the reference when VCFs were made.
singleFile=`ls *.vcf | head -1`
chromCount=`grep -v "#" $singleFile | awk '{print $1}' | sort | uniq -d | awk 'END {print NR}'`
chroms=`grep -v "#" $singleFile | awk '{print $1}' | sort -n | uniq -d`
echo "The number of chromosomes seen in VCF: $chromCount"

#################################################################################

# Remove selected files from comparison
# Use file:  /Volumes/TStuber/Results/_Mycobacterium/_TB-VCF/RemoveFromAnalysis.txt

function removeIsolates () {

cat ${RemoveFromAnalysis} | tr '\r' '\n' > /Volumes/TStuber/Results/_Mycobacterium/_TB-VCF/RemoveFromAnalysisUnixReady.txt

removeList=`cat /Volumes/TStuber/Results/_Mycobacterium/_TB-VCF/RemoveFromAnalysisUnixReady.txt`

for i in $removeList; do
    rm *${i}*
done

rm /Volumes/TStuber/Results/_Mycobacterium/_TB-VCF/RemoveFromAnalysisUnixReady.txt

}

#################################################################################

# If there are 2 vcf files with the same name one of the files might unknowingly
# get cut out of the analysis and keep the undesired vcf instead.  This will
# alert if 2 vcf with the same TB number are present.
# The regular expression used in sed should be changed based on vcf naming convention

function testDuplicates () {

echo "Checking for duplicate VCFs."

for i in *; do
    getbase=`basename "$i"`
    number=`echo $getbase | sed $tbNumberV | sed $tbNumberW`
    echo $number >> list
done
duplist=`sort list | uniq -d`
rm list
dupNumberSize=`echo $duplist | wc | awk '{print $3}'`
if [ $dupNumberSize -gt 4 ]
then
    echo "These are duplicated VCFs."
    echo "Must remove duplication, and restart script."
    echo $duplist
    exit 1 # Error status
else
    echo "Good! No duplicate VCFs present"
fi
}

#################################################################################

# Looks for defining positions in VCF files.
# If an AC=1 is found at a defined position it is flagged as a posible mixed infection.
# These defining positions must be SNPs found cluster's main branch

function AConeCallPosition () {

positionList=`awk ' { print $2 }' "${DefiningSNPs}" | awk ' NF > 0 '`

echo "AConeCallPosition.sh is running"
echo "*********************************************************************" >> log
echo "Possible Mixed Isolates" >> log
echo "Defining SNPs that are called as AC=1" >> log
echo "" >> log
for i in *.vcf; do
(echo "Finding possible AC1 matches for $i"; for pos in $positionList; do awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ "^"x"$" ) print FILENAME, "Pos:", $2, "QUAL:", $6, $8 }' $i; done | grep "AC=1;" | awk 'BEGIN {FS=";"} {print $1, $2}' >> log) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 2

echo "*********************************************************************" >> log
}

#################################################################################

# This function prepares the filter files.
# awk needs to see a number in the file, so if the file is blank 2 matching numbers are added.  2 numbers because duplicates are kept therefore allowing a number to be pasting into awk when comparing files.

function filterFilespreparation () {

python -u /Users/Shared/_programs/_my_scripts/python_scripts/inputXLS.py | sed 's/ u//g' | tr "," "\t" | sed "s/\'//g" | sed 's/\[//g' |sed 's/\]//g' |sed 's/ //g' | sed 's/^u//g' | sed 's/\.0//g' > ~/Desktop/filterFile.txt

echo "Waiting for filter file creation to complete"
#FilterFileCreations.sh
wait
curdr=`pwd`

cd "${FilterDirectory}"

echo "Preparing Filter Files"
for i in *.txt; do
    (getbase=`basename "$i"`
    number=`echo $getbase | sed 's/\(.*\)\..*/\1/'`
    echo $number
    cat $i | sort | uniq > "$number".num
    if [ $((chromCount)) -eq 1 ]; then
        echo "100000000" >> "$number.num"
        echo "100000000" >> "$number.num"
    elif [ $((chromCount)) -eq 2 ]; then
        echo "chrom1	100000000" >> "$number.num"
        echo "chrom1	100000000" >> "$number.num"
        echo "chrom2	100000000" >> "$number.num"
        echo "chrom2	100000000" >> "$number.num"
    else
        echo "Greater than 2 chromosomes present.  Exiting script."
        exit 1
    fi

        rm $i
        mv "$number.num" "$number.txt") &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 2

cd "${curdr}"
}

#################################################################################

# This function make sure there is a matching coverage file for each vcf.

function checkMatchingCoverageFile () {

echo "Checking for coverageFile for each vcf"

# Make list of tb numbers of vcfs.
for i in *; do
getbase=`basename "$i"`
number=`echo $getbase | sed $tbNumberV | sed $tbNumberW`
echo $number >> vcfFiles
done

# --> This is blocked out.  It needs to be added to script 1 to prepare the coverage just once and having it ready for script 2 <--
#for f in $coverageFiles/*; do
#    num=1
#  for i in $chroms; do
#    sed "s/$i/chrom${num}/g" $f > $coverageFiles/temp.vcf
#    mv $coverageFiles/temp.vcf $f
#     echo "$i was marked as chrom${num}"
#     num=$(( $num + 1 ))
# done
#done

# --> This is for hitting all coverage files in directory regardless of chromosome numbers in each
#for f in $coverageFiles/*; do
#    chroms=`awk '{a[$1]++} END{for (var in a) print var}' $f | sort -n`
#    num=1
#    echo "chroms are: $chroms"
#    read -p "$LINENO, ENTER"
#  for i in $chroms; do
#   sed "s/$i/chrom${num}/g" $f > $coverageFiles/temp.vcf
#   mv $coverageFiles/temp.vcf $f
#   echo "$i was marked as chrom${num}"
#   num=$(( $num + 1 ))
#  done
#done

# Make list of tb numbers of coverageFiles
for f in $coverageFiles/*; do
getbase=`basename "$f"`
number=`echo $getbase | sed $tbNumberV | sed $tbNumberW | sed 's/-coverage//g'`
echo $number >> coverageFiles
done

uniqlist=`cat vcfFiles coverageFiles | sort | uniq -u`

rm vcfFiles
rm coverageFiles

uniqNumberSize=`echo $uniqlist | wc | awk '{print $3}'`
if [ $uniqNumberSize -gt 4 ]
then
echo "No matching file was present"
echo "Either there is no coverageFile for a vcf"
echo "Or there is an extra coverageFile that is not needed"
echo "Must have matching files, and restart script."
echo "$uniqlist"

else
echo "coverageFiles are available for all vcfs"
echo "Perfect match, Good Job Chap!"
fi

}

#################################################################################

# Change SNPs with low QUAL values to N, based on parameter set above in variable settings

function changeLowCalls () {

for i in *.vcf; do
(echo "Changeing low calls"; base=`basename $i .vcf`; awk -v x=$lowEnd -v y=$highEnd 'BEGIN {OFS="\t"} { if ($6 >= x && $6 <= y) print $1, $2, $3, $4, "N", $6, $7, $8; else print $0 }' $i > ${base}.txt; rm $i; mv ${base}.txt ${base}.vcf) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 2

}

#################################################################################

#   Function: fasta and table creation
function fasta_table () {

# Loop through the directories
directories=`ls`
echo "$directories"

for d in $directories; do
    cd ./$d/
    dir=`basename $PWD`
    echo "Directory:  $dir"
    mkdir starting_files
    cp *.vcf ./starting_files
    echo "***Grabbing vcf file names"

    if [ $FilterGroups == yes ]; then
         if [ $((chromCount)) -eq 1 ]; then
            #Mark vcf allowing areas of the genome to be removed from the SNP analysis
            for i in *.vcf; do
                (m=`basename "$i"`; n=`echo $m | sed $dropEXT`
                echo "***Adding filter to $n***"
                grep -v "#" $i | awk '{print $2}' > $i.file
                cat "${FilterDirectory}/$d.txt" $i.file >> $i.catFile
                cat $i.catFile | sort | uniq -d > $i.txt
                pos=`cat $i.txt | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/'`
                echo $pos
                awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $i > $n.subfilter.vcf
                rm $i.file
                rm $i.catFile
                rm $i.txt
                grep -v "Not_Included" $n.subfilter.vcf > $n.nomisfits.vcf
                mv $n.nomisfits.vcf $i
                rm $n.nomisfits.vcf)  &
                let count+=1
                [[ $((count%NR_CPUS)) -eq 0 ]] && wait
            done
        elif [ $((chromCount)) -eq 2 ]; then
            #Mark vcf allowing areas of the genome to be removed from the SNP analysis
                for i in *.vcf; do m=`basename "$i"`; n=`echo $m | sed $dropEXT` # n is name with all right of "_" and "." removed.
                #Mark vcf allowing areas of the genome to be removed from the SNP analysis
                echo "***Adding filter to $n***"
                grep -v "#" $i | awk 'BEGIN {OFS="\t"} {print $1, $2}' > $i.file
                cat "${FilterDirectory}/$d.txt" $i.file >> $i.catFile
                cat $i.catFile | sort -k1,1 -k2,2 | uniq -d > $i.txt
                pos1=`cat $i.txt | grep "chrom1" | awk '{print $2}' | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/'`
                echo "pos1: $pos1"
                awk -v var1="chrom1" -v var2=$pos1 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $i > $n.subfilter1.vcf
                pos2=`cat $i.txt | grep "chrom2" | awk '{print $2}' | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/'`
                echo "pos2: $pos2"
                awk -v var1="chrom2" -v var2=$pos2 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $n.subfilter1.vcf > $n.subfilter.vcf
                rm $i.file
                rm $i.catFile
                rm $i.txt
                rm $n.subfilter1.vcf
                grep -v "Not_Included" $n.subfilter.vcf > $n.nomisfits.vcf
                mv $n.nomisfits.vcf $i
                rm $n.nomisfits.vcf
                done
        else
        echo "Greater than 2 chromosomes present.  Exiting script."
        exit 1
        fi
    wait
    sleep 2

    mkdir marked_files
    mv *.subfilter.vcf ./marked_files
    else
        echo "************ Group/Subgroup/Clade Filter not ran *************"
        echo "************ Group/Subgroup/Clade Filter not ran *************" >> ../../log
    fi

    # Make concatemer with the position and REF call.
    echo "***Making Concatemer"
    for i in *.vcf; do
#egrep -v "#" $i | egrep "AC=2;A" | awk -v Q="$QUAL" '$6 > Q' | awk '{if ($1 ~ /chrom1/) print "chrom1-" $2, $4; else if ($1 ~ /chrom2/) print "chrom2-" $2, $4; else print "awk script can not determine VCF input"}' >> concatemer
    egrep -v "#" $i | egrep "AC=2;A" | awk -v Q="$QUAL" '$6 > Q' | awk '{print $1 "-" $2, $4}' >> concatemer
    done

    # Get rid of duplicates in concatemer and list all the positions and REF calls
    echo "***Making total_pos"
    cat concatemer | sort -nk1 | uniq | sort -k1.6n -k1.8n > total_pos

    # Count the number of SNPs
    echo "***The total_pos of $d" >> ../../log
    grep -c ".*" total_pos >> ../../log

    echo "***Creating normalized vcf using AC2, QUAL > $QUAL"
        # Grab the name of the vcf file
        for i in *.vcf; do
            (m=`basename "$i"`; n=`echo $m | sed $dropEXT`
            # echo the name grabbed
            echo $n
            # Create .cut file that lists the positions and ALT calls
#egrep -v "#" $i | egrep "AC=2;A" | awk -v Q="$QUAL" '$6 > Q' | awk '{if ($1 ~ /chrom1/) print "chrom1-" $2, $5; else if ($1 ~ /chrom2/) print "chrom2-" $2, $5; else print "awk script can not determine VCF input"}' > $n.cut
egrep -v "#" $i | egrep "AC=2;A" | awk -v Q="$QUAL" '$6 > Q' | awk '{print $1 "-" $2, $5}' > $n.cut

            # Fill in the .cut file with REF calls at positions that were not called as SNPs
#cat $n.cut total_pos | awk '{ if (a[$1]++ == 0) print $0; }' |  sort -nk1 > $n.filledcutnoN
cat $n.cut total_pos | awk '{ if (a[$1]++ == 0) print $0; }' |  sort -k1.6n -k1.8n > $n.filledcutnoN

    # Run coverage function on .filledcut to replace zero to low coverage positions with N.
    echo "$n Replacing positions with coverage below $Ncov with an -"
    
    getbase=`basename "$n.filledcutnoN"`
    numberWithV=`echo $getbase | sed $tbNumberV | sed $tbNumberW`; number=${numberWithV%V} #Changed from $tbNumberOnly
    echo "numberWithV: $numberWithV"
    echo "number to grep list of coverage files: $number"
    # Get the zero coverage file matching the .filledcut file
    zeroFile=`ls "$coverageFiles" | grep "$number"`
    echo "This is the coverage file that was found: $zeroFile"
    # Get the position of the zeroCoverage file
    # Get the position of the .filledcut file
            if [ $((chromCount)) -eq 1 ]; then
                awk -v z=$Ncov '$3 <  z {print $2, $3}' $coverageFiles/$zeroFile | awk '{print $1}' > ${number}-zeroCoverage
                elif [ $((chromCount)) -eq 2 ]; then
                echo "*******************************************************************************"
                # Get the position of the zeroCoverage file at first chromosome
                awk -v z=$Ncov '$1 ~ /chrom1/ && $3 <  z {print chrom1","$2}' $coverageFiles/${number}-coverage > ${number}-zeroCoverage
                # Get the position of the zeroCoverage file at second chromosome
                awk -v z=$Ncov '$1 ~ /chrom2/ && $3 <  z {print chrom2","$2}' $coverageFiles/${number}-coverage >> ${number}-zeroCoverage
                echo "*******************************************************************************"
                else
                echo "Greater than 2 chromosomes present.  Exiting script."
                exit 1
            fi
            awk '{print $1}' $n.filledcutnoN > ${number}-filledcutNumbers
            # cat the two above together
            cat ${number}-zeroCoverage ${number}-filledcutNumbers > ${number}-NcatFile
            # Get duplicates of the cat file
            # This duplicats need to replace A,T,G or C to an N, Would like to make it a "." or "-" but Geneious does not recognize these characters.
            cat ${number}-NcatFile | sort | uniq -d > ${number}-duplicates
            echo "1000000000" >> ${number}-duplicates
            # Prepare duplicate regions to be used as variables in awk, to find and replace
            pos=`cat ${number}-duplicates | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/'`
            echo "Zero Coverage: $pos"
#awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ x ) print $1, "-"; else print $0}' $n.filledcutnoN > $n.filledcut
awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ x ) print $1, "-"; else print $0}' $n.filledcutnoN | sort -k1.6n -k1.8n > $n.filledcut

            rm ${number}-zeroCoverage
            rm ${number}-filledcutNumbers
            rm ${number}-NcatFile
            rm ${number}-duplicates
            rm $n.filledcutnoN
            rm $i) &
            let count+=1
            [[ $((count%NR_CPUS)) -eq 0 ]] && wait
            done
            echo "sleeping 10 seconds at line number: $LINENO"; sleep 10
            wait



        # Make a concatemer of the .filledcut files
        for i in *.filledcut; do
            cat $i >> cutConcatemer
            done
        echo "***Making the select file containing positions of interest"
        # Capture only positions that have more than one SNP type called at a position
        cat cutConcatemer | sort -k1.6n -k1.8n | uniq | awk '{print $1}' | uniq -d > select
        # Compare the positions in select with total_pos and output total_pos position that match select positions
        # but only with positions that are in the select file.
        # This getting rid of calls that are the same for all isolates being analyzed
        echo "***grepping the total_pos file"

        echo "***grepping the .filledcut files"

        grep -f select total_pos | sort -k1.6n -k1.8n > clean_total_pos
        # Begin the table
        awk '{print $1}' clean_total_pos | awk 'BEGIN{print "reference_pos"}1' | tr '\n' '\t' | sed 's/$//' >> $d.table.txt
        awk '{print $2}' clean_total_pos | awk 'BEGIN{print "reference_call"}1' | tr '\n' '\t' | sed 's/$//' >> $d.table.txt
                    # Make the fasta files:  Fill in positions with REF if not present in .clean file

        for i in *.filledcut; do
            (m=`basename "$i"`
            n=`echo $m | sed $dropEXT`
            # Compare the positions in select with "isolate".cut and output position for .cut that only matched select positions
            grep -f select $i | sort -k1.6n -k1.8n > $n.tod
            sed 's/chrom[0-9-]*//g' $n.tod | tr -d [:space:] | sed "s/^/>$n;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g'  > $n.fas
            # Add each isolate to the table
            awk '{print $2}' $n.tod | awk -v number="$n" 'BEGIN{print number}1' | tr '\n' '\t' | sed 's/$//' >> $d.table.txt) &
            let count+=1
            [[ $((count%NR_CPUS)) -eq 0 ]] && wait
            done
                wait
echo "sleeping 10 seconds at line number: $LINENO"; sleep 10
        #Make a reference fasta sequence
        awk '{print $2}' clean_total_pos > reference_seq
        cat reference_seq | tr -cd "[:print:]" | sed "s/^/>reference_seq;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > reference_seq.fas
        echo "***The clean_total_pos" >> ../../log
        grep -c ".*" clean_total_pos >> ../../log

        # Make a file containing all fasta files.
        cat *.fas > ${dir}_alignment.fasta

        #Clean-up
        echo "***Cleaning folder"
        rm *.cut
        rm *.filledcut
        rm concatemer
        rm cutConcatemer
        rm *.tod
        mkdir fasta
        mv *.fas ./fasta
        rm total_pos
        rm select
        rm reference_seq
        echo "***Done"
    cd ..
    done
}

#################################################################################
#################################################################################
#################################################################################
echo "****************************** START ******************************"
#################################################################################
#################################################################################
#################################################################################

# Clean the tag file that has been exported to Desktop
chmod 777 /Users/tstuberadmin/Desktop/Untitled.tab
cat /Users/tstuberadmin/Desktop/Untitled.tab | tr '\r' '\n' | awk -F '\t' 'BEGIN{OFS="\t";} {gsub("\"","",$5);print;}' | sed 's/\"##/##/' > /Users/tstuberadmin/Desktop/Untitled.txt

clean_tag.sh $genotypingcodes

# Test for duplicate VCFs
testDuplicates

#Prepare Filter files.
filterFilespreparation

#Test for match coverage file
checkMatchingCoverageFile

#copy the original vcfs to /starting_files
mkdir starting_files
for i in *.*; do
    cp $i ./starting_files
done

# Remove selected isolates from comparison
# This is optional, and should be turned on or off based on laboratories preference
removeIsolates

############################### Rename files ###############################

for i in *.txt; do
    mv $i ${i%.txt}.vcf
done

for i in *.vcf; do
    echo "******************** Naming convention ********************"
    echo "Original File: $i"
    base=`basename "$i"`
    searchName=`echo $base | sed $tbNumberV | sed $tbNumberW | sed 's/V//'`
    echo "searchName: $searchName"
    # Direct script to text file containing a list of the correct labels to use.
    # The file must be a txt file.
    p=`grep "$searchName" "$genotypingcodes"`
    echo "This is what was found in tag file: $p"
    newName=`echo $p | awk '{print $1}' | tr -d "[:space:]"` # Captured the new name
    n=`echo $base | sed $tbNumberV | sed $tbNumberW`
    noExtention=`echo $base | sed $dropEXT`
    VALtest=`echo $i | grep "VAL"`
    echo "VALtest: $VALtest"
    #Check if a name was found in the tag file.  If no name was found, keep original name, make note in log and cp file to unnamed folder.
    if [[ -z "$p" ]]; then # new name was NOT found
        if [[ -z "$VALtest" ]]; then
            name=$searchName
            echo "n is $n"
            echo "--> $name was not found in the FileMaker database" >> log
            mkdir -p FilesNotRenamed
            cp $i ./FilesNotRenamed
            mv $i ${name}.vcf
            echo "A"
        else
            name=${searchName}-Val
            mv $i ${name}.vcf
            echo "B"
        fi
    else # New name WAS found
        if [[ -z "$VALtest" ]]; then
            name=$newName
            mv $i ${name}.vcf
            echo "C"
        else
            name=${newName}-Val
            echo "newName is $name"
            mv $i ${name}.vcf
            echo "D"
        fi
    fi
done

##################### Start: Make Files Unix Compatiable #####################

#Fix validated (VAL) vcf files.  This is used in vcftofasta scripts to prepare validated vcf files opened and saved in Excel.
#Create list of isolates containing "VAL"
echo "Making Files Unix Compatiable"

for v in *.vcf; do
    dos2unix $v #Fixes files opened and saved in Excel
    cat $v | tr '\r' '\n' | awk -F '\t' 'BEGIN{OFS="\t";} {gsub("\"","",$5);print;}' | sed 's/\"##/##/' > temp
    mv temp $v
done

############## Capture the number of chromosomes and their name from a single VCF ##############

echo "The chromosome count is: $chromCount"
pwd

if [ $chromCount -gt 2 ]; then
    echo "More than 2 chromosomes detected."
    echo "Unable to handle more than 2."
    echo "Exiting script at line $LINENO"
    exit 1
fi
# Change chromosome identification to general chrom1 and/or chrom2

for f in *.vcf; do
echo "echoing f: $f"
    num=1
    for i in $chroms; do
        sed "s/$i/chrom${num}/g" $f > temp.vcf
        mv temp.vcf $f
        echo "$i was marked as chrom${num}"
        num=$(( $num + 1 ))
    done

done

########################################################################

date >> log
printf "%s\t%s\t%s\t%s\n" "TB Number" "Group" "Subgroup" "Clade" > FileMakerGroupImport.txt

AConeCallPosition
wait

# Change low QUAL SNPs to N, see set variables
changeLowCalls
wait

######################## Mark Files and Remove Marked Regions ########################

if [ $FilterAllVCFs == yes ]; then
echo "***Marking all VCFs and removing filtering regions"
echo "***Marking all VCFs and removing filtering regions was done." >> log
    # Label filter field for positions to be filtered in all VCFs
        if [ $((chromCount)) -eq 1 ]; then
        for i in *.vcf; do
            (m=`basename "$i"`; n=`echo $m | sed $dropEXT`; echo "********* $n **********"; grep -v "#" $i | awk '{print $2}' > $i.file; cat "${FilterDirectory}/FilterToAll.txt" $i.file >> $i.catFile; cat $i.catFile | sort | uniq -d > $i.txt; pos=`cat $i.txt | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/'`; awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $i > $n.filter.vcf; rm $i.file; rm $i.catFile; rm $i.txt; grep -v "Not_Included" $n.filter.vcf > $n.noPPE.vcf;  mv $n.noPPE.vcf $i) &
            let count+=1
            [[ $((count%NR_CPUS)) -eq 0 ]] && wait
        done

        elif [ $((chromCount)) -eq 2 ]; then
        for i in *.vcf; do
            m=`basename "$i"`
            n=`echo $m | sed $dropEXT` # n is name with all right of "_" and "." removed.
            #Mark vcf allowing areas of the genome to be removed from the SNP analysis
            echo "********* $n **********"
            grep -v "#" $i | awk 'BEGIN {OFS="\t"} {print $1, $2}' > $i.file
            cat "${FilterDirectory}/FilterToAll.txt" $i.file >> $i.catFile
            cat $i.catFile | sort -k1,1 -k2,2 | uniq -d > $i.txt
            pos1=`cat $i.txt | grep "chrom1" | awk '{print $2}' | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/'`
            echo "pos1: $pos1"
            awk -v var1="chrom1" -v var2=$pos1 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $i > $n.filter1.vcf
            pos2=`cat $i.txt | grep "chrom2" | awk '{print $2}' | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/'`
            echo "pos2: $pos2"
            awk -v var1="chrom2" -v var2=$pos2 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' $n.filter1.vcf > $n.filter.vcf
            rm $i.file
            rm $i.catFile
            rm $i.txt
            rm $n.filter1.vcf
            grep -v "Not_Included" $n.filter.vcf > $n.noPPE.vcf
            mv $n.noPPE.vcf $i
        done
        else
            echo "Greater than 2 chromosomes present.  Exiting script."
        exit 1
        fi

    wait
    sleep 2

        mkdir marked_files
        mv *.filter.vcf ./marked_files
    else
    echo "***Marking all VCFs and removing marked regions was NOT done."
    echo "***Marking all VCFs and removing marked regions was NOT done." >> log
fi

#################### Categorize VCFs into Groups, Subgroups and Clades #####################

# Print header
    printf "%10s %10s %10s %10s\n" "Isolate" "Group" "Subgroup" "Clade" >> log

for i in *.vcf; do
    # Get quality positions in VCF
    formatedpos=`grep -v "#" $i | egrep "AC=2;A" | awk -v Q="$QUAL" '$6 > Q' | awk '{print $2}' | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/'`
    
    # If a group number matches a quality position in the VCF (formatedpos) then print the position
    groupNumbers=`grep "Group" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}'`
    echo "This is the Group Numbers: $groupNumbers"

    # Typically a single group position is found, and the VCF will be placed into just one group.  It is posible that an isolate will need to go in more than one group because of were it falls on the tree.  In this case there may be 2 group, or more, group positions found.  The number of group positions found is captured in sizeGroup.
    sizeGroup=`echo $groupNumbers | wc | awk '{print $3}'`

    # Loop through the number of groups positions found
    loops=`grep "Group" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}' | awk 'END {print NR}'`
    if [ $sizeGroup -lt 2 ]; then # There was not a position found that places VCF into group
        echo "--> Group assignment for $i was not found.  Consider adding a new group position." >> log
        echo "$i was not assigned a Group"

        #If a name was found in the tag file the name is changed
        #And transferred to all 3 groups: All_vcfs, clade and subgroup
        #Filtered vcf is getting renamed (from Excel file) and shuttled to All_vcfs, Clades and subclades
        else
        counter=1
            while [ $counter -le $loops ]; do
                    if [ $loops -gt 1 ]; then
                        echo "========= $i was assigned to MORE than one Group =========" >> log
                    fi
                echo "The counter is at $counter"
                groupNumber=`echo "$groupNumbers" | tr "\n", "\t" | awk -v x=$counter '{print $x}'`
                echo $groupNumber
                mkdir -p Group-$groupNumber #Make groupNumber folder if one does not exist.
                cp $i ./Group-$groupNumber/ #Then copy to each folder
                let counter=counter+1
            done
        fi

    # Get the Subgroup number
    subgroupNumbers=`grep "Subgroup" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}'`
    echo "This is the Subgroup Numbers: $subgroupNumbers"
        # Check if a Group, subgroup or clade was called.
        sizeGroup=`echo $subgroupNumbers | wc | awk '{print $3}'`

    loops=`grep "Subgroup" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}' | awk 'END {print NR}'`

        if [ $sizeGroup -lt 2 ];
        then
        echo "$i was not assigned a Subgroup"

        #If a name was found in the tag file the name is changed
        #And transferred to all 3 groups: All_vcfs, clade and subgroup
        #Filtered vcf is getting renamed (from Excel file) and shuttled to All_vcfs, Clades and subclades
        else
        counter=1
            while [ $counter -le $loops ]; do
                    if [ $loops -gt 1 ]; then
                    echo "========= $i was assigned to MORE than one Subgroup =========" >> log
                    fi
                echo "The counter is at $counter"
                subgroupNumber=`echo "$subgroupNumbers" | tr "\n", "\t" | awk -v x=$counter '{print $x}'`
                echo $subgroupNumber
                mkdir -p Subgroup-$subgroupNumber #Make groupNumber folder if one does not exist.
                cp $i ./Subgroup-$subgroupNumber/ #Then copy to each folder
                let counter=counter+1
            done
        fi

    # Get the Clade number
    cladeNumbers=`grep "Clade" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}'`
    echo "This is the Clade Numbers: $cladeNumbers"

printf "%10s %10s %10s %10s\n" "${i%.vcf}" "$groupNumbers" "$subgroupNumbers" "$cladeNumbers" >> log
printf "%s\t%s\t%s\t%s\n" "${i%.vcf}" "$groupNumbers" "$subgroupNumbers" "$cladeNumbers" >> FileMakerGroupImport.txt
        # Check if a Group, subgroup or clade was called.
        sizeGroup=`echo $cladeNumbers | wc | awk '{print $3}'`
        loops=`grep "Clade" "${DefiningSNPs}" | awk -v x=$formatedpos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1}' | awk 'END {print NR}'`
        if [ $sizeGroup -lt 2 ];
        then

        echo "$i was not assigned a Clade"

        #If a name was found in the tag file the name is changed
        #And transferred to all 3 groups: All_vcfs, clade and subgroup
        #Filtered vcf is getting renamed (from Excel file) and shuttled to All_vcfs, Clades and subclades
        else
        counter=1
            while [ $counter -le $loops ]; do
                    if [ $loops -gt 1 ]; then
                    echo "========= $i was assigned to MORE than one Clade =========" >> log
                    fi
                echo "The counter is at $counter"
                cladeNumber=`echo "$cladeNumbers" | tr "\n", "\t" | awk -v x=$counter '{print $x}'`
                echo $cladeNumber
                mkdir -p Clade-$cladeNumber #Make groupNumber folder if one does not exist.
                cp $i ./Clade-$cladeNumber/ #Then copy to each folder
                let counter=counter+1
            done
        fi

    mkdir -p All_vcfs #Make All_vcfs folder if one does not exist.
    mv $i ./All_vcfs/

done

##################### Organize folders #####################

mkdir All_Groups
mv ./Group-*/ ./All_Groups
mkdir All_Subgroups
mv ./Subgroup*/ ./All_Subgroups/
mkdir All_Clades
mv ./Clade*/ ./All_Clades/

##################### Start: All vcf folder #####################
cd ./All_vcfs/
# Make concatemer with the position and REF call.
    # Factor in possible multiple chromosomes
    echo "***Making Concatemer"
    for i in *.vcf; do
#egrep -v "#" $i | egrep "AC=2;A" | awk -v Q="$QUAL" '$6 > Q' | awk '{if ($1 ~ /chrom1/) print "chrom1-" $2, $4; else if ($1 ~ /chrom2/) print "chrom2-" $2, $4; else print "awk script can not determine VCF input"}' >> concatemer
egrep -v "#" $i | egrep "AC=2;A" | awk -v Q="$QUAL" '$6 > Q' | awk '{print $1 "-" $2, $4}' >> concatemer
    done

# Get rid of duplicates in concatemer and list all the positions and REF calls
echo "***Making total_pos"
cat concatemer | sort -nk1 | uniq | sort -k1.6n -k1.8n > total_pos

# Count the number of SNPs
echo "***The total_pos of $d" >> ../log
grep -c ".*" total_pos >> ../log
echo "***Creating normalized vcf using AC2, QUAL > 150"
# Grab the name of the vcf file

for i in *.vcf; do
(m=`basename "$i"`; n=`echo $m | sed $dropEXT`; echo "$n"
#egrep -v "#" $i | egrep "AC=2;A" | awk -v Q="$QUAL" '$6 > Q' | awk '{if ($1 ~ /chrom1/) print "chrom1-" $2, $5; else if ($1 ~ /chrom2/) print "chrom2-" $2, $5; else print "awk script can not determine VCF input"}' > $n.cut
egrep -v "#" $i | egrep "AC=2;A" | awk -v Q="$QUAL" '$6 > Q' | awk '{print $1 "-" $2, $5}' > $n.cut

#cat $n.cut total_pos | awk '{ if (a[$1]++ == 0) print $0; }' |  sort -nk1 > $n.filledcutnoN
cat $n.cut total_pos | awk '{ if (a[$1]++ == 0) print $0; }' |  sort -k1.6n -k1.8n > $n.filledcutnoN

echo "$n Replacing positions with coverage below $Ncov with an -"
getbase=`basename "$n.filledcutnoN"`
numberWithV=`echo $getbase | sed $tbNumberV | sed $tbNumberW`; number=${numberWithV%V} #Changed from $tbNumberOnly
echo "numberWithV: $numberWithV"
echo "number: $number"
echo "number to grep list of coverage files: $number"
# Get the zero coverage file matching the .filledcut file
zeroFile=`ls "$coverageFiles" | grep "$number"`
echo "This is the coverage file that was found: $zeroFile"
echo "This is the number: $number"
echo "Zero coverage file found: $zeroFile"

# Get the position of the zeroCoverage file
if [ $((chromCount)) -eq 1 ]; then
        awk -v z=$Ncov '$3 <  z {print $2, $3}' $coverageFiles/$zeroFile | awk '{print "chrom1-" $1}' > ${number}-zeroCoverage
    elif [ $((chromCount)) -eq 2 ]; then
        echo "*******************************************************************************"
        # Get the position of the zeroCoverage file at first chromosome
        awk -v z=$Ncov '$1 ~ /chrom1/ && $3 <  z {print "chrom1-" $2}' $coverageFiles/${number}-coverage > ${number}-zeroCoverage
        # Get the position of the zeroCoverage file at second chromosome
        awk -v z=$Ncov '$1 ~ /chrom2/ && $3 <  z {print "chrom2-" $2}' $coverageFiles/${number}-coverage >> ${number}-zeroCoverage
        echo "*******************************************************************************"
    else
        echo "Greater than 2 chromosomes present.  Exiting script."
    exit 1
fi
# Get the position of the .filledcut file
awk '{print $1}' $n.filledcutnoN > ${number}-filledcutNumbers
# cat the two above together
cat ${number}-zeroCoverage ${number}-filledcutNumbers > ${number}-NcatFile
# Get duplicates of the cat file
# This duplicats need to replace A,T,G or C to an N, Would like to make it a "." or "-" but Geneious does not recognize these characters.
cat ${number}-NcatFile | sort | uniq -d > ${number}-duplicates
echo "1000000000" >> ${number}-duplicates
# Prepare duplicate regions to be used as variables in awk, to find and replace
pos=`cat ${number}-duplicates | tr "\n" "W" | sed 's/W/\$\|\^/g' | sed 's/\$\|\^$//' | sed 's/$/\$/' | sed 's/^/\^/'`
echo "Zero Coverage: $pos"
#awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ x ) print $1, "-"; else print $0}' $n.filledcutnoN | sort -nk1 > $n.filledcut
awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ x ) print $1, "-"; else print $0}' $n.filledcutnoN | sort -k1.6n -k1.8n > $n.filledcut

rm ${number}-zeroCoverage
rm ${number}-filledcutNumbers
rm ${number}-NcatFile
rm ${number}-duplicates
rm $n.filledcutnoN) &
let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait

done
wait
sleep 2

# Begin the table
awk '{print $1}' total_pos | awk 'BEGIN{print "reference_pos"}1' | tr '\n' '\t' | sed 's/$//' >> All_vcf.table.txt
awk '{print $2}' total_pos | awk 'BEGIN{print "reference_call"}1' | tr '\n' '\t' | sed 's/$//' >> All_vcf.table.txt
echo "***grepping the .filledcut files"
# Make the fasta files:  Fill in positions with REF if not present in .clean file

for i in *.filledcut; do
    echo " working on filled cut for $i"
    m=`basename "$i"`
    n=`echo $m | sed $dropEXT`
    # Compare the positions in select with "isolate".cut and output position for .cut that only matched select positions
    #grep -f select $i > $n.tod

    # Use this cat command to skip the time intensive grep
    # With All_vcfs this grep doesn't eliminate many snps.
    sed 's/chrom[0-9-]*//g' $i | tr -d [:space:] | sed "s/^/>$n;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g'  > $n.fas
    # Add each isolate to the table
    awk '{print $2}' $i | awk -v number="$n" 'BEGIN{print number}1' | tr '\n' '\t' | sed 's/$//' >> All_vcf.table.txt
    done
    wait

echo "grepping filledcut files is finished"
#Make a reference fasta sequence
awk '{print $2}' total_pos > reference_seq
cat reference_seq | tr -cd "[:print:]" | sed "s/^/>reference_seq;/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > reference_seq.fas
echo "***The clean_total_pos $d" >> ../log
grep -c ".*" total_pos >> ../log

        # Make a file containing all fasta files.
        cat *.fas > all_alignment.fasta

#Clean-up
mkdir starting_files
echo "***Cleaning folder"
mv *.vcf ./starting_files
rm *.cut
rm *.filledcut
rm *.filledcutnoN
rm concatemer
rm cutConcatemer
rm *.tod
mkdir fasta
mv *.fas ./fasta
#rm total_pos
rm reference_seq
echo "***Done"
echo "Full Directory: ${fulDir}"

##################### End: All vcf folder #####################

echo "***************************************************"
echo "***************** STARTING Groups *****************"
echo "***************************************************"
# Change directory to All_Groups
cd ${fulDir}/All_Groups
fasta_table &

echo "***************************************************"
echo "**************** STARTING SUBGROUPS ***************"
echo "***************************************************"
# Change directory to All_Subgroups
cd ${fulDir}/All_Subgroups
fasta_table &

echo "***************************************************"
echo "***************** STARTING CLADES *****************"
echo "***************************************************"
# Change directory to All_Clades
cd ${fulDir}/All_Clades
fasta_table &
wait
echo "At line $LINENO, sleeping 5 second"; sleep 5s
cd ${fulDir}
$PWD
cp /Volumes/Data_HD/Brucella/Table_Template.xlsx ./
cp "$0" "$PWD"
date >> log
#****************************************************************
function alignTable () {

# Beginning in fasta folder
echo "$d *********"
pwd

cat *.fas | sed '/reference/{N;d;}' >> fastaGroup.txt
clustalw2 -OUTFILE=alignment.txt -RANGE=1,2 -OUTPUT=FASTA -INFILE=fastaGroup.txt &&
grep ">" alignment.txt | sed 's/>//g' > cleanedAlignment.txt
awk 'NR==FNR{o[FNR]=$1; next} {t[$1]=$0} END{for(x=1; x<=FNR; x++){y=o[x]; print t[y]}}' cleanedAlignment.txt ../$d.table.txt > joined.txt
grep "reference" ../$d.table.txt > references.txt
cat references.txt joined.txt >> joined2.txt
mv joined2.txt ../$d.sortedTable.txt

rm alignment.txt
rm cleanedAlignment.txt
rm *.dnd
rm fastaGroup.txt
rm joined.txt
rm references.txt

echo "**** orgTable.sh Started ****"
cd ..
#get the number of columns
columnCount=`awk '$0 !~ /^$/ {print $0}' *sortedTable.txt | awk '{print NF-1; exit}'`
#echo "Column Count: $columnCount"

number=`jot - 1 $columnCount`
#echo "Numbers in list: $number"

#get the reference row
cat *sortedTable.txt | sed -n 2p | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' > referenceRow
cat referenceRow > out2.txt
cat *sortedTable.txt | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' > table
#remove first column from *sortedTable.txt

echo "countDif" > countOutput.txt
echo "countFirst" > firstOutput.txt
#iterate numbers up to number of columns
for n in $number; do
#use number to grab character in reference row i.e. C
letter=`awk -v x=$n '{print $x}' referenceRow`
#echo "Number: $n"
#echo "Letter: $letter"
awk -v var1="$n" '{print $var1}' table > column
grep "$letter" column | wc -l >> countOutput.txt
sed '1,2d' column > column2
cat column2 | awk -v var2=${letter} ' $0 !~ var2 {print NR; exit}' >> firstOutput.txt
done

#/var2/ {print NR; exit}' column
#Clear table2.txt
echo "" > table2.txt
#Add a \n to the end of file
sed -a '$a\' *sortedTable.txt >> table2.txt
#Prepare count line
cat countOutput.txt | sed 's/ //g' | tr "\n" "\t" > readyline.txt
#Create readytable.txt with count line.
cat table2.txt readyline.txt > orginizedTable2.txt
grep -v "^$" orginizedTable2.txt > orginizedTable3.txt

#Prepare firstOut line
cat firstOutput.txt | sed 's/ //g' | tr "\n" "\t" > readyFirstOut.txt
cat orginizedTable3.txt readyFirstOut.txt > orginizedTable4.txt

rm referenceRow
rm out2.txt
rm table
rm countOutput.txt
rm table2.txt
rm readyline.txt
rm orginizedTable2.txt

awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' orginizedTable4.txt > orginizedTable5.txt
awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' orginizedTable5.txt > orginizedTable6.txt

#Transpose
awk -f /Users/Shared/_programs/_my_scripts/awk_scripts/transpose.awk orginizedTable6.txt > orginizedTable7.txt

#Orgainize file based on 1st 2 columns
sort -n -k1 orginizedTable7.txt | sort -n -k2 > orginizedTable8.txt

#Convert spaces to tabs
awk -v OFS="\t" '$1=$1' orginizedTable8.txt > orginizedTable9.txt
awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' orginizedTable9.txt | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' > orginizedTable10.txt

#Transpose back
awk -f /Users/Shared/_programs/_my_scripts/awk_scripts/transpose.awk orginizedTable10.txt > orginizedTable11.txt

c=`basename $PWD`
#Convert spaces to tabs
awk -v OFS="\t" '$1=$1' orginizedTable11.txt > $c.organizedTable.txt

rm orginizedTable3.txt
rm orginizedTable4.txt
rm orginizedTable5.txt
rm orginizedTable6.txt
rm orginizedTable7.txt
rm orginizedTable8.txt
rm orginizedTable9.txt
rm orginizedTable10.txt
rm orginizedTable11.txt
rm column
rm column2
rm readyFirstOut.txt
rm firstOutput.txt
echo "**** orgTable.sh Finished ****"
}
#****************************************************************

echo "***************************************************"
echo "********** STARTING All_Clades Alignment **********"
echo "***************************************************"
cd ${fulDir}/All_Clades

workingdir=`basename $PWD`

if [ $workingdir == All_Clades ]
then
directories=`ls`
for d in $directories; do
    cd ${fulDir}/All_Clades/${d}/fasta
    echo "****************************************************"
    echo "************* Orginizing Table: $d *****************"
    echo "****************************************************"
    alignTable &

    pwd
done

else
echo "*** $workingdir not found ***"
fi

echo "***************************************************"
echo "********** STARTING All_Groups Alignment **********"
echo "***************************************************"
cd ${fulDir}/All_Groups

workingdir=`basename $PWD`

if [ $workingdir == All_Groups ]
then
directories=`ls`
for d in $directories; do
    cd ${fulDir}/All_Groups/${d}/fasta
    echo "****************************************************"
    echo "************* Orginizing Table: $d *****************"
    echo "****************************************************"
    alignTable &
pwd
done
else
echo "*** $workingdir not found ***"
fi

echo "***************************************************"
echo "******** STARTING All_SubGroups Alignment *********"
echo "***************************************************"
cd ${fulDir}/All_Subgroups

workingdir=`basename $PWD`

if [ $workingdir == All_Subgroups ]
then

directories=`ls`
for d in $directories; do
    cd ${fulDir}/All_Subgroups/$d/fasta
    echo "****************************************************"
    echo "************* Orginizing Table: $d *****************"
    echo "****************************************************"
    alignTable &
pwd
done
else
echo "*** $workingdir not found ***"
fi
wait
pwd
cd ${fulDir}
mv log log.txt
cp -r $PWD ${bioinfoVCF}
echo "******* $LINENO, $PWD"
fileName=`basename $0`

mail -s "$fileName $@ completed" tod.p.stuber@aphis.usda.gov < log.txt
mail -s "$fileName $@ completed" suelee.robbe-austerman@aphis.usda.gov < log.txt
mail -s "$fileName $@ completed" christine.r.quance@aphis.usda.gov < log.txt


echo "****************************** END ******************************"
pwd

#
#  Created by Stuber, Tod P - APHIS on 5/3/2014.
#
