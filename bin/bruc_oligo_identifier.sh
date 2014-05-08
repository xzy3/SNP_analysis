#!/bin/sh

#################################################################################

echo "**********************ID FILES START**********************"


pab1="AATTGTCGGATAGCCTGGCGATAACGACGC"
pab3="CACACGCGGGCCGGAACTGCCGCAAATGAC"
pab5="GCTGAAGCGGCAGACCGGCAGAACGAATAT"
pmel="TGTCGCGCGTCAAGCGGCGTGAAATCTCTG"
psuis1="TGCGTTGCCGTGAAGCTTAATTCGGCTGAT"
psuis2="GGCAATCATGCGCAGGGCTTTGCATTCGTC"
psuis3="CAAGGCAGATGCACATAATCCGGCGACCCG"
pceti1="GTGAATATAGGGTGAATTGATCTTCAGCCG"
pceti2="TTACAAGCAGGCCTATGAGCGCGGCGTGAA"
pcanis4="CTGCTACATAAAGCACCCGGCGACCGAGTT"
pcanis="ATCGTTTTGCGGCATATCGCTGACCACAGC"
povis="CACTCAATCTTCTCTACGGGCGTGGTATCC"

base=`basename $1`
forReads=`echo $1 | grep _R1`
echo "Forward Reads:  $forReads"

name=`echo $base | grep _R1`
n=`echo $name | sed 's/_.*//' | sed 's/\..*//'`
echo "working on: $n"

ab1=`grep $pab1 $forReads | wc -l`
ab3=`grep $pab3 $forReads | wc -l`
ab5=`grep $pab5 $forReads | wc -l`
mel=`grep $pmel $forReads | wc -l`
suis1=`grep $psuis1 $forReads | wc -l`
suis2=`grep $psuis2 $forReads | wc -l`
suis3=`grep $psuis3 $forReads | wc -l`
ceti1=`grep $pceti1 $forReads | wc -l`
ceti2=`grep $pceti2 $forReads | wc -l`
canis4=`grep $pcanis4 $forReads | wc -l`
canis=`grep $pcanis $forReads | wc -l`
ovis=`grep $povis $forReads | wc -l`

counts=`echo "$ab1 $ab3 $ab5 $mel $suis1 $suis2 $suis3 $ceti1 $ceti2 $canis4 $canis $ovis"`
echo $counts

binary=`echo $counts | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]"`
echo $binary

i=$binary

if [ $i == 111111111111 ]
then
	catch=`echo "*** Odd Isolate, Unexpected findings ***"`
    exit 1

elif [ $i == 011111111111 ]
then
	catch=`echo "Brucella abortus bv 1, 2 or 4"`

    `processZips.sh ab1 $catch | tee tee_processZips_out.txt` &
    echo "$catch" > tee_bruc_oligo_identifier_out2.txt

elif [ $i == 101111111111 ]
	then
	catch=`echo "Brucella abortus bv 3"`
    `processZips.sh ab1 $catch | tee tee_processZips_out.txt` &
    echo "$catch" > tee_bruc_oligo_identifier_out2.txt

elif [ $i == 110111111111 ]
    then
    catch=`echo "Brucella abortus bv 5, 6 or 9"`
    `processZips.sh ab1 $catch | tee tee_processZips_out.txt` &
echo "$catch" > tee_bruc_oligo_identifier_out2.txt

elif [ $i == 111011111111 ]
	then
    catch=`echo "Brucella melitensis"`
    `processZips.sh mel $catch | tee tee_processZips_out.txt` &
echo "$catch" > tee_bruc_oligo_identifier_out2.txt

elif [ $i == 111101111111 ]
	then
	catch=`echo "Brucella suis bv1"`
    `processZips.sh suis1 $catch | tee tee_processZips_out.txt` &
    echo "$catch" > tee_bruc_oligo_identifier_out2.txt

elif [ $i == 111110111111 ]
	then
	catch=`echo "Brucella suis bv2"`
    `processZips.sh suis2 $catch | tee tee_processZips_out.txt` &
    echo "$catch" > tee_bruc_oligo_identifier_out2.txt

elif [ $i == 111111011111 ]
	then
	catch=`echo "Brucella suis bv3"`
    `processZips.sh suis3 $catch | tee tee_processZips_out.txt` &
    echo "$catch" > tee_bruc_oligo_identifier_out2.txt

elif [ $i == 111111101111 ] || [ $i == 111111100111 ]
    then
    catch=`echo "Brucella ceti 1"`
    `processZips.sh ceti1 $catch | tee tee_processZips_out.txt` &
    echo "$catch" > tee_bruc_oligo_identifier_out2.txt

elif [ $i == 111111110111 ]
    then
    catch=`echo "Brucella ceti 2"`
    `processZips.sh ceti2 $catch | tee tee_processZips_out.txt` &
    echo "$catch" > tee_bruc_oligo_identifier_out2.txt

elif [ $i == 111111111011 ]
    then
    catch=`echo "Brucella suis bv4"`
    `processZips.sh suis4 $catch | tee tee_processZips_out.txt` &
    echo "$catch" > tee_bruc_oligo_identifier_out2.txt

elif [ $i == 111111111001 ]
    then
    catch=`echo "Brucella canis"`
    `processZips.sh canis $catch | tee tee_processZips_out.txt` &
    echo "$catch" > tee_bruc_oligo_identifier_out2.txt

elif [ $i == 111111111110 ]
    then
    catch=`echo "Brucella ovis"`
    `processZips.sh ovis $catch | tee tee_processZips_out.txt` &
    echo "$catch" > tee_bruc_oligo_identifier_out2.txt

else
	catch=`echo "*** Odd Isolate, Unexpected findings, See /Volumes/Data_HD/Brucella/bruc_oligo_identifier_output.txt ***"`
    echo "***bruc_oligo_identifier cannot find a pattern for $n, see line $LINENO of script***"
    exit 1
fi

tagname=`grep $n /Volumes/Data_HD/Brucella/bruc_tags.txt`
echo "$tagname"
echo "$tagname","$n","$binary","$catch" "$counts" >> /Volumes/Data_HD/Brucella/bruc_oligo_identifier_output.txt



#
#  Created by Stuber, Tod P - APHIS on 04/11/2014.
#
