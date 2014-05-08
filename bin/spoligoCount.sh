#!/bin/sh

#  spoligoCount.sh
#  Count the number of sequencing found.

publishedSpacer01="ATAGAGGGTCGCCGGTTCTGGATCA|TGATCCAGAACCGGCGACCCTCTAT"
publishedSpacer02="CCTCATAATTGGGCGACAGCTTTTG|CAAAAGCTGTCGCCCAATTATGAGG"

observedSpacer01="TGATCCAGAGCCGGCGACCCTCTAT|ATAGAGGGTCGCCGGCTCTGGATCA"
observedSpacer02="CAAAAGCTGTCGCCCAA|TTGGGCGACAGCTTTTG"

#workingDir="/Volumes/Data_HD/Mycobacterium/_TB-Data"

#cd $workingDir
# set working directory to correspond to "mylist" file
#mylist=`cat /Volumes/Data_HD/Mycobacterium/listFiles/${1}`
#echo "$mylist"

#echo "Isolate\tObserved Spacer01\tObserved Spacer02\tPublished Spacer01\tPublished Spacer" > /Volumes/Data_HD/Mycobacterium/spoligoCount.txt
for i in *; do
#cd ./${i}/Zips
#7za x *_R1*
    os1=`egrep -c "$observedSpacer01" $i`
    echo "$i first done"
    os2=`egrep -c "$observedSpacer02" $i`
    echo "$i second done"
    ps1=`egrep -c "$publishedSpacer01" $i`
    echo "$i third done"
    ps2=`egrep -c "$publishedSpacer02" $i`
    echo "$i fourth done"

echo "$i\t$os1\t$os2\t$ps1\t$ps2" >> /Volumes/Data_HD/Mycobacterium/spoligoCount.txt

#rm *.fastq
#cd $workingDir
done

#
#  Created by Stuber, Tod P - APHIS on 11/20/2013.
#
