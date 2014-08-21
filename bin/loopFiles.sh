#!/bin/bash

#  bovis_loopFiles.sh
#  Working directory must have Brucella or Bovis zip files
#  Brucelle files must begin with a "B"
#  Bovis files must NOT begin with a "B"

# grab the local configuration
INSTALL_ROOT=`dirname $0`
source "$INSTALL_ROOT/SNP_analysis_conf.sh"

if [ "$SEND_FILE_TO_NETWORK" = true ]; then
    date >> /Users/Shared/_WGS/coverageReport.txt
    date > /Users/Shared/_WGS/dailyReport.txt

    # Reset spoligo check file
    echo "" > /Users/Shared/_WGS/spoligoCheck.txt
    echo "WG Spoligo Check" >> /Users/Shared/_WGS/spoligoCheck.txt
fi

#Reset file
dateFile=`date "+%Y%m%d"`
printf "%s\t%s\n" "TB Number" "Octal Code" > "/Volumes/Data_HD/Mycobacterium/Analysis_new/${dateFile}_FileMakerSpoligoImport.txt"

if [ "$SEND_FILE_TO_NETWORK" = true ]; then
    echo "Isolate Total_Bases AveDep %>Q15" | awk '{ printf("%-30s %-12s %-10s %-10s\n", $1, $2, $3, $4) }' >> /Users/Shared/_WGS/dailyReport.txt

    if [ "$SEND_FILE_TO_NETWORK" = true ]; then
        date > /Users/Shared/_WGS/dailyStats.txt
    fi
fi

currentdir=`pwd`

for i in *.fastq.gz; do

    n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
    echo "n is : $n"

    mkdir -p $n
    mv $i $n/

done

for f in *; do
    cd $currentdir
    t=`echo $f | sed 's/\(^.\).*/\1/'`

    if [ $t == B ]; then
        #  Files beginning with a B are process as Brucella
        echo "***$f has been started as Brucella"

        ###
        cd ./$f
        mkdir ./temp
        cp *R1*.fastq.gz ./temp

        `gunzip ./temp/*R1*.fastq.gz && bruc_oligo_identifier.sh ./temp/*R1*.fastq | tee tee_bruc_oligo_identifier_out1.txt` &

        #`cd ./$f; /Users/Shared/_programs/_my_scripts/bruc_oligo_identifier.sh $R1; cd ..` &

        ###

        #`cd ./$f; /Users/Shared/_programs/_my_scripts/brucsuis_processZips.sh; cd ..` &
    elif [ $t == T ]; then
        #  Files beginning with a T are process as Taylorella
        echo "***$f has been started as Taylorella"
        `cd ./$f; /Users/Shared/_programs/_my_scripts/tayl_processZips.sh; cd ..` &

    elif [ $t == P ]; then
        #  Files beginning with a P are process as Pastuerella
        echo "***$f has been started as Pastuerella"
        `cd ./$f; /Users/Shared/_programs/_my_scripts/past_processZips.sh; cd ..` &


        #elif [ $t == C ]

        #  Files beginning with a P are process as Pastuerella
        #then

        #echo "***$f has been started as Pastuerella"

        #`cd ./$f; /Users/Shared/_programs/_my_scripts/mem_campy_processZips.sh; cd ..` &

        #  Files NOT beginning with a B are process as Bovis

    else
        echo "***$f has been started as Bovis"
        `cd ./$f; processZips.sh bovis; cd ..` &
        #`cd ./$f; /Users/Shared/_programs/_my_scripts/bovis_processZips.sh; cd..` &
    fi

done


#
#  Created by Tod Stuber on 11/05/12.
#
