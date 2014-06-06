#!/bin/bash

# Configuration file for SNP_analysis
# Author Seth Sims <xzy3@cdc.gov>

#   bwa, http://bio-bwa.sourceforge.net/bwa.shtml
BWA=`which bwa`
BWA_CORES='4'

#   samtools, http://samtools.sourceforge.net/samtools.shtml
SAMTOOLS=`which samtools`

#   picard, http://picard.sourceforge.net/command-line-overview.shtml
PICARD_TOOLS_ROOT='/Users/Shared/_programs/picard-tools-1.100/'

PICARD_BIG_JRE_SETTINGS='-Xmx4g'
PICARD_JRE_SETTINGS='-Xmx2g'

#   gatk, http://www.broadinstitute.org/gatk/
GENOME_ANALYSIS_TK='java -Xmx2g -jar /Users/Shared/_programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar'

#   bamtools
BAMTOOLS=`which bamtools`

# IGVTools
IGVTOOLS='java -Xmx2g -jar /Users/Shared/_programs/IGVTools/igvtools.jar'

# abyss-pe
ABYSS_PE='which abyss-pe'

# set defaults for variables that change the behavior of the script
SEND_FILES_TO_NETWORK=true

SNP_DATA_ROOT='/Volumes'


# load local configuration from user's home directory
if [ -e "$HOME/SNP_analysis_conf.sh" ] ; then
    source ~/SNP_analysis_conf.sh
fi

# the root and settings for picard should be ready now so make the variables for each individual
# tool.
PICARD_CreateSequenceDictonary="java $PICARD_BIG_JRE_SETTINGS -jar $PICARD_TOOLS_ROOT/CreateSequenceDictionary.jar"
PICARD_SamToFastq="java $PICARD_JRE_SETTINGS -jar $PICARD_TOOLS_ROOT/SamToFastq.jar"
PICARD_MarkDupicates="java $PICARD_JRE_SETTINGS -jar $PICARD_TOOLS_ROOT/MarkDuplicates.jar"
PICARD_QualityScoreDistribution="java $PICARD_JRE_SETTINGS -jar $PICARD_TOOLS_ROOT/QualityScoreDistribution.jar"
PICARD_CollectMultipleMetrics="java $PICARD_JRE_SETTINGS -jar $PICARD_TOOLS_ROOT/CollectMultipleMetrics.jar"
PICARD_CollectAlignmentSummaryMetrics="java $PICARD_JRE_SETTINGS -jar $PICARD_TOOLS_ROOT/CollectAlignmentSummaryMetrics.jar"
PICARD_CollectGcBiasMetrics="java $PICARD_JRE_SETTINGS -jar $PICARD_TOOLS_ROOT/CollectGcBiasMetrics.jar"
PICARD_CollectInsertSizeMetrics="java $PICARD_JRE_SETTINGS -jar $PICARD_TOOLS_ROOT/CollectInsertSizeMetrics.jar"
