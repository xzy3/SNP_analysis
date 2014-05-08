#!/bin/sh

#  email_loopFiles.sh

`loopFiles.sh` &&

echo "e-mailing files"

echo "" >> /Users/Shared/_WGS/dailyReport.txt
cat /Users/Shared/_WGS/spoligoCheck.txt >> /Users/Shared/_WGS/dailyReport.txt

cat /Users/Shared/_WGS/dailyStats.txt >> /Users/Shared/_WGS/dailyReport.txt
echo "" >> /Users/Shared/_WGS/dailyReport.txt

grep -v "*" /Users/Shared/_WGS/dailyReport.txt > /Users/Shared/_WGS/dailyReport2.txt

grep -v "Stats for BAM file" /Users/Shared/_WGS/dailyReport2.txt > /Users/Shared/_WGS/dailyReport3.txt

mail -s "WGS results" tod.p.stuber@aphis.usda.gov < /Users/Shared/_WGS/dailyReport3.txt

mail -s "WGS results" patrick.m.camp@aphis.usda.gov < /Users/Shared/_WGS/dailyReport3.txt

mail -s "WGS results" suelee.robbe-austerman@aphis.usda.gov < /Users/Shared/_WGS/dailyReport3.txt

mail -s "WGS results" Christine.R.Quance@aphis.usda.gov < /Users/Shared/_WGS/dailyReport3.txt

mail -s "WGS results" David.T.Farrell@aphis.usda.gov < /Users/Shared/_WGS/dailyReport3.txt

rm /Users/Shared/_WGS/dailyReport2.txt
rm /Users/Shared/_WGS/dailyReport3.txt

#
#  Created by Tod Stuber on 11/09/12.
#
