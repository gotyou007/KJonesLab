#!/bin/bash
awk 'BEGIN{OFS="\t";RS="\r\n"}{print $1" "$2}' samples.txt | while read old new; 
#add RS="\r\n" to suppress \r in shell; if you don't have \r, you can remove RS="\r\n"
do
   ls ${old}_* | while read file;
   do
      mv $file ${file/$old/$new}
   done
done


