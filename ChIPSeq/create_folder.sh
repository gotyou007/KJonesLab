#!/bin/bash
proj=20828X
for i in {2,3,6};
do
    mkdir -p $proj$i
    mv $proj${i}*.gz $proj$i # move all files starting with 20828X$i to 20828X$i folder
    sed "s/${proj}./$proj$i/" cmd.txt > 20828X$i/cmd.txt # replace current_name with new name
done

