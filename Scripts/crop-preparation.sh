#/usr/bin/bash

diameter=$1 # use 85 if no idea
prolatefile=$2 # use 85 if no idea

echo $diameter
echo $prolatefile

i=0
rm -rf temp* t-*

for file in *.txt
do
    let i=$i+1
    # name without extension
    name=${file%\.*}.mha
    outputname=${file%\.*}-c.mha
    
    echo $name
    cpstk crop -i $name -pr $prolatefile -d $diameter -o $outputname;
done
