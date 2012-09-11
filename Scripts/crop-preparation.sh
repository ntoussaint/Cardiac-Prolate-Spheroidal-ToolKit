#/usr/bin/bash

diameter=$1 # use 85 if no idea
prolatefile=$2 # use 'prolatetransform.tr' if no idea

echo $diameter
echo $prolatefile

for file in *.txt
do
    # name without extension
    name=${file%\.*}.mha
    outputname=${file%\.*}-c.mha
    
    echo $name
    cpstk crop -i $name -pr $prolatefile -d $diameter -o $outputname;
done
