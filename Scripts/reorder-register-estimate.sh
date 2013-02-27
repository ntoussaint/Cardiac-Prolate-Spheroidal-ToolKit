#/usr/bin/bash

numberoflevels=$1 # use 3   if no idea
factor=$2         # use 0.3 if no idea
anatomy=$3        # use 'anatomy.mha' if no idea
i=0
rm -rf temp* t-*

for file in *.txt
do
    let i=$i+1
    # name without extension
    name_without_ext=${file%\.*}
    name=${file%\.*}-c.mha
    gradname=${file%\.*}.txt
    outputname=${file%\.*}-tensors.mha
    
    echo $name
    echo $gradname
    
    # reorder the DWIs : this will average all DWIs present
    # and remove the actual B-0 image. The averaged DWI (mean-diffusivity)
    # will be put as the first image, as if it was the B-0, 
    # the rest of the DWIs will be untouched.
    # The constructer calculated mean-diffusivity image(s) will be removed
    cpstk reorder -i $name -g $gradname -o temp -t 2 -f $factor
    cp -f $name $name_without_ext-before.mha
    # Actually reconstruct the tensors.
    ttk estimate -i temp.mha -g temp.grad -b 0 -o $name_without_ext-tensors-initial.mha
    cp -f $name_without_ext-tensors-initial.mha $name_without_ext-tensors.mha
    # Unstack this 4D image to recover the mean-diffusivity image 
    # and store the DWIs as 'original'
    unstackimage -i $name -o original-dwi-
    
    # first iteration : put the arithmetic mean of all DWIs 
    # in front, it will be the reference for 1-to-1 registration
    cpstk reorder -i $name -g $gradname -o temp -t 1 -a 1.0
    
    # Multi-Level registration Process
    for level in $(seq 1 1 $numberoflevels)
    do	
	# Unstack this 4D image to recover the mean-diffusivity image
	unstackimage -i temp.mha -o t-
	# Take the mean-diffusivity image as reference for registration
	cp -f t-000.mha reference.mha
	j=0
	inputs=''
	for dwi in original-dwi-*.mha
	do
	    if [ "$level" -eq 1 ]
	    then
	        # Register each of the 'original' DWI to the j-est reference
		slicetosliceregistration -f reference.mha -m $dwi -o temp-$j -mt 0 -d 0
	    else
	        # Register each of the 'original' DWI to the j-est reference, using the previous translation as initialization
		slicetosliceregistration -f reference.mha -m $dwi -o temp-$j -t temp-$j.mat -mt 0 -d 0
	    fi
	    export inputs=${inputs}"\n"temp-$j.mha
	    let j=$j+1
	done
        
	# Stack back all registered DWIs images
	cp -f original-dwi-*00.mha temp-0.mha
	echo -e $j $inputs > input.blist
	stackimage -i input.blist -o temp.mha
	# This is now the 4D image to use for tensor estimation
	cp temp.mha $name_without_ext-after.mha

	# Put the new j-est mean-diffusivity in front of the stack
	cpstk reorder -i temp.mha -g temp.grad -o geom -t 2 -f $factor
	# Actually reconstruct the tensors.
	ttk estimate -i geom.mha -g geom.grad -b 0 -o $name_without_ext-tensors.mha

	# Put the new j-est mean-diffusivity in front of the stack
	cpstk reorder -i temp.mha -g temp.grad -o temp -t 1 -a 1.0
	# remove temporary files
	rm -rf t-* geom* temp-increased.mha input.blist
	
    done

    cpstk res-image -i $anatomy -r original-dwi-*00.mha -o reference.mha
    slicetosliceregistration -f reference.mha -m original-dwi-*00.mha -o $name_without_ext-transformed -mt 1 -d 1
    cpstk apply-tensors -i $name_without_ext-tensors.mha -t $name_without_ext-transformed.mat -o $name_without_ext-tensors-tr.mha
    #cp -f $name_without_ext-tensors.mha $name_without_ext-tensors-tr.mha
    rm -rf temp* original-dwi* input.blist reference.mha
done