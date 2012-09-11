#/usr/bin/bash

crop-preparation.sh 85 prolatetransform.tr
reorder-register-estimate.sh 3 0.3 domain.mha

for file in *-tensors-tr.mha
do
    name_without_ext=${file%\.*}
    ttk-utils apply_mask -t 1 -i $file -m domain.mha -o $name_without_ext-m.mha
done

n=0
inputs=''
for file in *-tensors-tr-m.mha
do
    export inputs=${inputs}"\n"$file
    let n=$n+1
done
echo -e $n $inputs > tensors.list

cpstk itk2vtk -i tensors.list -o tensors.vtk

# non-linear-registration.sh

# ln -s ~/shells/optimal-kernels.csv kernels.csv

cardiac-dti-pipeline-cleanup.sh

#cpstk extrapolate -p 1 -i tensors.vtk -d domain.mha -f1 1-to-2.mha -f2 2-to-1.mha -pr prolatetransform.tr -k 0 -o tensors.mha -u kernels.csv
#cpstk extract -i tensors.vtk -pr prolatetransform.tr -f1 1-to-2.mha -f2 2-to-1.mha -t helix -o helix-angles.csv
#cpstk extract -i tensors.vtk -pr prolatetransform.tr -f1 1-to-2.mha -f2 2-to-1.mha -t transverse -o transverse-angles.csv
#cpstk extract -i tensors.vtk -pr prolatetransform.tr -f1 1-to-2.mha -f2 2-to-1.mha -t sheet -o sheet-angles.csv
#cpstk extract -i tensors.vtk -pr prolatetransform.tr -f1 1-to-2.mha -f2 2-to-1.mha -t cl-cp-cs -o cl-cp-cs.csv
#ttk tractography -i tensors.mha -fa1 0.2 -fa2 0.2 -s 0.2 -fs 2.0 -n 2 -o fibres.fib -c 1
#cpstk colorify -i fibres.fib -pr prolatetransform.tr -f1 1-to-2.mha -f2 2-to-1.mha -t helix -o fibres-helix.vtk
