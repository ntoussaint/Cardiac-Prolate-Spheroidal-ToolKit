#/usr/bin/bash

mkdir processing
mv *after* *before* *tensors*.mha ./processing
mv *transformed* *tensors.mha coeff* *initial* ./processing
mv domain-ellipsoid.vtk coeff* *-landmarks* *-function* *-mesh* *-c.mha ./processing
mv domain-warped-to-ellipsoid* metricvalues.csv tensors.list ./processing
