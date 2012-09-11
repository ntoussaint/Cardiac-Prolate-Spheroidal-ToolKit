#/usr/bin/bash

mkdir processing
mv *after* *before* *tensors*.mha *transformed* ./processing
mv domain-ellipsoid* coeff* *-landmarks* *-function* *-mesh* *-c.mha ./processing
mv domain-warped-to-ellipsoid* metricvalues.csv tensors.list ./processing
