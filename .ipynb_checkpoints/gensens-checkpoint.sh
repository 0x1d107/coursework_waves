#!/bin/bash
rm -rf build
meson setup build -Dbuildtype=release
for RR in 32 64 128 
do
meson configure build -Dcpp_args=-DRR=1.0/$RR
bash plot.sh detector
mv sensor sensor_$RR
mv plot_detector.mp4 plot_$RR.mp4
done
