#!/bin/bash
rm -rf build
meson setup build -Dbuildtype=release
for RR in 8 16 32 64 100 128 150 200 250 300 
do
meson configure build -Dcpp_args=-DRR=1.0/$RR
bash plot.sh detector
mv sensor sensor_$RR
mv plot_detector.mp4 plot_$RR.mp4
done
