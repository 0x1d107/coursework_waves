#!/bin/bash
rm -rf build
meson setup build -Dbuildtype=release
# 8 16 32 64 100 128 150 200 250 300
for RR in $(seq 0 0.01 0.12)
do
meson configure build -Dcpp_args=-DRR=$RR
bash plot.sh detector
mv sensor sensor_$RR
mv plot_detector.mp4 plot_$RR.mp4
done
