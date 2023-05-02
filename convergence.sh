#!/bin/bash
ninja -C build
for M in $(seq 100 1000 11000)
do
    echo $M $(build/simple $M| cut -d= -f2)
done|tee conv.data
