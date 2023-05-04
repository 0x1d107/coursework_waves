#!/bin/bash
ninja -C build
for N in $(seq 10 10 300)
do
    echo $N $(build/simple $(( $N * $N )) $N| cut -d= -f2)
done|tee conv.data
