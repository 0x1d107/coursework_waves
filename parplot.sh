rm -rf frames
mkdir -p frames
#el=$(grep -c '^$' waves.data )
N=5000
seq 1 10 $N |awk '{print $0"\n"NR}' |parallel -n 2 --bar gnuplot -c parplot.gpl 
ffmpeg -y -f image2 -r 30 -i 'frames/frame_%03d.jpg' plot.mp4
