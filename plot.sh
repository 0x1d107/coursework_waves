EXEC=${1:-waves}
rm -f sensor_*.data
meson compile -C build || exit 1
build/$EXEC || exit 1
ffmpeg -i 'data/%04d.pgm' -r 24 -s 720x720 -y plot_$EXEC.mp4 
rm -r data
#bash parplot.sh
#gnuplot waves.gpl
#xdg-open waves.gif
