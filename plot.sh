meson compile -C build
build/waves
ffmpeg -i 'data/%04d.pgm' -r 24 -s 720x720 -y plot.mp4 
rm -r data
#bash parplot.sh
#gnuplot waves.gpl
#xdg-open waves.gif
