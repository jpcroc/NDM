#!/usr/bin/gnuplot

reset

# png
set terminal pngcairo size 640,480 enhanced font 'Verdana,10'

unset xtics
#unset xlabel
set xlabel "Reaction coordinate"
set ylabel "Energy (eV)"
set xrange [0:6]
set yrange [-0.05:0.7]

do for [ii=0:6] {
	outfile = sprintf('nrj%02.0f.png',ii)
   set output outfile
   if (ii == 6) {
      set arrow from 4,0 to 4,0.65 heads filled lw 2
      set key 5.2,0.3 title "E_m = 0.66 eV" font "Verdana-Bold,12" 
   }
   plot 'nrj.dat' u 1:($2-$3)-(-1137.97909-0.02914) index ii notitle w lp lw 2
}
