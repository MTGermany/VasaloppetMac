

##########################################################
# Line- und Punkttypen 
##########################################################

# MT 2020-03: Nun automatisch in ~/.gnuplot initfile

# set style line: linetype lt, point type pt

#gnuplot Ubuntu 18, Mai2019, Punkttypen
# pt 1:  +
# pt 2:  X
# pt 3:  hairs
# pt 4:  openSquare
# pt 5:  closedSquare
# pt 6:  ring
# pt 7:  bullet
# pt 8:  open upTriang
# pt 9:  closed upTriang
# pt 10: open downTriang
# pt 11: closed downTriang

#####################################################################
# dashed now explicit with dt (implicit over "dashed" + ls no longer works)
# BUG (2021, Ubuntu20): dt does not work with png (there always solid)
######################################################################

# post eps dashed no longer works but dashtype (dt) in ls
# specs! dt 1 = solid
 
set style line 1 lt 1 lw 3 pt 7 ps 1.9 dt 3 lc rgb "#000000" #schwarz, bullet
set style line 2 lt 1 lw 3 pt 9 ps 1.5 dt 3 lc rgb "#CC0022" #closedUpTriang
set style line 3 lt 1 lw 3 pt 4 ps 1.2 dt 3 lc rgb "#FF3300" #closedSquare
set style line 4 lt 1 lw 3 pt 4 ps 1.5 dt 3 lc rgb "#FFAA00" #gelb,
set style line 5 lt 1 lw 3 pt 5 ps 1.5 dt 3 lc rgb "#00DD22" #gruen,
set style line 6 lt 1 lw 3 pt 4 ps 1.5 dt 3 lc rgb "#00AAAA"
set style line 7 lt 1 lw 3 pt 4 ps 2.0 dt 3 lc rgb "#1100FF" #blau,
set style line 8 lt 1 lw 3 pt 8 ps 1.5 dt 3 lc rgb "#220088"
set style line 9 lt 1 lw 3 pt 9 ps 1.5 dt 3 lc rgb "#999999" #grau

set style line 11 lt 1 lw 6 pt 7 ps 1.9  lc rgb "#000000" 
set style line 12 lt 1 lw 6 pt 2 ps 1.5 dt 2 lc rgb "#CC0022" 
set style line 13 lt 8 lw 6 pt 4 ps 1.2  lc rgb "#FF3300"
set style line 14 lt 6 lw 6 pt 4 ps 1.5  lc rgb "#FFAA00"
set style line 15 lt 1 lw 6 pt 5 ps 1.5  lc rgb "#00DD22"
set style line 16 lt 5 lw 6 pt 7 ps 1.5  lc rgb "#00AAAA"
set style line 17 lt 1 lw 6 pt 7 ps 1.5  lc rgb "#1100FF"
set style line 18 lt 4 lw 6 pt 8 ps 1.5  lc rgb "#220088"
set style line 19 lt 7 lw 6 pt 9 ps 1.5  lc rgb "#999999"



##############################################################
# Beispiele fuer Funktionen 
##############################################################

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
filterRange(data,xmin,xmax)=((data>=xmin)&&(data<=xmax)) ? data : NaN
filterEndDigit(data,digit)=((floor(data))%digit==0) ? 1 : NaN



##############################################################
# Defined clippings/crop (to eliminate superfluous white space/make white
# space/make plotting area exactly rectangular,  to get identical margins for multiplot
# etc
##############################################################

#set lmargin at screen 0.15  # circumvent gnuplot clipping bugs
#set rmargin at screen 0.98
#set bmargin at screen 0.2
#set tmargin at screen 1.3



set term post eps enhanced color solid "Helvetica" 16


##############################################################
set out "sim2_traj.eps"; print "plotting sim2_traj.eps"
##############################################################


set key right bottom
set size 1,0.8

set xlabel "Time (after race start) [s]"
set ylabel "Distance [km]"

xmin=0
xmax=4000

set xrange [0:8000]
set yrange [0:4]

# update athlete indices "n=..." according to title in .traj file if changed

plot\
  "sim2.traj" u ($1):(0.001*filterRange($2,xmin,xmax)) t "n=100" w l ls 1,\
  "sim2.traj" u ($1):(0.001*filterRange($3,xmin,xmax)) t "n=500" w l ls 2,\
  "sim2.traj" u ($1):(0.001*filterRange($4,xmin,xmax)) t "n=1500" w l ls 3,\
  "sim2.traj" u ($1):(0.001*filterRange($5,xmin,xmax)) t "n=2500" w l ls 4,\
  "sim2.traj" u ($1):(0.001*filterRange($6,xmin,xmax)) t "n=4000" w l ls 5,\
  "sim2.traj" u ($1):(0.001*filterRange($7,xmin,xmax)) t "n=6000" w l ls 6,\
  "sim2.traj" u ($1):(0.001*filterRange($8,xmin,xmax)) t "n=8000" w l ls 7,\
  "sim2.traj" u ($1):(0.001*filterRange($9,xmin,xmax)) t"n=10000" w l ls 8,\
  "sim2.traj" u ($1):(0.001*filterRange($10,xmin,xmax))t"n=12000" w l ls 1,\
  "sim2.traj" u ($1):(0.001*filterRange($11,xmin,xmax))t"n=15000" w l ls 2

##############################################################
set out "sim2_HoegstaPunkten.eps"; print "plotting sim2_HoegstaPunkten.eps"
##############################################################

set key right top

set xlabel "Time [s]"
set ylabel "Flow [skiers/h]"
set auto y
plot\
  "sim2.HoegstaPunktenCalib" u 1:3 t "Data" w l ls 17,\
  "sim2.x2980" u 1:2 t "Simulation" w l ls 16

##############################################################
set out "sim2_HoegstaPunkten_ncum.eps"; print "plotting sim2_HoegstaPunkten_ncum.eps"
##############################################################

set key right bottom

set xlabel "Time [s]"
set ylabel "cumulated n-curve"
set auto y
plot\
  "sim2.HoegstaPunktenCalib" u 1:4 t "Data" w l ls 17,\
  "sim2.x2980" u 1:4 t "Simulation" w l ls 16



####################################################################
# contour prep
####################################################################

unset contour        # no contour lines
set contour surface  # Aktiviert Kontourlinien auf 3D-Flaeche
unset clabel         # dann lauter gleiche Kontourlinien; 
                     # Farbe und Typ mit "w l ls" beim splot-Kommando
                     # Achtung bug! explizit ("w l lt 3 lw 1" etc) gibt DOS!


set cntrparam bspline 

set cntrparam levels 12                    # n equidistant levels from min-max
#set cntrparam levels incr -10,1,-2        # min, increment, max
#set cntrparam levels discrete -8,-1.6,0,1 # freely set lines

set palette defined ( 0 "#dd00ff", 5 "blue", \
      10 "green", 15 "yellow", 30 "orange", 50 "red", 100 "#880000") 


set style line 99 lt 7 lw 0 linecolor rgb "#000000" # beliebige Farben:Schwarz
unset pm3d                                          # no color coding
set pm3d; set pm3d  map                             # color coded xy surface

unset surface        # pm3d off: kein Gitternetz gemacht
                     # pm3d on: Schnelles Plotten mit/ohne Netz
                     # (je nach pm3d Einstellung) mit Artefakten

set nogrid           # to avoid plotting an 2d grid on the xy surface


##############################################################
set out "sim2_rho.eps"; print "plotting sim2_rho.eps"
##############################################################

set size 1,0.8
set xlabel "Time (after race start) [s]"

set ylabel "Distance [km]"
set yrange [0:4]

set label 1 "Total Density [skiers/m]" at screen 0.75,0.78
set lmargin at screen 0.10  # circumvent gnuplot clipping bugs
set rmargin at screen 0.87
set bmargin at screen 0.12
set tmargin at screen 0.74


splot  "sim2.macro" u ($1):(0.001*$2):($3) t ""  w l ls 99




