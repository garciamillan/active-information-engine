set term epslatex standalone color size 4.8,3 header "\\usepackage{amsmath,amssymb}"
set out "fig_RTPtraj13_v2.tex"


FILE00 ="<grep POS 'test13_trajRTP.dat'"
tf = 15
set mxtics 5
a = 2

if (!exists("MP_LEFT"))   MP_LEFT = .15
if (!exists("MP_RIGHT"))  MP_RIGHT = 0.9
if (!exists("MP_BOTTOM")) MP_BOTTOM = .15
if (!exists("MP_TOP"))    MP_TOP = 0.95
if (!exists("MP_xGAP"))   MP_xGAP = 0.12
if (!exists("MP_yGAP"))   MP_yGAP = 0.01

set multiplot layout 3,1 margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_xGAP, MP_yGAP


set key samplen 2 
set key spacing 1.3
#self-propulsion
set ytics -0.5,0.5,1
set ylabel "self-propulsion" offset 1,0
set xtics ("" 0,"" 5,"" 10,"" 15,"" 20)
plot[0:tf][-0.8:0.8] \
	0. lc 'gray' t '', \
	FILE00 u ($3/10.):9 w l lc 3 lw a t "$v(t)$",\
	FILE00 u ($3/10.):9 every 3333::1 lc -1 pt 2 lw 3 t ''

#position
set ylabel 'position' offset 1,0
set key bottom right
set ytics (0, 0.5, 1)  #-0.5,0.25,0.75
unset yrange
plot[0:tf][-0.5:1.5] \
	x<11?0.:1/0 lc 'gray' t '', \
	FILE00 u ($3/10.):6 w l lw a t "$x(t)$", \
	FILE00 u ($3/10.):5 w l lc 'red' lw a t "$\\lambda(t)$"


FILE00 ="<grep WORK_C 'test13_trajRTP.dat'"

#accumulated work
set xlabel "$t/t_{\\text{f}}$" offset 0,0.75
set ylabel 'work' offset 1,0
set k top right
unset ytics
set ytics 0.1
set xtics 5
plot[0:tf][] \
	x<12?0.:1/0 lc 'gray' t '', \
	FILE00 u ($3/10.):($5+$6) w l lc 4 lw a t "$W(t)$", \
        -0.02005479986136146*x lc 'blue' t " $\\langle{W}\\rangle t$"
        


#unset multiplot

