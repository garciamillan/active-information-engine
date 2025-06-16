set term epslatex standalone color size 6.4,4 header "\\usepackage{amsmath,amssymb}"
set out "fig_AOUPtraj12.tex"

#set form y "$10^{%L}$"

FILE00 ="<grep POS 'test12_trajAOUP.dat'"
tf = 15
#set zeroaxis
#unset xlabel
#set xtics ("" 0, "" 5, "" 10, "" 15, "" 20)
set mxtics 5

if (!exists("MP_LEFT"))   MP_LEFT = .15
if (!exists("MP_RIGHT"))  MP_RIGHT = 0.9
if (!exists("MP_BOTTOM")) MP_BOTTOM = .15
if (!exists("MP_TOP"))    MP_TOP = 0.95
if (!exists("MP_xGAP"))   MP_xGAP = 0.12
if (!exists("MP_yGAP"))   MP_yGAP = 0.01

set multiplot layout 3,1 margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_xGAP, MP_yGAP
#set lmargin at screen 0.15


set key samplen 2 
set key spacing 1.3
#self-propulsion
set ytics -0.5,0.5,1
set ylabel "self-propulsion"
set xtics ("" 0,"" 5,"" 10,"" 15,"" 20)
plot[0:tf][-0.8:0.8] \
	0. lc 'gray' t '', \
	FILE00 u ($3/10.):9 w l lc 3 t "$v(t)$",\
	FILE00 u ($3/10.):9 every 3333::1 lc -1 pt 2 lw 2 t ''

#position
set ylabel 'position' offset 1,0
#set key at 19.5,0.37 #opaque
set key bottom right
set ytics 0.5  #-0.5,0.25,0.75
unset yrange
plot[0:tf][] \
	x<12?0.:1/0 lc 'gray' t '', \
	FILE00 u ($3/10.):6 w l t "$x(t)$", \
	FILE00 u ($3/10.):7 w l lc 1 t "$\\langle x(t) | \\{v(nt_{\\text{f}})\\}_{n\\in\\mathbb{N}}\\rangle$", \
	FILE00 u ($3/10.):5 w l t "$\\lambda(t)$" lc 'red'


FILE00 ="<grep WORK_C 'test12_trajAOUP.dat'"

#accumulated work
set xlabel "$t/t_{\\text{f}}$"
set ylabel 'work' offset 1,0
set k top right
#set k at 20,0.05 #opaque
#set ytics -0.1,0.05,0.05
unset ytics
set ytics 0.1
set xtics 5
plot[0:tf][] \
	x<12?0.:1/0 lc 'gray' t '', \
	FILE00 u ($3/10.):($5+$6) w l t "$W(t)$" lc 4, \
	FILE00 u ($3/10.):7 w l lc -1 t "$\\langle W(t) | \\{v(nt_{\\text{f}})\\}_{n\\in\\mathbb{N}} \\rangle$", \
        -0.02005479986136146*x lc 'blue' t " $\\mathcal{W} t$"
        
# old engine -0.004999999989694233*x lc 'blue' t " $\\mathcal{W} t$"


#unset multiplot

