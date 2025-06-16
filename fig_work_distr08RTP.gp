set term epslatex standalone color size 3.7,2.24
set out "fig_work_distr08RTP.tex"


ymax = 200
set arrow nohead dt '.' from -0.02005479986136146,0 to -0.02005479986136146,ymax

#WMOMEN  RTP  1000000  -20054.79986136146  1120.753374407114  -32.72739597626321  1.969129568901121  -0.05858811626075406  
#WMOMCENTRAL  RTP  1  0
#WMOMCENTRAL  RTP  2  0.0007185583769278499
#WMOMCENTRAL  RTP  3  1.857017751740248E-05
#WMOMCENTRAL  RTP  4  1.563050019326417E-06
#WMOMCENTRAL  RTP  5  8.465940492504699E-08

#set zeroaxis

set xlabel "$W$"
set ylabel "$\\mathcal{P}(W)$"

#set k top left
set xtics 0.01
set mxtics 2

plot [-0.05:-0.015][0:ymax]\
        "<grep WDISTR 'test08.dat'" u 4:($6/$5) w l lw 2 lc 2 t 'RTP'#, \
	'test03_04_MOMS.dat' u 3:2:4 w xerrorbar t ''





