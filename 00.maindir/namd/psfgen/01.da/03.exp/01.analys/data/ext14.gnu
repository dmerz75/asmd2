# Gnuplot script file for plotting data in file
# extension files "gnu.dat"

set terminal postscript eps color lw 15 "Helvetica" 20
set out 'ext14.eps'
set title "Chignolin-extension"
set xlabel "Timestep"
set ylabel "ext. (A)"

# set margin 1
set autoscale 
#set border lw 0.1

#set xr [1:7]
#set yr [3.5:8.2]
#set xtic .5
#set ytic .5


# set style data lines
set style line 1 lt 1 lw .1
set style line 2 lt 1 lw .1
set style line 3 lt 1 lw .1
set style line 4 lt 1 lw .1

unset log                              # remove any log-scaling
unset label                            # remove any previous labels

# EXTRAS
# set key 0.01,100
# set label "Yield Point" at 0.003,260
# set arrow from 0.0028,250 to 0.003,280

plot    "extCA1-4.dat" using 1:2 t 'atom 1-10', \
        "extCA1-4.dat" using 1:3 t 'atom 2-9', \
        "extCA1-4.dat" using 1:4 t 'atom 3-8', \
        "extCA1-4.dat" using 1:5 t 'atom 4-7'
