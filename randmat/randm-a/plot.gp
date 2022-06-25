# run: gnuplot plot.gp
#
#
set title "Randomly generated input data with matrix condition number 1e8"
#
set xlabel "Number of unknowns"
set ylabel "Time in seconds"
#
set terminal postscript eps enhanced color solid colortext 9
set output 'randm-a.eps'

# Now plot the data with lines and points
plot 'plot.data' using 2:4 w lp title 'GNU gama', \
     ''          using 2:3 w lp title 'ICGS', \
     ''          using 2:5 w lp title 'Fortran mtgso1'