set term png nocrop enhanced size 1280,720
set out './mach2d_output/mach2d_SIM_ID_residual.png'
set style data linespoints
set grid
set logscale y
set xlabel 'iteractions'
set ylabel 'norm L1 of the total residual (dimensionless)'
set time
set title 'Residual of the linear systems'

plot './mach2d_output/mach2d_SIM_ID_residual.dat' using 1:2 t''

