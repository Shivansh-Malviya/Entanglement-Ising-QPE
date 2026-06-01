set datafile separator whitespace
set xlabel "lambda = J/h"
set ylabel "Concurrence"
set grid

data_dir = "../results/sparse_trace/"
chain_sizes = "4 6 8 10 12 14 16 18 20"

plot for [n in chain_sizes] data_dir."sdt_con".n.".dat" using 1:2 with lines title "N=".n
