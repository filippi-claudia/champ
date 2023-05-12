set key center center
set xlabel "Optimization Step"
set ylabel "Total Energy (Ha)"
set yrange[-26.52:-26.30]

plot 'E_err_steps05_nblk1000_core_1' u :1:2 w errorbars t "Core=1, Steps=5, nblk=1000" lt 1 pt 7 ps 0.5,\
     'E_err_steps05_nblk250_core_1'  u :1:2 w errorbars t "Core=1, Steps=5, nblk= 250" lt 7 pt 7 ps 0.5,\
     'E_err_steps05_nblk250_core_2'  u :1:2 w errorbars t "Core=2, Steps=5, nblk= 250" lt 3 pt 7 ps 0.5,\
     'E_err_steps05_nblk500_core_1'  u :1:2 w errorbars t "Core=1, Steps=5, nblk= 500" lt 4 pt 7 ps 0.5,\
     'E_err_steps10_nblk1000_core_1' u :1:2 w errorbars t "Core=1, Steps=10,nblk=1000" lt 8 pt 7 ps 0.5
#    your data
#    'my_data' u :1:2 w errorbars t "my data" lt 2 pt 1 ps 1

