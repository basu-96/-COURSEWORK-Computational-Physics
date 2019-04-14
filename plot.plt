set terminal png
set output 'tb_2d_DOS_900_sites.png'
set xlabel 'Energy'
set ylabel 'DOS'
set title '2D TB with 900 sites'
#set obj 1 rect at 6,0.5 size 5,0.2
#set label 1 at 4,0.55 "Lattice size : 30x50 \n k_B : 1 \n Number of sweeps \n #for each temperature : 500" left offset 0,.5
plot 'tb_2d_dos.dat' notitle ps 0.01
