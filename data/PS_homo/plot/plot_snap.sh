#plot_snap.sh
gfortran -o separate.e separate.f num2char.f -g
./separate.e<<EOF
'../s0101'
'../s0201'
256 256 5
EOF
gnuplot gnu_snapplot
mv snap.eps snap_ps.eps

./separate.e<<EOF
'../s1101'
'../s1201'
256 256 5
EOF
gnuplot gnu_snapplot
mv snap.eps snap_p.eps

./separate.e<<EOF
'../s2101'
'../s2201'
256 256 5
EOF
gnuplot gnu_snapplot
mv snap.eps snap_s.eps

rm snap_h_* snap_v_* 
echo Snaps plotting, Done!
evince snap_ps.eps&
evince snap_p.eps&
evince snap_s.eps&


