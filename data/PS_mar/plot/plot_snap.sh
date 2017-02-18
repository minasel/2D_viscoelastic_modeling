#plot_snap.sh
input='../s0101'
output=s0101.eps
clip=0.01
gnuplot -e "inputname='${input}';outputname='${output}';titlename='P-wave stress';r='${clip}'" gnu_snap

input='../s0201'
output=s0201.eps
clip=0.01
gnuplot -e "inputname='${input}';outputname='${output}';titlename='P-wave stress';r='${clip}'" gnu_snap

input='../s1101'
output=s1101.eps
clip=0.01
gnuplot -e "inputname='${input}';outputname='${output}';titlename='P-wave stress';r='${clip}'" gnu_snap

input='../s1201'
output=s1201.eps
clip=0.01
gnuplot -e "inputname='${input}';outputname='${output}';titlename='P-wave stress';r='${clip}'" gnu_snap


input='../s2101'
output=s2101.eps
clip=0.01
gnuplot -e "inputname='${input}';outputname='${output}';titlename='P-wave stress';r='${clip}'" gnu_snap

input='../s2201'
output=s2201.eps
clip=0.01
gnuplot -e "inputname='${input}';outputname='${output}';titlename='P-wave stress';r='${clip}'" gnu_snap




