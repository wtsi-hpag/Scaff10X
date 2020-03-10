#!/bin/bash
function cleanFile
{
	FILE=$1
	shift
	cat $FILE | awk '{print $2 "\t" $3 }' | egrep -v "^0	"
}

#cleanFile human-hic-len.freq > bE378K21-screen.dat1.cleaned
#cleanFile hummingbird-hic.freq > bE378K21-screen.dat2.cleaned
#cleanFile fMasArm1-hic.freq > bE378K21-screen.dat3.cleaned
#cleanFile tdevil-hic.freq > bE378K21-screen.dat4.cleaned
#cleanFile human-mp15-hic.freq > bE378K21-screen.dat4.cleaned 

function plotcmd
{
        printf "set logscale x\n"
        printf "set logscale y\n"
	printf "set terminal svg\n"
        printf "set style line 1 lt 1 lw 3 pt 3 linecolor rgb \"red\"\n"
        printf "set style line 2 lt 1 lw 3 pt 3 linecolor rgb \"green\"\n"
        printf "set style line 3 lt 1 lw 3 pt 3 linecolor rgb \"blue\"\n"
        printf "set style line 4 lt 1 lw 3 pt 3 linecolor rgb \"violet\"\n"
        printf "set style line 5 lt 1 lw 3 pt 3 linecolor rgb \"cyan\"\n"
	printf "set xlabel \"Molecular Length (Minimum 100bp and 5 reads)\"\n"
	printf "set ylabel \"Frequency / 100\"\n"
        printf "plot [ 100 to 1000000 ] [ 0.001 to 10.0 ] \"human-bcl-5.freq\" title \"Human-10X\" with lines ls 1,\"hummingbird-bcl-5.freq\" title \"Hummingbird-10X\" with lines ls 2,\"fAnaTes1-bcl-5.freq\" title \"Fish fAnaTes1-10X\" with lines ls 3,\"fSimDai1-bcl-5.freq\" title \"Fish fSimDai1-10X\" with lines ls 4,\"sample-10x.freq\" title \"Test sample-10X\" with lines ls 5"
}

plotcmd | gnuplot > plot10x.svg
inkscape -z --export-text-to-path --export-pdf plot10x.pdf plot10x.svg
gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=plot10x.png plot10x.pdf

#rm -f bE378K21-raw.dat.cleaned bE378K21-screen.dat.cleaned data.svg
