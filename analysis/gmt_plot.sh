#!/bin/bash
set -e

# change parameters here
infile=test.nc
outfile=test.jpg
var=thk
cbartitle=Thickness
cbarunit=m
cbar=0/2000/100
#range=-124.8/-122.6/47/48.4
range=-124.5/-122.75/47.25/48.25

# set up the option parser
while [[ $# -gt 1 ]]
do
    key="$1"
    case $key in
        -i|--infile)
            infile="$2"
            shift
        ;;
        -o|--outfile)
            outfile="$2"
            shift
        ;;
        -v|--var)
            var="$2"
            shift
        ;;
        --cbar)
            cbar="$2"
            shift
        ;;
        --title)
            cbartitle="$2"
            shift
        ;;
        --unit)
            cbarunit="$2"
            shift
        ;;
    esac
shift
done

#python ../tools/nc2xyz.py -i $infile -o tmp.xyz -v $var
python ../tools/nc2gmt.py -i $infile -o tmp.xyz -v $var --srs "+init=epsg:26710" --interp True

#source activate gmt

gmt xyz2grd tmp.xyz -Gtmp.nc -R$range -I0.5m

gmt psbasemap -R$range -Jm6 -Ba0.5f0.25 -V -P -K -X1.5 -Y2 >> tmp.ps 
gmt makecpt -T$cbar -Cjet -Z >tmp.cpt
gmt grdimage tmp.nc -Ctmp.cpt -Jm -E300 -nb -Q -P -O -K >> tmp.ps
#gmt psscale -Dx15.5c/0.3c+w12c/0.5c -Ctmp.cpt -Baf -Bx+l$cbartitle -By+l$cbarunit -O >> tmp.ps
gmt psscale -Dx12.5c/0.11c+w8c/0.4c -Ctmp.cpt -Baf -Bx+l$cbartitle -By+l$cbarunit -O >> tmp.ps

#ps2pdf tmp.ps $outfile
gmt psconvert tmp.ps -A -Tj
mv tmp.jpg $outfile
rm gmt.history
rm tmp*

#source deactivate gmt

echo "GMT: plotting finished."
#exit

echo "moving to windows /jtlai/Work/Research/glacier/plot"
mv $outfile $WIN_HOME/Work/Research/glacier/plot
