#!/bin/bash

python ../tools/nc2xyz.py -i test.nc -o tmp.xyz -v topg

source activate gmt

outfile=test.jpg

xyz2grd tmp.xyz -Gtmp.nc -R-125/-122.6/47/48.5 -I1m

psbasemap -R-125/-122.6/47/48.5 -Jm6 -Ba1f0.5 -V -P -K -X1.5 -Y2 >> tmp.ps 
makecpt -T0/2200/200 -Cjet >tmp.cpt
grdimage tmp.nc -Jm6 -E300 -P -O -K >> tmp.ps
psscale -Dx16c/0.8c+w12c/0.5c -O -Ctmp.cpt >> tmp.ps

#ps2pdf tmp.ps $outfile
psconvert tmp.ps -A -Tj
mv tmp.jpg $outfile
rm gmt.history
rm tmp*

source deactivate gmt
