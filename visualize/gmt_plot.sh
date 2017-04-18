makecpt -Cjet -T5/7.5/.1 > CRcrust.cpt
grdimage CRcrust_mask.grd -R-106/-101/5/14 -Jm2 -CCRcrust.cpt -K -V -P > CRcrust_mask.ps
grdcontour CRcrust_mask.grd -Jm2 -C1 -R-106/-101/5/14 -Wblack -A+kblue+s8 -P -O -K >> CRcrust_mask.ps 
psscale -CCRcrust.cpt -D2.8/-1/6/0.2h -I -B.5:"Topo(m)": -P -O -K -V >> CRcrust_mask.ps
psbasemap -R-106/-101/5/14 -Jm2 -O -V -P -B1g1WeSn >> CRcrust_mask.ps
