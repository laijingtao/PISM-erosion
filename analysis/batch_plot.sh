#!/bin/bash
set -e

var=velsurf_mag
cbartitle='surface_velocity'
cbarunit=m/a
cbar=0/500

./gmt_plot.sh -i $PISM_DATA_DIR/elev_climate/state/olympics_g1000m_elev_sb_sia_ela_1500_mb_min_-1.0_mb_max_1.0_20000_30000a.nc -o ${var}_olympics_g1000m_elev_sb_sia_ela_1500_mb_min_-1.0_mb_max_1.0_30000a.jpg -v $var --title $cbartitle --unit $cbarunit --cbar $cbar

./gmt_plot.sh -i $PISM_DATA_DIR/elev_climate/state/olympics_g1000m_elev_sb_sia_ela_1500_mb_min_-3.0_mb_max_3.0_0_5000a.nc -o ${var}_olympics_g1000m_elev_sb_sia_ela_1500_mb_min_-3.0_mb_max_3.0_5000a.jpg -v $var --title $cbartitle --unit $cbarunit --cbar $cbar

./gmt_plot.sh -i $PISM_DATA_DIR/elev_climate/state/olympics_g1000m_elev_sb_sia_ela_1500_mb_min_-5.0_mb_max_5.0_0_5000a.nc -o ${var}_olympics_g1000m_elev_sb_sia_ela_1500_mb_min_-5.0_mb_max_5.0_5000a.jpg -v $var --title $cbartitle --unit $cbarunit --cbar $cbar
