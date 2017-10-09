#!/bin/bash
set -x -e

#ncap2 -O -s "delta_T=0.4*delta_T;" pism_dT.nc pism_scaled_dT.nc

python build_climate_file.py

/usr/bin/ncgen -b -o paleo_modifier.nc paleo_modifier.cdl
for P in 0.5 1.0 1.5; do
    for T in 1.5 0.0 -1.5; do
        ncap2 -O -s "delta_T(0)=${T};frac_P(0)=${P}" paleo_modifier.nc paleo_modifier_T_${T}_P_${P}.nc
    done
done

exit

for T in -7.0 -6.0 -5.0 -4.0; do
    ncap2 -O -s "delta_T(0)=${T};" paleo_modifier.nc paleo_modifier_T_${T}.nc
done

for P in 0.5 1.0 1.5; do
    ncap2 -O -s "frac_P(0)=${P};" paleo_modifier.nc paleo_modifier_P_${P}.nc
done

