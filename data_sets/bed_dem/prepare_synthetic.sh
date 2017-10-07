#!/bin/bash
set -x -e

version=99

mv test_dem.nc pism_Synthetic_v${version}.nc
cp pism_Synthetic_v${version}.nc ~/a/glacier/bed_dem/pism_Synthetic_1000m_v${version}.nc
