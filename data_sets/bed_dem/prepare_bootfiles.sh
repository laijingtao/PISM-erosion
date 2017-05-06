#!/bin/bash
set -x -e 

# Domain size
xmin=330000
ymin=5201000
xmax=530000
ymax=5381000

version=1

#cd $PISM_DATA_DIR
cd $HOME/glacier/data

infile=geotiff_olympics/Olympics_30m.tif

for GRID  in 50 100 200 250 500 1000; do
    outfile=pism_Olympics_${GRID}m_v${version}.nc
#    gdalwarp -overwrite -of netCDF -co "FORMAT=NC4" -r average -s_srs EPSG:26710 -t_srs EPSG:26710 -te $xmin $ymin $xmax $ymax -tr $GRID $GRID $infile $outfile
    gdalwarp -overwrite -of GTiff -r average -s_srs EPSG:26710 -t_srs EPSG:26710 -te $xmin $ymin $xmax $ymax -tr $GRID $GRID $infile tmp_${GRID}.tif
    gdal_translate -of netCDF -co "FORMAT=NC4" tmp_${GRID}.tif $outfile
    ncrename -v Band1,topg $outfile
    # ncatted -a standard_name,topg,o,c,"bedrock_altitude" -a units,topg,o,c,"m" -a proj4,global,o,c,"+init=EPSG:26710" $outfile
    ncatted -a standard_name,topg,o,c,"bedrock_altitude" -a units,topg,o,c,"m" -a _FillValue,topg,d,, $outfile
    ncap2 -O -s "where(topg<0) {topg=-0.;};"  $outfile  $outfile
done

rm tmp_*
#exit

# cut extra part
# Domain size
xmin=400000
ymin=5231000
xmax=523000
ymax=5341000

version=1


for GRID  in 50 100 200 250 500 1000; do
    outfile=pism_Olympics_mtns_${GRID}m_v${version}.nc
    gdalwarp -overwrite -of netCDF -co "FORMAT=NC4" -r average -s_srs EPSG:26710 -t_srs EPSG:26710 -te $xmin $ymin $xmax $ymax -tr $GRID $GRID $infile $outfile
    ncrename -v Band1,topg $outfile
    # ncatted -a standard_name,topg,o,c,"bedrock_altitude" -a units,topg,o,c,"m" -a proj4,global,o,c,"+init=EPSG:26710" $outfile
    ncatted -a standard_name,topg,o,c,"bedrock_altitude" -a units,topg,o,c,"m" $outfile
done


