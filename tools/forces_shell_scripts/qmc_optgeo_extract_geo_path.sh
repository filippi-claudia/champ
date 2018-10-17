maxiter=$1
maxiter_p1=$((maxiter+1))
numatoms=$2
labelfile=$3
echo " [MOLDEN FORMAT]"
echo " [N_GEO]"
echo "            $maxiter"
echo " [GEOCONV]"
echo " energy"
grep "current_energy " IterData.log | sed 's/^cur.* -/-/' | sed 's/+-.*$//'
echo " [GEOMETRIES] (XYZ)"
for ((i=1; i < maxiter_p1 ; i++))
do
  geofile=iter$i/geometry/NewGeoXYZ.geo
  echo " $numatoms "
  echo " "
  paste $labelfile $geofile | awk '{printf  " %2s     %11.7f   %11.7f   %11.7f \n",$1,$4*0.529177249,$5*0.529177249,$6*0.529177249}'
done
