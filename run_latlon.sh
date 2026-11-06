i=1
while [ $i -le 1 ] #!! <=== USER INPUT 1
do
lon_in=`sed "${i}q;d" stations_test.txt | cut -d ' ' -f1`
lat_in=`sed "${i}q;d" stations_test.txt | cut -d ' ' -f2`
#echo $lon_in $lat_in
DIR_IN=`echo "INPUT_"$lon_in"_"$lat_in`
FILE_OUT=`echo "OUT_"$lon_in"_"$lat_in".txt"`
LOG_OUT=`echo "LOG_"$lon_in"_"$lat_in".txt"`
ERR_OUT=`echo "ERR_"$lon_in"_"$lat_in".txt"`
echo "-------------------------------------------------------------------------"
echo "******************FORCING PARAMETRS EXTRACTION START********************"
echo "-------------------------------------------------------------------------"
mkdir OUTPUT
ferret -nojnl << EOF

let lon_in = $lon_in                 
let lat_in = $lat_in     
let clm_init = 1                   
set reg/t=01-jan-2014:31-dec-2018  !! <=== USER INPUT 2

sp mkdir $DIR_IN

sp rm -rf INPUT/
sp mkdir INPUT/

sp echo $lon_in > "./INPUT/details_argo11.txt"
sp echo $lat_in >> "./INPUT/details_argo11.txt"
use FORCING_DATA/T2M_2010_2022_1H_SHIFT_MEAN.nc ! 6-hr interval
use FORCING_DATA/TOT_PRE_2010_2022_1H_SHIFT_MEAN.nc ! 12-hr interval
use FORCING_DATA/SFX_2010_2022_1H_SHIFT_MEAN.nc ! 12-hr interval
use FORCING_DATA/SPE_HUM_2010_2022_1H_SHIFT_MEAN.nc ! 6-hr interval
use FORCING_DATA/U10_2010_2022_1H_SHIFT_MEAN.nc ! 6-hr interval
use FORCING_DATA/V10_2010_2022_1H_SHIFT_MEAN.nc ! 6-hr interval
use FORCING_DATA/SSS_2011_2022_3DAY_FILL.nc ! 3-day interval
use FORCING_DATA/WVEL_2013_2018_3DAY_250m_FILL_interp.nc ! 3-day interval

use FORCING_DATA/salt_5m_cli_fill_interp.nc ! Annual mean
use FORCING_DATA/temp_5m_cli_fill_interp.nc ! Annual mean
use FORCING_DATA/woa18_nitrate_01_12_MON_250m_fill_interp.nc !Seasonal mean
use FORCING_DATA/woa18_silicate_01_12_MON_250m_fill_interp.nc !Seasonal mean

use FORCING_DATA/TCO2_Climatological_Mean_interp.nc !Climatological Mean
use FORCING_DATA/TALK_Climatological_Mean_interp.nc !Climatological Mean

use FORCING_DATA/pCO2atm_2013_2018_daily.nc ! NEWWWW
use FORCING_DATA/Ontogenetic_interp.nc ! NEWWW

save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1/file="./INPUT/air_temp_data_for_argo11.txt" t2m[d=1] - 273.15  ! 6-hourly data
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1/file="./INPUT/ppt_data_for_argo11.txt" tp[d=2]/(60*60) ! for hourly total precipitation to rain in m/s ! 12-hourly data
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1/file="./INPUT/sfx_in_data_for_argo11.txt" ssr[d=3]/(60*60) ! for hourly total j/m2 to j/s/m2. 12-hourly data.
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1/file="./INPUT/specific_hum_data_for_argo11.txt" q[d=4] ! 6-hourly data in kg/kg
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1/file="./INPUT/uwnd_data_for_argo11.txt" u10[d=5] ! 6-hourly data in m/s
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1/file="./INPUT/vwnd_data_for_argo11.txt" v10[d=6] ! 6-hourly data in m/s
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1/file="./INPUT/salinity_clim_data_for_argo11.txt" salt_3dy_fill[d=7] ! 3 day data in p.s.u.
save/clobber/x=${lon_in}/y=${lat_in}/k=1:51/file=wvel_extracted_argo11.nc w_vel_fill[d=8] ! vertical velcoity in m/s

cancel reg
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/salt_in_data_for_argo11.txt" salt[d=9] ! salinity initial condition in p.s.u. (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/temp_in_data_for_argo11.txt" temp[d=10] ! temperature initial condition in oC (51 levels)

save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/NO3_in_data_for_argo11.DAT"/l=1 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/SIOH4_in_data_for_argo11.DAT"/l=1 SILICAT[d=12] ! sioh4 initial condition in microM/l (51 levels)

save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_jan_argo11.txt"/l=1 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_feb_argo11.txt"/l=2 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_mar_argo11.txt"/l=3 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_apr_argo11.txt"/l=4 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_may_argo11.txt"/l=5 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_jun_argo11.txt"/l=6 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_jul_argo11.txt"/l=7 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_aug_argo11.txt"/l=8 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_sep_argo11.txt"/l=9 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_oct_argo11.txt"/l=10 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_nov_argo11.txt"/l=11 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/no3_clim_extract_dec_argo11.txt"/l=12 NITRATE[d=11] ! no3 initial condition in microM/l (51 levels)

save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_jan_argo11.txt"/l=1 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels) 
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_feb_argo11.txt"/l=2 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels) 
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_mar_argo11.txt"/l=3 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels) 
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_apr_argo11.txt"/l=4 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels) 
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_may_argo11.txt"/l=5 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels) 
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_jun_argo11.txt"/l=6 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels) 
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_jul_argo11.txt"/l=7 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels) 
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_aug_argo11.txt"/l=8 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels) 
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_sep_argo11.txt"/l=9 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels) 
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_oct_argo11.txt"/l=10 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_nov_argo11.txt"/l=11 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/sioh4_clim_extract_dec_argo11.txt"/l=12 SILICAT[d=12] ! no3 initial condition in microM/l (51 levels)


save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/l=1/file="./INPUT/TCO2_in_data_for_argo11.DAT" TCO2[d=13]

save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/l=1/file="./INPUT/TALK_in_data_for_argo11.DAT" TALK[d=14]


save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/tco2_in_data_for_argo11.txt" tco2[d=13] ! tco2 initial condition in p.s.u. (51 levels)
save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1:51/file="./INPUT/talk_in_data_for_argo11.txt" talk[d=14] ! talk initial condition in p.s.u. (51 levels)


save/clobber/x=${lon_in}/y=${lat_in}/format=(f22.10)/nohead/k=1/file="./INPUT/pCO2_in_data_for_argo11.txt" co2[d=15] !NEWWW

save/clobber/format=(f22.10)/nohead/file="./INPUT/ontofact_data_for_argo11.txt" ONTORGR[d=16] ! 24-hourly data

sp rm "./INPUT/wvel_extracted_argo11.bin"
sp ./run_fort
sp ./read
q
EOF
echo "-------------------------------------------------------------------------"
echo "*****************FORCING PARAMETRS EXTRACTION END*********************"
echo "-------------------------------------------------------------------------"
argo_no='11'
#======================================================================
echo "-------------------------------------------------------------------------"
echo "*****************MODEL RUN START*********************"
echo "-------------------------------------------------------------------------"
sed  -e  "s/details_argo/details_argo$argo_no/g" main_ver7_argo_test.F  > main01.F
sed  -e  "s/air_temp_data/air_temp_data_for_argo$argo_no/g"  main01.F > main0.F
sed  -e  "s/uwnd_data/uwnd_data_for_argo$argo_no/g" main0.F > main1.F
sed  -e  "s/vwnd_data/vwnd_data_for_argo$argo_no/g" main1.F > main2.F
sed  -e  "s/sfx_in_data/sfx_in_data_for_argo$argo_no/g" main2.F > main3.F
sed  -e  "s/ppt_data/ppt_data_for_argo$argo_no/g" main3.F > main4.F
sed  -e  "s/wvel_data/wvel_extracted_argo$argo_no/g" main4.F > main5.F #Fhis file is actualy not needed for calculations.
sed  -e  "s/specific_hum_data/specific_hum_data_for_argo$argo_no/g" main5.F > main6.F
sed  -e  "s/salinity_clim_data/salinity_clim_data_for_argo$argo_no/g" main6.F > main7.F # in which location  salinity_clim extracted???
sed  -e  "s/temp_in_data/temp_in_data_for_argo$argo_no/g" main7.F > main8.F # in which location  salinity_clim extracted???
sed  -e  "s/salt_in_data/salt_in_data_for_argo$argo_no/g" main8.F > main9.F # in which location  salinity_clim extracted???
sed  -e  "s/NO3_in_data/NO3_in_data_for_argo$argo_no/g" main9.F > main10.F # in which location  salinity_clim extracted???
sed  -e  "s/SIOH4_in_data/SIOH4_in_data_for_argo$argo_no/g" main10.F > main11.F # in which location  salinity_clim extracted???
sed  -e  "s/temp_clim_data/temp_clim_data_for_argo$argo_no/g" main11.F > main12.F # in which location  salinity_clim extracted???

sed  -e  "s/tco2_in_data/tco2_in_data_for_argo$argo_no/g" main12.F > main13.F #NEW ADDITION
sed  -e  "s/talk_in_data/talk_in_data_for_argo$argo_no/g" main13.F > main14.F #NEW ADDITION

sed  -e  "s/pCO2_in_data/pCO2_in_data_for_argo$argo_no/g" main14.F > main15.F #NEWWW

sed  -e  "s/ontofact_data/ontofact_data_for_argo$argo_no/g" main15.F > main16.F #NEWWW



sed  -e  "s/argo/argo$argo_no/g" relaxation_nutrients_ver2_test.F > relax.F
#======================================================================


export TMPDIR=$HOME/tmp
mkdir -p $TMPDIR


rm ionmr.out
echo "ifort -Dsnap_tracer -I/opt/cray/pe/netcdf/4.4.1.1.3/INTEL/16.0/include -g -traceback -extend-source 132 -O3 -L/opt/cray/pe/netcdf/4.4.1.1.3/INTEL/16.0/lib main16.F nemuro_3D_ver2_test.F moto_v7_ver4_test.F implicit.F snaps.F stresstc.F  Irradiance_with_depth_test.F daily_average_test.F relax.F co2calc_SWS.F  co2flux.F  drtsafe.F ta_iter_SWS.F scco2.F -mkl -lnetcdff"
ifort -o ionmr.out -Dsnap_tracer -g -traceback -extend-source 132 -O3 \
$(nf-config --fflags) \
main16.F nemuro_3D_ver2_test.F moto_v7_ver4_test.F implicit.F snaps.F stresstc.F \
Irradiance_with_depth_test.F daily_average_test.F relax.F \
co2calc_SWS.F  co2flux.F  drtsafe.F ta_iter_SWS.F scco2.F \
$(nf-config --flibs) -mkl


echo "-------------------------------------------------------------------------"
echo "****************** MODEL IS RUNNING at $lon_in $lat_in ******************"
echo "-------------------------------------------------------------------------"
echo "                                                                         "

./ionmr.out > $FILE_OUT 2>$ERR_OUT
mv -v $FILE_OUT OUTPUT/
mv -v INPUT/* $DIR_IN
echo "                                                                         "
echo "-------------------------------------------------------------------------"
echo "****************** MODEL RUN END ************************"
echo "-------------------------------------------------------------------------"

######################################################################
i=$(( $i + 1 ))
done
