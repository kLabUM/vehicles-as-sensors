# vehicles-as-sensors
Code and data for "Vehicles as sensors: high accuracy precipitation maps from windshield wiper measurements"

## Data links

### Raw Data
- Video footage labeling (hi-res): 
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/camera_observations/hi_res_labeling_2014-08-11_raw.csv
- Raw NEXRAD data (2014-06-12): 
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/nexrad/NWS_NEXRAD_NXL3_KDTX_20140612000000_20140612235959.tar.gz
- Raw NEXRAD data (2014-06-28): 
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/nexrad/NWS_NEXRAD_NXL3_KDTX_20140628000000_20140628235959.tar.gz
- Raw NEXRAD data (2014-08-11): 
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/nexrad/NWS_NEXRAD_NXL3_KDTX_20140811000000_20140811235959.tar.gz
- NEXRAD DPR Netcdf (2014-06-12): 
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/nexrad/netcdf/DTX_DPR_20140612.nc
- NEXRAD DPR Netcdf (2014-06-28): 
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/nexrad/netcdf/DTX_DPR_20140628.nc
- NEXRAD DPR Netcdf (2014-08-11):
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/nexrad/netcdf/DTX_DPR_20140811.nc
- Rain gages:
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR13
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR22
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR24
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR26
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR33
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR34
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR38
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR40
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR41
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR44
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR47
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR49
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR5
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR51
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR55
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR56
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR57
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIANNAR59
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMICHELS1
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMICHELS10
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMICHELS6
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIDEXTE4
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMISALIN2
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMISALIN4
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMISALIN6
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMISALIN7
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMISALIN8
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMISALIN9
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIYPSIL10
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIYPSIL11
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIYPSIL13
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIYPSIL14
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/rain_gages/KMIYPSIL5
- Vehicle Data:
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/vehicle_data_r/20140612.csv
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/vehicle_data_r/20140628.csv
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/vehicle_data_r/20140629.csv
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/vehicle_data_r/20140701.csv
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/vehicle_data_r/20140805.csv
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/vehicle_data_r/20140811.csv
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/vehicle_data_r/20140819.csv
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/vehicle_data_r/20140820.csv
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/raw_data/vehicle_data_r/20140821.csv

### Processed Data
- Camera data:
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/processed_data/camera_observations_r/camera_observations_1m.csv
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/processed_data/camera_observations_r/camera_observations_comparison.csv
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/processed_data/camera_observations_r/camera_validation.csv
- Radar data:
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/processed_data/radar/data_20140612.nc
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/processed_data/radar/data_20140628.nc
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/processed_data/radar/data_20140811.nc
- Wiper-corrected product:
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/corrected_product/product_20140811_agg.nc
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/corrected_product/product_20140811_nocars.nc
  - https://s3.us-east-2.amazonaws.com/vehicles-as-sensors/corrected_product/roc_data.mat
