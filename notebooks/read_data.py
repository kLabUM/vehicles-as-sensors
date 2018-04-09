import os
import numpy as np
import pandas as pd
import netCDF4 as nc

homedir = os.path.expanduser('~')
aa_dir = homedir + '/data/prcp/aa_gauges'
wu_dir = homedir + '/data/prcp/wu_station'
vehicle_dir = homedir + '/data/vehicle_v2'
radar_dir = homedir + '/data/nexrad3/nc/inst_prcp_rate'

def read_ann_arbor_gages(aa_dir=os.path.join(homedir, aa_dir)):
    """
    Read ann arbor gages into a single table.

    Parameters:
    aa_dir: The directory containing the gage data

    Returns:
    aa_df: A single table containing all ann arbor city gage data
    aa_locs: A table containing lat and lon info for each gage

    """
    # Initialize empty dict to hold data
    d = {}
 
    # Loop through directory
    for fn in os.listdir(aa_dir):
        # Read files
        d[fn[0]] = pd.read_csv(aa_dir + '/' + fn)
        # Rename columns
        d[fn[0]].rename(columns={'Reading Date' : 'date',
                                'Rainfall (in.)' : 'prcp'},
                inplace=True)
        # Convert time column to datetime index
        d[fn[0]]['date'] = pd.to_datetime(d[fn[0]]['date'])
        d[fn[0]].set_index('date', inplace=True)
        d[fn[0]] = d[fn[0]]['prcp']

    # Join gages into a single dataframe
    aa_df = pd.concat([i for i in d.values()], axis=1)
    aa_df.columns = d.keys()
    # Put dates in UTC
    aa_df = aa_df.tz_localize('EST').tz_convert('UTC')
    del d

    aa_locs = {'C' : (42.294157, -83.754970),
                'J' : (42.284517, -83.795354),
                'S' : (42.253200, -83.733444),
                'N' : (42.294157, -83.710069),
                'B' : (42.306814, -83.754970)}

    aa_locs = pd.DataFrame.from_dict(aa_locs, orient='index').rename(columns={0:'lat', 1:'lon'})
    return aa_df, aa_locs


def read_wu_gages(wu_dir=os.path.join(homedir, wu_dir), var='HourlyPrecipIn', var_accum='dailyrainin'):
    """
    Read Weather Underground gages into a single table.

    Parameters:
    wu_dir: The directory containing the gage data

    Returns:
    wu_df: A single table containing all weather underground gage data
    wu_locs: A table containing lat and lon info for each gage

    """

    # Initialize empty dict to hold data
    wu_d = {}

    # Loop through directory
    for fn in os.listdir(wu_dir):
        station_id = fn
        if fn != 'station_locs':
            if fn in ['KMIANNAR33', 'KMISALIN8', 'KMIANNAR47', 'KMIANNAR49']:
                # Read files
                wu = pd.read_csv(os.path.join(wu_dir, fn), index_col=0).set_index('DateUTC')[var]
                wu.index = pd.to_datetime(pd.Series(wu.index)).values
                wu = wu.tz_localize('UTC')
                wu.index.name = 'time'
                wu.name = station_id
                wu_d[station_id] = wu
            else:
                df = pd.read_csv(os.path.join(wu_dir, fn), index_col=0)
                df['Time'] = pd.to_datetime(df['Time'])
                df['DateUTC'] = pd.to_datetime(df['DateUTC'])
                df.set_index('Time', inplace=True, drop=False)

                d = df.groupby(pd.TimeGrouper('d', closed='right'))[var_accum].diff()
                t = df.groupby(pd.TimeGrouper('d', closed='right'))['Time'].diff()
                d = d * (3600000000000 / t.astype(int))

                d[d < 0] = df.loc[d[d < 0].index, var_accum]
                d.name = fn
                wu = pd.concat([d, df['DateUTC']], axis=1).set_index('DateUTC')
                wu = wu.tz_localize('UTC')
                wu.index.name = 'time'
                wu_d[station_id] = wu

    # Concatenate separate gages
    wu_df = pd.concat([wu_d[i] for i in wu_d], axis=1)
    # Remove invalid entries
    wu_df[wu_df < 0] = np.nan
    # Read gage location file
    wu_locs = pd.read_csv(wu_dir + '/' + 'station_locs', index_col=0)
    return wu_df, wu_locs

def read_vehicle_data(veh_file):
    """
    Read vehicle data file for a single day

    Parameters:
    veh_file: Path of the vehicle data as a string.

    Returns:
    veh: A dataframe of the vehicle data for a given day.
    """

    # Read vehicle file
    veh = pd.read_csv(veh_file, header=None)
    # Fix column headings
    veh.columns = ['Device', 'Trip', 'Latitude', 'Longitude', 'Time', 'Wiper', 'GPS_Speed']
    # Convert time to datetime
    veh['Time'] = pd.to_datetime(veh['Time'])
    veh = veh.set_index('Time').tz_localize('UTC')
    return veh

def radar_to_panel(radar_path, var_name=None, dim_map={}, time_unit='s'):
    """
    Reads netcdf file and converts to a pandas Panel with sorted axes

    Parameters:
    -----------
    radar_path: Path to netCDF dataset
    var_name: Output variable name

    Returns:
    --------
    3-dimensional pandas Panel.
    """
    with nc.Dataset(radar_path, 'r') as d:
        dims = tuple(d.dimensions.keys())
        if dim_map:
            dims = tuple([dim_map[dim] for dim in dims])
        if not var_name:
            var_name = tuple(d.variables.keys())[0]

        p = pd.Panel(d.variables[var_name][:,:,:], items=d.variables[dims[0]][:], major_axis=d.variables[dims[1]][:], minor_axis=d.variables[dims[2]][:])
    p.items.name, p.major_axis.name, p.minor_axis.name = dims
    p.items = pd.to_datetime(p.items, unit=time_unit)
    p = p.sort_index(0).sort_index(1).sort_index(2)
    return p
