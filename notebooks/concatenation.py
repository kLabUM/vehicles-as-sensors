import os
import datetime
import numpy as np
import pandas as pd
from scipy import spatial
import netCDF4 as nc
import read_data

def geo_dist_approx(lon1, lat1, lon2, lat2, radians=True):
    """
    Calculate approximate distance (in km) between two geographic coordinates.
    """
    if not radians:
        lon1, lat1, lon2, lat2 = np.radians(lon1, lat1, lon2, lat2)
    R = 6371.  #radius of the earth in km
    x = (lon2 - lon1) * np.cos( 0.5*(lat2+lat1) )
    y = lat2 - lat1
    d = R * np.sqrt( x*x + y*y )
    return d

def join_gage_to_vehicle(gage_data, gage_loc, veh, time_error=5, return_dist=True):
    """
    Concatenate gage and vehicle datasets by space and time.

    Parameters:
    gage_data: A pandas dataframe containing the gage data
    gage_loc: A pandas dataframe containing the gage locations.
    date: The date in the form '20140101'
    vehicle_dir: The directory where vehicle data is stored.
    time error: Allowable difference in timestamps.

    Returns:
    Concatenated dataset containing both vehicle and gage data
    """

    # Get subset of gage data for a given day
    station_day = gage_data
    station_day.index.name = 'GageTime'

    # Compute distances between vehicles and gages for all timestamps
    space_distance = geo_dist_approx(np.radians(gage_loc[0]),
                                         np.radians(gage_loc[1]), 
                                         np.radians(veh['Longitude'].values),
                                         np.radians(veh['Latitude'].values))

    # Store vehicle and gage timestamps as integer arrays
    vehtime = veh.index.values.astype(np.int64) / 10**9
    gagetime = station_day.index.values.astype(np.int64) / 10**9

    # Return indices of the bins (gage) to which each vehicle timestamp belongs
    tq = np.digitize(vehtime, gagetime)

    if len(tq) > 0:
        # Pad the rightmost bin to prevent errors
        station_day[station_day.index.max()
                    + datetime.timedelta(minutes=1)] = np.nan
        # Concatenate vehicle and gage data
        cat_df = pd.concat([veh[['Device', 'Wiper',
                                 'Longitude', 'Latitude']].reset_index(),
                            station_day.iloc[tq].reset_index()], axis=1)
        if return_dist:
            cat_df['Dist_km'] = space_distance

        return cat_df

def time_join(axis_time, input_time):
    """
    Return the indices of axis_time nearest to input_time.

    Parameters:
    -----------
    axis_time:  Times to be indexed. Takes an array of
                datetime or datetime-like string/int.
    input_time: Times for which nearest index in axis_time will be found.
                Array of datetime or datetime-like string/int.

    Returns:
    --------
    dist: Pairwise smallest differences (in ns) between axis_time and input_time.
          Array-like with the same size as input_time.
    ix:   Indices of entries in axis_time that are nearest to each entry in
          input_time. Array-like with the same size as input_time.

    Example usage:
    --------------
    # dist, ix = time_join(t1, t2)
    # pairs = np.column_stack(t1[ix], t2)

    """
    if not 'datetime64' in str(np.asarray(axis_time).dtype):
        axis_time = pd.to_datetime(axis_time).values
        assert 'datetime64' in str(axis_time.dtype)
    if not 'datetime64' in str(np.asarray(input_time).dtype):
        input_time = pd.to_datetime(input_time).values
        assert 'datetime64' in str(input_time.dtype)

    axis_time = np.asarray(axis_time).astype(int).ravel()[:, None]
    input_time = np.asarray(input_time).astype(int).ravel()[:, None]

    tree = spatial.cKDTree(axis_time)
    dist, ix = tree.query(input_time)
    return dist, ix

def space_time_join(axis_time, axis_coords, input_time, input_coords):
    """
    Return the indices of axis_time and axis_coords nearest to input_time and
    input coords.

    Parameters:
    -----------
    axis_time:    Times to be indexed. Takes an array of
                  datetime or datetime-like string/int.
    axis_coords:  Spatial coordinates to be indexed. Takes an array of
                  shape (n, 2), where n is the number of points.
    input_time:   Times for which nearest index in axis_time will be found.
                  Array of datetime or datetime-like string/int.
    input_coords: Coordinates for which nearest index in axis_coords will be
                  found. Takes an array of shape (n, 2), where n is the number
                  of points.

    Returns:
    --------
    t_dist:  Pairwise smallest differences (in ns) between axis_time and
             input_time.
    xy_dist: Pairwise smallest differences between axis_coords and
             input_coords.
    i_t:     Indices of entries in axis_time that are nearest to each entry in
             input_time. Array-like with the same size as input_time.
    i_xy:    Indices of entries in axis_coords that are nearest to each entry in
             input_coords. Array-like with the same size as input_coords.

    """

    t_dist, i_t = time_join(axis_time, input_time)

    xy_tree = spatial.cKDTree(axis_coords)
    xy_dist, i_xy = xy_tree.query(input_coords)

    return t_dist, xy_dist, i_t, i_xy

def panel_spacetime_index(radar_panel, t_input, xy_input,
                          x_axis_n=2, y_axis_n=1, t_axis_n=0,
                          return_distance=False):
    """
    Get indices of panel at the t, x and y points given by t_input and xy_input.

    Parameters:
    -----------
    radar_panel: 3D pandas Panel of radar dataset
    xy_input: ndarray of x, y coordinates with shape (n, 2), where n is the
              number of points.
    t_input: 1D ndarray of time of length n, where n is the number of points.
    x_axis_n: Axis number corresponding to x coordinate indices in radar panel.
    y_axis_n: Axis number corresponding to y coordinate indices in radar panel.
    t_axis_n: Axis number corresponding to time indices in radar panel.
    return_distance: Return distances in time and space calculated by
                     nearest neighbors.

    Returns:
    --------
    i_0, i_1, i_2: Time and space indices of panel corresponding to input time
                   and space coordinates. Indices are given in the order of
                   [items, major_axis, minor_axis].
    t_dist, xy_dist: Distances in time and space between nearest neighbors.
                     (Optional)

    """
    x_axis = radar_panel.axes[x_axis_n]
    y_axis = radar_panel.axes[y_axis_n]
    t_axis = radar_panel.axes[t_axis_n]
    xy_axis = np.vstack(np.dstack(np.meshgrid(x_axis, y_axis)))
    t_dist, xy_dist, i_t, i_xy = space_time_join(t_axis, xy_axis,
                                                 t_input, xy_input)

    i_y, i_x = np.unravel_index(i_xy, (len(y_axis), len(x_axis)))
    z = [None, None, None]
    z[x_axis_n] = i_x
    z[y_axis_n] = i_y
    z[t_axis_n] = i_t
    if return_distance:
        return z, t_dist, xy_dist
    else:
        return z

def veh_gage_to_nc(date, outdir, radardir=read_data.radar_dir,
                   vehdir=read_data.vehicle_dir, return_err=True,
                   replace_nan=True):
    radar_panel = read_data.radar_to_panel(os.path.join(radardir,
                                                        'DTX_DPR_%s.nc' % date))
    veh = read_data.read_vehicle_data(os.path.join(vehdir, '%s.csv' % date))
    wu_df, wu_locs = read_data.read_wu_gages()
    wu_day = wu_df.loc[date]

    # Correct vehicle data
    ymin = radar_panel.major_axis.min()
    ymax = radar_panel.major_axis.max()
    xmin = radar_panel.minor_axis.min()
    xmax = radar_panel.minor_axis.max()

    veh = veh[(veh['Latitude'] > ymin) & (veh['Latitude'] < ymax) &
            (veh['Longitude'] > xmin) & (veh['Longitude'] < xmax)]
    veh = veh[veh['GPS_Speed'] > 1]
    veh = veh[~veh['Device'].isin([10139, 10589, 10615])]

    del wu_df

    wiper_panel = pd.Panel(np.repeat(np.nan, radar_panel.size).reshape(radar_panel.shape))
    wiper_panel.items = radar_panel.items
    wiper_panel.major_axis = radar_panel.major_axis
    wiper_panel.minor_axis = radar_panel.minor_axis

    # Need to remember to convert radar_panel to UTC?

    # Join vehicle to radar

    z0, v_tdist, v_xydist = panel_spacetime_index(radar_panel, veh.index, veh[['Longitude', 'Latitude']].values, return_distance=True)

    # Average wiper value at each t, y, x index.

    veh['t'] = z0[0]
    veh['y'] = z0[1]
    veh['x'] = z0[2]
    veh_w = veh.groupby(['t', 'y', 'x'])['Wiper'].mean().reset_index()

    wiper_panel.values[veh_w['t'].values, veh_w['y'].values, veh_w['x'].values] = veh_w['Wiper'].values

    # Vehicle to gages

    g_r, g_c = np.where(~np.isnan(wu_day))

    t = wu_day.index[g_r]
    xy = wu_locs.loc[wu_day.columns[g_c], ['lon', 'lat']].values 

    z1, g_tdist, g_xydist = panel_spacetime_index(radar_panel, t, xy, return_distance=True)

    gage_panel = pd.Panel(np.repeat(np.nan, radar_panel.size).reshape(radar_panel.shape))
    gage_panel.items = radar_panel.items
    gage_panel.major_axis = radar_panel.major_axis
    gage_panel.minor_axis = radar_panel.minor_axis

    gage_panel.values[z1[0], z1[1], z1[2]]  = wu_day.values[g_r, g_c]

    # Write to netCDF

    outfile = os.path.join(outdir, 'data_%s.nc' % date)

    if os.path.exists(outfile):
        os.remove(outfile)
    with nc.Dataset(outfile, 'w') as d:

        # Instantiate dimensions
        time = d.createDimension('time', radar_panel.items.size)
        lat = d.createDimension('lat', radar_panel.major_axis.size)
        lon = d.createDimension('lon', radar_panel.minor_axis.size)

        # Instantiate variables
        times = d.createVariable("time", "i8", ("time"))
        latitudes = d.createVariable("lat","f8",("lat",))
        longitudes = d.createVariable("lon","f8",("lon",))
        wiper = d.createVariable("wiper","f4",("time","lat","lon",))
        gage = d.createVariable("gage","f4",("time","lat","lon",))
        radar = d.createVariable("radar","f8",("time","lat","lon",))

        times[:] = radar_panel.items.values.astype(int)
        latitudes[:] = radar_panel.major_axis.values
        longitudes[:] = radar_panel.minor_axis.values
        if replace_nan:
            wiper[:, :, :] = np.where(np.isnan(wiper_panel.values), np.finfo("f4").min, wiper_panel.values)
            gage[:, :, :] = np.where(np.isnan(gage_panel.values), np.finfo("f4").min, gage_panel.values)
            radar[:, :, :] = np.where(np.isnan(radar_panel.values), np.finfo("f8").min, radar_panel.values)
        else:
            wiper[:, :, :] = wiper_panel.values
            gage[:, :, :] = gage_panel.values
            radar[:, :, :] = radar_panel.values

        del wiper_panel
        del gage_panel

    # Write to error netCDF file
    if return_err:
        outfile = os.path.join(outdir, 'error_%s.nc' % date)

        if os.path.exists(outfile):
            os.remove(outfile)
        with nc.Dataset(outfile, 'w') as e:

            # Instantiate dimensions
            time = e.createDimension('time', radar_panel.items.size)
            lat = e.createDimension('lat', radar_panel.major_axis.size)
            lon = e.createDimension('lon', radar_panel.minor_axis.size)

            # Instantiate variables
            times = e.createVariable("time", "i8", ("time"))
            latitudes = e.createVariable("lat","f8",("lat",))
            longitudes = e.createVariable("lon","f8",("lon",))
            wiper_t = e.createVariable("wiper_t","f8",("time","lat","lon",))
            wiper_xy = e.createVariable("wiper_xy","f4",("time","lat","lon",))
            gage_t = e.createVariable("gage_t","f8",("time","lat","lon",))
            gage_xy = e.createVariable("gage_xy","f4",("time","lat","lon",))

            times[:] = radar_panel.items.values.astype(int)
            latitudes[:] = radar_panel.major_axis.values
            longitudes[:] = radar_panel.minor_axis.values

            # Error arrays for vehicles

            err_df = pd.DataFrame(np.column_stack([z0[0], z0[1], z0[2], v_tdist, v_xydist])).groupby([0, 1, 2]).mean().reset_index()
            err_df[[0, 1, 2]] = err_df[[0, 1, 2]].astype(int)

            # Write to netCDF

            err_t = pd.Panel(np.repeat(np.nan, radar_panel.size).reshape(radar_panel.shape))
            err_t.items = radar_panel.items
            err_t.major_axis = radar_panel.major_axis
            err_t.minor_axis = radar_panel.minor_axis
            err_t.values[err_df[0].values, err_df[1].values, err_df[2].values] = err_df[3].values

            if replace_nan:
                wiper_t[:, :, :] = np.where(np.isnan(err_t.values), np.finfo("f8").min, err_t.values)
            else:
                wiper_t[:, :, :] = err_t.values

            del err_t

            err_xy = pd.Panel(np.repeat(np.nan, radar_panel.size).reshape(radar_panel.shape))
            err_xy.items = radar_panel.items
            err_xy.major_axis = radar_panel.major_axis
            err_xy.minor_axis = radar_panel.minor_axis

            err_xy.values[err_df[0].values, err_df[1].values, err_df[2].values] = err_df[4].values

            if replace_nan:
                wiper_xy[:, :, :] = np.where(np.isnan(err_xy.values), np.finfo("f8").min, err_xy.values)
            else:
                wiper_xy[:, :, :] = err_xy.values

            del err_xy

            # Error arrays for gages

            err_df = pd.DataFrame(np.column_stack([z1[0], z1[1], z1[2], g_tdist, g_xydist])).groupby([0, 1, 2]).mean().reset_index()
            err_df[[0, 1, 2]] = err_df[[0, 1, 2]].astype(int)

            err_t = pd.Panel(np.repeat(np.nan, radar_panel.size).reshape(radar_panel.shape))
            err_t.items = radar_panel.items
            err_t.major_axis = radar_panel.major_axis
            err_t.minor_axis = radar_panel.minor_axis
            err_t.values[err_df[0].values, err_df[1].values, err_df[2].values] = err_df[3].values

            if replace_nan:
                gage_t[:, :, :] = np.where(np.isnan(err_t.values), np.finfo("f8").min, err_t.values)
            else:
                gage_t[:, :, :] = err_t.values

            del err_t

            err_xy = pd.Panel(np.repeat(np.nan, radar_panel.size).reshape(radar_panel.shape))
            err_xy.items = radar_panel.items
            err_xy.major_axis = radar_panel.major_axis
            err_xy.minor_axis = radar_panel.minor_axis

            err_xy.values[err_df[0].values, err_df[1].values, err_df[2].values] = err_df[4].values
            if replace_nan:
                gage_xy[:, :, :] = np.where(np.isnan(err_xy.values), np.finfo("f8").min, err_xy.values)
            else:
                gage_xy[:, :, :] = err_xy.values

            del err_xy
    return 1

