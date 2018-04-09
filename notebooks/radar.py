import numpy as np
import pandas as pd
import netCDF4 as nc
from scipy import spatial

def radar_to_panel(radar_path, var_name=None):
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
        if not var_name:
            var_name = tuple(d.variables.keys())[0]

        p = pd.Panel(d.variables[var_name][:,:,:], items=d.variables[dims[0]][:], major_axis=d.variables[dims[1]][:], minor_axis=d.variables[dims[2]][:])
    p.items.name, p.major_axis.name, p.minor_axis.name = dims
    p.items = pd.to_datetime(p.items, unit='s')
    p = p.sort_index(0).sort_index(1).sort_index(2)
    return p

def get_radar_timeseries(radar_panel, xy_input, input_labels=None,
                         x_axis_n=2, y_axis_n=1,
                         t_axis_n=0):
    """
    Get a slice of the radar netCDF dataset along the time axis at the x, y
    points given by xy_input.

    Parameters:
    -----------
    radar_panel: 3D pandas Panel of radar dataset
    xy_input: ndarray of x, y coordinates with shape (n, 2), where n is the
              number of points.
    input_labels: Labels for x, y coordinate pairs.
    x_axis_n: Axis number corresponding to x coordinate indices in radar panel.
    y_axis_n: Axis number corresponding to y coordinate indices in radar panel.
    t_axis_n: Axis number corresponding to time indices in radar panel.

    Returns:
    --------
    Pandas DataFrame of time series. Each column is a time series
    corresponding to a given input x, y coordinate.

    """
    if input_labels is None:
        input_labels = np.arange(len(xy_input))
    x_axis = radar_panel.axes[x_axis_n]
    y_axis = radar_panel.axes[y_axis_n]
    t_axis = radar_panel.axes[t_axis_n]
    xy = np.vstack(np.dstack(np.meshgrid(x_axis, y_axis)))
    tree = spatial.cKDTree(xy)
    dist, ix = tree.query(xy_input)
    i_y, i_x = np.unravel_index(ix, (len(y_axis), len(x_axis)))
    i_t = [slice(None) for i in i_y]
    z = [None, None, None]
    z[x_axis_n] = i_x
    z[y_axis_n] = i_y
    z[t_axis_n] = i_t
    out_df = pd.concat([radar_panel.iloc[tuple(i)]
                        for i in zip(*z)], axis=1)
    out_df.columns = input_labels
    return out_df


#### Example Usage
# p = radar_to_panel('DTX_DPR_20140811.nc')
