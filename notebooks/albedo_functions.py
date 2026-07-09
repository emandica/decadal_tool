import xarray as xr
import pandas as pd
import numpy as np
import xskillscore as xs
import dask
from scipy import stats

import warnings

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import logging

def global_average(dset):     
    lat_dim = 'lat' if 'lat' in dset.dims else 'latitude'
    lon_dim = 'lon' if 'lon' in dset.dims else 'longitude'
    
    weights = np.cos(np.deg2rad(dset[lat_dim]))
    weights.name = "weights"
    
    return dset.weighted(weights).mean((lat_dim, lon_dim))
    
##################################################################################
def fix_longitude(ds):
    """Shift longitudes from 0–360 to -180–180 and sort them."""
    if 'lon' in ds.dims or 'longitude' in ds.dims:
        lon_name = 'lon' if 'lon' in ds.dims else 'longitude'
        ds = ds.assign_coords({lon_name: (((ds[lon_name] + 180) % 360) - 180)})
        ds = ds.sortby(lon_name)
    return ds
    
##################################################################################
def domain_selection(var, lat_min, lat_max, lon_min, lon_max):
#select domain
    var = var.sel(**{'lat': slice(lat_min,lat_max), 'lon': slice(lon_min,lon_max)})

#weighed mean
    weights = np.cos(np.deg2rad(var.lat))
    weights.name = "weights"

    var = var.weighted(weights)
    var = var.mean(("lon", "lat"))
    
    return var

##################################################################################
# Function to compute the slope and p-value for a single point (lat, lon)
def compute_slope_and_pvalue(albedo, temp):
    mask = ~np.isnan(albedo) & ~np.isnan(temp)
    if np.sum(mask) < 2 or np.all(x[mask] == x[mask][0]):  # If less than two valid data points, return NaN
        return np.nan, np.nan
    slope, intercept, r_value, p_value, std_err = stats.linregress(albedo[mask], temp[mask])
    return slope, p_value    
    
##################################################################################
def cdf_matching_single(obs, start, mid_1, mid_2, end, type='post'):
    obs_a = obs.sel(time=slice(start, mid_1))
    obs_a_mean = obs_a.mean('time')
    obs_a_std = obs_a.std('time')

    obs_b = obs.sel(time=slice(mid_2, end))
    obs_b_mean = obs_b.mean('time')
    obs_b_std = obs_b.std('time')

    if type == 'post':
        #alpha_parameter = obs_a_std / obs_b_std
        #obs_b = alpha_parameter * (obs_b - obs_b_mean) + obs_a_mean
        obs_b = (obs_b - obs_b_mean) + obs_a_mean
    elif type == 'pre':
        alpha_parameter = obs_b_std / obs_a_std
        obs_a = alpha_parameter * (obs_a - obs_a_mean) + obs_b_mean

    return xr.concat([obs_a, obs_b], dim='time')
    
##################################################################################
def detrend_dim(da, deg=1):
    """
    Detrend a dataset along a single dimension using polynomial regression.

    Parameters:
    da (xarray.DataArray or xarray.Dataset): Dataset to be detrended.
    deg (int, optional): Degree of the polynomial for detrending. Defaults to 1 (linear).

    Returns:
    xarray.DataArray or xarray.Dataset: Detrended dataset.
    """
    # Perform polynomial fit along the 'time' dimension
    p = da.polyfit('time', deg=deg)
    
    # Compute the fitted values
    fit = xr.polyval(da.time, p.polyfit_coefficients)

    # Detrend by subtracting the fitted values
    detrended_da = da - fit
    
    return detrended_da
    
##################################################################################
def DJF_seasonal_mean(dset):
    # Resample the dataset to compute seasonal means starting in December
    dset_resampled = dset.resample(time="QS-DEC").mean(dim="time")
    
    # Select only the times corresponding to the DJF season
    #djf_season = dset_resampled.sel(time=dset_resampled.time.dt.season=="DJF")
    djf_season = dset_resampled.sel(time=dset_resampled.time.dt.month.isin([12, 1, 2]))

    return djf_season

##################################################################################
def yearly_mean(dset):
    # Resample the dataset to compute yearly means
    dset_resampled = dset.resample(time="YE").mean(dim="time")
    
    return dset_resampled

##################################################################################
def time_values(dset):
    # Extract the time values from the numpy.datetime64 array
    time_values = dset['time'].values

    # Convert the time values to a pandas DateTimeIndex
    time_index = pd.DatetimeIndex(time_values)

    # Set the time to the first day of the month at midnight
    time_index = time_index.to_period('M').to_timestamp()

    # Update the time coordinate in the dataset
    dset['time'] = time_index
    return dset

##################################################################################
def albedo_data(obs, initial_date):
    # Select time slice
    obs = obs.sel(time=slice(str(initial_date) + '-11', '2014-02'))
    #obs['lon'] = obs['lon'] - 180

    # Calculate broadband albedo
    RSUN_UV = 0.459760
    RSUN_NI = 0.540240
    obs_bb_alb = RSUN_NI * obs['AL_DH_NI'] + RSUN_UV * obs['AL_DH_VI']
    obs_bb_alb = obs_bb_alb.to_dataset(name='bb_albedo')

    # Adjust longitude values
    obs_bb_alb['lon'] = (obs_bb_alb['lon'] + 180) % 360 - 180

    # Calculate DJF seasonal mean
    obs_bb_alb = DJF_seasonal_mean(obs_bb_alb)
    
    # Adjust time values
    #obs_bb_alb = time_values(obs_bb_alb) 

    # Handle invalid data
    obs_bb_alb_ = xr.where(obs_bb_alb < 0, np.nan, obs_bb_alb)
    obs_bb_alb = xr.where(obs_bb_alb_.std('time') < 0.001, np.nan, obs_bb_alb_)
    
    return obs_bb_alb

##################################################################################
def load_seasonal_experiment_data(DATA_PATH, exp, var, initial_date):
    '''
    Load and process seasonal experiment data.
    
    Parameters:
    DATA_PATH (str): Path to experiment folders
    exp (str): Experiment name
    var (str): Variable to retrieve
    initial_date (str): Initial date in 'YYYY' format
    
    Returns:
    xarray.Dataset: Processed dataset
    '''
    # Construct the file path
    path = f'{DATA_PATH}/{exp}/{var}/msmm_{exp}-fal_DJF_sm_1x1.nc'

    # Load the dataset
    dset = xr.open_dataset(path, chunks=-1).load()

    # Select the specified time slice
    dset=dset.sel(time=slice(f'{initial_date}-11', '2014-2'))

    # Adjust longitude values to the range [-180, 180]
    dset['lon'] = (dset['lon'] + 180) % 360 - 180
    dset = dset.sortby('lon')

    # Adjust time values
    dset = time_values(dset)

    # Handle invalid data    
    dset_ = xr.where(dset < 0, np.nan, dset)
    dset = xr.where(dset_.std('time') < 0.001, np.nan, dset_)
    
    return dset

##################################################################################
def alb_corr(obs, dset):
    """
    Calculate the Pearson correlation between observed data and the mean of an experimental dataset.

    Parameters:
    obs (xarray.DataArray): Observed data.
    dset (xarray.DataArray): Experimental dataset.

    Returns:
    tuple: Adjusted correlation values and p-values.
        - em_corr_ (xarray.DataArray): Correlation values with negative correlations set to 0.
        - em_p (xarray.DataArray): p-values for the correlations.
    """
    # Controllo su deviazioni standard nulle
    #if obs.std(dim='time').min() == 0:
    #    print("Attenzione: 'obs' ha regioni con deviazione standard pari a zero.")
    #if dset.mean('number').std(dim='time').min() == 0:
    #if dset.std(dim='time').min() == 0:
    #   print("Attenzione: 'dset' ha regioni con deviazione standard pari a zero.")

    # Pulizia dei dati
    #obs = obs.where(~np.isnan(obs))
    #dset = dset.where(~np.isnan(dset))

    # Calcolo della correlazione con gestione dei warning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        #em_corr = xs.pearson_r(obs, dset.mean('number'), dim='time', skipna=True)
        #em_p = xs.pearson_r_eff_p_value(obs, dset.mean('number'), dim='time', skipna=True)
        em_corr = xs.pearson_r(obs, dset.mean('member'), dim='time', skipna=True)
        em_p = xs.pearson_r_eff_p_value(obs, dset.mean('member'), dim='time', skipna=True)

    # Imposta a zero le correlazioni negative
    #em_corr_ = xr.where(em_corr < 0, 0, em_corr)

    # Mascheratura dei NaN per coerenza
    #em_corr_ = em_corr_.where(~np.isnan(obs) & ~np.isnan(dset.mean('number')))
    #em_p = em_p.where(~np.isnan(obs) & ~np.isnan(dset.mean('number')))

    return em_corr, em_p

##################################################################################
def land_seas_mask(dset, tipo = 'land'):

    dset = dset.assign_coords({
    "lon": dset.lon.assign_attrs({"standard_name": "longitude", "units": "degrees_east", "axis": "X"}),
    "lat": dset.lat.assign_attrs({"standard_name": "latitude", "units": "degrees_north", "axis": "Y"})
})
    # Load the land-sea mask dataset
    lsm = xr.open_dataset('/confess/dicarlo/00-dataset/lsm_nemo_ifs_r1x1.nc')

    # Adjust longitudes to the range [-180, 180]
    #lsm['lon'] = (lsm['lon'] + 180) % 360 - 180

    # Broadcast the mask to match the input dataset
    lsm = lsm.lsm.broadcast_like(dset)

    if tipo == 'land':
    # Apply the mask: retain values where lsm > 0.5 (land), otherwise set to NaN
        land = xr.where(lsm > 0.5, dset, np.nan)
    elif tipo == 'ocean':
        land = xr.where(lsm > 0.5, np.nan, dset)
    
    return land

##################################################################################
def bootstrap_quantile_correlation(sens, ctrl, ref, iterations=1000):
    """
    Perform bootstrap analysis to estimate uncertainty in correlation between sensitivity and control datasets with a reference dataset.

    Parameters:
    sens (xarray.DataArray or xarray.Dataset): Sensitivity dataset.
    ctrl (xarray.DataArray or xarray.Dataset): Control dataset.
    ref (xarray.DataArray or xarray.Dataset): Reference dataset.
    iterations (int): Number of bootstrap iterations.

    Returns:
    xarray.DataArray: Quantiles of correlation differences (delta_corr).
    """
    sens = sens.chunk({'member':-1,'time': -1, 'lon': -1, 'lat': -1})
    ctrl = ctrl.chunk({'member':-1,'time': -1, 'lon': -1, 'lat':-1})
    ref = ref.chunk({'time': -1, 'lon': -1, 'lat': -1})
    
    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        # Perform your slicing operation here    
        # Concatenate sensitivity and control datasets along 'number' dimension    
        combined = xr.concat([sens, ctrl], dim='member').chunk({'member':-1})
        
        # Perform bootstrap resampling once to get both f_ra and f_rb
        #bootstrapped = xs.resampling.resample_iterations(combined, iterations*2, 'member', replace=True).mean('member').squeeze()
        
        # Split resampled data
        #midpoint = combined.sizes['member'] // 2
        #f_ra = bootstrapped.isel(iteration=slice(0, midpoint))
        #f_rb = bootstrapped.isel(iteration=slice(midpoint, None))
        
        # Resample iterations with replacement and compute mean    
        f_ra = xs.resampling.resample_iterations(combined, iterations, 'member', replace=True).mean('member').squeeze().compute()
        f_rb = xs.resampling.resample_iterations(combined, iterations, 'member', replace=True).mean('member').squeeze().compute()

        # Calculate Pearson correlation with the reference dataset
        cor_f_ra = xs.pearson_r(f_ra, ref.chunk({'time':-1, 'lon':'auto', 'lat': 'auto'}), 'time', skipna=True)
        cor_f_rb = xs.pearson_r(f_rb, ref.chunk({'time':-1, 'lon':'auto', 'lat': 'auto'}), 'time', skipna=True)

        # Calculate difference in correlations    
        delta_corr = cor_f_ra - cor_f_rb

        # Calculate quantiles of delta_corr    
        sig_delta = delta_corr.chunk(dict(iteration=-1)).quantile([0.005,0.025,0.05,0.10,0.90,0.95,0.975,0.995], dim='iteration')
    
    return sig_delta

##################################################################################
def bootstrap_acc_significance(sens_anom, ctrl_anom, ref_anom, iterations=1000, alpha=0.05):
    """Significativita' della differenza di ACC (media d'ensemble) SENS vs CTRL,
    per serie temporali globali (dim member, time) — variante di
    bootstrap_quantile_correlation senza dimensioni spaziali.

    Costruisce la distribuzione nulla di Delta-ACC ricampionando i membri del
    pool combinato (SENS+CTRL) e confronta la differenza osservata con i quantili.

    Parametri:
        sens_anom, ctrl_anom (xarray.DataArray): anomalie con dim (member, time).
        ref_anom (xarray.DataArray): anomalie di riferimento con dim (time).
        iterations (int): iterazioni di bootstrap.
        alpha (float): livello (0.05 = 95% a due code).

    Ritorna:
        float: 1.0 se la differenza osservata cade fuori dai quantili nulli, altrimenti 0.0.
    """
    combined = xr.concat([sens_anom, ctrl_anom], dim='member')
    f_a = xs.resampling.resample_iterations(combined, iterations, 'member', replace=True).mean('member')
    f_b = xs.resampling.resample_iterations(combined, iterations, 'member', replace=True).mean('member')
    delta_null = (xs.pearson_r(f_a, ref_anom, 'time', skipna=True)
                  - xs.pearson_r(f_b, ref_anom, 'time', skipna=True))
    lo, hi = delta_null.quantile([alpha / 2, 1 - alpha / 2], dim='iteration').values
    obs = (float(xs.pearson_r(sens_anom.mean('member'), ref_anom, 'time', skipna=True))
           - float(xs.pearson_r(ctrl_anom.mean('member'), ref_anom, 'time', skipna=True)))
    return float((obs < lo) or (obs > hi))

##################################################################################
def bootstrap_single_correlation(ds, ref, iterations=1000):
    # Optimize chunking for ds and ref
    ds = ds.chunk({'time': -1, 'member': -1})
    ref = ref.chunk({'time': -1})
    
    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        # Resample iterations with replacement and compute mean    
        f_ra = xs.resampling.resample_iterations_idx(ds, iterations, 'time', replace=True, dim_max=len(ds.time))
        print(f_ra.iteration)
        f_ra = xs.resampling.resample_iterations_idx(ds, iterations, 'member', replace=True, dim_max=len(ds.member))
        f_ra = f_ra.mean('member').squeeze()
        
        # Calculate Pearson correlation with the reference dataset
        cor_f_ra = xs.pearson_r(f_ra, ref.chunk({'time':-1}), 'time', skipna=True)
        
        # Calculate quantiles of delta_corr    
        sig_delta = cor_f_ra.chunk(dict(iteration=-1)).quantile([0.005,0.025,0.05,0.10,0.90,0.95,0.975,0.995], dim='iteration')
    
    return sig_delta.compute()

import math

##################################################################################
def block_bootstrap(dataset, k):
    """
    Perform block bootstrap resampling on the dataset.

    Parameters:
    dataset (xarray.Dataset): The input dataset with dimensions (member, startdate).
    k (int): The number of bootstrap resamples.

    Returns:
    xarray.DataArray: The resampled dataset with an additional 'bootstrap' dimension.
    """
    dataset = dataset.chunk(-1)
    ancillary_a = []
    time_size = dataset.time.size
    member_size = dataset.member.size
    num_years_in_block = round(math.cbrt(time_size))
    block_size = round(time_size / num_years_in_block) + 1
    
    for _ in range(k):
        # Select random start dates with replacement
        r_time = np.random.choice(time_size, size=block_size, replace=True)
        
        # Generate blocks of 5 consecutive years, ensuring we handle edge cases
        ref_time = np.concatenate([
            np.arange(start, min(start + num_years_in_block, time_size))
            if start <= time_size - num_years_in_block
            else np.arange(time_size - num_years_in_block, time_size)
            for start in r_time
        ])[:time_size]

        # Resample the time dimension
        a = dataset.isel(time = ref_time)
        a.coords['time'] = dataset['time']
        
        # Resample members with replacement
        member_indices = np.random.choice(member_size, size=member_size, replace=True)
        a = a.isel(member = member_indices)

        # Compute the mean across the member dimension
        ancillary_a.append(a.mean(dim='member'))

    # Concatenate along a new 'bootstrap' dimension
    ancillary_a = xr.concat(ancillary_a, dim='bootstrap')
    
    return ancillary_a
    
##################################################################################
def map_plot(ds, ds_p, levels, title=None, antartica=True, type='parametric', sign=0.95, cmap='Reds', facecolor='gray'):
    """
    Create a geographical plot of a dataset with significance markers.

    Parameters:
    ds (xarray.DataArray): Dataset to be plotted.
    ds_p (xarray.DataArray or numpy.ndarray): Significance levels or quantiles.
    levels (list): Contour levels for the plot.
    title (str, optional): Title for the plot. Defaults to None.
    antartica (bool, optional): Whether to include Antarctica in the plot. Defaults to True.
    type (str, optional): Type of significance ('parametric' or 'quantile'). Defaults to 'parametric'.
    sign (float, optional): Significance level (e.g., 0.95, 0.90). Defaults to 0.95.
    cmap (str, optional): Colormap for the plot. Defaults to 'Reds'.

    Returns:
    None
    """

    # Create figure and axis
    fig = plt.figure(figsize=[16,9],dpi=900)
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.set_facecolor(facecolor)

    # Plot dataset    
    p = ds.plot.contourf(
        ax=ax,
        levels=levels,
        transform=ccrs.PlateCarree(),
        add_labels=False,
        add_colorbar=False,
        extend='both',
        cmap=cmap,
    )
    
    #plot significance markers
    try:
        significance_mask(ds, ds_p, type, sign, ax)
    except Exception as e:
        print(f"Error applying significance mask: {e}")
    
    #add coastlines
    ax.coastlines()

    # Optionally mask Antarctica    
    if not antartica:
        ax.set_ylim([-60,80])
    
    # add separate colorbar
    
    if antartica:
        cb = plt.colorbar(p, ticks=levels, shrink=0.6, extend='both')
    else:
        cb = plt.colorbar(p, ticks=levels, shrink=0.5, extend='both')
    cb.ax.tick_params(labelsize=20)

    #Drow gridlines and adjust labels
    gl = p.axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.ylocator = mticker.FixedLocator([-80, -60, -40, -20, 0, 20, 40, 60, 80])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 20, 'color': 'black'}
    gl.ylabel_style = {'size': 20, 'color': 'black'}
    
    # Set plot title    
    ax.set_title(title,fontsize=20)

##################################################################################
def significance_mask(ds, ds_sign, type, sign, ax):
    """
    Plot significance markers on a map based on significance levels.

    Parameters:
    ds (xarray.DataArray): Dataset to be checked for significance.
    ds_sign (xarray.DataArray or numpy.ndarray): Significance levels or quantiles.
    type (str): Type of significance ('parametric' or 'quantile').
    sign (float): Significance level (e.g., 0.95, 0.90).
    ax (matplotlib.axes._subplots.AxesSubplot): Matplotlib axis with Cartopy projection.

    Returns:
    None
    """
    try:
        # Define significance thresholds
        if type == 'parametric':
            thresholds = {0.99: (0.005, 0.995), 0.95: (0.025, 0.975), 0.90: (0.05, 0.95), 0.80: (0.10, 0.90)}
            if sign not in thresholds:
                raise ValueError(f"Unsupported significance level: {sign}")
            min_v, max_v = thresholds[sign]
            signific = np.where((ds_sign >= min_v) & (ds_sign <= max_v))

        elif type == 'quantile':
            indices = {0.99: (0, -1), 0.95: (1, -2), 0.90: (2, -3), 0.80: (3, -4)}
            if sign not in indices:
                raise ValueError(f"Unsupported significance level: {sign}")
            min_v, max_v = indices[sign]
            signific = np.where((ds >= ds_sign[min_v, :, :]) & (ds <= ds_sign[max_v, :, :]))
#        elif type == 'interval':
#            if sign not in indices:
#                raise ValueError(f"Unsupported significance level: {sign}")
            
        else:
            raise ValueError(f"Unsupported significance type: {type}")

        # Ensure we have longitude and latitude
        if 'lon' not in ds.coords or 'lat' not in ds.coords:
            raise ValueError("Dataset must contain 'lon' and 'lat' coordinates.")

        # Create a boolean mask for significant points
        signific_mask = np.zeros(ds.shape, dtype=bool)
        signific_mask[signific] = True

        # Create a meshgrid of longitude and latitude values
        lons, lats = np.meshgrid(ds.lon, ds.lat)

        # Subsample the data to reduce points for plotting
        lat_subsample_step = 2  # Example: keep every nth latitude
        lon_subsample_step = 2  # Example: keep every nrd longitude

        # Subsample indices for latitude and longitude
        lat_indices = np.arange(len(ds.lat))[::lat_subsample_step]
        lon_indices = np.arange(len(ds.lon))[::lon_subsample_step]

        # Get the subsampled meshgrid of lons and lats
        lons_subsampled, lats_subsampled = np.meshgrid(ds.lon[lon_indices], ds.lat[lat_indices])

        # Flatten the subsampled meshgrid for easier indexing
        lons_flat = lons_subsampled.flatten()
        lats_flat = lats_subsampled.flatten()

        # Identify the significant points within the subsampled grid
        significant_points = []
        for i in range(len(lons_flat)):
            lon_idx = np.where(ds.lon == lons_flat[i])[0][0]
            logging.info(f"lon_ids = {lon_idx}")
            lat_idx = np.where(ds.lat == lats_flat[i])[0][0]
            logging.info(f"lat_ids = {lat_idx}")
            if signific_mask[lat_idx, lon_idx]:
                significant_points.append((lons_flat[i], lats_flat[i]))
        
        # If there are significant points, unpack and plot them
        if significant_points:
            lons_significant, lats_significant = zip(*significant_points)
            ax.scatter(lons_significant, lats_significant, marker='x', 
                       s=7, c='k', alpha=0.3, transform=ccrs.PlateCarree(), linewidths=1)

    except Exception as e:
        print(f"Error in significance_mask: {e}")

###################################################################################

def lr_plot(ds_x, ds_y, y_pred, p, r, title=None):
    """
    Plots a scatter plot of observed data with a prediction line, 
    and optionally annotates the points with years if time data is available.
    """

    # Create figure
    fig = plt.figure(figsize=[12, 12], dpi=600)
    ax = fig.add_subplot(111)
    
    # Scatter plot
    ax.scatter(ds_x, ds_y, label="Observed Data", c='blue', alpha=0.7)
    
    # Prediction line
    ax.plot(ds_x, y_pred, color='red',
            label=f"Prediction Line (p={p:.2f}, r={r:.2f})",
            linewidth=2)
    
    # Annotate with years if 'time' is available
    if hasattr(ds_x, 'time'):
        try:
            time_values = ds_x.time.values
            time_years = [pd.to_datetime(t).year for t in time_values]
            for i, year in enumerate(time_years):
                ax.annotate(year, (ds_x[i], ds_y[i]),
                            fontsize=12, color='black')
        except Exception as e:
            print(f"Error annotating years: {e}")
    
    # Set axis limits
    ax.set_xlim(-1, 2)
    ax.set_ylim(-1, 2)
    ax.set_aspect('equal')
    
    # Labels and title
    ax.set_title(title if title else "Linear Regression Plot", fontsize=22)
    ax.set_xlabel('Albedo', fontsize=20)
    ax.set_ylabel('T2m', fontsize=20)
    
    # Increase tick label size
    ax.tick_params(axis='both', labelsize=16)
    
    # Add grid
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Add legend (fixed missing comma)
    ax.legend(fontsize=18, loc='lower right')
    
    # Show the plot
    plt.tight_layout()


#########################################################################################
from matplotlib.gridspec import GridSpec

def triple_figure(img1, img2, img3, title):
    # Crea la figura con dimensioni più grandi
    fig = plt.figure(figsize=(10, 8))

    # Definisci il layout: 2 righe, 2 colonne, con GridSpec
    gs = GridSpec(3, 2)
    gs.update(wspace=0.05, hspace=0.0)  # Riduci lo spazio orizzontale e verticale

    # Aggiungi i pannelli in alto
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Aggiungi il pannello in basso centrato (occupando entrambe le colonne)
    ax3 = fig.add_subplot(gs[1, :])

    # Mostra le immagini
    ax1.imshow(img1)
    ax1.axis('off')

    ax2.imshow(img2)
    ax2.axis('off')

    ax3.imshow(img3)
    ax3.axis('off')

    # Sistema il layout
    #plt.tight_layout()
    plt.savefig(title, dpi=600, bbox_inches = 'tight')