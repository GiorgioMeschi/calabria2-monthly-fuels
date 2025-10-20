
#%%

import json 
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import os
import rasterio as rio
from geospatial_tools import geotools as gt
from geospatial_tools import FF_tools

ft = FF_tools.FireTools()
Ras = gt.Raster()

from risico_operational.settings import HOME, DATAPATH, VS

#%%

# define function to plug into get risico point operational function. 
def plot_maps(year, month, outdir, hist_run):
    """
    Plot the fuel maps, susceptibility, dynamic and static inputs.
    """

    #drought variables
    # year, month = datetime.now().year, datetime.now().month
    if not hist_run:
        month = month - 1 # previos month wrt the actual for the avaialble data
    aggrs = [1, 3, 6]  # aggregation in months
    savepath = f'{HOME}/fuel_maps/calabria2-monthly-fuels/risico_operational/views/plot_maps/png'

    # SPI
    for aggr in aggrs:

        basep = f'/mnt/drought-ita-share/archive/Italy/SPI/MCM/maps/{year}/{month:02}'

        day = os.listdir(basep)[-1]
        name = f'SPI{aggr}-MCM_{year}{month:02}{day}.tif'
        path = f'{basep}/{day}/{name}'

        outfile = f'{savepath}/{year}-{month:02}-{day}/SPI{aggr}_{year}-{month:02}-{day}.png'
        output_folder_spi = os.path.dirname(outfile)
        os.makedirs(output_folder_spi, exist_ok=True)
        
        if not os.path.exists(outfile):

            title = f'SPI{aggr} - {year}/{month:02}/{day}'

            with rio.open(path) as src:
                arr = src.read(1)  # band 1
                arr_cropped = arr[710:980, 980:1180]

            fig, ax = Ras.plot_raster(arr_cropped, cmap='RdBu', title=title,
                                    dpi = 200)

            for img in ax.get_images(): 
                img.set_clim(vmin=-3, vmax=3)
            ax.figure.canvas.draw_idle()

            fig.savefig(outfile, dpi=200, bbox_inches='tight')
        

    #SPEI
    for aggr in aggrs:

        basep = f'/mnt/drought-ita-share/archive/Italy/SPEI/MCM-DROPS/maps/{year}/{month:02}'

        day = os.listdir(basep)[-1]
        name = f'SPEI{aggr}-MCM-DROPS_{year}{month:02}{day}.tif'
        path = f'{basep}/{day}/{name}'

        outfile = f'{savepath}/{year}-{month:02}-{day}/SPEI{aggr}_{year}-{month:02}-{day}.png'
        out_folder_spei = os.path.dirname(outfile)
        os.makedirs(out_folder_spei, exist_ok=True)
        
        if not os.path.exists(outfile):

            title = f'SPEI{aggr} - {year}/{month:02}/{day}'

            with rio.open(path) as src:
                arr = src.read(1)  # band 1
                arr_cropped = arr[710:980, 980:1180]

            fig, ax = Ras.plot_raster(arr_cropped, cmap='RdBu', title=title)

            for img in ax.get_images(): 
                img.set_clim(vmin=-5, vmax=5)
            ax.figure.canvas.draw_idle()

            fig.savefig(outfile, dpi=300, bbox_inches='tight')
        
    # check out folders are the same
    if output_folder_spi == out_folder_spei:
        print(f"Output folders are the same: {output_folder_spi}")

    ouput_folder = output_folder_spi


    # fuel map
    hazard_file = f'{outdir}/fuel12cl_wgs84.tif'
    crs = 'EPSG:4326'
    fires_file = f'{DATAPATH}/raw/burned_area/incendi_dpc_2007_2023_calabria_3857.shp'
    fire_col = 'date_iso'

    # readt filename in metaata.txt
    with open(f'{outdir}/metadata.txt', 'r') as f:
        meta = f.readlines()
        year = int(meta[0].strip().split('_')[1])
        month = int(meta[0].strip().split('_')[2])
        if hist_run:
            month = month + 1 # since the spi data are related to the month already in place.
    
    # year = 2025
    # month = 7

    settings = dict(
        fires_file=         fires_file,
        fires_col=          fire_col,
        crs=                crs,
        hazard_path=         hazard_file,
        xboxmin_hist=       0.1,
        yboxmin_hist=       0.1,
        xboxmin_pie=        0.1,
        yboxmin_pie=        0.7,
        out_folder=         ouput_folder,
        year=               year,
        month=              month,
        season=             False,
        haz_nodata=         0,
        pixel_to_ha_factor= 1,
        allow_hist=         False,
        allow_pie=          True,
        allow_fires=        False,
    )

    ft.plot_haz_with_bars(**settings)

    # susc

    susc_path = f'{DATAPATH}/susceptibility/{VS}/susc_calabria_{year}_{month}.tif'
    tr_path = f'{DATAPATH}/susceptibility/{VS}/thresholds/thresholds.json'
    thresholds = json.load(open(tr_path))
    tr1, tr2 = thresholds['lv1'], thresholds['lv2']

    # reproject the susc in epsg4326 to plot not distorted
    dem_file = f'{DATAPATH}/raw/dem/dem_calabria_20m_wgs84.tif'
    susc_path_wgs = f'{DATAPATH}/susceptibility/{VS}/susc_{year}_{month}_wgs84.tif'
    if os.path.exists(susc_path_wgs):
        os.remove(susc_path_wgs)
    Ras.reproject_raster_as_v2(in_file=susc_path,
                            out_file=susc_path_wgs,
                                reference_file=dem_file, 
                                input_crs = 'EPSG:3857', working_crs = 'EPSG:4326', interpolation = 'near')


    settings = dict(
        fires_file= fires_file,
        fires_col= 'date_iso', # 'finaldate',
        crs= 'epsg:4326',
        susc_path= susc_path_wgs,
        xboxmin_hist= 0.1,
        yboxmin_hist= 0.1,
        xboxmin_pie= 0.1,
        yboxmin_pie= 0.7,
        threshold1= tr1,
        threshold2= tr2,
        out_folder= ouput_folder,
        year= year,
        month= month,
        season= False,
        total_ba_period= 1,
        susc_nodata= -1,
        pixel_to_ha_factor= 1,
        allow_hist= False,
        allow_pie= True,
        allow_fires= False,
        normalize_over_y_axis= 10,
        limit_barperc_to_show= 1,
    )

    ft.plot_susc_with_bars(**settings)

    #dem slope northing easting
    slope_file = f'{DATAPATH}/raw/dem/slope_100m_wgs.tif'
    aspect_file = f'{DATAPATH}/raw/dem/aspect_100m_wgs.tif'
    veg_file = f'{DATAPATH}/raw/vegetation/vegetation_ml_wgs84.tif'
    static_outfolder = f'{savepath}/static'
    os.makedirs(static_outfolder, exist_ok=True)

    for file, title in zip([dem_file, slope_file, aspect_file], ['DEM', 'Slope', 'Aspect']):

        filename = f'{static_outfolder}/{title}.png'
        if not os.path.exists(filename):
            with rio.open(file) as src:
                arr = src.read(1)  # band 1
                # remove 0
                arr[np.isnan(arr)] = -9999
                arr = arr.astype('float32')
                arr[arr <= 0] = np.nan

                Ras.plot_raster(arr,
                                cmap='terrain',
                                title=title,
                                shrink_legend=0.6,
                                outpath=filename,
                                # figsize=(12, 10),
                                dpi = 200)


    # vegatation with discrete palette using matplotlib
    filename = f'{static_outfolder}/vegetation.png'

    if not os.path.exists(filename):

        with rio.open(veg_file) as src:
            arr = src.read(1)

        values = [-100000,  0.1, 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15]

        values_offset = [v + 0.5 for v in values]  # offset to center the values in the color bins
        color_dict = {
            "no data": "#ffffff",   # from negative and 0
            "not burnable": "#9e9e9e",
            "dune": "#565656",
            "corsi d acqua": "#323232",
            "Praterie": "#fff0b3",
            "Cespuglieti": "#e6c200",
            "Garighe": "#998000",
            "Macchia": "#4d3d00",
            "Conifere": "#e23b3b",
            "Latifoglie \ninfiammabili": "#ce09dc",
            "Latifoglie \nnon infiammabili": "#1aa31a",
            "Canneti": "#004700",
            "Rupi": "#002200",
            "Aree agricole": "#ffa64d",
            "Oliveti": "#34099a",
            "Vigneti": "#0a17a4",
            }
        

        Ras.plot_raster(arr,
                        array_classes = values_offset, array_colors = list(color_dict.values()), array_names = list(color_dict.keys()),
                        title='Vegetation',
                        shrink_legend=0.7,
                        labelsize = 8.5,
                        outpath = f'{static_outfolder}/vegetation.png',
                        # figsize=(12, 10),
                        dpi = 200)
        

    # return the list of output images
    plt.close('all')

    return ouput_folder, static_outfolder
                
    


#%%
