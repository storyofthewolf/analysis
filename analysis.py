#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# analysis.py                                                                        
#                                                                                    
# Author Eric Wolf                                                                   
# June 2023                                                                          
#                                                                                    
# Purpose:  Analysis of a single file or small batch of files
#                                                                                    
# Notes:  Currently this code only deals with global mean quantities.
#         Future development will include functions to extract vertical profiles, 
#         surface, zonal, and latitudinal contour slices.  The goal is to produce 
#         text file data that will be plotted elsewhere.  
#                                                                                    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import netCDF4 as nc
import numpy as np
import exocampy_tools as exo
import argparse
import analysis_utils

import argparse
import sys

# input arguments and options                                                                                                      
parser = argparse.ArgumentParser()
parser.add_argument('--quiet',      action='store_true', help='do not print to screen')
args = parser.parse_args()



root, num, filelist_short = analysis_utils.read_file_list()
filelist = np.empty(num, dtype='U200')
filelist[:] = root + '/' + filelist_short[:]

# define vector arrays of variables explored
nvars = 50
datacube = np.zeros((nvars, num), dtype=float)
varnames = np.empty(nvars, dtype='U100')

print(' ')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print(' Entering analysis.py ')
print(' files read in from files.in ')



# read in climate data from netcdf
for i in range(num):
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('~~~ ', filelist[i])
    ncid = nc.Dataset(filelist[i], 'r')
    lon = ncid.variables['lon'][:]
    nlon = lon.size
    lat = ncid.variables['lat'][:]
    nlat = lat.size
    lev = ncid.variables['lev'][:]
    nlev = lev.size

    # temperature variables
    TS     = ncid.variables['TS'][:]          ; TS = np.squeeze(TS)
    T      = ncid.variables['T'][:]           ; T  = np.squeeze(T)

    # ice and snow
    ICEFRAC = ncid.variables['ICEFRAC'][:]    ; ICEFRAC = np.squeeze(ICEFRAC)

    # water variables
    Q      = ncid.variables['Q'][:]           ; Q       = np.squeeze(Q)
    RELHUM = ncid.variables['RELHUM'][:]      ; RELHUM  = np.squeeze(RELHUM)
    TMQ    = ncid.variables['TMQ'][:]         ; TMQ     = np.squeeze(TMQ)

    # cloud 
    CLDLIQ   = ncid.variables['CLDLIQ'][:]    ; CLDLIQ    = np.squeeze(CLDLIQ)
    CLDICE   = ncid.variables['CLDICE'][:]    ; CLDICE    = np.squeeze(CLDICE)
    TGCLDLWP = ncid.variables['TGCLDLWP'][:]  ; TGCLDLWP  = np.squeeze(TGCLDLWP)
    TGCLDIWP = ncid.variables['TGCLDIWP'][:]  ; TGCLDIWP  = np.squeeze(TGCLDIWP)
    TGCLDCWP = ncid.variables['TGCLDCWP'][:]  ; TGCLDCWP  = np.squeeze(TGCLDCWP)


    # energy
    FLNT   = ncid.variables['FLNT'][:]    ; FLNT      = np.squeeze(FLNT)
    FSNT   = ncid.variables['FSNT'][:]    ; FSNT      = np.squeeze(FSNT)

    if 'FSDTOA' in ncid.variables:
        # Variable exists, you can proceed to read it
        FSDTOA = ncid.variables['FSDTOA'][:]  ; FSDTOA    = np.squeeze(FSDTOA)

    FLNS    = ncid.variables['FLNS'][:]   ; FLNS      = np.squeeze(FLNS)
    FSNS    = ncid.variables['FSNS'][:]   ; FSNS      = np.squeeze(FSNS)
    LHFLX   = ncid.variables['LHFLX'][:]  ; LHFLX     = np.squeeze(LHFLX)
    SHFLX   = ncid.variables['SHFLX'][:]  ; SHFLX     = np.squeeze(SHFLX)

    # fluxes
    FUL   = ncid.variables['FUL'][:]    ; FUL    = np.squeeze(FUL)
    FDL   = ncid.variables['FDL'][:]    ; FDL    = np.squeeze(FDL)
    FUS   = ncid.variables['FUS'][:]    ; FUS    = np.squeeze(FUS)
    FDS   = ncid.variables['FDS'][:]    ; FDS    = np.squeeze(FDS)
    ncid.close()


    ############################################################################
    ########  compute area weighted averages and derived quantities  ###########
    TS_gmean             = exo.area_weighted_avg(lon, lat, TS)
    ICEFRAC_gmean        = exo.area_weighted_avg(lon, lat, ICEFRAC)
    TMQ_gmean            = exo.area_weighted_avg(lon, lat, TMQ)
    TGCLDLWP_gmean       = exo.area_weighted_avg(lon, lat, TGCLDLWP)
    TGCLDIWP_gmean       = exo.area_weighted_avg(lon, lat, TGCLDIWP)
    TGCLDCWP_gmean       = exo.area_weighted_avg(lon, lat, TGCLDCWP)
    FLNT_gmean           = exo.area_weighted_avg(lon, lat, FLNT)
    FSNT_gmean           = exo.area_weighted_avg(lon, lat, FSNT)
    FLNS_gmean           = exo.area_weighted_avg(lon, lat, FLNS)
    FSNS_gmean           = exo.area_weighted_avg(lon, lat, FSNS)
    SHFLX_gmean          = exo.area_weighted_avg(lon, lat, SHFLX)
    LHFLX_gmean          = exo.area_weighted_avg(lon, lat, LHFLX)
    if 'FSDTOA' in ncid.variables:
        FSDTOA_gmean     = exo.area_weighted_avg(lon, lat, FSDTOA)
    temp                 = FUL[1,:,:] ; temp = np.squeeze(temp)
    FULTOA_gmean         = exo.area_weighted_avg(lon, lat, temp)
    temp                 = FDL[1,:,:] ; temp = np.squeeze(temp)
    FDLTOA_gmean         = exo.area_weighted_avg(lon, lat, temp)
    temp                 = FUS[1,:,:] ; temp = np.squeeze(temp)
    FUSTOA_gmean         = exo.area_weighted_avg(lon, lat, temp)
    temp                 = FDS[1,:,:] ; temp = np.squeeze(temp)
    FDSTOA_gmean         = exo.area_weighted_avg(lon, lat, temp)
    toa_albedo_gmean     = FUSTOA_gmean/FDSTOA_gmean
    toa_balance_gmean    = FSNT_gmean - FLNT_gmean
    srf_balance_gmean    = FSNS_gmean - FLNS_gmean - SHFLX_gmean - LHFLX_gmean
    energy               = FSNT[:,:] - FLNT[:,:] 
    energy_gmean         = exo.area_weighted_avg(lon, lat, energy)

    # top layer temperature, water vapor and clouds
    TTOP = T[1,:,:] 
    TTOP_gmean =  exo.area_weighted_avg(lon, lat, TTOP)
    QTOP = Q[1,:,:] 
    QTOP_gmean =  exo.area_weighted_avg(lon, lat, QTOP)
    CLDICE_TOP = CLDICE[1,:,:] 
    CLDICE_TOP_gmean =  exo.area_weighted_avg(lon, lat, CLDICE_TOP)
    CLDLIQ_TOP = CLDLIQ[1,:,:] 
    CLDLIQ_TOP_gmean =  exo.area_weighted_avg(lon, lat, CLDLIQ_TOP)

    if args.quiet == True:
        print("quiet mode....")
    else:
        ########  print global mean quantities  ###########    
        # These are a set of outputs of common interest for
        # print to screen applications
        print("TS mean ", TS_gmean)
        print("ICEFRAC", ICEFRAC_gmean)
        print("toa albedo ", toa_albedo_gmean)
        print("olr ", FULTOA_gmean)
        print("TMQ TGCLDLWP TGCLDIWP TGCLDCWP ", TMQ_gmean, TGCLDIWP_gmean, TGCLDLWP_gmean, TGCLDCWP_gmean)
        print("TOA ", toa_balance_gmean, energy_gmean)
        print("SRF ", srf_balance_gmean)
        print("FLNT FSNT", FLNT_gmean, FSNT_gmean)
        if 'FSDTOA' in ncid.variables: print("FSDTOA", FSDTOA_gmean)
        print("LW FLUXES ", FULTOA_gmean, FDLTOA_gmean, FULTOA_gmean - FDLTOA_gmean)
        print("SW FLUXES ", FUSTOA_gmean, FDSTOA_gmean)
        print("TOP ", TTOP_gmean, QTOP_gmean, CLDICE_TOP_gmean, CLDLIQ_TOP_gmean)


    # Presently, the data sent to print to file routines are user specified here
    # Later I might create a namelist around these instead
    datacube[0,i] = TS_gmean          ; varnames[0] = 'TS'   
    datacube[1,i] = ICEFRAC_gmean     ; varnames[1] = 'ICEFRAC'   
    datacube[2,i] = toa_albedo_gmean  ; varnames[2] = 'TOAALB'   
    datacube[3,i] = FULTOA_gmean      ; varnames[3] = 'FULTOA'   
    datacube[4,i] = FLNT_gmean        ; varnames[4] = 'FLNT'   
    datacube[5,i] = toa_balance_gmean ; varnames[5] = 'toaEBAL'   
    datacube[6,i] = srf_balance_gmean ; varnames[6] = 'srfEBAL'   
    datacube[7,i] = TMQ_gmean         ; varnames[7] = 'TMQ'   
    datacube[8,i] = TGCLDLWP_gmean    ; varnames[8] = 'TGCLDLWP'   
    datacube[9,i] = TGCLDIWP_gmean    ; varnames[9] = 'TGCLDIWP'   
    

analysis_utils.print_data_to_file(num, filelist_short, datacube, varnames)



print(' Exiting analysis.py ... ')
print(' ... i hope you found the answers to all your questions ')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print(' ')
