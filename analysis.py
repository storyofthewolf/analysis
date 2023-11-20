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
parser.add_argument('--quiet',            action='store_true', help='do not print to screen')
parser.add_argument('--vertical',         action='store_true', help='calculate vertical profiles')
parser.add_argument('--synchronous',      action='store_true', help='calculate substellar/antistellar means')
args = parser.parse_args()

root, num, filelist_short, grav, mwdry = analysis_utils.read_file_list()
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

if args.vertical == True:
    print("If using Z vertical coordinates, make sure to set gravity and mwdry in files.in")


# read in climate data from netcdf
for i in range(num):

    if args.quiet == False: 
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('~~~ ', filelist[i])
    if args.vertical == True:
        print('~~~ gravity = ', grav[i])
        print('~~~ mwdry = ', mwdry[i])

    ncid = nc.Dataset(filelist[i], 'r')

    lon  = ncid.variables['lon'][:]
    lat  = ncid.variables['lat'][:]
    lev  = ncid.variables['lev'][:]
    hyai = ncid.variables['hyai'][:]
    hybi = ncid.variables['hybi'][:]
    hyam = ncid.variables['hyam'][:]
    hybm = ncid.variables['hybm'][:]

    nlat = lat.size
    nlon = lon.size
    nlev = lev.size


    # pressures
    PS     = ncid.variables['PS'][:]          ; PS = np.squeeze(PS)
    P0     = ncid.variables['P0'][:]          ; P0 = np.squeeze(P0)
   

    # temperature variables
    TS     = ncid.variables['TS'][:]          ; TS = np.squeeze(TS)
    T      = ncid.variables['T'][:]           ; T  = np.squeeze(T)

    # ice and snow
    ICEFRAC = ncid.variables['ICEFRAC'][:]    ; ICEFRAC = np.squeeze(ICEFRAC)

    # water variables
    Q      = ncid.variables['Q'][:]           ; Q       = np.squeeze(Q)
    RELHUM = ncid.variables['RELHUM'][:]      ; RELHUM  = np.squeeze(RELHUM)
    TMQ    = ncid.variables['TMQ'][:]         ; TMQ     = np.squeeze(TMQ)

    # cloud water
    CLDLIQ   = ncid.variables['CLDLIQ'][:]    ; CLDLIQ    = np.squeeze(CLDLIQ)
    CLDICE   = ncid.variables['CLDICE'][:]    ; CLDICE    = np.squeeze(CLDICE)
    TGCLDLWP = ncid.variables['TGCLDLWP'][:]  ; TGCLDLWP  = np.squeeze(TGCLDLWP)
    TGCLDIWP = ncid.variables['TGCLDIWP'][:]  ; TGCLDIWP  = np.squeeze(TGCLDIWP)
   # TGCLDCWP = ncid.variables['TGCLDCWP'][:]  ; TGCLDCWP  = np.squeeze(TGCLDCWP)

    # cloud fractions
    CLOUD    = ncid.variables['CLOUD'][:]      ; CLOUD    = np.squeeze(CLOUD)
    CLDTOT   = ncid.variables['CLDTOT'][:]     ; CLDTOT    = np.squeeze(CLDTOT)
    CLDLOW   = ncid.variables['CLDLOW'][:]     ; CLDLOW    = np.squeeze(CLDLOW)
    CLDMED   = ncid.variables['CLDMED'][:]     ; CLDMED    = np.squeeze(CLDMED)
    CLDHGH   = ncid.variables['CLDHGH'][:]     ; CLDHGH    = np.squeeze(CLDHGH)

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

    # heating/cooling rates
    QRS   = ncid.variables['QRS'][:]    ; QRS    = np.squeeze(QRS)
    QRL   = ncid.variables['QRL'][:]    ; QRL    = np.squeeze(QRL)

    ncid.close()


    ############################################################################
    ########  compute area weighted averages and derived quantities  ###########
    TS_gmean             = exo.area_weighted_avg(lon, lat, TS)
    ICEFRAC_gmean        = exo.area_weighted_avg(lon, lat, ICEFRAC)
    TMQ_gmean            = exo.area_weighted_avg(lon, lat, TMQ)
    TGCLDLWP_gmean       = exo.area_weighted_avg(lon, lat, TGCLDLWP)
    TGCLDIWP_gmean       = exo.area_weighted_avg(lon, lat, TGCLDIWP)
    CLDTOT_gmean         = exo.area_weighted_avg(lon, lat, CLDTOT)
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

    # define global pressure coordinate arrays
    G = grav[i]
    R = 8.314462/(mwdry[i]/1000.)
    lev_P, ilev_P = exo.hybrid2pressure(nlon, nlat, nlev, PS, P0, hyam, hybm, hyai, hybi)
    lev_Z, ilev_Z = exo.hybrid2height(nlon, nlat, nlev, PS, P0, hyam, hybm, hyai, hybi, T, G, R)
    
    # do vertical profiles
    # this is slow, so only do when requested
    if args.vertical == True:
        # define global mean profiles
        Pmid_profile = analysis_utils.calc_gmean_profiles(lon, lat, lev_P)
        Pint_profile = analysis_utils.calc_gmean_profiles(lon, lat, ilev_P)
        Tmid_profile = analysis_utils.calc_gmean_profiles(lon, lat, T)
        Qmid_profile = analysis_utils.calc_gmean_profiles(lon, lat, Q)
        # create Tint_profile
        Tint_profile = np.zeros((nlev+1), dtype=float)
        Tint_profile[nlev] = TS_gmean    
        Tint_profile[0] = Tmid_profile[0]    
        for z in range(nlev-1):
            Tint_profile[z+1] = (Tmid_profile[z] + Tmid_profile[z+1])/2.
        Zmid_profile = analysis_utils.calc_gmean_profiles(lon, lat, lev_Z)
        Zint_profile = analysis_utils.calc_gmean_profiles(lon, lat, ilev_Z)
        # run temperature profile diagnostics for lapse rate, 
        # stratosphere max temperature and tropopause
        lapse_rate, itropo, istrat = analysis_utils.tprofile_diags(Pmid_profile, Tmid_profile, Zint_profile, Tint_profile)
        T_STRAT_gmean = Tmid_profile[istrat]
        T_TROPO_gmean = Tmid_profile[itropo]
        Q_STRAT_gmean = Qmid_profile[itropo]
        if args.quiet == False:
            print("-------------- midlayer profile ----------------")
            for g in range(nlev):
                print(g, Pmid_profile[g], Zmid_profile[g], Tmid_profile[g], lapse_rate[g])
            print("-------------- interface profile ----------------")
            for g in range(nlev+1):
                print(g, Pint_profile[g], Zint_profile[g], Tint_profile[g])

        # function to print profile information to a text file
        # analysis_utils.print_vertical_to_file(num, filelist_short, data)


    # top layer temperature, water vapor and clouds
    PTOP = lev_P[1,:,:] 
    PTOP_gmean =  exo.area_weighted_avg(lon, lat, PTOP)
    TTOP = T[1,:,:] 
    TTOP_gmean =  exo.area_weighted_avg(lon, lat, TTOP)
    QTOP = Q[1,:,:] 
    QTOP_gmean =  exo.area_weighted_avg(lon, lat, QTOP)
    CLDICE_TOP = CLDICE[1,:,:] 
    CLDICE_TOP_gmean =  exo.area_weighted_avg(lon, lat, CLDICE_TOP)
    CLDLIQ_TOP = CLDLIQ[1,:,:] 
    CLDLIQ_TOP_gmean =  exo.area_weighted_avg(lon, lat, CLDLIQ_TOP)

    # compute substellar and antistellar means
    if args.synchronous == True:
        TS_SS = np.zeros((nlat, nlon), dtype=float)
        TS_AS = np.zeros((nlat, nlon), dtype=float)
        for x in range(nlon):
            for y in range(nlat):
                if (FDS[0,y,x] >  0.0):
                    TS_SS[y,x] = TS[y,x]
                    TS_AS[y,x] = -999.0
                else:
                    TS_SS[y,x] = -999.0
                    TS_AS[y,x] = TS[y,x]
        TS_SS_gmean = exo.area_weighted_avg(lon, lat, TS_SS)
        TS_AS_gmean = exo.area_weighted_avg(lon, lat, TS_AS)
        


    if args.quiet == False:
        ########  print global mean quantities  ###########    
        # These are a set of outputs of common interest for
        # print to screen applications
        print("------------------ global mean ------------------")
        print("TS mean ", TS_gmean)
        print("ICEFRAC", ICEFRAC_gmean)
        print("toa albedo ", toa_albedo_gmean)
        print("olr ", FULTOA_gmean)
        print("TMQ TGCLDLWP TGCLDIWP ", TMQ_gmean, TGCLDIWP_gmean, TGCLDLWP_gmean)
        print("TOA ", toa_balance_gmean, energy_gmean)
        print("SRF ", srf_balance_gmean)
        print("FLNT FSNT", FLNT_gmean, FSNT_gmean)
        if 'FSDTOA' in ncid.variables: print("FSDTOA", FSDTOA_gmean)
        print("LW FLUXES ", FULTOA_gmean, FDLTOA_gmean, FULTOA_gmean - FDLTOA_gmean)
        print("SW FLUXES ", FUSTOA_gmean, FDSTOA_gmean)
        print("TOP ", PTOP_gmean, TTOP_gmean, QTOP_gmean)
        if args.synchronous == True:
            print("TS_SS, TS_AS ", TS_SS_gmean, TS_AS_gmean)
        
    # Presently, the data sent to print to file routines are user specified here
    # Later I might create a namelist around these instead
    x=0
    datacube[x,i] = TS_gmean          ; varnames[x] = 'TS'        ; x=x+1
    if (args.vertical == True):
        datacube[x,i] = T_STRAT_gmean ; varnames[x] = 'T_STRAT'   ; x=x+1
        datacube[x,i] = T_TROPO_gmean ; varnames[x] = 'T_TROPO'   ; x=x+1
    datacube[x,i] = ICEFRAC_gmean     ; varnames[x] = 'ICEFRAC'   ; x=x+1
    datacube[x,i] = toa_albedo_gmean  ; varnames[x] = 'TOAALB'    ; x=x+1
    datacube[x,i] = FULTOA_gmean      ; varnames[x] = 'OLR'       ; x=x+1
    datacube[x,i] = toa_balance_gmean ; varnames[x] = 'toaEBAL'   ; x=x+1
    datacube[x,i] = srf_balance_gmean ; varnames[x] = 'srfEBAL'   ; x=x+1
    datacube[x,i] = TMQ_gmean         ; varnames[x] = 'TMQ'       ; x=x+1
    datacube[x,i] = TGCLDLWP_gmean    ; varnames[x] = 'TGCLDLWP'  ; x=x+1
    datacube[x,i] = TGCLDIWP_gmean    ; varnames[x] = 'TGCLDIWP'  ; x=x+1
    datacube[x,i] = CLDTOT_gmean      ; varnames[x] = 'CLDTOT'    ; x=x+1
    if (args.vertical == True):
        datacube[x,i] = Q_STRAT_gmean ; varnames[x] = 'Q_STRAT'   ; x=x+1
    
# output global mean quantities a text file
analysis_utils.print_data_to_file(num, filelist_short, datacube, varnames)


print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print(' Exiting analysis.py ... ')
print(' ... i hope you found the answers that you seek ')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print(' ')
