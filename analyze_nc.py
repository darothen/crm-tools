import numpy as np
import netCDF4 as nc

import os
import time

data_src = "squall_master.nc"
data_out = "squall_analyzed.nc"

model_initial_time = "2011-04-25 18:00:00"

nx, ny, nz = 500, 101, 80
dx, dy, dz = 1000., 1000., 250.

## Some physical constants

Rd = 287.058
Cp = 1004.
kappa = Rd/Cp
P0 = 1e5 # Pa

def make_file(filename, t, x, y, z, netcdf_format="NETCDF4"):

    output = nc.Dataset(filename, "w", format=netcdf_format)

    nt = len(t)
    #nz = len(z)
    #ny = len(y)
    #nx = len(x)
    
    output.createDimension("xc", nx)
    output.createDimension("yc", ny)
    output.createDimension("zc", nz)
    output.createDimension("time", nt)    
    
    times = output.createVariable("time", "i4", ("time", ))
    xs = output.createVariable("xc", "f4", ("xc", ))
    ys = output.createVariable("yc", "f4", ("yc", ))
    zs = output.createVariable("zc", "f4", ("zc", )) 

    zs[:] = z[:]
    zs.units = "meter"
    zs.long_name = "z-coordinate in Cartesian system"
    zs.axis = "z"
    zs.positive = "up"

    xs[:] = x[:]
    xs.units = "meter"
    xs.axis = "x"
    xs.long_name = "x-coordinate in Cartesian system"
    
    ys[:] = y[:]
    ys.units = "meter"
    ys.axis = "y"
    ys.long_name = "y-coordinate in Cartesian system"
    
    times[:] = t[:]
    times.units = "seconds since %s" % model_initial_time
    times.calendar = 'gregorian'
    times.long_name = "time"

    return output

def make_var(ncfile, data, varname, dim_shape, units, long_name):
    print "%s..." % varname
    print "   writing to master file..."
    dtw = ncfile.createVariable(varname, "f4", dim_shape, fill_value=-9999.)
    dtw[:] = data[:]
    dtw.units = units
    dtw.long_name = long_name
    ncfile.sync()
    print "   done"

def cal_ze(q, n, species, rho):
    qmin = 1e-5 # kg/kg
    nmin = 1e-1 # 1/kg
    zemin = 1e-3

    q = q*1e-3 # g/kg -> kg/kg
    n = n*1e3  #  1/L -> 1/m3
    n_kg = n/rho # 1/m3 -> 1/kg

    coeffs = { "rain"   : 7.295e13,
               "cloud"  : 2.043e13, 
               "ice"    : 508543.02e18,
               "bullet" : 2191334.8e18,
               "plate"  : 508543.03e18,
               "snow"   : 364583.33e18,
               "graupel": 4.74e15, }

    if species not in coeffs:
        return np.zeros_like(q)

    mask = (q > qmin) & (n_kg > nmin)

    ze = np.zeros_like(q)

    if species in ["rain", "graupel" ]:
        ze[mask] += coeffs[species]*((q[mask]*rho[mask])**2.)/n[mask]
    elif species in ["ice", "plate", "bullet", "snow"]:
        ze[mask] += coeffs[species]*((q[mask]*rho[mask])**3.)/(n[mask]**2.)

    return ze

if __name__ == "__main__":

    print "%s -> %s" % (data_src, data_out)

    src = nc.Dataset(data_src, "r")
    data = {}
    for key in src.variables.keys():
        data[key] = src.variables[key][:]
    ref_var = data["P"]

    times = src.variables["time"]
    xs = src.variables["xc"]
    xsu = src.variables["xcu"]
    ys = src.variables["yc"]
    ysv = src.variables["ycv"]
    zs = src.variables["zc"]
    zsw = src.variables["zcw"]

    out = make_file(data_out, times[:], xs[:], ys[:], zs[:])
    dim_shape_3D = ("time", "zc", "yc", "xc")
    dim_shape_prof = ("zc", )

    #################################################################
    ## 1) Map the winds to the cell-centered grid
    print "Mapping winds to cell-centered grid"

    u = data["U"]
    u_new = np.zeros_like(ref_var)
    u_new[:] = 0.5*(u[:,:,:,:-1] + u[:,:,:,1:])
    make_var(out, u_new, "U", dim_shape_3D, "m/s", "x-velocity")

    w = data["W"]
    w_new = np.zeros_like(ref_var)
    w_new[:] = 0.5*(w[:,:-1,...] + w[:,1:,...])
    make_var(out, w_new, "W", dim_shape_3D, "m/s", "z-velocity")

    #################################################################
    ## 2) Thermodynamic state

    ## a) Pressure in understandable units
    pressure = P0*(data['P']**(1./kappa))
    pressure = pressure*1e-2 # Pa -> hPa
    make_var(out, pressure, "PRES", dim_shape_3D, "hPa", "atmospheric pressure")

    ## b) Pressure perturbation
    pressure_pert = np.swapaxes(pressure, 1, 3) - data["P0"]
    pressure_pert = np.swapaxes(pressure_pert, 1, 3)
    make_var(out, pressure_pert, "PRES_PERT", dim_shape_3D, "hPa", "perturbation atmospheric pressure")

    ## c) Base state potential temperature
    pt0 = (data["T0"] + 273.15)*((100.*data["P0"]/P0)**(-kappa))
    make_var(out, pt0, "PT0", dim_shape_prof, "K", "base state potential temperature")

    ## d) Potential temperature
    theta = data["T"]/data["P"]
    make_var(out, theta, "THETA", dim_shape_3D, "K", "potential temperature")

    ## e) Potential temperature perturbation
    theta_pert = np.swapaxes(theta, 1, 3) - pt0
    theta_pert = np.swapaxes(theta_pert, 1, 3)
    make_var(out, theta_pert, "THETA_PERT", dim_shape_3D, "K", "perturbation to potential temperature")

    ## f) Virtual tempeature
    tv = data["T"]*(1. + 0.61*data["QV"])
    make_var(out, tv, "TV", dim_shape_3D, "K", "virtual tempeature")

    ## g) Density
    rho = 100.*pressure/Rd/tv
    make_var(out, rho, "RHO", dim_shape_3D, "kg/m^3", "air density")

    ## h) Density Potential Temperature
    theta_rho = theta*(1. + 0.61*data["QV"] - data["QTT"]*1e-3)
    make_var(out, theta_rho, "THETA_RHO", dim_shape_3D, "K", "density potential temperature")

    ## i) Density potential temperature perturbation
    theta_rho_pert = theta_pert*(1. + 0.61*data["QV"] - data["QTT"]*1e-3)
    make_var(out,theta_rho_pert,"THETA_RHO_PERT", dim_shape_3D, "K", "density potential temperature perturbation")

    #################################################################
    ## 3) Radar Reflectivity

    ## a) rain
    ze_rain = cal_ze(data["QR"], data["NR"], "rain", rho)
    make_var(out, ze_rain, "ZER", dim_shape_3D, "mm^6/m^3", "rain - radar reflectivity")

    ## b) snow
    ze_snow = cal_ze(data["QS"], data["NS"], "snow", rho)
    make_var(out, ze_snow, "ZES", dim_shape_3D, "mm^6/m^3", "snow - radar reflectivity")

    ## c) graupel
    ze_graup = cal_ze(data["QG"], data["NG"], "graupel", rho)
    make_var(out, ze_graup, "ZEG", dim_shape_3D, "mm^6/m^3", "graupel - radar reflectivity")

    ## d) ice
    ze_ice = cal_ze(data["QI"], data["NI"], "ice", rho)
    make_var(out, ze_ice, "ZEI", dim_shape_3D, "mm^6/m^3", "ice - radar reflectivity")

    ## e) plate
    ze_plate = cal_ze(data["QP"], data["NP"], "plate", rho)
    make_var(out, ze_plate, "ZEP", dim_shape_3D, "mm^6/m^3", "plate - radar reflectivity")

    ## f) bullet
    ze_bullet = cal_ze(data["QB"], data["NB"], "bullet", rho)
    make_var(out, ze_bullet, "ZEB", dim_shape_3D, "mm^6/m^3", "bullet - radar reflectivity")

    ## h) cloud
    ze_cloud = cal_ze(data["QC"], data["NC"], "cloud", rho)
    make_var(out, ze_cloud, "ZEC", dim_shape_3D, "mm^6/m^3", "cloud - radar reflectivity")

    ## last) dbz
    ze_sum = ze_rain + ze_snow + ze_graup + ze_ice + ze_plate + ze_bullet + ze_cloud
    mask = ze_sum > 1e-15
    ze_sum[mask] = 10.*np.log10(ze_sum[mask])
    ze_sum[~mask] = -99.
    make_var(out, ze_sum, "REFLEC", dim_shape_3D, "dBZ", "radar return")

    #################################################################
    ## 4) Other data it would be useful to have
    make_var(out, data["QTT"], "QTT", dim_shape_3D, "g/kg", "total liquid/ice water")

    out.title = "CRM71 - analyzed output"
    out.history = "Created " + time.ctime(time.time())
    out.source = "crm_tools/analyze_nc.py"
    out.Conventions = "CF-1.5"
    
    out.close()

print "\nfinished."