"""
ncdfvdfcreate -timedims t -vars QC:QR:QS:QG test3d_master_nofilter.nc mydata.vdf
ncdf2vdf -timedims t -vars QC:QR:QS:QG test3d_master_nofilter.nc mydata.vdf
"""

import numpy as np
import pandas as pd
from model_util import crm_io
import vis.soundings as s
import netCDF4 as nc

import os
import time

## Record in master files? Or write a separate file for each variable?
MASTER_FILE  = True

## Which output sets to write?
DYNAMICS     = True
THERMO       = True
MICRO        = True
GAS_CHEM     = False
DIAG         = True
DIAG_3D      = False

suffix = None
#data_src = "/Users/daniel/workspace/crm71_squall_2d/OUTPUT/"
data_src = "/Volumes/legion_home/models/crm71_squall/OUTPUT/"
#data_src = "/Volumes/legion_storage02/crm_testing/OUTPUT/"

#suffix = "500ccn_1600in"
#data_src = "/Volumes/legion_storage02/crm_testing/kshv_2d_largedomain/kshv_%s/" % suffix

model_initial_time = "2011-04-25 18:00:00"

nx, ny, nz = 500, 101, 80
dx, dy, dz = 1000., 1000., 250.
dt         = 5.0 
iax        = -2

t_end = 1.*(60*60) \
      + 56.*(60)  \
      + 0.

if iax < 0:
    dt_out = -1.*iax*60
else:
    dt_out = iax*dt
nt = int(t_end/dt_out)
print nt

## Only edit below here if you are commenting out variables to ignore

dynamics_vars = {
    "U": {"long": "x-velocity", "units": "m/s", "3D": True},
    "V": {"long": "y-velocity", "units": "m/s", "3D": True},
    "W": {"long": "updraft velocity", "units": "m/s", "3D": True},
    #"K": {"long": "TKE", "units": "m^2/s^2", "3D": True},
}

thermo_vars = {
    "T":  {"long": "temperature", "units": "K", "3D": True},
    "PT": {"long": "potential temperature", "units": "K", "3D": True},
    "P":  {"long": "exner function deviation", "units": "unitless", "3D": True},
    "QV": {"long": "water vapor mix rat", "units": "kg/kg", "3D": True},
}

micro_vars = {
    "CCN": {"long": "cloud condensation nuclei", "units": "1/cc", "3D": True},
    "IN":  {"long": "ice nuclei", "units": "1/L", "3D": True},

    "QC":  {"long": "cloud drop mixing ratio", "units": "g/kg", "3D": True},
    "NC":  {"long": "cloud drop number", "units": "1/cm^3", "3D": True},
     
    "QR":  {"long": "rain mixing ratio", "units": "g/kg", "3D": True},
    "NR":  {"long": "rain number", "units": "1/L", "3D": True},
     
    "QS":  {"long": "snow mixing ratio", "units": "g/kg", "3D": True},
    "NS":  {"long": "snow number", "units": "1/L", "3D": True},
     
    "QG":  {"long": "graupel mixing ratio", "units": "g/kg", "3D": True},
    "NG":  {"long": "graupel number", "units": "1/L", "3D": True},
     
    "QI":  {"long": "crystal mixing ratio", "units": "g/kg", "3D": True},
    "NI":  {"long": "crystal number", "units": "1/L", "3D": True},
     
    "QB":  {"long": "bullet mixing ratio", "units": "g/kg", "3D": True},
    "NB":  {"long": "bullet number", "units": "1/L", "3D": True},
    
    "QP":  {"long": "plate mixing ratio", "units": "g/kg", "3D": True},
    "NP":  {"long": "plate number", "units": "1/L", "3D": True},

    "QTT": {"long": "total vapor mix rat", "units": "g/kg", "3D": True},
     
    "PRECIP": {"long": "surface precipitation", "units": "kg/m^2", "3D": False,
                "valid_range": [0., 0.02]},
}

gas_chem_vars = {
    "CH2O": {"long": "formaldehyde gas", "units": "pptv", "3D": True},
    "CH4": {"long": "methane gas", "units": "ppbv", "3D": True},
    "CH3O2H": {"long": "acetic acid gas", "units": "pptv", "3D": True},
    "DMS": {"long": "dimethyl sulfide gas", "units": "0.1 pptv", "3D": True},
    "HO": {"long": "hydroxide gas", "units": "10^3 1/cm^3", "3D": True},
    "HO2": {"long": "hydroperoxyl gas", "units": "ppmv", "3D": True},
    "HNO3": {"long": "nitric acid gas", "units": "pptv", "3D": True},
    "H2O2": {"long": "hydrogen peroxide gas", "units": "pptv", "3D": True},
    "H2SO4": {"long": "sulfuric acid gas", "units": "pptv", "3D": True},
    "O3": {"long": "ozone gas", "units": "ppbv", "3D": True},
    "NO": {"long": "nitric oxide gas", "units": "ppmv", "3D": True},
    "NO2": {"long": "nitrogen dioxide gas", "units": "ppmv", "3D": True},
    "SO2": {"long": "sulfur dioxide gas", "units": "pptv", "3D": True},
}
diag_vars = [
    ("nuci1", "nucleated ice number concentration", "1/kg", True), # tendency of nucleated ice crystals, accounting for existing number of crystals and number of available IN (either from prognostic aerosol or precribed IN concentrations)
    ("nuci2", "unused ice number concentration", "", False),
    ("frozenq", "freezing of cloud droplet mass", "kg/kg", True), # tendency of QC to freeze into bullets
    ("frozenn", "freezing of cloud droplet number", "1/kg", True), # tendency of QC to freeze into bullets
    ("pr", "rain mass flux", "kg/m^2/s", True), 
    ("pnr", "rain number flux", "1/m^2/s", True),
    ("pi", "snow/ice/plate/bullet mass flux", "kg/m^2/s", True),
    ("pni", "snow/ice/plate/bullet number flux", "kg/m^2/s", True),
    ("pg", "graupel mass flux", "kg/m^2/s", True),
    ("png", "graupel number flux", "kg/m^2/s", True),
    ("auq", "autoconversion of cloud to rain mass rate", "kg/kg/s", True),
    ("aun", "autoconversion of cloud to rain number rate", "1/kg/s", True),
    ("aig", "autoconversion of ice/bullet/plate/snow to graupel mass rate", "kg/kg/s", True),
    ("aign", "autoconversion ice/bullet/plate/snow to graupel number rate", "1/kg/s", True),
    ("frg", "rainwater freezing mass rate", "kg/kg/s", True),
    ("frgn", "rainwater freezing number rate", "1/kg/s", True), 
    ("mir", "melting of ice/bullet/plate/snow to rain mass rate", "kg/kg/s", True), 
    ("mirn", "melting of ice/bullet/plate/snow to rain number rate", "1/kg/s", True), 
    ("mgr", "melting of graupel to rain mass rate", "kg/kg/s", True),
    ("mgrn", "melting of graupel to rain number rate", "1/kg/s", True),
    ("clr", "collection of cloud by rain mass rate", "kg/kg/s", True),
    ("clrn", "collection of cloud by rain number rate", "1/kg/s", True),
    ("cli", "collection of cloud by ice/bullet/plate/snow mass rate", "kg/kg/s", True),
    ("clin",  "collection of cloud by ice/bullet/plate/snow mass rate", "1/kg/s", True),
    ("clg", "collection of cloud by graupel mass rate", "kg/kg/s", True),
    ("clgn", "collection of cloud by graupel number rate", "1/kg/s", True),
    ("cri", "collection of rain by ice/bullet/plate/snow mass rate", "kg/kg/s", True),
    ("crin", "collection of rain by ice/bullet/plate/snow number rate", "1/kg/s", True),
    ("crg", "collection of rain by graupel mass rate", "kg/kg/s", True),
    ("crgn", "collection of rain by graupel number rate", "1/kg/s", True),
    ("crrn", "rain self-collection/breakup number", "1/kg", True),
    ("cir", "collection of ice by rain mass rate", "kg/kg/s", True), # Currently hard-disabled in model
    ("cirn", "collection of ice by rain number rate", "1/kg/s", True),
    ("cig", "collection of ice/bullet/plate/snow by graupel mass rate", "kg/kg/s", True),
    ("cign", "collection of ice/bullet/plate/snow by graupel number rate", "1/kg/s", True),
    ("ciin", "autoconversion of ice/bullet/plate to snow number rate", "1/kg/s", True),
    ("cd", "tendency in cloudwater mass due to condensation", "kg/kg/s", True),
    ("cdn", "tendency in cloudwater number due to condensation", "1/kg/s", True),
    ("cdr", "tendency in rainwater mass due to condensation", "kg/kg/s", True),
    ("cdrn", "tendency in rainwater number due to condensation", "1/kg/s", True),
    ("dep1", "graupel evaporation due to riming number rate", "1/kg/s", True),
    ("dep2", "tendency in ice/plate/bullet/snow mass due to condensation", "kg/kg/s", True),
    ("depn", "tendency in ice/plate/bullet/snow number due to condensation", "1/kg/s", True),
    ("depg1", "graupel evaporation due to riming mass rate", "kg/kg/s", True),
    ("depg2", "tendency in graupel mass due to condensation", "kg/kg/s", True),
    ("depgn", "tendency in graupel number due to condensation", "1/kg/s", True),
]

diag_3d_vars = {
    "AUTO_CR": {"long": "autoconversion of cloud to rain", "units": "kg/kg/s", "3D": True}, 
    "AUTO_IG": {"long": "autoconversion of ice to graupel", "units": "kg/kg/s", "3D": True}, 
}

all_vars = {}
if DYNAMICS: all_vars.update(dynamics_vars)
if THERMO: all_vars.update(thermo_vars)
if MICRO: all_vars.update(micro_vars)
if GAS_CHEM: all_vars.update(gas_chem_vars)
if DIAG_3D: all_vars.update(diag_3d_vars)

print sorted(all_vars.keys())

if __name__ == "__main__":

    def make_file(filename, heights, netcdf_format="NETCDF4"):

        output = nc.Dataset(filename, "w", format=netcdf_format)
        
        output.createDimension("xc", nx)
        output.createDimension("xcu", nx+1)

        output.createDimension("yc", ny)
        output.createDimension("ycv", ny+1)

        output.createDimension("zc", nz)
        output.createDimension("zcw", nz+1)
        output.createDimension("time", nt)    
        
        times = output.createVariable("time", "i4", ("time", ))
        xs = output.createVariable("xc", "f4", ("xc", ))
        ys = output.createVariable("yc", "f4", ("yc", ))
        zs = output.createVariable("zc", "f4", ("zc", )) 

        xus = output.createVariable("xcu", "f4", ("xcu", ))
        yvs = output.createVariable("ycv", "f4", ("ycv", ))
        zws = output.createVariable("zcw", "f4", ("zcw", ))

        zs[:] = (dz/2)+heights[:]
        zs.units = "meter"
        zs.long_name = "z-coordinate in Cartesian system"
        zs.axis = "z"
        zs.positive = "up"

        #print len(heights[:]), len([heights[:], ] + [heights[-1]+dz, ])
        all_heights = list(heights[:])
        all_heights.append(heights[-1]+dz)
        zws[:] = np.array(all_heights)
        zws.units = "meter"
        zws.long_name = "z-coordinate in staggered Cartesian system for z-velocity"
        zws.axis = "z"
        zws.positive = "up"


        xs[:] = (dx/2)+np.arange(nx)*dx
        xs.units = "meter"
        xs.axis = "x"
        xs.long_name = "x-coordinate in Cartesian system"

        xus[:] = np.arange(nx+1)*dx 
        xus.units = "meter"
        xus.axis = "x"
        xus.long_name = "x-coordinate in staggered Cartesian system for x-velocity"
        
        ys[:] = (dy/2)+np.arange(ny)*dy
        ys.units = "meter"
        ys.axis = "y"
        ys.long_name = "y-coordinate in Cartesian system"

        yvs[:] = np.arange(ny+1)*dy
        yvs.units = "meter"
        yvs.axis = "y"
        yvs.long_name = "y-coordinate in staggered Cartesian system for y-velocity"
        
        times[:] = np.arange(nt)*dt_out + dt_out
        print "dt_out = %.1f seconds" % dt_out 
        print "times", times[:]
        times.units = "seconds since %s" % model_initial_time
        times.calendar = 'gregorian'
        times.long_name = "time"

        return output

    ## Read initial profile
    with open(data_src+"cm.prt", "r") as f:
        lines = f.readlines()

        case_line = lines[1]
        case = case_line.strip().split()[5]
        
        ## Fast forward to line beginning with "Z(M)"
        for i, line in enumerate(lines):
            if line.strip().startswith("Z(M)"): break
        
        model_meta = lines[:i]

        header = line.strip()
        initial_profile_lines = map(lambda x: x.strip(), lines[i+2:i+2+nz])
        
        ## Convert to floats, numpy array
        initial_profiles = []
        for line in initial_profile_lines:
            print line
            bits = line.strip().split()[1:]
            bits = map(float, bits)
            initial_profiles.append(bits)
        initial_profiles = np.array(initial_profiles)
        
        ## Build into pandas DataFrame
        
        initial = pd.DataFrame(initial_profiles, reversed(range(1, nz+1)),
                               header.strip().split())
        print initial
        
        temperature = initial["T0(K)"] - 273.15 # C
        pressure = initial["P(mb)"]             # mb
        exner = initial["PAI0"]                 # unitless
        wv = initial["QV0(G/KG)"]/1000.         # to g/g
        heights = initial["Z(M)"].values[::-1]
        
    ## From the initial profiles, calculate dewpoint temperature at
    ## each model level in order to plot sounding. Use the Magnus formula
    epsilon = 0.622
    es_T = s.calc_es(temperature)
    p = pressure*100.
    e_env = p*wv/epsilon
    RH = (e_env/es_T)*100.

    a = 6.112 # mb
    b = 17.67 #
    c = 243.5 # C
    gamma = lambda T, RH: np.log(RH/100.) + b*T/(c+T)
    dewpoint = c*gamma(temperature, RH)/(b-gamma(temperature, RH))

    sounding_data = { 
        "T0": {'data': temperature, 'long': "basic state temperature", 'units': "K"},
        "P0": {'data': pressure, 'long': "basic state pressure", 'units': "hPa"},
        "QV0": {'data': wv, 'long': "basic state water vapor mixing ratio", 'units': "kg/kg"}, 
        "PI0": {'data': exner, 'long': "basic state exner function", 'units': "unitless"},
    }

    #######################################################
    
    if MASTER_FILE:
        if not suffix:
            master_output = make_file(case+"_master.nc", heights)
        else:
            master_output = make_file(case+"_%s.nc" % suffix, heights)
        print "WRITING MASTER OUTPUT"

        ## Add the reference profiles/data to the master data
        for var, d in sounding_data.iteritems():
            dtw = master_output.createVariable(var, "f4", ("zc",))
            dtw.units = d['units']
            dtw.long_name = d['long']
            dtw[:] = d['data'].values[::-1]
        ps = master_output.createVariable("PS", "f4", ())
        ps.units = "Pa"
        ps.long_name = "reference ps in exner function"
        ps[:] = 1e5

    for var, d in all_vars.iteritems():
        print var, d['long']        

        if not os.path.exists(data_src+var):
            print "no output found; skipping"
            continue
        
        print "reading file... ",
    
        ## Set up the dimension tags and transposition necessary to convert the
        ## CRM binary output to valid CF-1.5 shapes/formats
        if d['3D']:
            print "(3D) ",
            if var == "U":
                fortran_shape = (nx+1, ny, nz)
                dim_shape = ("time", "zc", "yc", "xcu")
            elif var == "V": 
                fortran_shape = (nx, ny+1, nz)
                dim_shape = ("time", "zc", "ycv", "xc")
            elif var == "W":
                fortran_shape = (nx, ny, nz+1)
                dim_shape = ("time", "zcw", "yc", "xc")
            else:
                fortran_shape = (nx, ny, nz)
                dim_shape = ("time", "zc", "yc", "xc")
        else:
            fortran_shape = (nx, ny, 1)
            dim_shape = ("time", "yc", "xc")
        print dim_shape, fortran_shape

        ## Choose the method to read the CRM binary output
        reader = crm_io.read

        fortran_data = reader(data_src+var, nt, *fortran_shape)
        data = np.squeeze(np.transpose(fortran_data, [0, 3, 2, 1]))

        print "success!"
        print data.shape
        print data.min(), data.max()

        if MASTER_FILE:
            print "   writing var to master file... ",
            dtw = master_output.createVariable(var, "f4", dim_shape, fill_value=-9999.)        
            
            ## Now, need to do some re-ordering. It seems as if the '0' index
            ## on the time dimension corresponds to the terminal model output,
            ## and the first few indices afterwards are spin-up. For now, simply 
            ## re-order things and let the model spin-up worry about itself.
            dtw[:-1] = data[1:]
            dtw[-1] = data[0]
            dtw.units = d['units']
            dtw.long_name = d['long']
            if "valid_range" in d:
                dtw.valid_range = d['valid_range']
            master_output.sync()
            print "done"
        else:
            print "   writing var to individual file... ",
            var_output = make_file(case+"_%s.nc" % var, heights)
            dtw = var_output.createVariable(var, "f4", dim_shape, fill_value=-9999.)

            dtw[:-1] = data[1:]
            dtw[-1] = data[0]
            dtw.units = d['units']
            dtw.long_name = d['long']
            if "valid_range" in d:
                dtw.valid_range = d['valid_range']

            var_output.close()
            print "done"

    if MASTER_FILE:
        ## File attributes
        master_output.title = "CRM71 - %s" % case # replace with case name
        master_output.history = "Created " + time.ctime(time.time())
        master_output.source = "crm_tools/proc_out.py"
        master_output.Conventions = "CF-1.5"
        master_output.meta = model_meta
        
        master_output.close()

    if DIAG: 

        nt = crm_io.read_diag(data_src+"DIAG", nz, True)
        all_time, all_tdiag = crm_io.save_diag(data_src+"DIAG", nz, nt, True)

        if not suffix:
            print case+"-diag.nc"
            diag = nc.Dataset(case+"_diag.nc", "w", format="NETCDF4")
        else:
            print case + "_" + suffix
            diag = nc.Dataset(case+"_%s_diag.nc" % suffix, "w", format="NETCDF4")
        diag.createDimension("zc", nz)
        diag.createDimension("time", nt)
        
        times = diag.createVariable("time", "i4", ("time", ))
        zs = diag.createVariable("zc", "f4", ("zc", )) 

        zs[:] = (dz/2)+heights[:]
        zs.units = "meter"
        zs.long_name = "z-coordinate in Cartesian system"
        zs.positive = "up"

        times[:] = all_time
        times.units = "seconds since %s" % model_initial_time
        times.calendar = 'gregorian'
        times.long_name = "time"

        for i, (var, longname, units, use) in enumerate(diag_vars):
            print i, var, longname, units
            if use:
                dtw = diag.createVariable(var, "f4", ("zc", "time",))
                dtw[:] = all_tdiag[:,:,i].T
                dtw.units = units
                dtw.long_name = longname

        diag.title = "CRM71 - %s - Microphys Diag" % case # replace with case name
        diag.history = "Created " + time.ctime(time.time())
        diag.close()


