""" Port of calze.pro

Calculate radar reflectivity factor from CRM output

"""

import iris
import numpy as np

filename = "kshv_master.nc"
qmin = 1e-5 # kg/kg
nmin = 1e-1 # 1/kg
zemin = -30.
ps = 1e5 # for calculating exner deviations

cubes = iris.load(filename)
data = {}
for c in cubes:
    if c.var_name in ["P", "T", "QV", "QR", "NR", "QI", "NI", "QG", "NG", "PI0"]:
        data[c.var_name] = c
ze = np.ones_like(c.data)*zemin

## Conversions following the calze.pro script

## Compute pressure:
Rd, Cp, ps = 287., 1004., 1e5
# 1) Take P; it's the *total exner function* (basic state plus deviation)
pi = data["P"].data
p = ps*(pi**(Cp/Rd)) # pressure in Pa

# 2) Take temperature
t = data["T"].data # temperature in K

# 3) Esimate density. 
Tv = t*(1. + 0.61*data['QV'].data)
rho = p/(Rd*Tv) # density in kg/m^3

## Read in reflective. We want mass mixing ratios in kg/kg, number conc in 1/kg
## which we can later convert to volumes with the calculated rho
qr = data["QR"].data*1e-3 # g/kg -> kg/kg
nr = data["NR"].data*1e3  # 1/L  -> 1/m^3 

mask = (qr > qmin) & (nr > nmin)
ze[mask] += 7.295e13*rho[mask]*(qr[mask]**2.0)/nr[mask]

qi = data["QI"].data*1e-3 # same
ni = data["NI"].data*1e3

mask = (qi > qmin) & (ni > nmin)
ze[mask] += 4.4327*rho[mask]*(qi[mask]**3.0)/(ni[mask]**2.0)

qg = data["QG"].data*1e-3 # same
ng = data["NG"].data*1e3

mask = (qg > qmin) & (ng > nmin)
ze[mask] += 8.95845e14*rho[mask]*(qg[mask]**2.0)/ng[mask]

## Final correction
mask = ze > 1e-3
ze[mask] = 10.*np.log10(ze[mask])

## Just for fun, make a test plot of radar reflectivity
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Slider

def radar_colormap():
    nws_reflectivity_colors = [
    "#646464", # ND
    "#ccffff", # -30
    "#cc99cc", # -25
    "#996699", # -20
    "#663366", # -15
    "#cccc99", # -10
    "#999966", # -5
    "#646464", # 0
    "#04e9e7", # 5
    "#019ff4", # 10
    "#0300f4", # 15
    "#02fd02", # 20
    "#01c501", # 25
    "#008e00", # 30
    "#fdf802", # 35
    "#e5bc00", # 40
    "#fd9500", # 45
    "#fd0000", # 50
    "#d40000", # 55
    "#bc0000", # 60
    "#f800fd", # 65
    "#9854c6", # 70
    "#fdfdfd" # 75
    ]

    cmap = mpl.colors.ListedColormap(nws_reflectivity_colors)
    cmap.set_over(alpha=0.)
    cmap.set_under(alpha=0.)

    return cmap

fig = plt.figure(figsize=(15, 4))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2, left=0.1)

cmap = radar_colormap()
norm = mpl.colors.Normalize(vmin=-29.9, vmax=80)

def plot_data(t):
    t = int(t)
    assert 0 <= t < ze.shape[0]
    c = ax.contourf(ze[t, :, 0, :], levels=np.arange(-35, 80.1, 5),
                    cmap=cmap, norm=norm)

cf = ax.contourf(ze[30, :, 0, :], levels=np.arange(-35, 80.1, 5),
                 cmap=cmap, norm=norm)
cb = plt.colorbar(cf, orientation='vertical')

ax_t = plt.axes([0.1, 0.1, 0.8, 0.02])
t_slider = Slider(ax_t, 'time', 0, t.shape[0], 5, "%d")
def on_change(val):
    t = t_slider.val
    ax.cla()
    plot_data(t)
t_slider.on_changed(on_change)


