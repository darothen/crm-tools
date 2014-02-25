"""
Make an interactive figure which displays QC, W, and the horizontal wind vectors
for user-controllable heights and times.
"""

import iris
import iris.quickplot as qplt
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.widgets import Slider

## Read in the data
filename = "supercell_master.nc"

cubes = iris.load(filename)
data = {}
for c in cubes: 
    data[c.var_name] = c
ex_scalar = data['QC']
time = ex_scalar.coord(var_name="time")

## Re-map the horizontal winds to the cell-centered axis

xc  = ex_scalar.coord(var_name='xc')
xcu = data['U'].coord(var_name='xcu')
yc  = ex_scalar.coord(var_name='yc')
ycv = data['V'].coord(var_name='ycv')
zc  = ex_scalar.coord(var_name='zc')
zcw = data['W'].coord(var_name='zcw')

print "W"
w_new = iris.analysis.interpolate.linear(data['W'], 
                                         [(zcw.long_name, zc.points), ])
w_new.remove_coord(zcw.long_name)
w_new.add_dim_coord(zc, 1)
data['W'] = w_new

print "V" 
v_new = iris.analysis.interpolate.linear(data['V'], 
                                         [(ycv.long_name, yc.points), ])
v_new.remove_coord(ycv.long_name)
v_new.add_dim_coord(yc, 2)
data['V'] = v_new

print "U"
u_new = iris.analysis.interpolate.linear(data['U'], 
                                         [(xcu.long_name, xc.points), ])
u_new.remove_coord(xcu.long_name)
u_new.add_dim_coord(xc, 3)
data['U'] = u_new

print "done"

## Make the interactive plot.

fig = plt.figure(figsize=(12., 7.))
ax = fig.add_subplot(111, axisbg='w')
fig.subplots_adjust(bottom=0.25, left=0.1, right=0.9)

ax_sl_alt = plt.axes([0.1, 0.1, 0.75, 0.02])
ax_sl_time = plt.axes([0.1, 0.14, 0.75, 0.02])
cax = plt.axes ([0.92, 0.25, 0.02, .65])
alt_slider = Slider(ax_sl_alt, 'Alt', 0, 15e3, 6875, "%5.0f m")
time_slider = Slider(ax_sl_time, 'Time', 0, 3600*3, 3600, "%5d s")

def plot_data(ax, z, t=360):
    ax = plt.subplot(111)
    
    assert z in zc.points
    assert t in time.points
    
    alt = iris.Constraint(**{zc.name(): z})
    ts = iris.Constraint(time=t)

    w = data['W'].extract(alt & ts)
    v = data['V'].extract(alt & ts)
    u = data['U'].extract(alt & ts)
    qc = data['QC'].extract(alt & ts)
        
    c = qplt.contour(qc, coords=[xc.name(), yc.name()],
                     colors='grey', levels=[1e-3, ],
                     alpha=0.8, linewidths=2.5,
                     title="")
    ax.set_title("")
    
    levels = np.linspace(-20, 20, 201)
    #cw = qplt.contourf(w, coords=[xc.name(), yc.name()],
    #                  cmap=plt.cm.RdBu, vmin=-20., vmax=20.)
    cw = ax.contourf(xc.points, yc.points, w.data,
                     cmap=plt.cm.RdBu_r, levels=levels, extend="both")
    cb = plt.colorbar(cw, cax=cax, orientation='vertical')
    
    ax.set_title("")
    
    sl = 10
    ax.quiver(xc.points[::sl], yc.points[::sl], 
              u.data[::sl, ::sl], v.data[::sl, ::sl],
              units='inches', scale=75, headwidth=2)
    
    hour = t // 3600
    minute = (t % 3600) / 60
    ax.set_title("Time = %02d:%02d | Alt = %5d m" % (hour, minute, z),
                 loc='left', fontsize=11)
plot_data(ax, 6875., 3600.)

def on_change(val):
    z = alt_slider.val
    t = time_slider.val
    
    ## Find closest value in the coord dimenstions
    def _find_nearest(array, val):
        idx = (np.abs(array-val)).argmin()
        return array[idx]
    z = _find_nearest(zc.points, z)
    t = _find_nearest(time.points, t)
    assert z > 0 and t > 0
    #print z, t
    
    ax = plt.subplot(111)
    ax.cla()
    plot_data(ax, z, t)
    fig.canvas.draw()
    
alt_slider.on_changed(on_change)
time_slider.on_changed(on_change)