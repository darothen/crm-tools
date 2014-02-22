
import datetime
import os.path

import brewer2mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import iris
import iris.coords as coords

coord_map = {
    "y": "y-coordinate in Cartesian system",
    "x": "x-coordinate in Cartesian system",
}

runs = { 12: { 'dptil': 3, 'dqv': 0.9 }, 
         16: { 'dptil': 3, 'dqv': 0.7 }, 
         17: { 'dptil': 3, 'dqv': 0.8 },
}

bmap = brewer2mpl.get_map('Dark2', 'Qualitative', 8)
colors = {}
for i, rn in enumerate(runs.keys()):
    colors[rn] = bmap.hex_colors[i]

init = datetime.datetime(2011, 4, 25, 18)

def run_name_callback(cube, field, filename):
    run_name = os.path.basename(filename)
    run_coord = coords.AuxCoord(run_name, long_name="Run Number", units="no_unit")

    cube.add_aux_coord(run_coord)

constrained_y = iris.Constraint(coord_values={coord_map["y"]:500})
for rn, rd in runs.iteritems():
    cubes = iris.load("supercell_run%d.nc" % rn, 
                         constraints=constrained_y,
                         callback=run_name_callback)
    for c in cubes:
        if c.var_name == "PRECIP": precip = c
    domain_max_precip = precip.collapsed(coord_map['x'], iris.analysis.MAX)
    domain_tot_precip = precip.collapsed(coord_map['x'], iris.analysis.SUM)

    time = domain_max_precip.coords('time')[0].points
    conv_times = [init + datetime.timedelta(seconds=int(t)) for t in time]

    rd['max'] = pd.Series(domain_max_precip.data, index=conv_times)
    rd['tot'] = pd.Series(domain_tot_precip.data, index=conv_times)

## Re-shape to make pandas DataFrame
max_precips, tot_precips = {}, {}
for rn, rd in runs.iteritems():
    max_precips["run %d" % rn] = rd['max']
    tot_precips["run %d" % rn] = rd['tot']
max_precips = pd.DataFrame(max_precips)
tot_precips = pd.DataFrame(tot_precips)

fig, [ax_max, ax_tot] = plt.subplots(2, 1, sharex=True, figsize=(10, 8))
plt.subplots_adjust(hspace=0.1)

plot_kwargs = { 'lw': 2, }
max_precips.plot(ax=ax_max,
                 colormap=bmap.mpl_colormap, **plot_kwargs)
ax_max.set_ylabel("domain max precip (mm)")

tot_precips.plot(ax=ax_tot, style="--", 
                 colormap=bmap.mpl_colormap, legend=False, **plot_kwargs)
ax_tot.set_ylabel("domain total precip (mm)")
ax_tot.set_xlabel("simulation time")






