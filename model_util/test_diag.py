from crm_io import read_diag, save_diag

#filename = "/Volumes/legion_home/models/crm71_2d/OUTPUT/DIAG"
filename = "/Volumes/legion_storage02/crm_testing/kshv_2d_largedomain/kshv_500ccn_100in/DIAG"
nz = 65
spmd = True
ts = 20

nt = read_diag(filename, nz, spmd)
all_time, all_tdiag = save_diag(filename, nz, nt, spmd)
print all_time
print all_tdiag[ts,:,36]

from pylab import *
import numpy as np
fig = plt.figure(1)
plt.clf()
d = np.ma.masked_invalid(all_tdiag[ts,:,36])
d = np.ma.masked_outside(d, -10, 10)
d = np.ma.filled(d, 0.)
print d.shape
plt.plot(d, range(nz), "o")