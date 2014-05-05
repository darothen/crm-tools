import numpy as np
from wrf_io import get_sounding

filename = "input_sounding"
n = 1000
dry = False

#nl = get_sounding(filename,  n)
zk, p, pd, theta, rho, u, v, qv, nl = get_sounding(filename, dry, n)

zk = zk[:nl]
p = p[:nl]
pd = pd[:nl]
theta = theta[:nl]
rho = rho[:nl]
u = u[:nl]
v = v[:nl]
qv = qv[:nl]

"""
print "input levels ", nl
print " sounding" 
print "  k  height(m)  press (Pa) pd(Pa)     theta (K)   dden(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) "
for k in xrange(nl):
    print " %3d" % k,
    print " %10.3e" % zk[k],
    print " %10.3e" % p[k],
    print " %10.3e" % pd[k],
    print " %10.3e" % theta[k],
    print " %10.3e" % rho[k],
    print " %10.3e" % u[k],
    print " %10.3e" % v[k],
    print " %10.3e" % qv[k]
"""

