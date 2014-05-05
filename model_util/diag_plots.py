"""
Quick diagnostic plots for looking at properties of soundings and atmospheric profiles.

Author: Daniel Rothenberg (darothen@mit.edu)

"""

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

def hydrostatic_plot(zs, pressure, temperature, qv=None,
                     dp_bnds=None, resid_bnds=None):
    """
    Compute the deviation from hydrostatic balance from a given
    moist or dry atmospheric profile of temperature and pressure 
    (and optionally qv), and use it to plot (1) dP/dz and -rho*g,
    and (2) the difference between the two terms.
    
    pressure: Pa
    temperature: K
    qv: kg/kg
    """
    
    nz = len(zs)
    dP_dz = np.zeros(nz)
    
    if qv is None:
        rho = pressure/(287.*temperature)
    else: # use virtual temperature
        tv = temperature*(1. + 0.61*qv)
        rho = pressure/(287.*tv)
    
    for i in xrange(1, nz-1):
        dP_dz[i] = (pressure[i+1] - pressure[i-1])/(zs[i+1] - zs[i-1])
    dP_dz = np.ma.masked_equal(dP_dz, 0)
    
    fig, axes = plt.subplots(1, 2, figsize=(10,4), sharey=True)
    plt.subplots_adjust(wspace=0.1)
    ax_dp, ax_resid = axes
    
    ax_dp.plot(dP_dz, zs/1e3, "k", -rho*9.81, zs/1e3, "--rx")
    ax_dp.set_ylabel("height (km)")
    ax_dp.set_title("dP/dz (K/m)", loc='right', color='k')
    ax_dp.set_title("-$\\rho$g (K/m)", loc='left', color='r')
    ax_dp.set_ylim(0, 15)
    
    ax_resid.plot(1e3*(dP_dz - (-rho*9.81)), zs/1e3, "k")
    ax_resid.set_title("dP/dz + $\\rho$g (10$^{-3}$ K/m)", loc="right")

def gravity_wave_plot(zs, pressure, temperature, 
                      th_bnds=None, dth_bnds=None, c_bnds=None, 
                      nsq_bnds=None):
    """
    Compute the Brunt-Vaisala frequency from a given dry 
    atmosphere temperature and pressure profile, and use it to
    plot (1) $\theta$(z), (2) d$\theta$/dz, (3) gravity wave phase speed

    pressure: Pa
    temperature: K
    """
    nz = len(zs)
    nsq = np.zeros(nz)
    d_theta_dz = np.zeros_like(nsq)
    
    # exner function
    pi = (pressure/1e5)**(287./1004.)
    # potential temperature
    theta = temperature/pi

    fig, axes = plt.subplots(2, 2, figsize=(10,8), sharey=True)
    plt.subplots_adjust(wspace=0.1, hspace=0.2)
    ax_theta, ax_dth_dz, ax_cstar, ax_Nsq = axes.ravel()
    
    ax_theta.plot(theta, zs/1e3, "k")
    ax_theta.set_ylabel("height (km)")
    #ax_theta.set_xlabel("K")
    ax_theta.set_ylim(0, 15)
    ax_theta.set_title("$\\theta$ (K)", loc="right")
    if th_bnds: ax_theta.set_xlim(*th_bnds)
    
    for i in xrange(1, nz-1):
        d_theta_dz[i] = (theta[i+1] - theta[i-1])/(zs[i+1] - zs[i-1])
        nsq[i] = 9.81/theta[i]*d_theta_dz[i]

    d_theta_dz = np.ma.masked_equal(d_theta_dz, 0)
    nsq = np.ma.masked_equal(nsq, 0)
    
    ax_dth_dz.plot(d_theta_dz, zs/1e3, "k")
    #ax_dth_dz.set_xlabel("K/m")
    ax_dth_dz.set_title("d$\\theta$/dz (K/m)", loc="right")
    if dth_bnds: ax_dth_dz.set_xlim(*dth_bnds)

    ax_Nsq.plot(nsq*1e3, zs/1e3, "k")
    ax_Nsq.set_title("N$^2$ (10$^{-3}$ s$^{-1}$)", loc="right")
    if nsq_bnds: ax_Nsq.set_xlim(*nsq_bnds)
        
    cstar = np.sqrt(nsq)*zs[-1]/np.pi
    ax_cstar.plot(cstar, zs/1e3, "k")
    ax_cstar.set_ylabel("height (km)")
    ax_cstar.set_title("c$^*$ (m/s)", loc="right")
    if c_bnds: ax_cstar.set_xlim(*c_bnds)