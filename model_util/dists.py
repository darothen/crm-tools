""" Collection of functions and utilities for building and analzying 
hydrometeor size distributions in the CRM.

"""

def cal_lambda(b, beta, am0, q, n):
    gamma = b + beta
    for m in xrange(b+beta-1, b, -1):
        gamma *= m
    return (am0*float(gamma)*n/q)**(1.0/float(beta))