#!/usr/bin/env python
# coding: utf-8

import numpy as np
from math import sqrt
import planet_data as pd
from OrbitPropagator import orbitPropagator as OP
import tools as T

tspan = 3600*24*1.0
dt = 100

if __name__ == '__main__':
    r_mag = pd.earth['radius']+1000
    v_mag = sqrt(pd.earth['mu']/r_mag)
    r0 = [r_mag, 0 ,0]
    v0 = [0, v_mag, 0]
    
    r_mag = pd.earth['radius']+1500
    v_mag = sqrt(pd.earth['mu']/r_mag)*1.3
    r00 = [r_mag, 0 ,0]
    v00 = [0, v_mag, 0]
    
    op0 = OP(r0,v0, tspan, dt)
    op00 = OP(r00,v00, tspan, dt)
    
    op0.propagate_orbit()
    op00.propagate_orbit()
    
    T.plot_n_orbits([op0.rs, op00.rs], labels = ['0', '1'], show_plot = True)

