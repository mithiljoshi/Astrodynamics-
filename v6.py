#!/usr/bin/env python
# coding: utf-8

import numpy as np
import tools as T
import matplotlib.pyplot as plt
from OrbitPropagator import orbitPropagator as OP
import planet_data as pd

tspan= 3600*24*1.0
dt = 10.0

cb =pd.earth

if __name__ == "__main__":
    #ISS
    c0 =[cb['radius']+1000.0, 0.0006189, 51.6393,0.0,234.1955,105.6372]
    #GEO
    c1 = [cb['radius']+35000, 0.0,0.0,0.0,0.0,0.0]

    op0= OP(c0, tspan, dt, coes = True)
    op1 = OP(c1, tspan, dt,coes=True)

    op0.propagate_orbit()
    op1.propagate_orbit()
    T.plot_n_orbits([op0.rs, op1.rs], labels=['ISS', 'GEO'], show_plot=True)
    



