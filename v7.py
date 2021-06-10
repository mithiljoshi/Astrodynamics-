from sys import path
import planet_data as pd
import tools as T
from OrbitPropagator import orbitPropagator as OP

tspan = 3600*24*1.0
dt = 10.0

cb = pd.earth

if __name__ == '__main__':
    op0 = OP(T.tle2coes('iss.txt'), tspan, dt, coes = True, deg = False)
    op1 = OP(T.tle2coes('COSMOS 2251 DEB.txt'), tspan, dt, coes = True, deg= False)

    T.plot_n_orbits([op0.rs, op1.rs], labels=['iss', 'debri'], show_plot=True)
