from sys import path
import planet_data as pd
import tools as T
from OrbitPropagator import orbitPropagator as OP
from OrbitPropagator import null_perts

tspan = 3600*24*1.0
dt = 10.0

cb = pd.earth

if __name__=='__main__':
    perts = null_perts()
    perts['J2']= True
    op = OP(T.tle2coes('iss.txt'), tspan, 100.0, coes = True, perts=perts)
    # op.plot3d(show_plot = True)
    op.calculate_coes()
    