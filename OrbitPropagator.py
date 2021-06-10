#!/usr/bin/env python
# coding: utf-8
from math import degrees
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import planet_data as pd
import tools as T

def null_perts():
    return {
        'J2': False,
        'aero': False,
        'moon_grav': False,
        'solar_grav':False
    }

class orbitPropagator:
    def __init__(self,state0,tspan, dt,coes= False, deg = True, cb = pd.earth, perts = null_perts()):
        if coes:
            self.r0, self.v0 = T.coes2rv(state0, deg = deg, mu = cb['mu'])
        else:
            self.r0 = state0[:3]
            self.v0 = state0[3:]
        
        self.y0 = self.r0.tolist()+self.v0.tolist()
        self.tspan = tspan
        self.dt = dt
        self.cb = cb
        #total number of steps
        self.n_steps = int(np.ceil(self.tspan/self.dt))
        #initialize arrays
        self.ys = np.zeros((self.n_steps, 6))
        self.ts = np.zeros((self.n_steps, 1))

        #initial conditions
        self.y0  = [self.r0[0],self.r0[1], self.r0[2], self.v0[0], self.v0[1],self.v0[2]]
        self.ys[0] = np.array(self.y0)
        self.step=1

        #initialize solver
        self.solver = ode(self.diff_eqn)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.y0, 0)
        
        # define perturbations dictionary
        self.perts = perts
        
        self.propagate_orbit()

    def propagate_orbit(self):
        #propogate orbit
        while self.solver.successful() and self.step< self.n_steps:
            self.solver.integrate(self.solver.t+self.dt)
            self.ts[self.step] =  self.solver.t
            self.ys[self.step] = self.solver.y
            self.step+=1
        self.rs = self.ys[:,:3]
        self.vs = self.ys[:,3:]
    
    def diff_eqn(self, t, y):  
        # unpack the state
        rx, ry, rz, vx, vy, vz = y
        r = np.array([rx, ry, rz])

        #norm of radius vector
        norm_r = np.linalg.norm(r)

        #two body acceleration
        a = -r*self.cb['mu']/norm_r**3

        if self.perts['J2']:
            z2 = r[2]**2
            r2 = norm_r**2
            tx = r[0]/norm_r*(5*z2/r2-1)
            ty = r[1]/norm_r*(5*z2/r2-1)
            tz = r[2]/norm_r*(5*z2/r2-1)

            a_j2 = 1.5*self.cb['J2']*self.cb['mu']*self.cb['radius']**2/norm_r**4*np.array([tx,ty,tz])
            a+=a_j2
        return [vx,vy,vz, a[0], a[1], a[2]]
    
    def plot3d (self, show_plot = False ,save_plot = False ):
        fig=plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111, projection='3d')
        # plot trajectory
        ax.plot(self.rs[:,0], self.rs[:,1], self.rs[:,2], 'k', label = 'trajectory')
        ax.plot([self.rs[0,0]], [self.rs[0,1]], [self.rs[0,2]], 'k*', label = 'initial position')

        #plot central body
        _u, _v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        _x = self.cb['radius']*np.cos(_u)*np.sin(_v)
        _y = self.cb['radius']*np.sin(_u)*np.sin(_v)
        _z = self.cb['radius']*np.cos(_v)
        ax.plot_wireframe(_x,_y,_z, cmap = "Blues")

        #plot x, y, z, vectors
        l = self.cb['radius']*2
        x,y,z = [[0,0,0], [0,0,0], [0,0,0]]
        u,v,w = [[1,0,0], [0,1,0], [0,0,1]]
        ax.quiver(x,y,z,u,v,w,color='k')

        max_val = np.max(np.abs(self.rs))

        ax.set_xlim([-max_val, max_val])
        ax.set_ylim([-max_val, max_val])
        ax.set_zlim([-max_val, max_val])

        ax.set_xlabel(['X (km)'])
        ax.set_ylabel(['Y (km)'])
        ax.set_zlabel(['Z (km)'])

        #ax.set_aspect('equal')

        ax.set_title('Two Body Visual')
        plt.legend(['Trajectory', 'Initial Position'])
        for ii in range(0,360,1):
            ax.view_init(elev=270,azim=ii)
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig('Two body'+'.png', dpi = 300)

    def calculate_coes(self):
        print('Calculating COEs...')

        self.coes = np.zeros((self.n_steps,6))

        for n in range(self.n_steps):
            self.coes[n, :] = T.rv2coes(self.rs[n,:], self.vs[n,:], mu = self.cb['mu'], degrees = degrees)

    def plot_coes(self, hours = False, days = False, show_plot = False, save_plot = False, title = 'COEs', figsize = (16,6)):
        print("Plotting COEs...")

        fig, axs = plt.subplots(nrows=2, ncols=3, figsize = figsize)
        fig.suptitle(title, fontsize = 20)

        # x axis
        if hours:
            ts = self.ts/3600.0
            xlabel = 'Time Lapsed (hours)'
        elif self.days:
            ts = self.ts/3600.0/24.0
            xlabel = 'Time Lapsed (days)'
        else:
            ts = self.ts
            xlabel = 'Time Lapsed (seconds)'



        #plot true anomaly
        axs[0,0].plot(self.ts, self.coes[:,3])
        axs[0,0].set_title('True Anomaly vs. Time')
        axs[0,0].grid(True)
        axs[0,0].set_ylabel('Angle (degrees)')

        #plot semi major axis
        axs[1,0].plot(self.ts, self.coes[:,0])
        axs[1,0].set_title('Semi Major Axis vs. Time')
        axs[1,0].grid(True)
        axs[1,0].set_ylabel('Semi Major Axis (km')
        axs[1,0].set_xlabel(xlabel)

        #plot eccentricity 
        axs[0,1].plot(self.ts, self.coes[:,1])
        axs[0,1].set_title('Eccentricity vs. Time')
        axs[0,1].grid(True)

        #plot argument of periapse
        axs[0,2].plot(self.ts, self.coes[:,4])
        axs[0,2].set_title('Argument of Periapse vs. Time')
        axs[0,2].grid(True)

        #plot inclination
        axs[1,1].plot(self.ts, self.coes[:,2])
        axs[1,1].set_title('Inclination vs. Time')
        axs[1,1].grid(True)
        axs[1,1].set_ylabel('Angle (degrees)')
        axs[1,1].set_xlabel(xlabel)

        #plot RAAN
        axs[1,2].plot(self.ts, self.coes[:,5])
        axs[1,2].set_title('RAAN vs. Time')
        axs[1,2].grid(True)
        axs[1,2].set_xlabel(xlabel)

        if show_plot:
            plt.show()
        