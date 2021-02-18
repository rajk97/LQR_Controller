
from BuggySimulator import *
import numpy as np
import scipy
import math
from scipy.ndimage import gaussian_filter1d
from util import *
from scipy import signal
import matplotlib.pyplot as plt
    
class controller():
    
    def __init__(self,traj,vehicle, e=0, e111=0):
        self.vehicle=vehicle
        self.traj=traj
        self.prev_vx_error=0
        self.integral_vx_error=0        
        self.e = e
        self.e111 = e111
        # self.curv=self.compute_curvature()

    # def compute_curvature(self):

    #     sigma_gaus = 10
    #     traj=self.traj
    #     xp = scipy.ndimage.filters.gaussian_filter1d(input=traj[:,0], sigma=sigma_gaus,order=1)
    #     xpp = scipy.ndimage.filters.gaussian_filter1d(input=traj[:,0], sigma=sigma_gaus,order=2)
    #     yp = scipy.ndimage.filters.gaussian_filter1d(input=traj[:,1], sigma=sigma_gaus,order=1)
    #     ypp = scipy.ndimage.filters.gaussian_filter1d(input=traj[:,1], sigma=sigma_gaus,order=2)
    #     curv=np.zeros(len(traj))
    #     for i in range(len(xp)):
    #         curv[i] = (xp[i]*ypp[i] - yp[i]*xpp[i])/(xp[i]**2 + yp[i]**2)**1.5
    #     return curv
        
        

    
    
    def control_update(self):

        traj=self.traj
        vehicle=self.vehicle 
        
        lr = vehicle.lr
        lf = vehicle.lf
        Ca = vehicle.Ca
        Iz = vehicle.Iz
        f = vehicle.f
        m = vehicle.m
        g = vehicle.g

        delT = 0.05

        #reading current vehicle states
        X = vehicle.state.X
        Y = vehicle.state.Y
        xdot = vehicle.state.xd
        ydot = vehicle.state.yd
        phi = vehicle.state.phi
        phidot = vehicle.state.phid
        delta = vehicle.state.delta
        vx = xdot


        # ---------------|Lateral Controller|-------------------------
        
        dist, index1 = closest_node(X, Y, traj)

        n = 180

        index = index1 + n

        if(index>=8203):
            index = 8202

        X_new = traj[index][0]
        Y_new = traj[index][1]

        K_p = 2500
        K_i = 0.0001
        K_d = 1.5



        dist1 = math.sqrt((Y_new - Y)**2 + (X_new - X)**2)

        v_d = dist1/((n-80)*delT)
        v_d = 6

        v = math.sqrt(xdot**2 + ydot**2)

        # e22 = wrap2pi(phi - np.arctan2(Y_new-Y, X_new-X))

        


        
        # curv=self.curv
        # ------------|Lateral Controller|---------------------------------
       
        A = np.array([[0,1,0,0], [0,-(4*Ca)/(m*xdot),(4*Ca)/m,(-(2*Ca*lf)+(2*Ca*lr))/(m*xdot)], [0,0,0,1], [0,(-(2*Ca*lf)+(2*Ca*lr))/(Iz*xdot),((2*Ca*lf)-(2*Ca*lr))/(Iz),(-(2*Ca*lf*lf)-(2*Ca*lr*lr))/(Iz*xdot)]])

        # print(np.shape(A))
        B = np.array([[0], [(2*Ca)/m], [0], [(2*Ca*lf)/Iz]])
        C = np.identity(4)
        D = [[0],[0],[0],[0]]
        
        syscont = signal.StateSpace(A,B,C,D)
        
        sysdisc = syscont.to_discrete(delT)
        Ad = sysdisc.A
        Bd = sysdisc.B

        R = np.zeros((1, 1))
        R[0][0] = 0.06

        Q = 0.8*np.identity(4)
        
        S = np.matrix(scipy.linalg.solve_discrete_are(Ad, Bd, Q, R))

        K = -np.matrix(scipy.linalg.inv(Bd.T*S*Bd+R)*(Bd.T*S*Ad))

        

        

        
        r  = wrap2pi(phi - np.arctan2(Y_new-Y, X_new-X))
        e21 = 0.1
        e22 = ydot - (vx*r)
        e23 = r
        e24 = phidot

        e = np.array([[e21], [e22], [e23], [e24]])

        deltad = float(-K*(e))

        # delta = np.matmul(K, e)

        # deltad = -delta[0][0]
        












        
        #--------|Longitudinal Controller|------------------------------
        e11 = (v_d - v)
        
        F = K_p*(e11) + K_i*(self.e + e11) + K_d*(e11 - self.e111)
        self.e += e11
        self.e111 = e11 
        # -----------------------------------------------------------------

        # Communicating the control commands with the BuggySimulator
        controlinp = vehicle.command(F,deltad)
        # F: Force
        # deltad: desired rate of steering command
       

        return controlinp



