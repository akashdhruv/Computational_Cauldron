# Importing Libraries
import INS
import POISSON
import numpy as np
from math import *
import matplotlib.pyplot as plt

# Domain Length and Limits

Dx_min = -0.5
Dx_max =  0.5

Dy_min = -0.5
Dy_max =  0.5

Lx = Dx_max - Dx_min
Ly = Dy_max - Dy_min

# Block size

Nxb = 100
Nyb = 100

dx = Lx/Nxb
dy = Ly/Nyb

# physical variables

x = np.linspace(Dx_min,Dx_max,Nxb+1)
y = np.linspace(Dy_min,Dy_max,Nyb+1)

p = np.zeros((Nxb+2,Nyb+2),dtype=float)

u = np.zeros((Nxb+2,Nyb+2),dtype=float)
v = np.zeros((Nxb+2,Nyb+2),dtype=float)

ut = np.zeros((Nxb+2,Nyb+2),dtype=float)
vt = np.zeros((Nxb+2,Nyb+2),dtype=float)

[X,Y] = np.meshgrid(x,y)

u_old = np.zeros((Nxb+2,Nyb+2),dtype=float)
v_old = np.zeros((Nxb+2,Nyb+2),dtype=float)

G1 = np.zeros((Nxb,Nyb),dtype=float)
G2 = np.zeros((Nxb,Nyb),dtype=float)

G1_new = np.zeros((Nxb,Nyb),dtype=float)
G2_new = np.zeros((Nxb,Nyb),dtype=float)

p_RHS = np.zeros((Nxb,Nyb),dtype=float)

p_new = np.zeros((Nxb+2,Nyb+2),dtype=float)

# ins parameters

ins_inRe = 0.001  
ins_sig  = 0.01
ins_cfl  = 0.15

# driver parameters

dt_sig = ins_sig*(min(dx,dy)**2)/ins_inRe
dt_cfl = ins_cfl*min(dx,dy)

dt = min(dt_sig,dt_cfl)

t = 40.0

nt = int(t/dt)

Maxit = 1500

# Physics Squence
tstep = 0

while(tstep<=nt):

	u_old = u.copy()
        v_old = v.copy()
	
	#_____________________________Predictor_____________________________#

	G1_new,G2_new,ut,vt = INS.predictor(u,v,G1,G2,ins_inRe,dx,dy,dt,tstep,Nxb,Nyb)

        G1 = G1_new
        G2 = G2_new

        #__________________Predictor Boundary Conditions_____________________#

        # LOW X
        ut[0,:]  =  0.0
        vt[0,:]  = -vt[1,:]

        # HIGH X
        ut[-2,:] =  0.0
        ut[-1,:] =  0.0	
        vt[-1,:] = -vt[-2,:]

        # LOW Y
        vt[:,0]  =  0.0
        ut[:,0]  = -ut[:,1]

        # HIGH Y
        vt[:,-1] =  0.0
        vt[:,-2] =  0.0
        ut[:,-1] = 2.0 -ut[:,-2]

        #_____________________________Poisson Solver________________________#

        p_RHS = -((1/(dy*dt))*(vt[1:-1,1:-1]-vt[1:-1,:-2]))-((1/(dx*dt))*(ut[1:-1,1:-1]-ut[:-2,1:-1]))

        p_counter = 0

        while(p_counter < Maxit):

        	p_new = POISSON.solver(p,p_RHS,dx,dy,Nxb,Nyb)

		#____________________Poisson Boundary Conditions____________#

                # LOW X
		p_new[0,:]  =  p_new[1,:]
                       
		# HIGH X
		p_new[-1,:] =  p_new[-2,:]

		# LOW Y
		p_new[:,0]  =  p_new[:,1] 

		# HIGH Y
		p_new[:,-1] =  p_new[:,-2]  

                p_res = sqrt(np.sum((p_new-p)**2)/np.size(p))

                p = p_new

                p_counter += 1

                if(p_res<10**-7 and p_res != 0.):
			break

	#________________________________Corrector____________________________#

	u,v = INS.corrector(ut,vt,p,dt,dx,dy,Nxb,Nyb)

        #__________________Corrector Boundary Conditions_____________________#

        # LOW X
	u[0,:]  =  0.0
	v[0,:]  = -v[1,:]

        # HIGH X
	u[-2,:] =  0.0
	u[-1,:] =  0.0
	v[-1,:] = -v[-2,:]

        # LOW Y
	v[:,0]  =  0.0
	u[:,0]  = -u[:,1]

        # HIGH Y
	v[:,-1] =  0.0
	v[:,-2] =  0.0
	u[:,-1] = 2.0 -u[:,-2]

	u_res = sqrt(np.sum((u_old-u)**2)/np.size(u))
	v_res = sqrt(np.sum((v_old-v)**2)/np.size(v))

	print "---------------------------PARAMETER DISPLAY-----------------------"
	print "Simulation Time     : ",tstep*dt," s"
	print "U velocity Residual : ",u_res
	print "V velocity Residual : ",v_res
	print "Pressure Residual   : ",p_res
	print "Poisson Counter     : ",p_counter

	tstep += 1

	if(u_res<10**-8 and u_res != 0. and v_res<10**-8 and v_res != 0.):
		break

uu = 0.5*(u[:-1,:-1] + u[:-1,1:])
vv = 0.5*(v[:-1,:-1] + v[1:,:-1])
pp = 0.25*(p[:-1,:-1] + p[1:,:-1] + p[:-1,1:] + p[1:,1:])

X  = X.T
Y  = Y.T

plt.figure()
plt.contourf(X,Y,np.sqrt(uu**2+vv**2))
plt.axis('equal')

plt.figure()
plt.contourf(X,Y,pp)
plt.axis('equal')

plt.show()
