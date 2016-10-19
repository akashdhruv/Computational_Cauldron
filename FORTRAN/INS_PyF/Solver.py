# Importing Libraries
import INS
import POISSON
import HEAT
import numpy as np
from math import *
from mpi4py import MPI
import matplotlib.pyplot as plt


# Defining MPI communication function

def MPI_applyBC(u,x_id,y_id,x_procs,y_procs,x_comm,y_comm):

	if(x_procs > 1):

		if(x_id%2 == 0):

			if(x_id == 0):

				x_comm.send(u[-2,:],dest=(x_id+1)%x_procs,tag=1)
				u[-1,:] = x_comm.recv(source=(x_id+1)%x_procs,tag=2)

			elif(x_id == nblockx - 1):

				x_comm.send(u[1,:],dest=(x_id-1+x_procs)%x_procs,tag=3)
				u[0,:] = x_comm.recv(source=(x_id-1+x_procs)%x_procs,tag=4)

			else:
				x_comm.send(u[-2,:],dest=(x_id+1)%x_procs,tag=1)
				u[-1,:] = x_comm.recv(source=(x_id+1)%x_procs,tag=2)

				x_comm.send(u[1,:],dest=(x_id-1+x_procs)%x_procs,tag=3)
				u[0,:] = x_comm.recv(source=(x_id-1+x_procs)%x_procs,tag=4)

		elif(x_id%2 == 1):

			if(x_id == nblockx - 1):
				
				x_comm.send(u[1,:],dest=(x_id-1+x_procs)%x_procs,tag=2)
				u[0,:] = x_comm.recv(source=(x_id-1+x_procs)%x_procs,tag=1)

			else:

				x_comm.send(u[1,:],dest=(x_id-1+x_procs)%x_procs,tag=2)
				u[0,:] = x_comm.recv(source=(x_id-1+x_procs)%x_procs,tag=1)
				
				x_comm.send(u[-2,:],dest=(x_id+1)%x_procs,tag=4)	
				u[-1,:] = x_comm.recv(source=(x_id+1)%x_procs,tag=3)

	if(y_procs > 1):

		if(y_id%2 == 0):

			if(y_id == 0):

				y_comm.send(u[:,-2],dest=(y_id+1)%y_procs,tag=5)
				u[:,-1] = y_comm.recv(source=(y_id+1)%y_procs,tag=6)

			elif(y_id == nblocky - 1):

				y_comm.send(u[:,1],dest=(y_id-1+y_procs)%y_procs,tag=7)
				u[:,0] = y_comm.recv(source=(y_id-1+y_procs)%y_procs,tag=8)

			else:

				y_comm.send(u[:,-2],dest=(y_id+1)%y_procs,tag=5)
				u[:,-1] = y_comm.recv(source=(y_id+1)%y_procs,tag=6)

				y_comm.send(u[:,1],dest=(y_id-1+y_procs)%y_procs,tag=7)
				u[:,0] = y_comm.recv(source=(y_id-1+y_procs)%y_procs,tag=8)


		elif(y_id%2 == 1):

			if(y_id == nblocky - 1):

				y_comm.send(u[:,1],dest=(y_id-1+y_procs)%y_procs,tag=6)
				u[:,0] = y_comm.recv(source=(y_id-1+y_procs)%y_procs,tag=5)
		
			else:

				y_comm.send(u[:,1],dest=(y_id-1+y_procs)%y_procs,tag=6)
				u[:,0] = y_comm.recv(source=(y_id-1+y_procs)%y_procs,tag=5)


				y_comm.send(u[:,-2],dest=(y_id+1)%y_procs,tag=8)
				u[:,-1] = y_comm.recv(source=(y_id+1)%y_procs, tag=7)

	return u

# Initializing MPI environment
nblockx = 2
nblocky = 2

comm = MPI.COMM_WORLD
myid = comm.Get_rank()
procs = comm.Get_size()

x_comm = comm.Split(myid/nblockx,myid%nblockx)
y_comm = comm.Split(myid%nblockx,myid/nblockx)

x_id = x_comm.Get_rank()
x_procs = x_comm.Get_size()

y_id = y_comm.Get_rank()
y_procs = y_comm.Get_size()

# Domain Length and Limits

Dx_min = -0.5
Dx_max =  0.5

Dy_min = -0.5
Dy_max =  0.5

Lx = Dx_max - Dx_min
Ly = Dy_max - Dy_min

gr_Lx = Lx/nblockx
gr_Ly = Ly/nblocky

# Block size

Nxb = 50
Nyb = 40

dx = gr_Lx/Nxb
dy = gr_Ly/Nyb

# physical variables

x = Dx_min + (myid%nblockx)*gr_Lx + dx*np.linspace(0,Nxb,Nxb+1)

y = Dy_min + (myid/nblockx)*gr_Ly + dy*np.linspace(0,Nyb,Nyb+1)

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

T = np.zeros((Nxb+2,Nyb+2),dtype=float)
T_new = np.zeros((Nxb+2,Nyb+2),dtype=float)

T[:,:] = 313.0

# ins parameters

ins_inRe = 0.001  
ins_sig  = 0.01
ins_cfl  = 0.15

# heat parameters

ht_Pr = 0.7

# driver parameters

dt_sig = ins_sig*(min(dx,dy)**2)/ins_inRe
dt_cfl = ins_cfl*min(dx,dy)

dt = min(dt_sig,dt_cfl)

t = 40.0

nt = int(t/dt)

Maxit = 1500

p_res = 0.
u_res = 0.
v_res = 0.
T_res = 0.

ins_p_res = 0.
ins_v_res = 0.
ins_v_res = 0.
ins_T_res = 0.

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

	ut = MPI_applyBC(ut,x_id,y_id,x_procs,y_procs,x_comm,y_comm)
	vt = MPI_applyBC(vt,x_id,y_id,x_procs,y_procs,x_comm,y_comm)

        # LOW X
	if(x_id == 0):
		ut[0,:]  =  0.0
        	vt[0,:]  = -vt[1,:]

        # HIGH X
	if(x_id == nblockx-1):
        	ut[-2,:] =  0.0
        	ut[-1,:] =  0.0	
        	vt[-1,:] = -vt[-2,:]

        # LOW Y
	if(y_id == 0):
        	vt[:,0]  =  0.0
        	ut[:,0]  = -ut[:,1]

        # HIGH Y
	if(y_id == nblocky-1):
        	vt[:,-1] =  0.0
        	vt[:,-2] =  0.0
        	ut[:,-1] = 2.0 -ut[:,-2]

        #_____________________________Poisson Solver________________________#

        p_RHS = -((1/(dy*dt))*(vt[1:-1,1:-1]-vt[1:-1,:-2]))-((1/(dx*dt))*(ut[1:-1,1:-1]-ut[:-2,1:-1]))

        p_counter = 0

        while(p_counter < Maxit):

        	p_new = POISSON.solver(p,p_RHS,dx,dy,Nxb,Nyb)

		#___________________Pressure Boundary Conditions____________#

		p_new = MPI_applyBC(p_new,x_id,y_id,x_procs,y_procs,x_comm,y_comm)

                # LOW X
		if(x_id == 0):	
			p_new[0,:]  =  p_new[1,:]
                       
		# HIGH X
		if(x_id == nblockx-1):
			p_new[-1,:] =  p_new[-2,:]

		# LOW Y
		if(y_id == 0):
			p_new[:,0]  =  p_new[:,1] 

		# HIGH Y
		if(y_id == nblocky-1):
			p_new[:,-1] =  p_new[:,-2]  
                 
                #_________________Residuals and Convergence Check__________#

                p_res = np.sum((p_new-p)**2)

                p = p_new

                p_counter += 1
              
		ins_p_res = comm.allreduce(p_res, op=MPI.SUM)
		ins_p_res = sqrt(ins_p_res/(np.size(p)*procs))

                if(ins_p_res<10**-6 and ins_p_res != 0.):
			break

	#________________________________Corrector____________________________#

	u,v = INS.corrector(ut,vt,p,dt,dx,dy,Nxb,Nyb)

        #__________________Corrector Boundary Conditions_____________________#

	u = MPI_applyBC(u,x_id,y_id,x_procs,y_procs,x_comm,y_comm)
	v = MPI_applyBC(v,x_id,y_id,x_procs,y_procs,x_comm,y_comm)

        # LOW X
	if(x_id == 0):
		u[0,:]  =  0.0
		v[0,:]  = -v[1,:]

        # HIGH X
	if(x_id == nblockx - 1):
		u[-2,:] =  0.0
		u[-1,:] =  0.0
		v[-1,:] = -v[-2,:]

        # LOW Y
	if(y_id == 0):
		v[:,0]  =  0.0
		u[:,0]  = -u[:,1]

        # HIGH Y
	if(y_id == nblocky - 1):
		v[:,-1] =  0.0
		v[:,-2] =  0.0
		u[:,-1] = 2.0 -u[:,-2]

        #___________________________Residuals_______________________________#

	u_res = np.sum((u_old-u)**2)
	v_res = np.sum((v_old-v)**2)

	ins_u_res = comm.allreduce(u_res, op=MPI.SUM)
	ins_u_res = sqrt(ins_u_res/(np.size(p)*procs))

	ins_v_res = comm.allreduce(v_res, op=MPI.SUM)
	ins_v_res = sqrt(ins_v_res/(np.size(p)*procs))

        #_______________________Heat Advection Diffusion____________________#

        T_new = HEAT.tempsolver(T,u,v,dx,dy,dt,ins_inRe,ht_Pr,Nxb,Nyb)

        #____________________Temperature Boundary Conditions________________#

	T_new = MPI_applyBC(T_new,x_id,y_id,x_procs,y_procs,x_comm,y_comm)

	# LOW X
	if(x_id == 0):
		T_new[0,:]  =  T_new[1,:]

	# HIGH X
	if(x_id == nblockx - 1):
		T_new[-1,:] =  T_new[-2,:]

	# LOW Y
	if(y_id == 0):
		T_new[:,0]  =  T_new[:,1]

	# HIGH Y
	if(y_id == nblocky - 1):
		T_new[:,-1] = 2*383.15 - T_new[:,-2]

        #___________________________Residuals_______________________________#

        T_res = np.sum((T_new-T)**2)

	ins_T_res = comm.allreduce(T_res, op=MPI.SUM)
	ins_T_res = sqrt(ins_T_res/(np.size(T)*procs))
	
	T = T_new

	if(myid == 0 and tstep%5 == 0):

		print "---------------------------PARAMETER DISPLAY-----------------------"
		print "Simulation Time     : ",tstep*dt," s"
		print "U velocity Residual : ",ins_u_res
		print "V velocity Residual : ",ins_v_res
        	print "Temperature Residual: ",ins_T_res
		print "Pressure Residual   : ",ins_p_res
		print "Poisson Counter     : ",p_counter

	tstep += 1

	if(ins_u_res<10**-7 and ins_u_res != 0. and ins_v_res<10**-7 and ins_v_res != 0.):
		break

uu = 0.5*(u[:-1,:-1] + u[:-1,1:])
vv = 0.5*(v[:-1,:-1] + v[1:,:-1])
pp = 0.25*(p[:-1,:-1] + p[1:,:-1] + p[:-1,1:] + p[1:,1:])
tt = 0.25*(T[:-1,:-1] + T[1:,:-1] + T[:-1,1:] + T[1:,1:])

X  = X.T
Y  = Y.T

"""
plt.figure()
plt.contourf(X,Y,np.sqrt(uu**2+vv**2))
plt.axis('equal')

plt.figure()
plt.contourf(X,Y,pp)
plt.axis('equal')

plt.figure()
plt.contourf(X,Y,tt)
plt.axis('equal')

plt.show()

"""
