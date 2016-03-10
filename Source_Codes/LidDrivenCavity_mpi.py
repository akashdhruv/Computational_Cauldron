import numpy as np
from math import *
from mpi4py import MPI
import time

## Parallel Navier Stokes Solver

# Initialization MPI environment

comm=MPI.COMM_WORLD
rank=MPI.COMM_WORLD.Get_rank()
size=MPI.COMM_WORLD.Get_size()
name=MPI.Get_processor_name()

if rank==0:
	t1=time.time()
comm.Barrier()

# Grid Parameters

Lx=1.
Ly=1.  

Nxb=128
Nyb=32

iProcs=1
jProcs=4

tProcs=iProcs*jProcs

Nx=Nxb*iProcs
Ny=Nyb*jProcs

dx=Lx/(Nx)
dy=Ly/(Ny)

rms_error=0.
norm_u=0.
norm_v=0.

# Global Constants

inRe=0.01
U=1.0

w=1. # Relaxation factor

# Staggered Grid

x=np.linspace(0.,Lx,Nxb+1)

#y=((Nyb)*rank)*dy+np.linspace(0,Nyb*dy,Nyb+1)
y= np.linspace(0.,1.,Nx)[(Nyb*rank):(Nyb*rank)+Nyb]

#p=np.loadtxt('Re100P%d.dat' % rank)
p=np.zeros((Nxb+2,Nyb+2),dtype=float)
p_old=np.empty_like(p)

#u=np.loadtxt('Re100U%d.dat' % rank)
u=np.zeros((Nxb+1,Nyb+2),dtype=float)
ut=np.zeros((Nxb+1,Nyb+2),dtype=float)

#v=np.loadtxt('Re100V%d.dat' % rank)
v=np.zeros((Nxb+2,Nyb+2),dtype=float)
vt=np.zeros((Nxb+2,Nyb+2),dtype=float)

[X,Y]=np.meshgrid(x,y)

# Differential equation functions

def F1C(ue,uw,us,un,vs,vn,dx,dy):
    	F1C=-((ue**2)-(uw**2))/dx-((un*vn)-(us*vs))/dy
    	return F1C

def FV(uP,uE,uW,uN,uS,dx,dy,inRe):
    	FV=(inRe/dx)*(((uE-uP)/dx)-((uP-uW)/dx))+(inRe/dy)*(((uN-uP)/dy)-((uP-uS)/dy))
    	return FV

def F2C(vn,vs,ve,vw,ue,uw,dx,dy):
    	F2C=-((ue*ve)-(uw*vw))/dx-((vn**2)-(vs**2))/dy
    	return F2C

# Simulation Parameters   

t=0
tmax=300.

dt=0.5*(dx**2)*(dy**2)/(inRe*(dx**2+dy**2))
#dt=0.5*(min(dx,dy))**2/inRe
#dt=0.5*min(dx,dy)/U

nt=int((tmax-t)/dt)

maxit=1500

comm.Barrier()

tstep=0

# Simulation Start
while(tstep<nt):

	u_old=u.copy()
	v_old=v.copy()
	
	# Predictor 
        ue1=(u[1:-1,1:-1]+u[2:,1:-1])/2
        uw1=(u[1:-1,1:-1]+u[:-2,1:-1])/2
        us1=(u[1:-1,1:-1]+u[1:-1,:-2])/2
        un1=(u[1:-1,1:-1]+u[1:-1,2:])/2
        vs1=(v[1:-2,:-2]+v[2:-1,:-2])/2
        vn1=(v[1:-2,1:-1]+v[2:-1,1:-1])/2
        
	if tstep==0 or tstep:
        	G1=(F1C(ue1,uw1,us1,un1,vs1,vn1,dx,dy)+FV(u[1:-1,1:-1],u[2:,1:-1],u[:-2,1:-1],u[1:-1,2:],u[1:-1,:-2],dx,dy,inRe))
        	ut[1:-1,1:-1]=u[1:-1,1:-1]+(dt/1)*G1
        	G1_old=G1.copy()
	else:
		G1=(F1C(ue1,uw1,us1,un1,vs1,vn1,dx,dy)+FV(u[1:-1,1:-1],u[2:,1:-1],u[:-2,1:-1],u[1:-1,2:],u[1:-1,:-2],dx,dy,inRe))
                ut[1:-1,1:-1]=u[1:-1,1:-1]+(dt/2)*(3*G1-G1_old)
                G1_old=G1.copy()	
        
        ve12=(v[1:-1,1:-1]+v[2:,1:-1])/2
        vw12=(v[1:-1,1:-1]+v[:-2,1:-1])/2
        vs12=(v[1:-1,1:-1]+v[1:-1,:-2])/2
        vn12=(v[1:-1,1:-1]+v[1:-1,2:])/2
        us12=(u[:-1,1:-1]+u[:-1,2:])/2
        un12=(u[1:,1:-1]+u[1:,2:])/2
        
	if tstep==0 or tstep:
        	G2=(F2C(ve12,vw12,vs12,vn12,us12,un12,dx,dy)+FV(v[1:-1,1:-1],v[2:,1:-1],v[:-2,1:-1],v[1:-1,2:],v[1:-1,:-2],dx,dy,inRe))
        	vt[1:-1,1:-1]=v[1:-1,1:-1]+(dt/1)*G2
		G2_old=G2.copy()
	else:
		G2=(F2C(ve12,vw12,vs12,vn12,us12,un12,dx,dy)+FV(v[1:-1,1:-1],v[2:,1:-1],v[:-2,1:-1],v[1:-1,2:],v[1:-1,:-2],dx,dy,inRe))
                vt[1:-1,1:-1]=v[1:-1,1:-1]+(dt/2)*(3*G2-G2_old)
                G2_old=G2.copy()
	
	# MPI and Domain Boundary Conditions
                                                                        
        vt[0,:]=0-vt[1,:] # Domain BC
        vt[-1,:]=0-vt[-2,:] # Domain BC
	ut[0,:]=0.
	ut[-1,:]=0.

	comm.Barrier()
                                                                                    
	if(rank%2==0):
                if rank==0:
			ut[:,0]=2*0-ut[:,1] # Domain BC
                        comm.send(ut[:,-2],dest=rank+1,tag=1)
                        ut[:,-1]=comm.recv(source=rank+1,tag=2)
			vt[:,0]=0.
             
                elif rank==size-1:
			ut[:,-1]=2*U-ut[:,-2] # Domain BC
			comm.send(ut[:,1],dest=rank-1,tag=7)
                        ut[:,0]=comm.recv(source=rank-1,tag=8)
			vt[:,-2]=0.
			vt[:,-1]=vt[:,-2]
                else:
                        comm.send(ut[:,-2],dest=rank+1,tag=1)
                        ut[:,-1]=comm.recv(source=rank+1,tag=2)
                     
                        comm.send(ut[:,1],dest=rank-1,tag=7)
                        ut[:,0]=comm.recv(source=rank-1,tag=8)

	elif(rank%2==1):
                if rank==size-1:
			ut[:,-1]=2*U-ut[:,-2] # Domain BC
                        ut[:,0]=comm.recv(source=rank-1,tag=1)
                        comm.send(ut[:,1],dest=rank-1,tag=2)
			vt[:,-2]=0.
			vt[:,-1]=vt[:,-2]
                else:
                        ut[:,0]=comm.recv(source=rank-1,tag=1)
                        comm.send(ut[:,1],dest=rank-1,tag=2)
            
                        ut[:,-1]=comm.recv(source=rank+1,tag=7)
                        comm.send(ut[:,-2],dest=rank+1,tag=8)

        comm.Barrier()

        # Data exchange between processes

        if(rank%2==0):
                if rank==0:
                        comm.send(vt[:,-2],dest=rank+1,tag=5)
                        vt[:,-1]=comm.recv(source=rank+1,tag=6)
                elif rank==size-1:
                        comm.send(vt[:,1],dest=rank-1,tag=11)
                        vt[:,0]=comm.recv(source=rank-1,tag=12)
                else:
                        comm.send(vt[:,-2],dest=rank+1,tag=5)
                        vt[:,-1]=comm.recv(source=rank+1,tag=6)
                        comm.send(vt[:,1],dest=rank-1,tag=11)
                        vt[:,0]=comm.recv(source=rank-1,tag=12)

        elif(rank%2==1):

                if rank==size-1:
                        vt[:,0]=comm.recv(source=rank-1,tag=5)
                        comm.send(vt[:,1],dest=rank-1,tag=6)
                else:
                        vt[:,0]=comm.recv(source=rank-1,tag=5)
                        comm.send(vt[:,1],dest=rank-1,tag=6)
                        vt[:,-1]=comm.recv(source=rank+1,tag=11)
                        comm.send(vt[:,-2],dest=rank+1,tag=12)

        comm.Barrier()

	# Poisson Solver

	p_counter=0

	while (p_counter<maxit):
	
		p_old=p.copy()

		p[1:Nx+1,1:Nyb+1]=w*(1/((2/(dx*dx))+(2/(dy*dy))))*((p_old[1:Nxb+1,2:Nyb+2]+p_old[1:Nxb+1,:Nyb])/(dy*dy)+\
				  (p_old[2:Nxb+2,1:Nyb+1]+p_old[:Nxb,1:Nyb+1])/(dx*dx)-(1/(dy*dt))*(vt[1:Nxb+1,1:Nyb+1]-\
				vt[1:Nxb+1,:Nyb])-(1/(dx*dt))*(ut[1:Nxb+1,1:Nyb+1]-ut[:Nxb,1:Nyb+1]))+(1-w)*p_old[1:Nxb+1,1:Nyb+1]
		
		comm.Barrier()

		# Pressure BC
		if(rank%2==0):

			if(rank==0):
				p[:,0]=p[:,1]	
				comm.send(p[:,-2],dest=rank+1,tag=21)
				p[:,-1]=comm.recv(source=rank+1,tag=22)
				
			elif(rank==size-1):
				p[:,-1]=p[:,-2]	
				comm.send(p[:,1],dest=rank-1,tag=23)
				p[:,0]=comm.recv(source=rank-1,tag=24)

			else:
				comm.send(p[:,-2],dest=rank+1,tag=21)
				p[:,-1]=comm.recv(source=rank+1,tag=22)
				comm.send(p[:,1],dest=rank-1,tag=23)
				p[:,0]=comm.recv(source=rank-1,tag=24)	


		elif(rank%2==1):
		
			if(rank==size-1):
				p[:,-1]=p[:,-2]
				p[:,0]=comm.recv(source=rank-1,tag=21)
				comm.send(p[:,1],dest=rank-1,tag=22)

			else:
				p[:,0]=comm.recv(source=rank-1,tag=21)
				comm.send(p[:,1],dest=rank-1,tag=22)
				p[:,-1]=comm.recv(source=rank+1,tag=23)
				comm.send(p[:,-2],dest=rank+1,tag=24)        
                
		comm.Barrier()

                                                   
                p[0,:]=p[1,:]
                p[-1,:]=p[-2,:]

		p_counter+=1

		rms_error = sqrt(sum(sum((p-p_old)**2))/np.size(p))
	
		comm.Barrier()
		v_rms_error=comm.gather(rms_error, root = 0)
	
		if(rank==0):
			if(max(v_rms_error)<10**-6  and max(v_rms_error)!=0):
				exit_val=1
			else:
				exit_val=0

		else: exit_val = None
		
		comm.Barrier()
		exit_val=comm.bcast(exit_val, root=0)

		if(exit_val==1):
			break

	# Corrector

	u[1:-1,1:-1]=ut[1:-1,1:-1]-(dt/(dx))*(p[2:-1,1:-1]-p[1:-2,1:-1])
	v[1:-1,1:-1]=vt[1:-1,1:-1]-(dt/(dy))*(p[1:-1,2:]-p[1:-1,1:-1])
	
 	# MPI and Domain Boundary Conditions

	comm.Barrier()

	if(rank%2==0):
		if rank==0:
			u[:,0]=2*0-u[:,1] # Domain BC
			comm.send(u[:,-2],dest=rank+1,tag=1)
			u[:,-1]=comm.recv(source=rank+1,tag=2)
			v[:,0]=0. 		
	
		elif rank==size-1:
			u[:,-1]=2*U-u[:,-2] # Domain BC
			comm.send(u[:,1],dest=rank-1,tag=7)
			u[:,0]=comm.recv(source=rank-1,tag=8)
			v[:,-2]=0.
			v[:,-1]=v[:,-2]
		else:
			comm.send(u[:,-2],dest=rank+1,tag=1)
			u[:,-1]=comm.recv(source=rank+1,tag=2)

			comm.send(u[:,1],dest=rank-1,tag=7)
			u[:,0]=comm.recv(source=rank-1,tag=8)
 
         
	elif(rank%2==1):
		if rank==size-1: 
			u[:,-1]=2*U-u[:,-2] # Domain BC
			u[:,0]=comm.recv(source=rank-1,tag=1)
			comm.send(u[:,1],dest=rank-1,tag=2)
			v[:,-2]=0.
			v[:,-1]=v[:,-2]
		else:
			u[:,0]=comm.recv(source=rank-1,tag=1)
			comm.send(u[:,1],dest=rank-1,tag=2)

			u[:,-1]=comm.recv(source=rank+1,tag=7)
			comm.send(u[:,-2],dest=rank+1,tag=8)
	
	comm.Barrier()                                              
	
	v[0,:]=0-v[1,:]
        v[-1,:]=0-v[-2,:]
	u[0,:]=0.
	v[-1,:]=0.

	comm.Barrier()

        if(rank%2==0):
                if rank==0:
                        comm.send(v[:,-2],dest=rank+1,tag=21)
                        v[:,-1]=comm.recv(source=rank+1,tag=22)
                elif rank==size-1:
                        comm.send(v[:,1],dest=rank-1,tag=23)
                        v[:,0]=comm.recv(source=rank-1,tag=24)
                else:
                        comm.send(v[:,-2],dest=rank+1,tag=21)
                        v[:,-1]=comm.recv(source=rank+1,tag=22)
                        comm.send(v[:,1],dest=rank-1,tag=23)
                        v[:,0]=comm.recv(source=rank-1,tag=24)

        elif(rank%2==1):

                if rank==size-1:
                        v[:,0]=comm.recv(source=rank-1,tag=21)
                        comm.send(v[:,1],dest=rank-1,tag=22)
                else:
                        v[:,0]=comm.recv(source=rank-1,tag=21)
                        comm.send(v[:,1],dest=rank-1,tag=22)
                        v[:,-1]=comm.recv(source=rank+1,tag=23)
                        comm.send(v[:,-2],dest=rank+1,tag=24)

        comm.Barrier()

	tstep+=1
	
	#rms_u=sum(sum((u-u_old)**2))/np.size(u)

	norm_u=np.linalg.norm(u-u_old)
	norm_v=np.linalg.norm(v-v_old)

	comm.Barrier()
	u_norm_u=comm.gather(norm_u,root = 0)
	v_norm_v=comm.gather(norm_v,root = 0)

	if rank==0:
		if tstep%5==0:

			print "---------------------------PARAMETER DISPLAY-----------------------"
			print "Simulation Time: ",tstep*dt," s"      	
			print "U velocity Residual: ", max(u_norm_u)
			print "V velocity Residual: ", max(v_norm_v)
			print "Pressure Residual: ",max(v_rms_error)
			print "Poisson Counter: ", p_counter

		if(max(v_norm_v)<10**-8 and max(v_norm_v)!=0 and max(u_norm_u)<10**-8 and max(u_norm_u)!=0):
			v_exit_val=1
		else: v_exit_val=0


	else:
		v_exit_val= None

	comm.Barrier()

	v_exit_val=comm.bcast(v_exit_val,root=0)

	if(v_exit_val==1):
		break

#	if tstep%10==0:

#		uu=(u[:,:-1]+u[:,1:])*0.5
#		vv=(v[1:,1:]+v[:-1,1:])*0.5

#		uu=uu.T
#		vv=vv.T

#		uu=np.reshape(uu,np.size(uu))
#		vv=np.reshape(vv,np.size(vv))

#		DataOut=np.column_stack((X.T,Y.T,uu.T,vv.T))
		#filename='/home/pi/Parallel_NS/IOData/data'+`ind`+'_0'+`rank`+'.dat'
#		np.savetxt('/home/pi/Parallel_NS/IOData/data0%d.dat' % rank,DataOut) 
		
#		ind=ind+1

	#print u_old-u

uu=(u[:,:-2]+u[:,1:-1])*0.5
vv=(v[1:,1:-1]+v[:-1,1:-1])*0.5

uu=uu.T
vv=vv.T

X=np.reshape(X,np.size(X))
Y=np.reshape(Y,np.size(Y))
uu=np.reshape(uu,np.size(uu))
vv=np.reshape(vv,np.size(vv))

DataOut=np.column_stack((X.T,Y.T,uu.T,vv.T))
np.savetxt('LidData20%d.dat' % rank,DataOut)
#np.savetxt('Re100U%d.dat' % rank,u)
#np.savetxt('Re100V%d.dat' % rank,v)
#np.savetxt('Re100P%d.dat' % rank,p)

if rank==0:
	t2=time.time()
	print t2-t1
comm.Barrier()

