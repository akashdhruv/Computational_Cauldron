import matplotlib.pyplot as plt
import numpy as np

x=np.loadtxt('X0.dat')
y=np.loadtxt('Y0.dat')
u=np.loadtxt('U0.dat')
v=np.loadtxt('V0.dat')

M=128+1
N=128+1

x=np.reshape(x,(N,M))
y=np.reshape(y,(N,M))
u=np.reshape(u,(N,M))
v=np.reshape(v,(N,M))

plt.figure()
plt.title('Re = 1000, Grid = 128 x 128, Solver = FORTRAN')
#plt.contourf(X,Y,U,density=5)
plt.streamplot(x,y,u,v,density=5)
plt.plot(x[:,0],y[:,0],'k')
plt.plot(x[:,-1],y[:,-1],'k')
plt.plot(x[0,:],y[0,:],'k')
plt.plot(x[-1,:],y[-1,:],'k')
#plt.ylim(-0.2,1.2)
#plt.xlim(-0.2,2.2)
plt.xlabel('X')
plt.ylabel('Y')

plt.axis('equal')
plt.show()
