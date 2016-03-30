import matplotlib.pyplot as plt
import numpy as np

k=3
d=2

M=20+1
N=20+1

X=np.zeros((N*k,M*k),dtype=float)
Y=np.zeros((N*k,M*k),dtype=float)
U=np.zeros((N*k,M*k),dtype=float)
V=np.zeros((N*k,M*k),dtype=float)

for i in range(0,k**d):

	x=np.loadtxt('X%d.dat' % i)
	y=np.loadtxt('Y%d.dat' % i)
	u=np.loadtxt('U%d.dat' % i)
	v=np.loadtxt('V%d.dat' % i)
        
	x=np.reshape(x,[N,M])
	y=np.reshape(y,[N,M])
	u=np.reshape(u,[N,M])
	v=np.reshape(v,[N,M])
		
	X[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=x
	Y[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=y
	U[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=u
	V[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=v
 

plt.figure()
#plt.title('Re = 1000, Grid = 128 x 128, Solver = FORTRAN')
#plt.contourf(X,Y,V,density=5)
plt.streamplot(X,Y,U,V,density=5)
plt.plot(X[:,0],Y[:,0],'k')
plt.plot(X[:,-1],Y[:,-1],'k')
plt.plot(X[0,:],Y[0,:],'k')
plt.plot(X[-1,:],Y[-1,:],'k')
#plt.ylim(-0.2,1.2)
#plt.xlim(-0.2,2.2)
plt.xlabel('X')
plt.ylabel('Y')

plt.axis('equal')
plt.show()
