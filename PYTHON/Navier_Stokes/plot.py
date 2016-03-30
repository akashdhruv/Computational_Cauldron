import matplotlib.pyplot as plt
import numpy as np

data1=np.loadtxt('LidData200.dat')
data2=np.loadtxt('LidData201.dat')
data3=np.loadtxt('LidData202.dat')
data4=np.loadtxt('LidData203.dat')
data5=np.loadtxt('LidData204.dat')
data6=np.loadtxt('LidData205.dat')
data7=np.loadtxt('LidData206.dat')
data8=np.loadtxt('LidData207.dat')

M=32+1
N=4

x1=np.zeros((N,M),dtype=float)
y1=np.zeros((N,M),dtype=float)
u1=np.zeros((N,M),dtype=float)
v1=np.zeros((N,M),dtype=float)

x2=np.zeros((N,M),dtype=float)
y2=np.zeros((N,M),dtype=float)
u2=np.zeros((N,M),dtype=float)
v2=np.zeros((N,M),dtype=float)

x3=np.zeros((N,M),dtype=float)
y3=np.zeros((N,M),dtype=float)
u3=np.zeros((N,M),dtype=float)
v3=np.zeros((N,M),dtype=float)

x4=np.zeros((N,M),dtype=float)
y4=np.zeros((N,M),dtype=float)
u4=np.zeros((N,M),dtype=float)
v4=np.zeros((N,M),dtype=float)

X=np.zeros((M,M),dtype=float)
Y=np.zeros((M,M),dtype=float)
U=np.zeros((M,M),dtype=float)
V=np.zeros((M,M),dtype=float)

x5=np.zeros((N,M),dtype=float)
y5=np.zeros((N,M),dtype=float)
u5=np.zeros((N,M),dtype=float)
v5=np.zeros((N,M),dtype=float)

x6=np.zeros((N,M),dtype=float)
y6=np.zeros((N,M),dtype=float)
u6=np.zeros((N,M),dtype=float)
v6=np.zeros((N,M),dtype=float)

x7=np.zeros((N,M),dtype=float)
y7=np.zeros((N,M),dtype=float)
u7=np.zeros((N,M),dtype=float)
v7=np.zeros((N,M),dtype=float)

x8=np.zeros((N,M),dtype=float)
y8=np.zeros((N,M),dtype=float)
u8=np.zeros((N,M),dtype=float)
v8=np.zeros((N,M),dtype=float)

k=0

for i in range(N):
	for j in range(M):
	
		x1[i,j]=data1[k,0]
		y1[i,j]=data1[k,1]
		u1[i,j]=data1[k,2]
		v1[i,j]=data1[k,3]

		x2[i,j]=data2[k,0]
                y2[i,j]=data2[k,1]
                u2[i,j]=data2[k,2]
                v2[i,j]=data2[k,3]

		x3[i,j]=data3[k,0]
                y3[i,j]=data3[k,1]
                u3[i,j]=data3[k,2]
                v3[i,j]=data3[k,3]

		x4[i,j]=data4[k,0]
                y4[i,j]=data4[k,1]
                u4[i,j]=data4[k,2]
                v4[i,j]=data4[k,3] 

		x5[i,j]=data5[k,0]
                y5[i,j]=data5[k,1]
                u5[i,j]=data5[k,2]
                v5[i,j]=data5[k,3]

                x6[i,j]=data6[k,0]
                y6[i,j]=data6[k,1]
                u6[i,j]=data6[k,2]
                v6[i,j]=data6[k,3]

                x7[i,j]=data7[k,0]
                y7[i,j]=data7[k,1]
                u7[i,j]=data7[k,2]
                v7[i,j]=data7[k,3]

                x8[i,j]=data8[k,0]
                y8[i,j]=data8[k,1]
                u8[i,j]=data8[k,2]
                v8[i,j]=data8[k,3]
		
		k=k+1



X=np.concatenate((x1,x2,x3,x4,x5,x6,x7,x8))
Y=np.concatenate((y1,y2,y3,y4,y5,y6,y7,y8))
U=np.concatenate((u1,u2,u3,u4,u5,u6,u7,u8))
V=np.concatenate((v1,v2,v3,v4,v5,v6,v7,v8))

plt.figure()
plt.title('Re = 1000, Grid = 100 x 100')
#plt.contour(X,Y,np.sqrt(U**2+V**2),density=5)
plt.streamplot(X,Y,U,V,density=5)
plt.plot(X[:,0],Y[:,0],'k')
plt.plot(X[:,-1],Y[:,-1],'k')
plt.plot(X[0,:],Y[0,:],'k')
plt.plot(X[-1,:],Y[-1,:],'k')
#plt.ylim(-0.2,1.2)
#plt.xlim(-0.2,2.2)
plt.axis('equal')
plt.show()
