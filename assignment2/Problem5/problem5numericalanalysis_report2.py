# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 12:02:28 2020

@author: Cristopher Morales Ubal
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp

a=0
b=16
c=0
d=8
alpha=1
ti=0
tf=40

# definition function k

def k(x,y,alpha):
    if((x>=1 and x<=2 and y>=1 and y<=2) or (x>=1 and x<=3 and y>=3 and y<=5) or (x>=4 and x<=7 and y>=4 and y<=7) or (x>=9 and x<=12 and y>=4 and y<=6) or (x>=13 and x<=15 and y>=1 and y<=3)):
        return alpha
    else:
        return 0
    
# sparse matrix
        
def sparsematrix(a,b,hx):
    N=int((b-a)/hx)
    updiag=np.ones(N-1)
    downdiag=-1*np.ones(N-1)
    dim=(N,N-1)
    D=1/hx*(sp.diags([updiag,downdiag],[0,-1],shape=dim))
    return D

# Laplace matrix
    
def Laplace2dmatrix(a,b,c,d,hx,hy):
    Dx=sparsematrix(a,b,hx)
    Dy=sparsematrix(c,d,hy)
    Lxx=Dx.transpose().dot(Dx)
    Lyy=Dy.transpose().dot(Dy)
    Ix=sp.eye(int(b-a)/hx-1)
    Iy=sp.eye(int(d-c)/hy-1)
    L=sp.kron(Iy,Lxx)+sp.kron(Lyy,Ix)
    return L

# defining grid steps
hs=[0.1,0.08,0.04]

# defining time steps

#ts=[0,1,2,3,5,10,20,40]
t1=0
import time

for h in hs:
    dt=0.99*h**2/4 #linear criterion
    Nt=int((tf-ti)/dt)
    hx=h
    hy=h

    # meshgrid

    Nx=int((b-a)/hx)
    Ny=int((d-c)/hy)
    x,y=h*np.mgrid[0:Nx+1,0:Ny+1]
    
    L=Laplace2dmatrix(a,b,c,d,hx,hy)  
    
    k1=[]
    for j in range(1,Ny):
        for i in range(1,Nx):
            k10=k(i*hx,j*hy,alpha)
            k1.append(k10)
    k1 = np.array(k1)
    u0=np.exp(-2*(x-1.5)**2-2*(y-1.5)**2)
    u0=np.array(np.reshape(u0[1:Nx,1:Ny],newshape=(Nx-1)*(Ny-1),order='F'))
    ts=np.linspace(ti,Nt*dt,Nt+1)
    fig=plt.figure()
    u1sol=np.reshape(u0,newshape=(Ny-1,Nx-1))
    plt.imshow(u1sol,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
    plt.title('Fisher equation FE solution t=0.000' )
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.colorbar() #orientation="horizontal"
    plt.show()
    t1=time.time()
    for t in ts[1:]:
        u1=u0+dt*(-L*u0+k1 * u0 *(np.ones((Nx-1)*(Ny-1))-u0))   
        if (np.abs(t-1)<dt/2 or np.abs(t-2)<dt/2 or np.abs(t-3)<dt/2 or np.abs(t-5)<dt/2 or np.abs(t-10)<dt/2 or np.abs(t-20)<dt/2 or np.abs(t-40)<dt):
            u1sol=np.reshape(u1,newshape=(Ny-1,Nx-1))
            plt.imshow(u1sol,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
            plt.title('Fisher equation FE solution t= %4.3f s' %t)
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.colorbar() #orientation="horizontal"
            plt.show()
        u0=u1
    t2 = time.time()
    print(t2-t1)



