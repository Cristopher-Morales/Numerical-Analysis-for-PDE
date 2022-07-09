# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 14:39:59 2020

@author: Cristopher Morales Ubal
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la

a=0
b=4
c=0
d=4
dt=0.015
hx=0.08
hy=0.08
alpha=5

# meshgrid

Nx=int((b-a)/hx)
Ny=int((d-c)/hy)
x,y=hx*np.mgrid[0:Nx+1,0:Ny+1]

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

L=Laplace2dmatrix(a,b,c,d,hx,hy)

# PART C)

# forward euler solution

#initial condition
u0 = np.exp(-alpha*(x-2)**2-alpha*(y-2)**2)

# reshape
u0=np.reshape(u0[1:Nx,1:Ny],newshape=(Nx-1)*(Ny-1),order='F')

t0=0
tf=0.15
Nt=int(tf/dt)
t1=0
j=0
i=0

fig1=plt.figure(figsize=(10,10))

for t in range(1,Nt+2):
    u=u0+(dt*-L*u0)
    usol=np.reshape(u,newshape=(Ny-1,Nx-1))
    if (j==int(0/dt) or j==int(0.045/dt) or j==int(0.09/dt) or j==int(0.15/dt)):
        ax=fig1.add_subplot(2, 4, i+1)
        nb=ax.imshow(usol,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
        plt.title('FDM-FE t= %4.3f s' %t1)
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.colorbar(nb,orientation="horizontal")
        i=i+1
    t1=t1+dt
    u0=u
    j=j+1  

# Backward euler solution

#initial condition

u0 = np.exp(-alpha*(x-2)**2-alpha*(y-2)**2)

# reshape
u0=np.reshape(u0[1:Nx,1:Ny],newshape=(Nx-1)*(Ny-1),order='F')

t0=0
tf=0.15
Nt=int(tf/dt)
t1=0  
j=0

for t in range(1,Nt+2):
    A=sp.eye((Nx-1)*(Ny-1))+L*dt
    u=la.spsolve(A,u0)
    usol=np.reshape(u,newshape=(Ny-1,Nx-1))
    if (j==int(0/dt) or j==int(0.045/dt) or j==int(0.09/dt) or j==int(0.15/dt)):
        ax=fig1.add_subplot(2, 4, i+1)
        nb=ax.imshow(usol,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
        plt.title('FDM-BE t= %4.3f s' %t1)
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.colorbar(nb,orientation="horizontal")
        i=i+1
    t1=t1+dt
    u0=u
    j=j+1
plt.show()

# PART D

dts=[0.005,0.003,0.0015]

for dt in dts:
    # forward euler solution
    #initial condition
    u0 = np.exp(-alpha*(x-2)**2-alpha*(y-2)**2)
    # reshape
    u0=np.reshape(u0[1:Nx,1:Ny],newshape=(Nx-1)*(Ny-1),order='F')
    t0=0
    tf=0.15
    Nt=int(tf/dt)
    t1=0
    j=0
    i=0
    fig1=plt.figure(figsize=(10,10))
    for t in range(1,Nt+2):
        u=u0+(dt*-L*u0)
        usol=np.reshape(u,newshape=(Ny-1,Nx-1))
        if (j==int(0/dt) or j==int(0.045/dt) or j==int(0.09/dt) or j==int(0.15/dt)):
            ax=fig1.add_subplot(2, 4, i+1)
            nb=ax.imshow(usol,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
            plt.title('FDM-FE t= %4.3f s' %t1)
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.colorbar(nb,orientation="horizontal")
            i=i+1
        t1=t1+dt
        u0=u
        j=j+1  
    # Backward euler solution
    #initial condition
    u0 = np.exp(-alpha*(x-2)**2-alpha*(y-2)**2)
    # reshape
    u0=np.reshape(u0[1:Nx,1:Ny],newshape=(Nx-1)*(Ny-1),order='F')
    dt1=0.015
    t0=0
    tf=0.15
    Nt=int(tf/dt1)
    t1=0  
    j=0
    for t in range(1,Nt+2):
        A=sp.eye((Nx-1)*(Ny-1))+L*dt1
        u=la.spsolve(A,u0)
        usol=np.reshape(u,newshape=(Ny-1,Nx-1))
        if (j==int(0/dt1) or j==int(0.045/dt1) or j==int(0.09/dt1) or j==int(0.15/dt1)):
            ax=fig1.add_subplot(2, 4, i+1)
            nb=ax.imshow(usol,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
            plt.title('FDM-BE t= %4.3f s' %t1)
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.colorbar(nb,orientation="horizontal")
            i=i+1
        t1=t1+dt1
        u0=u
        j=j+1
    plt.show()

# PART F)
    


# meshgrid

hxs=[0.08,0.04,0.02,0.01]
a1=0
import time    
for hx in hxs:
    dt=(hx**2)/4
    hy=hx
    Nx=int((b-a)/hx)
    Ny=int((d-c)/hy)
    x,y=np.mgrid[0:(b-a):hx,0:(d-c):hy]
    L=Laplace2dmatrix(a,b,c,d,hx,hy)
    # forward euler solution
    #initial condition
    u0 = np.exp(-alpha*(x-2)**2-alpha*(y-2)**2)
    # reshape
    u0=np.reshape(u0[1:Nx,1:Ny],newshape=(Nx-1)*(Ny-1),order='F')
    t0=0
    tf=0.15
    Nt=int(tf/dt)
    t1=0
    j=0
    i=0
    fig1=plt.figure(figsize=(10,10))
    a1 = time.time()
    for t in range(1,Nt+2):
        u=u0+(dt*-L*u0)
        usol=np.reshape(u,newshape=(Ny-1,Nx-1))
        if (j==int(0/dt) or j==int(0.045/dt) or j==int(0.09/dt) or j==int(0.15/dt)):
            ax=fig1.add_subplot(2, 4, i+1)
            nb=ax.imshow(usol,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
            plt.title('FDM-FE t= %4.3f s' %t1)
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.colorbar(nb,orientation="horizontal")
            i=i+1
        t1=t1+dt
        u0=u
        j=j+1  
    a2 = time.time()
    print(a2-a1)
    # Backward euler solution
    #initial condition
    u0 = np.exp(-alpha*(x-2)**2-alpha*(y-2)**2)
    # reshape
    u0=np.reshape(u0[1:Nx,1:Ny],newshape=(Nx-1)*(Ny-1),order='F')
    dt1=0.015
    t0=0
    tf=0.15
    Nt=int(tf/dt1)
    t1=0  
    j=0
    a3 = time.time()
    for t in range(1,Nt+2):
        A=sp.eye((Nx-1)*(Ny-1))+L*dt1
        u=la.spsolve(A,u0)
        usol=np.reshape(u,newshape=(Ny-1,Nx-1))
        if (j==int(0/dt1) or j==int(0.045/dt1) or j==int(0.09/dt1) or j==int(0.15/dt1)):
            ax=fig1.add_subplot(2, 4, i+1)
            nb=ax.imshow(usol,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
            plt.title('FDM-BE t= %4.3f s' %t1)
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.colorbar(nb,orientation="horizontal")
            i=i+1
        t1=t1+dt1
        u0=u
        j=j+1
    plt.show()
    a4 = time.time()
    print(a4-a3)




    

    
    




    








