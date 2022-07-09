# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 11:43:13 2020

@author: Cristopher Morales Ubal
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp

a=0
b=4
c=0
d=12
ti=0
tf=8
h=0.02
v=4
alpha=-50

cwave=2




hx=h
hy=h
dt=np.sqrt(1/2)*0.99*h/cwave
Nt=int((tf-ti)/dt)
print(Nt)
# meshgrid

Nx=int((b-a)/hx)
Ny=int((d-c)/hy)
x,y=h*np.mgrid[0:Nx+1,0:Ny+1]

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



u0=np.zeros((Ny-1)*(Nx-1))

u1=0.5*((cwave*dt)**2)*np.sin(2*np.pi*ti*v)*np.exp(alpha*(x-2)**2+alpha*(y-2)**2)
u1=np.reshape(u1[1:Nx,1:Ny],newshape=(Nx-1)*(Ny-1),order='F')

ts=np.linspace(ti,Nt*dt,Nt+1)

fig=plt.figure()
for t in ts[1:]:
    f=np.sin(2*np.pi*t*v)*np.exp(alpha*(x-2)**2+alpha*(y-2)**2)
    f1=np.reshape(f[1:Nx,1:Ny],newshape=(Nx-1)*(Ny-1),order='F')
    u2=2*u1-u0+(cwave*dt)**2*(-L*u1+f1)   
    if (np.abs(t-1)<dt/2 or np.abs(t-2)<dt/2 or np.abs(t-3)<dt/2 or np.abs(t-4)<dt/2 or np.abs(t-5)<dt/2 or np.abs(t-6)<dt/2 or np.abs(t-7)<dt/2 or np.abs(t-8)<dt):
        u2sol=np.reshape(u2,newshape=(Ny-1,Nx-1))
        plt.imshow(u2sol,extent = [a , b, c , d],origin='lower',interpolation='none')
        plt.title('wave equation solution t= %1.0f s' %t)
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.colorbar(orientation="vertical")
        plt.show()
    u0=u1
    u1=u2


