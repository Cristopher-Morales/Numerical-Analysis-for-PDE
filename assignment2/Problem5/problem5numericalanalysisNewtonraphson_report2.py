# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 09:59:51 2020

@author: Cristopher Morales Ubal
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la
from utils import SparseMatrices as spm

a=0
b=16
c=0
d=8
alpha=1
ti=0
tf=40
tol=1/1000

# definition function k

def k(x,y,alpha):
    if((x>=1 and x<=2 and y>=1 and y<=2) or (x>=1 and x<=3 and y>=3 and y<=5) or (x>=4 and x<=7 and y>=4 and y<=7) or (x>=9 and x<=12 and y>=4 and y<=6) or (x>=13 and x<=15 and y>=1 and y<=3)):
        return alpha
    else:
        return 0

# defining grid steps
hs=[0.1,0.08,0.04]

# defining time steps

#ts=[0,1,2,3,5,10,20,40]
t1=0
import time
error=[]
i=0
positions=[]
for h in hs:
    dt=0.15#0.99*h**2/4 #linear criterion
    Nt=int((tf-ti)/dt)
    hx=h
    hy=h

    # meshgrid

    Nx=int((b-a)/hx)
    Ny=int((d-c)/hy)
    x,y=h*np.mgrid[0:Nx+1,0:Ny+1]
    
    L=spm.Laplace2dmatrix(a,b,c,d,hx,hy)  
    
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
        ui1=u0
        while True:
            fui2=-L*ui1+(k1 * ui1 *(np.ones((Ny-1)*(Nx-1))-ui1))
            diag=sp.diags([ui1 * k1],[0],shape=[(Ny-1)*(Nx-1),(Ny-1)*(Nx-1)])
            diag2=sp.diags([k1],[0],shape=[(Ny-1)*(Nx-1),(Ny-1)*(Nx-1)])
            A=sp.eye((Nx-1)*(Ny-1))-(diag2-L-2*diag)*dt
            vi=la.spsolve(A,ui1-u0-dt*fui2)
            ui2=ui1-vi
            error.append(np.linalg.norm(ui2-ui1))
            i=i+1
            if (np.linalg.norm(ui2-ui1)<tol):
                positions.append(i)
                break
            else:
                ui1=ui2
            
        if (np.abs(t-1)<dt/2 or np.abs(t-2)<dt/2 or np.abs(t-3)<dt/2 or np.abs(t-5)<dt/2 or np.abs(t-10)<dt/2 or np.abs(t-20)<dt/2 or np.abs(t-40)<dt):
            u1sol=np.reshape(ui2,newshape=(Ny-1,Nx-1))
            plt.imshow(u1sol,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
            plt.title('Fisher equation FE solution t= %4.3f s' %t)
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.colorbar() #orientation="horizontal"
            plt.show()
        u0=ui2
    t2 = time.time()
    print(t2-t1)