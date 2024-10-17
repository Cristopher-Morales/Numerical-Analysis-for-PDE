# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 07:18:31 2019

@author: Cristopher Morales Ubal
"""

# problem 4

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as la
from utils.meshGrid import meshgrid as mg
import utils.SparseMatrices as spm

a=0
b=2
c=0
d=1
hx=0.02
hy=0.02
Nx=int((b-a)/hx)
Ny=int((d-c)/hy)

Dx=spm.sparsematrix(a,b,hx)
Dy=spm.sparsematrix(c,d,hy)
plt.spy(Dx,marker='o')
plt.show()
L=spm.Laplace2dmatrix(a,b,c,d,hx,hy)
plt.spy(L,marker='o')
plt.show()
x,y=mg(a,b,c,d,hx,hy)
print(x,y)
f = 20*np.sin(1.5*x*np.pi+np.pi)*np.sin(np.pi*y)
plt.imshow(f,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
plt.title('Source Function $f(x,y)$')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.colorbar()
plt.show()

gridpointsx=np.linspace(a,b,Nx+1)
gridpointsy=np.linspace(c,d,Ny+1)

# boundary conditions

uleft=np.sin(2*np.pi*gridpointsy)
uright=np.sin(2*np.pi*gridpointsy)
ubottom=np.sin(0.5*np.pi*gridpointsx)
utop=np.zeros(Nx+1)

# function f modified at grid point (1,1),(1,Ny-1)(Nx-1,1) and (Nx-1,Ny-1) adding boundary conditions
f[1,1]=20*np.sin(1.5*hx*np.pi+np.pi)*np.sin(np.pi*hy)+(ubottom[1]+uleft[1])/hx**2
f[Nx-1,1]=20*np.sin(1.5*(Nx-1)*hx*np.pi+np.pi)*np.sin(np.pi*hy)+(ubottom[Nx-1]+uright[1])/hx**2
f[Nx-1,Ny-1]=20*np.sin(1.5*(Nx-1)*hx*np.pi+np.pi)*np.sin(np.pi*(Ny-1)*hy)+(utop[Nx-1]+uright[Ny-1])/hx**2
f[1,Ny-1]=20*np.sin(1.5*hx*np.pi+np.pi)*np.sin(np.pi*(Ny-1)*hy)+(utop[1]+uleft[Ny-1])/hx**2

# function f modified at grid values of the form (i,1) ,(i,Ny-1), (1,j) and (Nx-1,j) adding boundary conditions

f[2:Nx-2,1]=20*np.sin(1.5*gridpointsx[2:Nx-2]*hx*np.pi+np.pi)*np.sin(np.pi*hy)+ubottom[2:Nx-2]/hx**2
f[2:Nx-2,Ny-1]=20*np.sin(1.5*gridpointsx[2:Nx-2]*hx*np.pi+np.pi)*np.sin(np.pi*(Ny-1)*hy)+utop[2:Nx-2]/hx**2
f[1,2:Ny-2]=20*np.sin(1.5*hx*np.pi+np.pi)*np.sin(np.pi*gridpointsy[2:Ny-2]*hy)+uleft[2:Ny-2]/hx**2
f[Nx-1,2:Ny-2]=20*np.sin(1.5*(Nx-1)*hx*np.pi+np.pi)*np.sin(np.pi*gridpointsy[2:Ny-2]*hy)+uright[2:Ny-2]/hx**2

print(np.size(f))
flexiord=np.reshape(f[1:Nx,1:Ny],newshape=(Nx-1)*(Ny-1),order='F')

# solving the equation

u=la.spsolve(L,flexiord)
usol=np.reshape(u,newshape=(Ny-1,Nx-1))
plt.imshow(usol,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
plt.title('Solution  2D-Poisson')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.colorbar()
plt.show()
