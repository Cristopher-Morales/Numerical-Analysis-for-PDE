# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:46:39 2019

@author: Cristopher Morales Ubal
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la

a=0
b=8
c=0
d=4
alpha=-5
hx=0.02
hy=0.02
Nx=int((b-a)/hx)
Ny=int((d-c)/hy)
gridx=np.linspace(a,b,Nx+1)
gridy=np.linspace(c,d,Ny+1)
x,y=hx*np.mgrid[0:Nx+1,0:Ny+1]

# function k(x,y)

def k(x,y):
    return 1+4*x+6*y

# Source Function f(x,y)

def f(x,y,alpha):
    return np.exp(alpha*(x-1)**2+alpha*(y-1)**2)+np.exp(alpha*(x-3)**2+alpha*(y-1)**2)+np.exp(alpha*(x-5)**2+alpha*(y-1)**2)+np.exp(alpha*(x-7)**2+alpha*(y-1)**2)+np.exp(alpha*(x-1)**2+alpha*(y-3)**2)+np.exp(alpha*(x-3)**2+alpha*(y-3)**2)+np.exp(alpha*(x-5)**2+alpha*(y-3)**2)+np.exp(alpha*(x-7)**2+alpha*(y-3)**2)
    

#print(gridpointsx,gridpointsy)

maindiagonal=[]
updiagonal=[]
downdiagonal=[]
downerdiagonal=[]
upperdiagonal=[]

coeffifunction=[]
#print(np.size(gridpointsx[1:Nx]))
#print(np.size(gridpointsx[1:Ny]))

#for j in gridpointsy[1:Ny]:
   # for i in gridpointsx[1:Nx]:
       # source=f(i,j,alpha)
        # sourcefunction.append(source)
        
#sourcefunction
        
sourcefunction1=f(x,y,alpha)
sourcefunction=np.reshape(sourcefunction1[1:Nx,1:Ny],newshape=(Nx-1)*(Ny-1),order='F')

for j in range(np.size(gridy[1:])+1):
    for i in range(np.size(gridx[1:])+1):
        coeffifunction1=k(i*hx,j*hy)
        coeffifunction.append(coeffifunction1)

for j in range(np.size(gridy[1:Ny])):
    for i in range(np.size(gridx[1:Nx])):
        maindiagonal1=(k((i+1+0.5)*hx,(j+1)*hy)+k((i+1-0.5)*hx,(j+1)*hy))/(hx**2)+(k((i+1)*hx,(j+1+0.5)*hy)+k((i+1)*hx,(j+1-0.5)*hy))/(hy**2)
        maindiagonal.append(maindiagonal1)
        
for j in range(np.size(gridy[1:Ny])):
    for i in range(np.size(gridx[1:Nx])):
        updiagonal1=-(k((i+1+0.5)*hx,(j+1)*hy))/(hx**2)
        updiagonal.append(updiagonal1)
        
for j in range(np.size(gridy[1:Ny])):
    for i in range(np.size(gridx[1:Nx])):
        downdiagonal1=-(k((i+1-0.5)*hx,(j+1)*hy))/(hx**2)
        downdiagonal.append(downdiagonal1)

for j in range(Ny-2):
    downdiagonal[(j+1)*(Nx-1)]=0


for j in range(np.size(gridy[1:Ny])):
    for i in range(np.size(gridx[1:Nx])):
        upperdiagonal1=-k((i+1)*hx,(j+1+0.5)*hy)/(hy**2)
        upperdiagonal.append(upperdiagonal1)

for j in range(np.size(gridy[1:Ny])):
    for i in range(np.size(gridx[1:Nx])):
        downerdiagonal1=-k((i+1)*hx,(j+1-0.5)*hy)/(hy**2)
        downerdiagonal.append(downerdiagonal1)


updiagonal=updiagonal[:(Nx-1)*(Ny-1)-1]
downdiagonal=downdiagonal[1:]
upperdiagonal=upperdiagonal[:(Nx-1)*(Ny-1)-Nx+1]
downerdiagonal=downerdiagonal[Nx-1:]

for j in range(Ny-2):
    updiagonal[(j+1)*(Nx-1)-1]=0

#print(diagonal)
#print(updiagonal)
#print(downdiagonal)
#print(upperdiagonal)
#print(downerdiagonal)

#print(np.size(sourcefunction))
L=sp.diags([maindiagonal, updiagonal,downdiagonal, upperdiagonal, downerdiagonal],[0,1,-1,Nx-1,-(Nx-1)],format='csc')#.toarray()
plt.spy(L,marker='o')
plt.show()

# solution
u=la.spsolve(L,sourcefunction)
print(min(u))

sourcefunction=np.reshape(sourcefunction,newshape=(Ny-1,Nx-1))
plt.imshow(sourcefunction,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
plt.title('Source Function $f(x,y)$')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.colorbar()
plt.show()

coeffifunction=np.reshape(coeffifunction,newshape=(Ny+1,Nx+1))
plt.imshow(coeffifunction,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
plt.title('Coefficient Function $k(x,y)$')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.colorbar()
plt.show()

usol=np.reshape(u,newshape=(Ny-1,Nx-1))
#usol=u.reshape(Nx-1,Ny-1)
plt.imshow(usol,extent = [np.min(x) , np.max(x), np.min(y) , np.max(y)],origin='lower',interpolation='none')
plt.title('Solution 2D-Poisson ')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.colorbar()
plt.show()


