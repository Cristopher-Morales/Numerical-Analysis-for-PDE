# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 09:19:52 2019

@author: Cristopher Morales Ubal
"""

# problema 3
import numpy as np
import math
import matplotlib.pyplot as plt 
import scipy.optimize as optimize

# part b)

h=0.2
a=0
b=1
u0=1
uN=2

def gridpoints(a,b,h):
    N=int(((b-a)/h)+1)
    gridpoints=np.linspace(a,b,N)
    return [gridpoints,N]

def Laplacematrix(a,b,h):
    N=int(((b-a)/h)+1)
    diagonal = 2.0*np.ones(N-2) 
    upanddown = -1.0*np.ones(N-3)
    A=np.diag(diagonal,0)+np.diag(upanddown,-1)+np.diag(upanddown,1)
    Laplacematrix=(1/h**2)*A
    return Laplacematrix

def GlobalError(k):
    hk=0.2/2**k
    gridpoint1= gridpoints(a,b,hk)[0]
    N=gridpoints(a,b,hk)[1]
    Lmatrix1=Laplacematrix(a,b,hk)
    f2=np.exp(gridpoint1[1:N-1])
    f2[0]=np.exp(gridpoint1[0])+u0/hk**2
    f2[N-3]=np.exp(gridpoint1[N-3])+uN/hk**2
    u2=np.zeros(N)
    u2[1:N-1]=np.linalg.solve(Lmatrix1,f2)
    u2[0]=u0
    u2[N-1]=uN
    u2exact=-np.exp(gridpoint1)+math.e*gridpoint1+2
    error=np.linalg.norm(u2exact-u2)
    return error

def func(h,c,a):
    return c*h**a

def GlobalErrorNonUniGrid(k):
    
    hk=0.2/2**k
    N=int((b-a)/hk)
    gridpoint2=np.zeros(N)
    gridpoint2[1]=2.0*hk
    gridpoint2[2:]=np.linspace((2.0+1.0)*hk,b,N-2)
    Lmatrix2=np.zeros((N-2,N-2))
    Lmatrix2[0,0]=3/hk**2
    Lmatrix2[0,1]=-2/hk**2
    Lmatrix2[1:,0:]=Laplacematrix(a,b,hk)[1:N-2,0:N-2]
    f2=np.exp(gridpoint2[1:N-1])
    f2[0]=4*np.exp(gridpoint2[1])+u0/hk**2
    f2[N-3]=np.exp(gridpoint2[N-2])+uN/hk**2
    u2=np.zeros(N)
    u2[1:N-1]=np.linalg.solve(Lmatrix2,f2)
    u2[0]=u0
    u2[N-1]=uN
    u2exact=-np.exp(gridpoint2)+math.e*gridpoint2+2
    error=np.linalg.norm(u2exact-u2)
    return error

gridpoint1= gridpoints(a,b,h)[0]
print(gridpoint1)
N=gridpoints(a,b,h)[1]
print(N)
Lmatrix1=Laplacematrix(a,b,h)
print(Lmatrix1)

# function f1=1

f1=np.ones(N-2)
f1[0]=1+u0/h**2
f1[N-3]=1+uN/h**2
print(f1)

# function f2=exp(x)

print(gridpoint1[1:N-1])
f2=np.exp(gridpoint1[1:N-1])
f2[0]=np.exp(gridpoint1[0])+u0/h**2
f2[N-3]=np.exp(gridpoint1[N-3])+uN/h**2
print(f2)

# part c)

u1=np.zeros(N)
u1[1:N-1]=np.linalg.solve(Lmatrix1,f1)
u1[0]=u0
u1[N-1]=uN
print(u1)
u1exact=-0.5*gridpoint1**2+3*gridpoint1/2+1
print(u1exact)
u2=np.zeros(N)
u2[1:N-1]=np.linalg.solve(Lmatrix1,f2)
u2[0]=u0
u2[N-1]=uN
print(u2)
u2exact=-np.exp(gridpoint1)+math.e*gridpoint1+2
print(u2exact)

# plot solutions 
plt.figure(1)
plt.subplot(121)
plt.plot(gridpoint1,u1,'ro-',label='Numerical Solution')
plt.plot(gridpoint1,u1exact,'bx-',label='Exact Solution')
plt.title('$f_{1}(x)=1$')
plt.xlabel('$x$')
plt.ylabel('$u_{1}(x)$')
plt.legend()
plt.show( )

plt.subplot(122)
plt.plot(gridpoint1,u2,'ro-',label='Numerical Solution')
plt.plot(gridpoint1,u2exact,'bx-',label='Exact Solution')
plt.title('$f_{2}(x)=exp(x)$')
plt.xlabel('$x$')
plt.ylabel('$u_{2}(x)$')
plt.legend()
plt.show( )

# Compute error

error1=np.linalg.norm(u1exact-u1)
print(error1)
error2=np.linalg.norm(u2exact-u2)
print(error2)

# part d

k=1
hk=0.2/2**k
gridpoint1= gridpoints(a,b,hk)[0]
N=gridpoints(a,b,hk)[1]
Lmatrix1=Laplacematrix(a,b,hk)

print(gridpoint1[1:N-1])
f2=np.exp(gridpoint1[1:N-1])
f2[0]=np.exp(gridpoint1[0])+u0/hk**2
f2[N-3]=np.exp(gridpoint1[N-3])+uN/hk**2
print(f2)
u2=np.zeros(N)
u2[1:N-1]=np.linalg.solve(Lmatrix1,f2)
u2[0]=u0
u2[N-1]=uN
print(u2)
u2exact=-np.exp(gridpoint1)+math.e*gridpoint1+2
print(u2exact)

plt.subplot(122)
plt.plot(gridpoint1,u2,'ro-',label='Numerical Solution')
plt.plot(gridpoint1,u2exact,'bx-',label='Exact Solution')
plt.title('$h=0.1$')
plt.xlabel('$x$')
plt.ylabel('$u_{2}(x)$')
plt.legend()
plt.show( )

error2_1=np.linalg.norm(u2exact-u2)
print(error2_1)

k=np.linspace(0,4,5)
error2=np.zeros(5)
print(k)

for i in range(5):
  error2[i]=GlobalError(i)

print(error2)

h=0.2/2**np.linspace(0,4,5)
print(h)
plt.plot(h,error2,'bo-')
plt.title('Behaviour of Global Error')
plt.xlabel('$h$')
plt.ylabel('Error')
parameters,cov=optimize.curve_fit(func,h,error2)
print(parameters)
plt.show()

# part e) non uniform grid

h1=0.2
N=int((b-a)/h1)
gridpoint2=np.zeros(N)
gridpoint2[1]=2.0*h1
gridpoint2[2:]=np.linspace((2.0+1.0)*h1,b,N-2)
print(gridpoint2)
Lmatrix2=np.zeros((N-2,N-2))
Lmatrix2[0,0]=3/h1**2
Lmatrix2[0,1]=-2/h1**2
Lmatrix2[1:,0:]=Laplacematrix(a,b,h1)[1:N-2,0:N-2]
print(Lmatrix2)

# function f1=1

f1=np.ones(N-2)
f1[0]=4+u0/h1**2
f1[N-3]=1+uN/h1**2
print(f1)

# function f2=exp(x)

print(gridpoint2[1:N-1])
f2=np.exp(gridpoint2[1:N-1])
f2[0]=4*np.exp(gridpoint2[1])+u0/h1**2
f2[N-3]=np.exp(gridpoint2[N-2])+uN/h1**2
print(f2)

# solve the lineal algebraic problem

u1=np.zeros(N)
u1[1:N-1]=np.linalg.solve(Lmatrix2,f1)
u1[0]=u0
u1[N-1]=uN
print(u1)
u1exact=-0.5*gridpoint2**2+3*gridpoint2/2+1
print(u1exact)
u2=np.zeros(N)
u2[1:N-1]=np.linalg.solve(Lmatrix2,f2)
u2[0]=u0
u2[N-1]=uN
print(u2)
u2exact=-np.exp(gridpoint2)+math.e*gridpoint2+2
print(u2exact)

# plot solutions 

plt.figure(1)
plt.subplot(121)
plt.plot(gridpoint2,u1,'ro-',label='Numerical Solution')
plt.plot(gridpoint2,u1exact,'bx-',label='Exact Solution')
plt.title('Non Uniform Grid $f_{1}(x)=1$')
plt.xlabel('$x$')
plt.ylabel('$u_{1}(x)$')
plt.legend()
plt.show( )

plt.subplot(122)
plt.plot(gridpoint2,u2,'ro-',label='Numerical Solution')
plt.plot(gridpoint2,u2exact,'bx-',label='Exact Solution')
plt.title('Non Uniform Grid $f_{2}(x)=exp(x)$')
plt.xlabel('$x$')
plt.ylabel('$u_{2}(x)$')
plt.legend()
plt.show( )

# Compute error

error1=np.linalg.norm(u1exact-u1)
print(error1)
error2=np.linalg.norm(u2exact-u2)
print(error2)

# behaviour error with respect to step


k=np.linspace(0,4,5)
error2=np.zeros(5)
print(k)

for i in range(5):
  error2[i]=GlobalErrorNonUniGrid(i)

print(error2)

h=0.2/2**np.linspace(0,4,5)
plt.plot(h,error2,'bo-')
plt.title('Behaviour of Global Error Non Uniform Grid')
plt.xlabel('$h$')
plt.ylabel('Error')
parameters,cov=optimize.curve_fit(func,h,error2)
print(parameters)
plt.show()
