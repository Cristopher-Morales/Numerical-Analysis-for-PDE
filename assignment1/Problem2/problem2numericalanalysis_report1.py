# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 06:09:52 2019

@author: Cristopher Morales Ubal
"""

import numpy as np
import math
import matplotlib.pyplot as plt 

# part a) 

L = 1.0
h = 0.1 
N = int((L/h)+1)
u0 = 0.0 
uN = 0.0 

gridpoints = np.linspace(0.0,1.0,N) 
print(gridpoints)

# part d) 

diagonal = 2.0*np.ones(N-2) 
upanddown = -1.0*np.ones(N-3)
A=np.diag(diagonal,0)+np.diag(upanddown,-1)+np.diag(upanddown,1)
Laplacematrix=(1/h**2)*A
print(A)
print(Laplacematrix)
print(Laplacematrix[0,:])
print(Laplacematrix[1,:])
print(Laplacematrix[8,:])
plt.spy(Laplacematrix,precision=0,marker='o',markerfacecolor='g')
plt.show( )

# part e)

i=np.linspace(1,9,9)
eigvalueoriginal=(np.pi*i/L)**2
eigvalueLmatrix=(4/h**2)*(np.sin(np.pi*i*h/(2*L)))**2
[eigenvalues,eigenvectors]=np.linalg.eig(Laplacematrix)

print(eigvalueLmatrix)
print(eigvalueoriginal)
print(eigenvalues)
print(eigenvectors)

# table 
 
# plot of eigenvalues

ima = np.zeros(9)
plt.plot(eigenvalues,ima,'bx',label='numerical')
plt.plot(eigvalueoriginal,ima,'ro',label='Original Problem')
plt.title('Eigenvalues')
plt.xlabel(r'Re($λ_{i}$)')
plt.ylabel(r'Im($λ_{i}$)')
plt.legend()
plt.show( )

# part f)
N=10
j=np.linspace(0,N,N+1)
sampleeigvector=np.zeros((N-2,N+1))
k=0
while k<N-2:
    sampleeigvector[k,:]=np.sin(np.pi*(k+1)*j/N)
    plt.plot(j,sampleeigvector[k,:])
    plt.show()
    k=k+1
    
    


    
    
    