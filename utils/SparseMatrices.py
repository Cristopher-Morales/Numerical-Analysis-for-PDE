#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 11:48:49 2024

@author: cristopher Morales Ubal
"""
import numpy as np
import scipy.sparse as sp
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

def Laplacematrix(a,b,h):
    N=int(((b-a)/h)+1)
    diagonal = 2.0*np.ones(N-2) 
    upanddown = -1.0*np.ones(N-3)
    A=np.diag(diagonal,0)+np.diag(upanddown,-1)+np.diag(upanddown,1)
    Laplacematrix=(1/h**2)*A
    return Laplacematrix
