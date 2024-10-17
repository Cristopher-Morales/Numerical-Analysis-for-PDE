#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 18:03:49 2024

@author: Cristopher Morales Ubal
"""
import numpy as np

def meshgrid(a,b,c,d,hx,hy):
    Nx=int((b-a)/hx)
    Ny=int((d-c)/hy)
    A=np.mgrid[0:Nx+1,0:Ny+1]   
    return hx*A

def gridpoints(a,b,h):
    N=int(((b-a)/h)+1)
    gridpoints=np.linspace(a,b,N)
    return [gridpoints,N]
