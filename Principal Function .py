#!/usr/bin/env python
# coding: utf-8

# In[4]:



import numpy as np
import random as random


# In[ ]:


# Script to calculate numerical simulation of the number of cells on the markovian and no markovian models 

# Authors: Natalia G. Lavalle 

#           nglavalle@gmail.com
          
#           Instituto de Física de Líquidos y Sistemas Biológicos (IFLySiB) — Universidad Nacional de La Plata and             CONICET, Calle 59 n. 789, B1900BTE La Plata, Argentina
          
#           Osvaldo Chara

#           School of Biosciences, University of Nottingham, Sutton Bonington Campus,Nottingham LE12 5RD, UK                   Instituto de Tecnología, Universidad Argentina de la Empresa, Buenos Aires C1073AAO, Argentina

#           Tomás S.Grigera
          
#           CCT CONICET La Plata, Consejo Nacional de Investigaciones Científicas y Técnicas, Argentina
#           Departamento de Física, Facultad de Ciencias Exactas, Universidad Nacional de La Plata, Argentina
#          Istituto dei Sistemi Complessi, Consiglio Nazionale delle Ricerche, via dei Taurini 19, 00185 Rome,Italy
         
# Last update: 23 May 2023


# In[ ]:


# nt is a function that returns an array of M elements, where each element is a realization of the number of cells over time
# M is the number of samples of tissues with the same characteristics
# T is the value of refractory period 
# tmax is the time 
# σ is the rate of proliferation
# κ is the rate of apoptosis 
# ni is the inicial number of cells 


def nt(M,T,tmax,σ,κ,ni):

    N=np.zeros((tmax,T+1),dtype=int)
    l=0
    Nm=[]
    H= σ + κ
    for k in range (0,M):
        N=np.zeros((tmax,T+1),dtype=int) 
        N[0,T]=ni
        Nt=[ni]
        for t in range(1,tmax):
            for i in range(0,T):
                N[t,i]=N[t-1,i+1]
            N[t,0]+=N[t-1,0]
            for n in range(0,N[t,0]):
                x=random.random()
                if x<=κ:
                    N[t,0]-=1
                elif x>κ and x<H:
                    N[t,T]+=2
                    N[t,0]-=1
                            
            for j in range(1,T):
                for k in range(0,N[t,j]):
                    y=random.random()
                    if y < κ:
                        N[t,j]-=1

            l=sum(N[t])
            Nt.append(l)
        Nm.append(Nt)
    Nm=np.asarray(Nm)
    return Nm

