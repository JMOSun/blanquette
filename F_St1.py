# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 12:36:47 2015

@author: user
"""

import numpy as np
#import xlrd
#import CoolProp.CoolProp as CP
#import matplotlib.pylab as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Ei", help="Le stock initial en pourcentage du stock maximal")
parser.add_argument("Power", help="La puissance électrique exigée")
args=parser.parse_args()

Ei=float(args.Ei)
Power=float(args.Power)



global tf,V0,P0,T0,fE,N,fP,T_coeffs

#N2011=np.load('N_2011.npy')



# Fonction reliant la pression dans le réservoir et la quantité d'énergie stockée
fP=np.poly1d([ 0.31849889,  0.27168558,  0.4072201 ])
         

T_coeffs=np.load('T_coeffs.npy')

class S(object): #Classe représentant l'état du sysème définie par sa pression, son volume et l'énergie stockée
    """
    Classe représentant l'état du système par son volume, sa pression et l'énergie emmagasiné"""
    def __init__(self,E): #Etat définie à partir de la valeur sa pression
        if E>1:
            E=1
        if  E<0:
            E=0
        self.E=E
        self.P=fP(E)*250e5
        #rho=CP.PropsSI('D','T',T0,'P',self.P,'Air')
        #self.V=rho0*V0/rho
        
    def Latt(Liste,Attribut): #Liste contenant la valeur de l'attribut
        T=[]
        for k in Liste:
            T.append(k.__getattribute__(Attribut))
            return T

    def ListToTab(L):
        n=len(L)
        T=np.zeros(n)
        i=0
        for k in L:
            T[i]=k
            i+=1
        return T
    
    def Tab(Liste,Attribut): #Tableau contenant la valeur de l'attribut
        return ListToTab(Latt(Liste,Attribut))
    
    def eta(self,Pu):
        P=self.P/250e5
        Pw=abs(Pu)/100e3
        rmt=0
        k=np.zeros(5)
        som=[]
        tmp=0
        for i in np.arange(5):
            tmp=0
            for j in np.arange(2,7):
                tmp=tmp+T_coeffs[j,i]*Pw**(j-2)
            som.append(tmp)
        k[0]=1/(T_coeffs[0,0]/(Pw**(T_coeffs[1,0]))+som[0])
        k[1]=(T_coeffs[0,1]/(Pw**(T_coeffs[1,1]))+som[1])
        k[2]=-(T_coeffs[0,2]/(Pw**(T_coeffs[1,2]))+som[2])   
        k[3]=(T_coeffs[0,3]/(Pw**(T_coeffs[1,3]))+som[3])        
        k[4]=-(T_coeffs[0,4]/(Pw**(T_coeffs[1,4]))+som[4]) 
        i=0
        for i in range(5):
            rmt=rmt+k[i]*P**i
        if Pw <0.1:
            rmt=0.00000000001
        return rmt
    
    def S10(self,Pu):
        if Pu<0:
            c=-1
        else :
            c=1
        Ec=Pu*self.eta(Pu)**c*10
        Eadim=Ec/1744009560.0
        Sn=S(self.E+Eadim)
        return Sn
        
    
   
        
        
def Latt(Liste,Attribut):
    T=[]
    for k in Liste:
        T.append(k.__getattribute__(Attribut))
    return T

def ListToTab(L):
    n=len(L)
    T=np.zeros(n)
    i=0
    for k in L:
        T[i]=k
        i+=1
    return T

def Tab(Liste,Attribut):
    return ListToTab(Latt(Liste,Attribut))

    
S0=S(Ei)

#print(S0.P)
#print(S0.eta(40e3))

#print(Power)

S1=S0.S10(Power)
print(S1.E)


