# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 12:36:47 2015

@author: user
"""

import numpy as np
import xlrd
import CoolProp.CoolProp as CP
import matplotlib.pylab as plt
import scipy.optimize as optimization
import random as rdm
import numpy.random as alea
import lmfit
import pickle
import os.path
import time as temps

from scipy.integrate import simps

debut=temps.time()

global tf,V0,P0,T0,fE,N

#N2011=np.load('N_2011.npy')

compt_plot=0

# Fonction reliant la pression dans le réservoir et la quantité d'énergie stockée
fE=np.poly1d([   3.50954367,  -25.82030936,   84.77714828, -163.75296739,
        206.48972107, -178.40534573,  108.13744199,  -46.87890618,
         15.55363657,   -2.64106872])
         
#Définition de l'état où l'énergie est nulle         
V0=200
Vf=V0/2.5
V1=200
d0=0.1 #Diamètre
P0=100e5
T0=293.15
rho0=CP.PropsSI('D','T',T0,'P',P0,'Air')
cyl=242
rpm0=1500

tf=6*3600

path1=r'./eta_elec_100kW.xlsx'

book1 = xlrd.open_workbook(path1)
first_sheet1 = book1.sheet_by_index(0)


path2=r'./eta_elec_100kW.xlsx'

book2 = xlrd.open_workbook(path2)
first_sheet2 = book2.sheet_by_index(0)
    


def rdm_el(Pu): #Rendement mécanique -> électrique
    x=[k.value for k in first_sheet1.col(0)]
    y=[k.value for k in first_sheet1.col(1)]
    return np.interp(abs(Pu),x,y)


def var_el(Pu): #Rendement de l'inverter
    x=[k.value for k in first_sheet2.col(0)]
    y=[k.value for k in first_sheet2.col(1)]
    return np.interp(abs(Pu),x,y)    
    
    
def eta_hydrau(rpm,P): #Rendement volumètrique de la pompe
    S=1/(-2.987894025497691342e+00*abs(rpm)/P-6.333293557613091096e-01)+1-1.758296915625696477e-03
    S=np.clip(S,0.01,1)
    return S

def Cf(rpm): #Coefficient de frottement dans la pompe
    a=0.000000014
    b=0.000000
    c=0.23
    return a*rpm**2+b*abs(rpm)+c

def dQ(rho): #Chaleur à évacuer
    return 9.021557E+01*np.log(rho) - 2.368101E+01
    
def f_eta_stock(Pstar): #Rendement du stockage 
    z=[-0.85423949,  1.85148906]
    f = np.poly1d(z)
    return f(Pstar)

def f_Pstar(Psurf): #Influence de la puissance thermique à évacuer sur l'écart entre la pression et la pression isotherme
    z=[ -4.31007918e-08,   1.57353400e-04,   1.00956153e+00]
    f = np.poly1d(z)
    return f(Psurf)

class S(object): #Classe représentant l'état du sysème définie par sa pression, son volume et l'énergie stockée
    """
    Classe représentant l'état du système par son volume, sa pression et l'énergie emmagasiné"""
    def __init__(self,P): #Etat définie à partir de la valeur sa pression
        if P>250e5:
            P=250e5
        if  P<100e5:
            P=100e5
        self.P=P
        rho=CP.PropsSI('D','T',T0,'P',P,'Air')
        self.V=rho0*V0/rho
        self.E=fE(P/250e5)*(500*3600*1e3)
        
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
    
    def f_ch(self,rpm,tf): #Fonction calculant l'état S1 à partir de S0, de la vitesse de rotation et du temps d'éxécution
        
        class output(object):
            Pu_meca = None
            Pu_elec = None
            Pu_hydr = None
            eta_tot = None
            P       = None
            V       = None
            Pn      = None
            
    
        Nt=2
        #time=np.linspace(0,tf,Nt)
        Dt=tf
        V=np.zeros(Nt)
        Pn=np.zeros(Nt)
        rho=np.zeros(Nt)
        
        
        V[0]=self.V
        rho[0]=CP.PropsSI('D','T',T0,'P',self.P,'Air')
        Pn[0]=self.P
        
        
        Q=np.zeros(Nt-1)
        Psurf=np.zeros(Nt-1)
        Pstar=np.zeros(Nt-1)
        eta_stock=np.zeros(Nt-1)
        #P=np.zeros(Nt-1)
        Pu_meca=np.zeros(Nt-1)
        Pu_elec=np.zeros(Nt-1)
        Pu_elec_res=np.zeros(Nt-1)
        Pu_hydr=np.zeros(Nt-1)
        eta_tot=np.zeros(Nt-1)
        
        
            
        Q_ideal=cyl*rpm/1.0e6/60.
        V_ideal=V[0]-Q_ideal*Dt
        rho_ideal=rho[0]*V[0]/V_ideal
        P_ideal=CP.PropsSI('P','T',T0,'D',rho_ideal,'Air')
        
        err=1.    
        Pp=P_ideal   
    
        if rpm >0:
            sign=-1
        else :
            sign=1    
        
        while err>1e-3:
            Qq=cyl*rpm/1.0e6/60.*eta_hydrau(rpm,Pp/1.0e5)**(-sign)
            Vv=V[0]-Qq*Dt
            rho_r=rho[0]*V[0]/Vv
            Pp2=CP.PropsSI('P','T',T0,'D',rho_r,'Air')
            err=abs(Pp2-Pp)/Pp
            Pp=Pp2
        V[1]=Vv
        Pn[1]=CP.PropsSI('P','T',T0,'D',rho[0]*V[0]/Vv,'Air')
        rho[1]=rho_r
        #P[i]=Pp
        Q=Qq
        Psurf=(dQ(rho[1])-dQ(rho[0]))*1000*rho[0]*V[0]/(2*(V[0]+V[1])/d0)/(Dt)
        Pstar=f_Pstar(Psurf)
        eta_stock=f_eta_stock(Pstar)
        Couple=(cyl*(Pn[1]+Pn[0])/2/1.0e5/63+Cf(rpm)*cyl/4.9)*eta_stock**(sign*0.5)
        Pu_meca=Couple*rpm*2*np.pi/60
        Pu_elec=Pu_meca*rdm_el(Pu_meca)**(sign)
        Pu_elec_res=Pu_elec*var_el(Pu_elec)**(sign)
        Pu_hydr=(Pn[1]+Pn[0])/2*Q
        eta_tot=(Pu_hydr/Pu_elec_res)**(-sign)
        
        out=output()
        out.eta_tot=eta_tot
        out.Pu_hydr=Pu_hydr
        out.Pu_elec=Pu_elec_res      
        out.Pu_meca=Pu_meca
        
        if Pn[0]==250e5 and rpm >0:
            out.eta_tot=0
            out.Pu_hydr=0
            out.Pu_elec=0
            out.Pu_meca=0
            Pn[1]=Pn[0]
        
        if Pn[0]==100e5 and rpm <0:
            out.eta_tot=0
            out.Pu_hydr=0
            out.Pu_elec=0
            out.Pu_meca=0
            Pn[1]=Pn[0]
            
        S1=S(Pn[1])
   
        return S1,out
        
    def f_ch_T(self,rpm,tf): 
        Nt=len(rpm)
        time=np.linspace(0,tf,Nt+1)
        Dt=time[1]-time[0]
        TS=[]
        TS.append(self)
        Tpow=[]
        for k in rpm:
            Si,Pu = TS[-1].f_ch(k,Dt)
            TS.append(Si)
            Tpow.append(Pu)
        return TS,Tpow
    
    

    def f_ch2(self,Pelec): #Effectue la même chose que f_ch mais utilise la puissance électrique au lieu de la vitesse de rotation     
        #Puissance Pelec pendant 10 secondes (intervalle de temps choisi de façon à avoir une vitesse de rotation constante sur la période)
        def residual(xdat,ydat):
            S,P=self.f_ch(xdat,10)
            r=ydat-P.Pu_elec
            return r
        
        if Pelec >0:
            x0=1500 
        else :
            x0=-1500
            
        
           
        Pu_cible=Pelec
            
        p1, success = optimization.leastsq(residual, x0, args=Pu_cible)
        #print(residual(p1,Pelec))
        S1,Dat_out=self.f_ch(p1,10)
        return S1,Dat_out, p1
    
    def f_ch2_S(self,Pelec):
        a,b,c=self.f_ch2(Pelec)
        return a
        
    
    def Tr(self,q,t,r): #Fonction similaire à la précédente où les bornes de la puissance sont à indiquer
        global N2011
        
        def sigN(i):
            T=np.zeros(i)
            for k in range(i):
                x=alea.randn()
                if (x>-1) and (x<1) :
                    T[k]=x
                elif x<0:
                    T[k]=-1
                else :
                    T[k]=1 
            return T
        
        
        P0=q/t
        N=t/10
        time=np.linspace(0,t,num=N+1)
        TS=[]
        Tpow=[]
        sigN1=sigN(N)
        TS.append(self)
        i=0
        for k in time[1:]:
            P=P0+sigN1[i]*r*(100-abs(P0))
            TS.append(TS[-1].f_ch2(P)[0])
            Tpow.append(TS[-1].f_ch2(P)[1])
            i+=1
        return TS,Tpow,time    
    
    def T(self,q,t): #Renvoie les états successifs et les puissances mises en jeu à partir d'une quantité d'énergie à emmagasiner q sur un temps t
        TS,Tpow,time=self.Tr(q,t,0.5)
        return TS,Tpow,time
    

    def TSn(self,q,t):
        TS=self.T(q,t)[0]
        return TS
    
    def Tpow(self,q,t):
        TP=self.T(q,t)[1]
        return TP
    """
    # Fonction qui recherche les bornes permettant de maximiser le rendement    
    def Bestr(self,q,t):
        compt=0
        def y(r):
            global compt
            tab=np.zeros(10)
            for k in range(10):
                Pu=self.Tr(q,t,r)[1]
                eta=ListToTab(Latt(Pu,'eta_tot'))
                tab[k]=np.mean(eta)
            m=np.mean(tab)
            compt+=1
            print(compt)
            return 1-m
        res = optimization.minimize_scalar(y,bounds=(0, 0.9), method='bounded')
        #M=lmfit.minimize(y,p,method='nelder')
        return res.x
    """
    
    # Fonction qui recherche les bornes permettant de maximiser l'énergie stocker    
    def Bestr(self,q,t):
        compt=0
        def y(r):
            global compt
            tab=np.zeros(10)
            for k in range(10):
                TSf=self.Tr(q,t,r)[0]
                E=ListToTab(Latt(TSf,'E'))
                tab[k]=E[-1]
            m=np.mean(tab)
            compt+=1
            print(compt)
            return 1-m
        res = optimization.minimize_scalar(y,bounds=(0, 1), method='bounded')
        #M=lmfit.minimize(y,p,method='nelder')
        return res.x
    """
    # Fonction qui recherche les bornes permettant de minimiser les pertes    
    def Bestr(self,q,t):
        compt=0
        def y(r):
            global compt
            tab=np.zeros(10)
            for k in range(10):
                TPw=self.Tr(q,t,r)[1]
                E=Tab(TPw,'Pu_hydr')-Tab(TPw,'Pu_elec')
                tab[k]=np.linalg.norm(E,2)
            m=np.mean(tab)
            compt+=1
            print(compt)
            return m
        res = optimization.minimize_scalar(y,bounds=(0, 1), method='bounded')
        #M=lmfit.minimize(y,p,method='nelder')
        return res.x
    """
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

def DrawTab(Liste,Attribut,new=True,r=None,mrkr='+',ls=''):
    global compt_plot
    #if new:
    #    plt.show()
    units={'P':'Pa', 'V':'m³', 'E':'kWh', 'eta_tot':'%', 'Pu_hydr':'kW','Pu_elec':'kW','Pu_meca':'kW',}
    time=Liste[-1]
    if Attribut in ['P','V','E']:
        T=Tab(Liste[0],Attribut)
        boo=False
    else:
        T=Tab(Liste[1],Attribut)
        boo=True
    if boo :
        time=time[1:]
    if new:
        plt.figure(compt_plot)
    else:
        compt_plot+=1
    if r<> None:
        lbl1='r = '+str(r)
    else :
        lbl1=Attribut
    plt.plot(time,T,label=lbl1,marker=mrkr,linestyle=ls)
    plt.xlabel('[s]')
    plt.ylabel('['+units[Attribut]+']')
    plt.legend(loc='upper left')
    


"""

"""
U=['Pa','m³','kWh','kW','%']

"""

#S1,power=S0.f_ch(-1500,10)
#print(S1.__getattribute__('P'))

Pavg=[]
Eta_avg=[]


#for k in range (10):


r=[]
"""




compt=0.0

resT=[]


S0=S(15567190.3125)
print(S0.E/1744009560.0)
S1=S0.f_ch2_S(100e3)

fin=temps.time()

#print(fin-debut)

print(S1.P)
print(S1.E/1744009560.0)

"""
i=0

Pu0=60
#for k in [60, 40, -40]: 
k=60
if os.path.isfile('r'+str(k)): 
    with open('r'+str(k), 'rb') as fichier:
        mon_depickler = pickle.Unpickler(fichier)
        rlist=mon_depickler.load()
        #T=mon_depickler.load()

boo2=True
for j in rlist[:2]:
    T=S0.Tr(60*3600*1000,3600,j)
    DrawTab(T,'Pu_elec',new=boo2,r=j)
    boo2=False
    
    
boo2=True
for j in rlist[:2]:
    T=S0.Tr(60*3600*1000,3600,j)
    DrawTab(T,'E',new=boo2,r=j)
    boo2=False
    
"""
"""
for j in range(2):
    debut=temps.time()
    r=S0.Bestr(k*3600*1000,3600)
    fin=temps.time()
    print(r)
    print(fin-debut)
    #T=S0.Tr(Pu0*3600*1000,3600,r)
    #T2=S0.Tr(Pu0*3600*1000,3600,0.2)
    rlist.append(r)
    
with open('r'+str(k), 'wb') as fichier:
    mon_pickler = pickle.Pickler(fichier)
    mon_pickler.dump(rlist)
"""

"""
DrawTab(T,'Pu_elec',new=True,r=0.045261488015355414,mrkr='+')
DrawTab(T2,'Pu_elec',new=False,r=0.2,mrkr='+')  
"""
"""
for k in range(5):
    if k==0:
        boo2=True
    else :
        boo2=False
    T=S0.Tr(100*3600*1000,3600,r)
    T2=S0.Tr(100*3600*1000,3600,0.2)  
    DrawTab(T,'E',new=boo2,r=r,mrkr='+')
    DrawTab(T2,'E',new=False,r=0.2,mrkr='+') 
"""    
"""
with open('donnees', 'wb') as fichier:
    mon_pickler = pickle.Pickler(fichier)
    mon_pickler.dump(r)
    mon_pickler.dump(resT)


Sf=S(250e5)
"""

"""
T3600=S0.T(100*3600*1000,3600)

S3600=S0.TSn(100*3600*1000,3600)

Pu3600=S0.Tpow(100*3600*1000,3600)

plt.figure(1)
plt.plot(np.linspace(0,3600,num=360),Latt(Pu3600,'Pu_elec'),'+')

plt.figure(2)
plt.plot(np.linspace(0,3600,num=361),ListToTab(Latt(S3600,'E'))/3600/1000,'+')
"""
"""
TPu=ListToTab(Tab(Pu3600,'Pu_elec'))

Eta=ListToTab(Tab(Pu3600,'eta_tot'))

t=T3600[2]


Pu_moy=simps(TPu,t[1:])/(t[-1]-t[1])
Pavg.append(Pu_moy)    

if abs(Pu_moy-100e3)/100e3 <0.01:
    compt+=1

eta_moy=simps(Eta,t[1:])/(t[-1]-t[1])
Eta_avg.append(eta_moy)
n=k+1

print('Part des puissances moyennes considérés égales à 100 kW')
print(str(compt/n*100.00)+' %')

#T0=T3600[-1].TSn(-100*3600*1000,3600)
"""
"""
rpm1=1500*np.ones(50)

TS50,Tpow50=S0.f_ch_T(rpm1,500)

S2=S0.f_ch2_S(-100e3)

S3=S2.f_ch2_S(-100e3)
"""
