# -*- coding: utf-8 -*-
"""
Created on Thu Sep 03 10:54:15 2015

@author: Alexandre Malon
"""
import numpy as np
import xlrd
import CoolProp.CoolProp as CP
import matplotlib.pylab as plt
import scipy.optimize as optimization

V0=200
Vf=V0/2.5
V1=200
d0=0.1
P0=100e5
T0=293.15
rho0=CP.PropsSI('D','T',T0,'P',P0,'Air')
cyl=242

Pu_test=np.linspace(10.0e3,100.0e3,10)

plt.ion()



path1=r'./eta_elec_100kW.xlsx'

book1 = xlrd.open_workbook(path1)
first_sheet1 = book1.sheet_by_index(0)


path2=r'./eta_elec_100kW.xlsx'

book2 = xlrd.open_workbook(path2)
first_sheet2 = book2.sheet_by_index(0)

def rdm_el(Pu):
    x=[k.value for k in first_sheet1.col(0)]
    y=[k.value for k in first_sheet1.col(1)]
    return np.interp(Pu,x,y)


def var_el(Pu):
    x=[k.value for k in first_sheet2.col(0)]
    y=[k.value for k in first_sheet2.col(1)]
    return np.interp(Pu,x,y)    
    
    
def eta_hydrau(rpm,P):
    S=1/(-2.987894025497691342e+00*rpm/P-6.333293557613091096e-01)+1-1.758296915625696477e-03
    S=np.clip(S,0.01,1)
    return S

def Cf(rpm):
    a=0.000000014
    b=0.000000
    c=0.23
    return a*rpm**2+b*rpm+c

def dQ(rho):
    return 9.021557E+01*np.log(rho) - 2.368101E+01
    
def f_eta_stock(Pstar):
    z=[-0.85423949,  1.85148906]
    f = np.poly1d(z)
    return f(Pstar)

def f_Pstar(Psurf):
    z=[ -4.31007918e-08,   1.57353400e-04,   1.00956153e+00]
    f = np.poly1d(z)
    return f(Psurf)    
    
z=[]

for Pu0 in Pu_test:

    tf=6*3600*(100e3/Pu0)**1.2
    
    rpm0=2000*Pu0/100e3
    
    def f_ch(rpm,tf):
        
        class output(object):
            Pu_meca = None
            Pu_elec = None
            Pu_hydr = None
            eta_tot = None
            P       = None
            V       = None
            Pn      = None
        
        Nt=len(rpm)+1
        time=np.linspace(0,tf,Nt)
        V=np.zeros(Nt)
        Pn=np.zeros(Nt)
        rho=np.zeros(Nt)
        
        
        V[0]=V1
        rho[0]=rho0*V0/V1
        Pn[0]=P0
        
        
        Q=np.zeros(Nt-1)
        Psurf=np.zeros(Nt-1)
        Pstar=np.zeros(Nt-1)
        eta_stock=np.zeros(Nt-1)
        P=np.zeros(Nt-1)
        Pu_meca=np.zeros(Nt-1)
        Pu_elec=np.zeros(Nt-1)
        Pu_elec_res=np.zeros(Nt-1)
        Pu_hydr=np.zeros(Nt-1)
        eta_tot=np.zeros(Nt-1)
        
        for i in range(len(rpm)):
            
            Q_ideal=cyl*rpm[i]/1.0e6/60.
            V_ideal=V[i]-Q_ideal*(time[i+1]-time[i])
            rho_ideal=rho0*V0/V_ideal
            P_ideal=CP.PropsSI('P','T',T0,'D',rho_ideal,'Air')
            
            err=1.    
            Pp=P_ideal   
            
            while err>1e-3:
                Qq=cyl*rpm[i]/1.0e6/60.*eta_hydrau(rpm[i],Pp/1.0e5)
                Vv=V[i]-Qq*(time[i+1]-time[i])
                rho_r=rho0*V0/Vv
                Pp2=CP.PropsSI('P','T',T0,'D',rho_r,'Air')
                err=abs(Pp2-Pp)/Pp
                Pp=Pp2
            V[i+1]=Vv
            Pn[i+1]=CP.PropsSI('P','T',T0,'D',rho0*V0/Vv,'Air')
            rho[i+1]=rho_r
            P[i]=Pp
            Q[i]=Qq
            Psurf[i]=(dQ(rho[i+1])-dQ(rho[i]))*1000*rho0*V0/(2*(V[i]+V[i+1])/d0)/(time[i+1]-time[i])
            Pstar[i]=f_Pstar(Psurf[i])
            eta_stock[i]=f_eta_stock(Pstar[i])
            Couple=(cyl*P[i]/1.0e5/63+Cf(rpm[i])*cyl/4.9)/np.sqrt(eta_stock[i])
            Pu_meca[i]=Couple*rpm[i]*2*np.pi/60
            Pu_elec[i]=Pu_meca[i]/rdm_el(Pu_meca[i])
            Pu_elec_res[i]=Pu_elec[i]/var_el(Pu_elec[i])
            Pu_hydr[i]=P[i]*Q[i]
            eta_tot[i]=Pu_hydr[i]/Pu_elec_res[i]
            
            out=output()
            out.eta_tot=eta_tot
            out.Pu_hydr=Pu_hydr
            out.Pu_elec=Pu_elec_res
            out.Pu_meca=Pu_meca
            out.P=P
            out.Pn=Pn
            out.V=V      
            
        return out
    
    def residual(xdat,ydat):
        r=ydat-f_ch(xdat,tf).Pu_elec
        return r
    
    
    
    
    x0=rpm0*np.ones(50)
    
    
    
    Pu_cible=Pu0*np.ones(len(x0))
    
    
    
    
    p1, success = optimization.leastsq(residual, x0, args=Pu_cible)
       
    Dat_out=f_ch(p1,tf)
    
    time=np.linspace(0,tf,len(Dat_out.Pu_elec)+1)
    
    W_in=np.zeros(len(time))
    
    for i in range(len(Dat_out.Pu_elec)):
        W_in[i+1]=W_in[i]+Dat_out.Pu_hydr[i]*(time[i+1]-time[i])
    
    
    """
    plt.figure(0)
    x=W_in/(5*3600*100e3)
    y=Dat_out.Pn/250e5
    
    plt.plot(x,y,'+')
    
    zf=np.polyfit(x, y, 3)
    fz = np.poly1d(zf)
    
    x_new = np.linspace(min(x), max(x), 50)
    y_new = fz(x_new)
    
    plt.plot(x_new,y_new)    
    plt.xlabel('Energie stockée [J]')
    plt.ylabel('Pression [Pa]')
    
    plt.figure(1)
    x=Dat_out.P
    y=Dat_out.eta_tot
    
    plt.plot(x,y,'+')
    plt.xlim([100e5, 250e5])
    
    # calculate polynomial
    z.append(np.polyfit(x, y, 4))
    f = np.poly1d(z[-1])
    
    # calculate new x's and y's
    x_new = np.linspace(min(x), max(x), 50)
    y_new = f(x_new)
    
    plt.plot(x_new,y_new)
    plt.xlim([100e5, 250e5])
    
    
    plt.figure(2)
    plt.plot(time[1:],Dat_out.P)
    
    plt.figure(3)
    plt.plot(time,Dat_out.V)
    """

