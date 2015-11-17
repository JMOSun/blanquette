# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 10:58:07 2015

@author: Alexandre Malon
"""

import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import xlrd

path=r'\\solar-lyon\Data\SunRZone-InterneLyon\SunR_SmE\SCD-i\02_Maquette\Pour ASEO\coeffs.xlsx'

book1 = xlrd.open_workbook(path)
first_sheet1 = book1.sheet_by_index(0)

kc=[[k.value for k in first_sheet1.col(i)] for i in range(5)]


def f_rdm(P,Pu):
    k0=lambda Pu,coeff: [[1/(coeff[0]/(kk/100.0e3)**coeff[1]+coeff[2]+coeff[3]*(kk/100.0e3)+coeff[4]*(kk/100.0e3)**2+coeff[5]*(kk/100.0e3)**3+coeff[6]*(kk/100.0e3)**4) for kk in kP] for kP in Pu]
    k1=lambda Pu,coeff: [[(coeff[0]/(kk/100e3)**coeff[1]+coeff[2]+coeff[3]*(kk/100e3)+coeff[4]*(kk/100e3)**2+coeff[5]*(kk/100e3)**3+coeff[6]*(kk/100e3)**4) for kk in kP] for kP in Pu]
    k2=lambda Pu,coeff: [[-(coeff[0]/(kk/100e3)**coeff[1]+coeff[2]+coeff[3]*(kk/100e3)+coeff[4]*(kk/100e3)**2+coeff[5]*(kk/100e3)**3+coeff[6]*(kk/100e3)**4) for kk in kP] for kP in Pu]
    k3=lambda Pu,coeff: [[(coeff[0]/(kk/100e3)**coeff[1]+coeff[2]+coeff[3]*(kk/100e3)+coeff[4]*(kk/100e3)**2+coeff[5]*(kk/100e3)**3+coeff[6]*(kk/100e3)**4) for kk in kP] for kP in Pu]
    k4=lambda Pu,coeff: [[-(coeff[0]/(kk/100e3)**coeff[1]+coeff[2]+coeff[3]*(kk/100e3)+coeff[4]*(kk/100e3)**2+coeff[5]*(kk/100e3)**3+coeff[6]*(kk/100e3)**4) for kk in kP] for kP in Pu]
    
    
    rdm=k0(Pu,kc[0]) + k1(Pu,kc[1])*P+ k2(Pu,kc[2])*P**2+ k3(Pu,kc[3])*P**3+ k4(Pu,kc[4])*P**4
    return rdm

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'




##################################################
# Define our surface

x = np.arange(100.0e5, 250.0e5, 1.0e5)
z = np.arange(10.0e3, 110.0e3, 100.0)
X, Z = np.meshgrid(x, z)

Y = f_rdm(X,Z)



##################################################
# Make contour labels using creative float classes
# Follows suggestion of Manuel Metz
##################################################
plt.figure()

# Basic contour plot
levels = np.linspace(10.0, 110.0, 11)
CS = plt.contour(X, Y, Z/1000.,colors='black',levels=levels)

# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'

class nf(float):
     def __repr__(self):
         str = '%.1f' % (self.__float__(),)
         if str[-1]=='0':
             return '%.0f' % self.__float__()
         else:
             return '%.1f' % self.__float__()

# Recast levels to new class
CS.levels = [nf(val) for val in CS.levels ]

# Label levels with specially formatted floats
if plt.rcParams["text.usetex"]:
     fmt = r'%r \kW'
else:
     fmt = '%r kW'
plt.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=10)

plt.xlim([1.0e7, 2.5e7])
plt.ylim([0.1, 1.0])

plt.xlabel('Pression [Pa]')
plt.ylabel('Rendement [-]')

plt.grid()
plt.show()