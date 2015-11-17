# -*- coding: utf-8 -*-
"""
Created on Mon Sep 07 17:00:59 2015

@author: Alexandre Malon
"""

k4=lambda Pu,coeff: [-(coeff[0]/(kk/100e3)**coeff[1]+coeff[2]+coeff[3]*(kk/100e3)+coeff[4]*(kk/100e3)**2+coeff[5]*(kk/100e3)**3+coeff[6]*(kk/100e3)**4) for kk in Pu]

def residual_k4(coeff,y_k):
    return [y_k[i]-k4(Pu_test,coeff)[i] for i in range(len(y_k))]

coeff_4_0=[5.00E-33, 2.4, -5.0e-31, 4.0e-35, -8.0e-40, 8.0e-45, -3.0e-50]


coeff_4_opt, success = optimization.leastsq(residual_k4, coeff_4_0, args=[k[0] for k in z])

x=Pu_test
y=[k[0] for k in z]

plt.figure()
plt.plot(x,y,'+')
    
    
x_new = np.linspace(min(x), max(x), 50)
y_new = k4(x_new,coeff_4_opt)
plt.plot(x_new,y_new)

""""""
k3=lambda Pu,coeff: [(coeff[0]/(kk/100e3)**coeff[1]+coeff[2]+coeff[3]*(kk/100e3)+coeff[4]*(kk/100e3)**2+coeff[5]*(kk/100e3)**3+coeff[6]*(kk/100e3)**4) for kk in Pu]

def residual_k3(coeff,y_k):
    return [y_k[i]-k3(Pu_test,coeff)[i] for i in range(len(y_k))]

coeff_3_0=[5.00E-25, 2.4, 0, 0,0,0,0]


coeff_3_opt, success = optimization.leastsq(residual_k3, coeff_3_0, args=[k[1] for k in z])

x=Pu_test
y=[k[1] for k in z]

plt.figure()
plt.plot(x,y,'+')
    
    
x_new = np.linspace(min(x), max(x), 50)
y_new = k3(x_new,coeff_3_opt)
plt.plot(x_new,y_new)

""""""
k2=lambda Pu,coeff: [-(coeff[0]/(kk/100e3)**coeff[1]+coeff[2]+coeff[3]*(kk/100e3)+coeff[4]*(kk/100e3)**2+coeff[5]*(kk/100e3)**3+coeff[6]*(kk/100e3)**4) for kk in Pu]

def residual_k2(coeff,y_k):
    return [y_k[i]-k2(Pu_test,coeff)[i] for i in range(len(y_k))]

coeff_2_0=[5.00E-17, 2.0, 0,0,0,0,0]


coeff_2_opt, success = optimization.leastsq(residual_k2, coeff_2_0, args=[k[2] for k in z])

x=Pu_test
y=[k[2] for k in z]

plt.figure()
plt.plot(x,y,'+')
    
    
x_new = np.linspace(min(x), max(x), 50)
y_new = k2(x_new,coeff_2_opt)
plt.plot(x_new,y_new)

""""""
k1=lambda Pu,coeff: [(coeff[0]/(kk/100e3)**coeff[1]+coeff[2]+coeff[3]*(kk/100e3)+coeff[4]*(kk/100e3)**2+coeff[5]*(kk/100e3)**3+coeff[6]*(kk/100e3)**4) for kk in Pu]

def residual_k1(coeff,y_k):
    return [y_k[i]-k1(Pu_test,coeff)[i] for i in range(len(y_k))]

coeff_1_0=[0.0, 0.8, 0,0,0,0,0]


coeff_1_opt, success = optimization.leastsq(residual_k1, coeff_1_0, args=[k[3] for k in z])

x=Pu_test
y=[k[3] for k in z]

plt.figure()
plt.plot(x,y,'+')
    
    
x_new = np.linspace(min(x), max(x), 50)
y_new = k1(x_new,coeff_1_opt)
plt.plot(x_new,y_new)

""""""
k0=lambda Pu,coeff: [1/(coeff[0]/(kk/100.0e3)**coeff[1]+coeff[2]+coeff[3]*(kk/100.0e3)+coeff[4]*(kk/100.0e3)**2+coeff[5]*(kk/100.0e3)**3+coeff[6]*(kk/100.0e3)**4) for kk in Pu]

def residual_k0(coeff,y_k):
    return [y_k[i]-k0(Pu_test,coeff)[i] for i in range(len(y_k))]

coeff_0_0=[1, 0.6, 0.0,0.0,0.0,0.0,0.0]


coeff_0_opt, success = optimization.leastsq(residual_k0, coeff_0_0, args=[k[4] for k in z])



x=Pu_test
y=[k[4] for k in z]

plt.figure()
plt.plot(x,y,'+')
    
    
x_new = np.linspace(min(x), max(x), 50)
y_new = k0(x_new,coeff_0_opt)
plt.plot(x_new,y_new)


def eta_corr(Pu,P):
    k_0=k0(Pu,coeff_0_opt)
    k_1=k1(Pu,coeff_1_opt)
    k_2=k2(Pu,coeff_2_opt)
    k_3=k3(Pu,coeff_3_opt)
    k_4=k4(Pu,coeff_4_opt)
    
    eta= k_0+k_1*P+k_2*P**2+k_3*P**3+k_4*P**4
    return eta


plt.figure(1)
for Pu in Pu_test:
    Pu=[Pu]
    x=np.linspace(100e5,250e5,50)
    y=eta_corr(Pu,x)
    plt.plot(x,y)

