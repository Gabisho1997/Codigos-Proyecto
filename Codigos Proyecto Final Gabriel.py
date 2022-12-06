# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 18:09:04 2022

@author: Maquina
"""

import scipy.integrate as spi
import numpy as np
import pylab as pl
import scipy.optimize as optimize

import os
import sys
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#%matplotlib inline
#%precision 4
plt.style.use('ggplot')
np.random.seed(12345)
# Modelo en términos de ODEs
def ode_SIR1(INP, t, gammah, sigmah, betah, betap ,gammap,sigmap,p ):  
    Y = np.zeros((11))
    V = INP    
    
    
    
    Y[0] = -betap*V[0]*V[2]                        #dS_p/dt
    Y[1] = betap*V[0]*V[2]-gammap*V[1]             #dE_p/dt
    Y[2] =  gammap*V[1] -sigmap*V[2]               #dI_p/dt 
    Y[3]=  sigmap*V[2]                             #dD_p/dt
    Y[4]= -betah*V[4]*V[3]                         #dS_h/dt
    Y[5]= betah*V[4]*V[3]-gammah*V[5]              #dE_h/dt
    Y[6]=  gammah*V[5]*(1-p)-sigmah*V[6]           #dI_h/dt
    Y[7]= p*gammah*V[6]                            #dV_h/dt
    Y[8]= sigmah*V[6]                              #dD_h/dt                     
    Y[9]= gammap*V[1]                        #Incidencia acumulada perros
    Y[10]= gammah*V[6]                     #incidencia acumulada humanos
    return Y  

# parámetros y condiciones iniciales
betah= 0.2
gammah= 0.35
sigmah = 0.17
betap= 0.11
gammap = 0.2
sigmap=1.35
p= 0.83

Nh=1000
Np=500
Ip0 = 1
Dp0= 0
Ep0=0
Sp0 = Np-Ip0-Ep0
Ih0=0
Dh0=0
Eh0=0
Vh0=0
Sh0 = Nh-Ih0-Eh0-Vh0

INPUT = (Sp0, Ip0, Ep0, Dp0, Sh0,Eh0,Ih0,Vh0, Dh0, 0.0,0.0)

t_start = 0.0; 
t_end = 30; 
t_inc = 0.5
t_range = np.arange(t_start, t_end+t_inc, t_inc)


SOL = spi.odeint(ode_SIR1,INPUT,t_range,args=(gammah, sigmah, betah, betap ,gammap,sigmap,p))
incp=np.diff(SOL[:,9]) #incidencia perros
pl.plot(incp, '-g', label='Incidencias perros')

inch=np.diff(SOL[:,10]) #incidencias humanos
pl.plot(inch, '-r', label='Incidencias humanos')
#pl.plot(SOL[:,2])
pl.legend(loc=0)
#pl.title('Simulación de la población de perros')
#pl.xlabel('Time')
pl.ylabel('Incidencias')
pl.show()

ruip=np.random.poisson(lam=incp)
ruih=np.random.poisson(lam=inch)

pl.plot(incp, label='Incidencias perros  por semana')
pl.plot(ruip, label='Incidencias perros con ruido por semana')



pl.plot(inch, label='Incidencias humanos  por semana')
pl.plot(ruih,  label='Incidencias humanos con ruido por semana')
#pl.plot(SOL[:,2])
pl.legend(loc=0)
#pl.title('Simulación de la población de perros')
#pl.xlabel('Time')
pl.ylabel('Incidencias con ruido Poisson')
pl.show()








pl.figure(figsize=(15, 15))
pl.subplot(211)
pl.plot(SOL[:,0], '-g', label='Susc perros')
pl.plot(SOL[:,1], '-r', label='Exp perros')
pl.plot(SOL[:,2], '-k', label='infec perros')
pl.plot(SOL[:,3], '-c', label='death perros')
pl.legend(loc=0)
pl.title('Simulación de la población de perros')
#pl.xlabel('Time')
pl.ylabel('Sp, Ep,Ip and Dp')
#pl.subplot(212)
#pl.plot(SOL[:,1], '-r', label='Infec indiv')
pl.legend(loc=0)
pl.xlabel('Time')
#pl.ylabel('Infectious')
pl.show()




pl.figure(figsize=(15, 15))
pl.subplot(211)
pl.plot(SOL[:,4], '-g', label='Susc humanos')
pl.plot(SOL[:,5], '-r', label='Exp humanos')
pl.plot(SOL[:,6], '-k', label='infec humanos')
pl.plot(SOL[:,7], '-b', label='Vacc humanos')
pl.plot(SOL[:,8], '-c', label='death humanos')
pl.legend(loc=0)
pl.title('Simulación de la población de humanos')
#pl.xlabel('Time')
pl.ylabel('Sh, Eh,Ih,Vh and Dh')
#pl.subplot(212)
#pl.plot(SOL[:,1], '-r', label='Infec indiv')
pl.legend(loc=0)
pl.xlabel('Time')
#pl.ylabel('Infectious')
pl.show()










def ode_SIR(INP,t,ps): 
    try:
        betap = ps['betap'].value
        gammap = ps['gammap'].value
        sigmap = ps['sigmap'].value
        betah = ps['betah'].value
        gammah = ps['gammah'].value
        sigmah = ps['sigmah'].value
        p = ps['p'].value
        Nh = ps['Nh'].value
        Np = ps['Np'].value
    except:
        betah, gammah, sigmah, betap,gammap,sigmap,p,Nh = ps,Np=ps
    Y = np.zeros((11))
    V = INP    
    Y[0] = -betap*V[0]*V[2]                        #dS_p/dt
    Y[1] = betap*V[0]*V[2]-gammap*V[1]             #dE_p/dt
    Y[2] =  gammap*V[1] -sigmap*V[2]               #dI_p/dt 
    Y[3]=  sigmap*V[2]                             #dD_p/dt
    Y[4]= -betah*V[4]*V[3]                         #dS_h/dt
    Y[5]= betah*V[4]*V[3]-gammah*V[5]              #dE_h/dt
    Y[6]=  gammah*V[5]*(1-p)-sigmah*V[6]           #dI_h/dt
    Y[7]= p*gammah*V[6]                            #dV_h/dt
    Y[8]= sigmah*V[6]                              #dD_h/dt                     
    Y[9]= gammap*V[1]                        #Incidencia acumulada perros
    Y[10]= gammah*V[6]                     #incidencia acumulada humanos
    return Y   
    
       

def g(t,INP,ps):
    SOL = spi.odeint(ode_SIR,INPUT,t_range,args=(ps,))
    return SOL





###Ajuste del modelo 
def loss_funcnorm(params):
  betah,betap,gammah,gammap,sigmah,sigmap,p= params
  Solc=spi.odeint(ode_SIR1,INPUT,t_range,args=( betah,betap,gammah,gammap,sigmah,sigmap,p))
  incp2=np.diff(Solc[:,9])
  inch2=np.diff(Solc[:,10])
  loss1=np.log(np.sum((incp2-ruip)**2))
  loss2=np.log(np.sum((inch2-ruih)**2))
  return np.array([loss1+loss2])

initial_guess=np.array([0.08,0.08, 0.7, 1.0, 0.36, 1.2,0.6])
bnds= ((0.0001, 2), (0.1, 2), (0.1, 2), (0.1, 2), (0.1, 2),(0, 2), (0, 1))
result = optimize.minimize(loss_funcnorm,initial_guess,bounds=bnds)

result.x


SOLec = spi.odeint(ode_SIR1,INPUT,t_range,args=(0.1806, 0.2123, 0.8902, 1.0008, 0.1035, 1.32  , 0.595))


pl.figure(figsize=(15, 15))
pl.subplot(211)

pl.plot(SOL[:,0], '-g', label='Susc perros')
pl.plot(SOL[:,1], '-r', label='Exp perros')
pl.plot(SOL[:,2], '-k', label='infec perros')
pl.plot(SOL[:,3], '-c', label='death perros')

pl.plot(SOLec[:,0], '--g', label='Susc perros')
pl.plot(SOLec[:,1], '--r', label='Exp perros')
pl.plot(SOLec[:,2], '--k', label='infec perros')
pl.plot(SOLec[:,3], '--c', label='death perros')
pl.legend(loc=0)
pl.title('Simulación de la población de perros')
#pl.xlabel('Time')
pl.ylabel('Sp, Ep,Ip and Dp')
#pl.subplot(212)
#pl.plot(SOL[:,1], '-r', label='Infec indiv')
pl.legend(loc=0)
pl.xlabel('Time')
#pl.ylabel('Infectious')
pl.show()




pl.figure(figsize=(15, 15))
pl.subplot(211)
pl.plot(SOL[:,4], '-g', label='Susc humanos')
pl.plot(SOL[:,5], '-r', label='Exp humanos')
pl.plot(SOL[:,6], '-k', label='infec humanos')
pl.plot(SOL[:,7], '-b', label='Vacc humanos')
pl.plot(SOL[:,8], '-c', label='death humanos')



pl.plot(SOLec[:,4], '--g', label='Susc humanos')
pl.plot(SOLec[:,5], '--r', label='Exp humanos')
pl.plot(SOLec[:,6], '--k', label='infec humanos')
pl.plot(SOLec[:,7], '--b', label='Vacc humanos')
pl.plot(SOLec[:,8], '--c', label='death humanos')



pl.legend(loc=0)
pl.title('Simulación de la población de humanos')
#pl.xlabel('Time')
pl.ylabel('Sh, Eh,Ih,Vh and Dh')
#pl.subplot(212)
#pl.plot(SOL[:,1], '-r', label='Infec indiv')
pl.legend(loc=0)
pl.xlabel('Time')
#pl.ylabel('Infectious')
pl.show()



incpec=np.diff(SOLec[:,9]) #incidencia perros
inchec=np.diff(SOL[:,10]) #incidencias humanos


pl.plot(ruip, label='Incidencias perros con ruido por semana')
pl.plot(incpec,"--", label='Incidencias perros con ruido por semana')

pl.plot(ruih,  label='Incidencias humanos con ruido por semana')
pl.plot(inchec,"--",  label='Incidencias humanos con ruido por semana')
#pl.plot(SOL[:,2])
pl.legend(loc=0)
#pl.title('Simulación de la población de perros')
#pl.xlabel('Time')
pl.title('Incidencias con ruido Poisson y estimadas con funciòn de pèrdida')
pl.show()








