# -*- coding: utf-8 -*-
"""Created on Sat Mar 14 00:19:00 2020
@author: uribe"""

from newton_raphson import newton_raphson
import sympy as sym
x= sym.Symbol('x',real=True)
import math as math
#Parameter known
row=1000 #density of water Kg/m3
g=9.81 # gravity m/s2
#Headrace tunnel: (i)
ri=2.20
#Concrete lining: (a)
ra=2.40; ka=1E-8;
#Gouted zone: (g)
rg=2.50; kg=1E-7
#Pressure internal
Pi=4.2
##############SOLVE Pa: Concrete lining ###################
#Fuction (Bouvard 1975)
F_sim=x-2*sym.pi*ka/sym.log(ra/ri)*(1E+5*Pi/(row*g)-x/(2*sym.pi*kg)*sym.log(x/(sym.pi*kg*ra))-3/4*ra)
#Derivative
DF_sim=sym.diff(F_sim,x)
#lambdify for use fuction
F= sym.lambdify([x],F_sim,'numpy')
DF = sym.lambdify([x],DF_sim,'numpy')
#Eval
f = lambda x: F(x)
Df = lambda x: DF(x)
qa=newton_raphson(f,Df,1,1E-8,10) #Solve
print('qa=',"%.2E" % qa,'m3/s')
qa_km=qa*1E+6 #l/s/km
print('qa_km=',"%.2f" % qa_km,'l/s/km')
Pa=Pi-qa*(1E-5*row*g*math.log(ra/ri)/(2*math.pi*ka))
print('Pa=',"%.2f" % Pa,'bares')

