#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Libraries standard
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly as py
import plotly.graph_objs as go
import plotly.express as px
from scipy.optimize import minimize
import scipy.optimize as spo
from plotly.subplots import make_subplots
import streamlit as st

#image = Image.open('201103_INEOS_POETS_V_CMYK.jpg')
#st.sidebar.image(image)



FR = st.slider('Flow (kg/h): ',1, 10000,1000)
st.write('Flow :', FR,'(kg/h)')

Elem=st.slider('Elements: ',1, 100,20)
st.write('Numbers of elements: ', FR)

col1, col2 = st.columns(2)

with col1: 
    st.header('Process inputs:')
    Temp= st.number_input('Temperature (°C)')
    Den= st.number_input('Density (kg/m3)')
    P_in= st.number_input('Pressure inlet (bara): ')
    P_out= st.number_input('Pressure inlet (bara): ')
    Visc= st.number_input('Viscosity (cP): ')
    MW= st.number_input('Molecular weight (g/mol): ')
    k=st.number_input('Compressibility factor: ')

with col2:
    st.header('Pipe characteristics:')
    rou=st.number_input('Roughness (mm): ')
    l=st.number_input('Equiv. length (m): ')
    Di=st.number_input('Diameter (mm): ')

critical_speed=pow(k*8.3143*(Temp+273.15)/(MW/1000),0.5)
    
l_s=[]
for i in range(Elem-1):
    l_s.append(((Elem-(i+1))*2-1)/(pow(Elem-1,2))*l)

def Df2(sp,ls,vis,di,rou):
    
    def diff(pc):
        diff=100000*((0.0625/((math.log10(rou/di/3.7+5.74/(Re**0.9))))**2)-(pc*(di/1000)*9.81/(2*ls*(sp**2))))
        return diff
    
    pc = spo.root(diff,0)
    DP_f=pc.x
    
    return float(DP_f)

DP_g=[]
sp_g=[]
DP_B=[]
PF=[]

Pf=P_inlet

for i in range(len(l_s)):

    rho=Dens*(Pf/P_inlet)**2
    sp=4*FR/(rho*np.pi*(0.001*di)**2)/3600
    Re=rho*sp*(diameter/1000)/(Visc/1000)
     
    res=Df2(sp,l_s[0][i],Visc,diameter,roughness)  
    Pf=(Pf-res*rho*9.81*0.00001)
    
    DP_g.append(res)
    DP_B.append(res*rho*9.81*0.00001)
    sp_g.append(sp)
    PF.append(Pf)

df=pd.DataFrame()

df['DP(Bar)']=DP_B
df['DP_g(m)']=DP_g
df['Velocity(m/s)']=sp_g
df['P(bar)']=PF
df['Elem DP (m)']=l_s

var1=st.selectbox('Select the variables:', ['DP(Bar)','DP_g','Velocity(m/s)','P(bar)','Elem DP (m)'])
var2=st.selectbox('Select the variables:', ['DP(Bar)','DP_g','Velocity(m/s)','P(bar)','Elem DP (m)'])

fig = make_subplots(specs=[[{"secondary_y": True}]])

fig.add_trace(go.Scatter(y=df[var1],x=df.index),secondary_y=False)
fig.add_trace(go.Scatter(y=df[var2],x=df.index),secondary_y=True)

st.plotly_chart(fig, use_container_width=True)

if df['Velocity'].max()>=critical_speed:
    st.warning('Velocity maximum: ',df['Velocity'].max(), '> Critical Speed :',critical_speed, icon="⚠️")
else:
    st.info('Velocity maximum: ',df['Velocity'].max(),' < Critical Speed :',critical_speed, icon="ℹ️")

st.caption('Application developed by Adilton Lopes da Silva (INEOS Polymers Engineering & Technology Support)')

