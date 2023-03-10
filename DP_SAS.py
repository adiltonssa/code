#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# In[ ]:
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
from PIL import Image

password = "streamlit123"

def check_password():
    """Returns `True` if the user had the correct password."""

    def password_entered():
        """Checks whether a password entered by the user is correct."""
        if st.session_state["password"] == st.secrets["password"]:
            st.session_state["password_correct"] = True
            del st.session_state["password"]  # don't store password
        else:
            st.session_state["password_correct"] = False

    if "password_correct" not in st.session_state:
        # First run, show input for password.
        st.text_input(
            "Password", type="password", on_change=password_entered, key="password"
        )
        return False
    elif not st.session_state["password_correct"]:
        # Password not correct, show input + error.
        st.text_input(
            "Password", type="password", on_change=password_entered, key="password"
        )
        st.error("😕 Password incorrect")
        return False
    else:
        # Password correct.
        return True

def flow(P_inlet,Pf,FR,rou,Visc,di,l_s,Dens):
        Pf=P_inlet
        #DP_g=[]
        #sp_g=[]
        DP_B=[]
        #PF=[]
        rho=Dens*(Pf/P_inlet)**2
        
        for i in range(len(l_s)):

            sp=4*FR/(rho*np.pi*(0.001*di)**2)/3600
            Re=rho*sp*(di/1000)/(Visc/1000)
           
            res=Df2(sp,l_s[0][i],Visc,di,rou,Re)
            Pf=(Pf-res*rho*9.81*0.00001)
            rho=Dens*(Pf/P_inlet)**2
            
            DP_B.append(res*rho*9.81*0.00001)
        
        return sum(DP_B)
    #########################################################
    
    #Find root two equations:
def Df2(sp,ls,vis,di,rou,Re):
    
    def diff(pc):
        diff=100000*((0.0625/((np.log10(rou/di/3.7+5.74/(Re**0.9))))**2)-(pc*(di/1000)*9.81/(2*ls*(sp**2))))
        return diff
    
    pc = spo.root(diff,0)
    DP_f=pc.x
    
    return float(DP_f)
    #########################################################    
       
    
if check_password():

    image = Image.open('201103_INEOS_POETS_V_CMYK.jpg')
    st.sidebar.image(image)

    st.title('Delta Pressure Multisteps - Slurry Process INEOS')


    l=st.slider('Equiv. length (m): ',1,50,11)
    Elem=st.slider('Elements: ',1, 100,25)
   
    st.sidebar.header('Process inputs:')
    Temp= st.sidebar.number_input('Temperature (°C)',value=80,min_value=10, max_value=200)
    Dens= st.sidebar.number_input('Density (kg/m3)',value=11.0,min_value=1.0, max_value=1000.0,step=0.1)
    P_in= st.sidebar.number_input('Pressure inlet (bara):',value=3.5,min_value=1.0,step=0.1, max_value=100.0) 
    P_out= st.sidebar.number_input('Pressure out (bara):',value=3.0,min_value=0.1,step=0.1, max_value=100.0)
    Visc= st.sidebar.number_input('Viscosity (cP): ',value=0.0010,min_value=0.0001, max_value=5.000,step=0.001,format="%.4f")
    MW= st.sidebar.number_input('Molecular weight (g/mol): ',value=58.5,min_value=10.0, max_value=500.0,step=0.1)
    k=st.sidebar.number_input('Compressibility factor: ',value=1.15)

    st.sidebar.header('Pipe characteristics:')
    rou=st.sidebar.number_input('Roughness (mm): ',value=0.045,min_value=0.010, max_value=5.000,step=0.001,format="%.3f")
    di=st.sidebar.number_input('Diameter (mm): ',value=100.0,min_value=1.0, max_value=200.0,step=0.1)

    critical_speed=pow(k*8.3143*(Temp+273.15)/(MW/1000),0.5)
    Dpr=P_in-P_out
    
    #########################################################
    l_s=[]
    for i in range(Elem-1):
        l_s.append(((Elem-(i+1))*2-1)/(pow(Elem-1,2))*l)
    #########################################################
    
    l_s=pd.DataFrame(l_s) 

    #########################################################
    FR=1000
        
    crf=8
    inter=0
    tol=1
    Fl1=flow(P_in,P_out,FR,rou,Visc,di,l_s,Dens)
    dift=Dpr-Fl1

    while dift>0:
        FR=2*FR
        dift=Dpr-flow(P_in,P_out,2*FR,rou,Visc,di,l_s,Dens)
        
        
    FR0=0.5*FR
    FR1=2*FR
    
    while  crf>=tol | inter<100:

        FRm=(FR1+FR0)*0.5
    
        Flm=flow(P_in,P_out,FRm,rou,Visc,di,l_s,Dens)
        crm=Dpr-Flm
   
        cr=Dpr-Fl1

        if crm<0:
            FR1=FRm
            FR0=FR0
        else :
            FR1=FR1
            FR0=FRm
    
        inter=inter+1
        crf=1000*(FR1-FR0)/FR1
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric(label="Flow rate calculaded (Kg/h) =", value=round(FRm,2))
    with col2:
        st.metric(label="Differencial Pressure (bar) =", value=round(-Dpr,3))
             
    ########################################################    
        
    #st.write('Flow rate calculaded =:', FR1, 'kg/h.')
    #st.write('Differencial Pressure =:', -Dpr, 'Bar.')
    
    Ft=FRm

    DP_g=[]
    sp_g=[]
    DP_B=[]
    PF=[]

    Pf=P_in
    rho=Dens*(Pf/P_in)**2
    ########################################################
    for i in range(len(l_s)):

        
        sp=4*Ft/(rho*np.pi*(0.001*di)**2)/3600
        Re=rho*sp*(di/1000)/(Visc/1000)
        
        res=Df2(sp,l_s[0][i],Visc,di,rou,Re)  
        Pf=(Pf-res*rho*9.81*0.00001)
        rho=Dens*(Pf/P_in)**2
    
        DP_g.append(res)
        DP_B.append(res*rho*9.81*0.00001)
        sp_g.append(sp)
        PF.append(Pf)
        ########################################################
    df=pd.DataFrame()

    df['DP(Bar)']=DP_B
    df['DP(m)']=DP_g
    df['Velocity(m/s)']=sp_g
    df['P(bar)']=PF
    df['Elem DP (m)']=l_s

    col1, col2 = st.columns(2)
    with col1:
        var1=st.selectbox('Select the variable 1:', ['P(bar)','DP(Bar)','DP(m)','Velocity(m/s)','Elem DP (m)'])
    with col2:
        var2=st.selectbox('Select the variable 2:', ['Velocity(m/s)','DP(Bar)','DP(m)','P(bar)','Elem DP (m)'])

        if var1=='P(bar)':
            nam1='P(bar)'
        elif var1=='DP(Bar)':
            nam1='DP(Bar)'
        elif var1=='DP(m)':
            nam1='DP(m)'
        elif var1=='Velocity(m/s)':
            nam1='Velocity(m/s)'
        else:
            nam1='Elem DP (m)'

        if var2=='P(bar)':
            nam2='P(bar)'
        elif var2=='DP(Bar)':
            nam2='DP(Bar)'
        elif var2=='DP(m)':
            nam2='DP(m)'
        elif var2=='Velocity(m/s)':
            nam2='Velocity(m/s)'
        else:
            nam2='Elem DP (m)'   
    
        
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(go.Scatter(y=df[var1],x=df.index,name=nam1),secondary_y=False)
    fig.add_trace(go.Scatter(y=df[var2],x=df.index,name=nam2),secondary_y=True)
    fig.update_layout(height=600, width=800, title_text="Delta Pressure Multisteps - Slurry INEOS")
    fig.update_xaxes(title_text="Elements")

    st.plotly_chart(fig, use_container_width=True)
    
    vm=df['Velocity(m/s)'].max()

    if vm>=critical_speed:
        st.warning('Velocity maximum (m/s)> Critical Speed (m/s)', icon="⚠️")
    else:
        st.info('Velocity maximum (m/s) < Critical Speed(m/s)', icon="ℹ️")

    st.write('Critical speed(m/s): ', round(critical_speed,2))
    st.write('Maximum speed(m/s): ', round(vm,2))           

               
st.caption('Application developed by Adilton Lopes da Silva (INEOS Polymers Engineering & Technology Support)')

