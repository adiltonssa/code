#!/usr/bin/env python
# coding: utf-8

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

def flow(P_inlet,Pf,FR,rou,Visc,di,l_s):
        Pf=P_inlet
        DP_B=[]
        for i in range(len(l_s)):

            rho=Dens*(Pf/P_inlet)**2
            sp=4*FR/(rho*np.pi*(0.001*di)**2)/3600
            Re=rho*sp*(di/1000)/(Visc/1000)
              
            res=Df2(sp,l_s[0][i],Visc,di,rou,Re)
        
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


    l=st.slider('Equiv. length (m): ',1,100,11)
    Elem=st.slider('Elements: ',1, 100,25)
   
    st.sidebar.header('Process inputs:')
    FR= st.number_input('Flow rate - Estimated (kg/h)',value=4000,min_value=10, max_value=20000)
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
    Dpr=P_out-P_in
    
               
st.caption('Application developed by Adilton Lopes da Silva (INEOS Polymers Engineering & Technology Support)')


