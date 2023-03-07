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


st.set_page_config(
        page_title="Equlibrium",
        page_icon="🛰️",
        )

image = Image.open('201103_INEOS_POETS_V_CMYK.jpg')
st.sidebar.image(image)

st.title('Flow Rate - Slurry Process INEOS')


l=st.slider('Equiv. length (m): ',1,50,11)
Elem=st.slider('Elements: ',1, 100,25)
V=st.slider('Volume (m3): ',1, 50,5)
   
st.sidebar.header('Process inputs:')
Temp= st.sidebar.number_input('Temperature (°C)',value=80,min_value=10, max_value=200)
Dens= st.sidebar.number_input('Density (kg/m3)',value=20.0,min_value=1.0, max_value=1000.0,step=0.1)
Pi= st.sidebar.number_input('High Pressure(bara):',value=9.0,min_value=1.0,step=0.1, max_value=100.0) 
pi= st.sidebar.number_input('Low Pressure (bara):',value=1.35,min_value=0.1,step=0.1, max_value=100.0)
Visc= st.sidebar.number_input('Viscosity (cP): ',value=0.0010,min_value=0.0001, max_value=5.000,step=0.001,format="%.4f")
MW= st.sidebar.number_input('Molecular weight (g/mol): ',value=58.5,min_value=10.0, max_value=500.0,step=0.1)
k=st.sidebar.number_input('Compressibility factor: ',value=1.15)
DT=st.sidebar.number_input('Tempo (s): ',value=0.6,min_value=0.1, max_value=10.0)

st.sidebar.header('Pipe characteristics:')
rou=st.sidebar.number_input('Roughness (mm): ',value=0.045,min_value=0.010, max_value=5.000,step=0.001,format="%.3f")
di=st.sidebar.number_input('Diameter (mm): ',value=100.0,min_value=1.0, max_value=200.0,step=0.1)

#cst_H = st.sidebar.checkbox('HP constant')
#cst_L = st.sidebar.checkbox('LP constant')

#if cst_H:
#    C_H=0
#else:
#     C_H=1

#if cst_L:
#    C_L=0
#else:
#     C_L=1


cs=pow(k*8.3143*(Temp+273.15)/(MW/1000),0.5)

    #########################################################
l_s=[]
for i in range(Elem-1):
    l_s.append(((Elem-(i+1))*2-1)/(pow(Elem-1,2))*l)
    #########################################################
    
l_s=pd.DataFrame(l_s) 


#Functions

def flow(P_inlet,FR,rou,Visc,di,l_s,Dens):
    Pf=P_inlet
    DP_g=[]
    sp_g=[]
    DP_B=[]
    PF=[]
    rho=Dens*(Pf/P_inlet)**2
    
    for i in range(len(l_s)):
        
        sp=4*FR/(rho*np.pi*(0.001*di)**2)/3600
        Re=rho*sp*(di/1000)/(Visc/1000)  
        res=Df2(sp,l_s[0][i],Visc,di,rou,Re)
        Pf=(Pf-res*rho*9.81*0.00001)
        rho=Dens*(Pf/P_inlet)**2
        #DP_g.append(res)
        DP_B.append(res*rho*9.81*0.00001)
        #print(Pf)
        #print(rho)
        #print(sp)
        #sp_g.append(sp)
        #PF.append(Pf)
        
    return sum(DP_B)
def Df2(sp,ls,vis,di,rou,Re):
    
    def diff(pc):
        diff=100000*((0.0625/((np.log10(rou/di/3.7+5.74/(Re**0.9))))**2)-(pc*(di/1000)*9.81/(2*ls*(sp**2))))
        return diff
    
    pc = spo.root(diff,0)
    DP_f=pc.x
    
    return float(DP_f)
    #########################################################    
      
def flowrate(P_i,P_o,rou,Visc,di,l_s,Dens,cs):
    
    Dpr=P_i-P_o
    FR=1000
    Fl1=flow(P_i,FR,rou,Visc,di,l_s,20.6)
    dift=Dpr-Fl1
    crf=10
    inter=0
    tol=1
    
    #Find the speed to reach the DP defined: 
    while dift>0:
        FR=2*FR
        dift=Dpr-flow(P_i,FR,rou,Visc,di,l_s,Dens)
    
    #Speed> Critical speed - Find flow to garantee the seep lower than Critical speed.
    if DP_f(P_i,P_o,FR,rou,Visc,di,l_s,Dens)>cs:      
    
        #First guess is flow using critical speed/2:
        A=np.pi*(((di*0.5)/1000)**2)
        FR=A*cs*3600*Dens
        dif=0
        FF=FR/2
    
        # First guess less 500 kg to find a max speed lower than critical speed.
    
        while dif<0.1:
            FR=FR-500
            vm=DP_f(P_i,P_o,FR,rou,Visc,di,l_s,Dens)
            dif=cs-vm             
    
    #Speed> Critical speed: bissection method to find the root.
    else:
        
        FR0=0.5*FR
        FR1=2*FR
        crf=10


        while  crf>=tol | inter<100:

            FR=(FR1+FR0)*0.5
    
            Flm=flow(P_i,FR,rou,Visc,di,l_s,Dens)
            crm=Dpr-Flm

            if crm<0:
                FR1=FR
                FR0=FR0
            else :
                FR1=FR1
                FR0=FR
    
            inter=inter+1
            crf=10000*(FR1-FR0)/FR1            
            
        
    return FR

def DP_f(P_inlet,Pf,FR,rou,Visc,di,l_s,Dens):
    
    Pf=P_inlet
    DP_g=[]
    sp_g=[]
    DP_B=[]
    PF=[]
    for i in range(len(l_s)):
        rho=Dens*(Pf/P_inlet)**2
        sp=4*FR/(rho*np.pi*(0.001*di)**2)/3600
        Re=rho*sp*(di/1000)/(Visc/1000)              
        res=Df2(sp,l_s[0][i],Visc,di,rou,Re)
        
        Pf=(Pf-res*rho*9.81*0.00001)

        sp_g.append(sp)
       
    v_m=max(sp_g)     
    
    return v_m

def equil_H(HP,MI,DT,V,k,pi,Pi,Dens,FR_c):

    M_AD=DT/3600*FR_c

    N_M=MI-M_AD

    if HP==1:
        RHO_n=N_M/V
    else:
        RHO_n=Dens

    P_N=(RHO_n/Dens)**(1/k)*Pi

    return  RHO_n,MI,M_AD,N_M,P_N

def equil_L(LP,mi,DT,V,k,pi,Pi,rho_i,FR_c):

    m_ad=DT/3600*FR_c

    N_M=m_ad+mi

    if LP==1:
        rho_n=N_M/V
    else:
        rho_n=rho_i

    P_n=(rho_n/rho_i)**(1/k)*pi

    return  rho_n,mi,m_ad,N_M,P_n


RHO_n=[]
M_I=[]
M_AD=[]
N_M=[]
P_N=[]
rho_n=[]
m_i=[]
m_ad=[]
n_M=[]
p_n=[]
FL=[]

dif_p=1
inter=0
D=Dens
rho_i=((pi/Pi)**k)*Dens
mi=rho_i*V
MI=Dens*V

    
while dif_p>-1 | inter<100:
    
    
    F_out=flowrate(Pi,pi,rou,Visc,di,l_s,D,cs)

    O=equil_H(1,MI,DT,V,k,pi,Pi,D,F_out)
    o=equil_L(1,mi,DT,V,k,pi,Pi,rho_i,F_out)

    Pi=O[4]
    pi=o[4]
    MI=O[3]
    mi=o[3]
    D=O[0]
    rho_i=o[0]
   
    
    dif_p=100*(Pi-pi)
    inter=inter+1
    
    FL.append(F_out)
    RHO_n.append(O[0])
    M_I.append(O[1])
    M_AD.append(O[2])
    N_M.append(O[3])
    P_N.append(Pi)
    rho_n.append(o[0])
    m_i.append(o[1])
    m_ad.append(o[2])
    n_M.append(o[3])
    p_n.append(pi)

    
df_n=pd.DataFrame()
df_n['Flow Rate (kg/h)']=FL
df_n['LP (bar)']=p_n
df_n['LP - rho (kg/m3)']=rho_n
df_n['LP - New P (bar)']=p_n
df_n['LP - Mass initial (kg)']=m_i
df_n['LP - Add mass (kg)']=m_ad
df_n['LP - Mass Final (kg)']=n_M
df_n['HP (bar)']=P_N
df_n['HP - rho (kg/m3)']=RHO_n
df_n['HP - New P (bar)']=P_N
df_n['HP - Mass initial (kg)']=M_I
df_n['HP - Add mass (kg)']=M_AD
df_n['HP - Mass Final (kg)']=N_M

df_n.index=df_n.index*DT

col1, col2 = st.columns(2)
with col1:
    VAR1=st.selectbox('Select the variable 1 - Fig 1:', ['Flow Rate (kg/h)','LP (bar)','LP - rho (kg/m3)','LP - New P (bar)'
        'LP - Mass initial (kg)', 'LP - Add mass (kg)','LP - Mass Final (kg)','HP (bar)','HP - rho (kg/m3)',
        'HP - New P (bar)','HP - Mass initial (kg)','HP - Add mass (kg)','HP - Mass Final (kg)'])
with col2:
    VAR2=st.selectbox('Select the variable 2 - Fig 1:', ['LP (bar)','Flow Rate (kg/h)','LP - rho (kg/m3)','LP - New P (bar)'
        'LP - Mass initial (kg)', 'LP - Add mass (kg)','LP - Mass Final (kg)','HP (bar)','HP - rho (kg/m3)',
        'HP - New P (bar)','HP - Mass initial (kg)','HP - Add mass (kg)','HP - Mass Final (kg)'])

    
    if VAR1=='LP (bar)':
        NAM1='LP (bar)'
    elif VAR1=='Flow Rate (kg/h)':
        NAM1='Flow Rate (kg/h)'
    elif VAR1=='LP - rho (kg/m3)':
        NAM1='LP - rho (kg/m3)'
    elif VAR1=='LP - New P (bar)':
        NAM1='LP - New P (bar)'
    elif VAR1=='LP - Mass initial (kg)':
        NAM1='LP - Mass initial (kg)'
    elif VAR1=='LP - Add mass (kg)':
        NAM1='LP - Add mass (kg)'
    elif VAR1=='LP - Mass Final (kg)':
        NAM1='LP - Mass Final (kg)'
    elif VAR1=='HP (bar)':
        NAM1='HP (bar)'
    elif VAR1=='HP - rho (kg/m3)':
        NAM1='HP - rho (kg/m3)'
    elif VAR1=='HP - New P (bar)':
        NAM1='HP - New P (bar)'
    elif VAR1=='HP - Mass initial (kg)':
        NAM1='HP - Mass initial (kg)'
    elif VAR1=='HP - Add mass (kg)':
        NAM1='HP - Add mass (kg)'       
    else:
        NAM1='HP - Mass Final (kg)'

    if VAR2=='LP (bar)':
        NAM2='LP (bar)'
    elif VAR2=='Flow Rate (kg/h)':
        NAM2='Flow Rate (kg/h)'
    elif VAR2=='LP - rho (kg/m3)':
        NAM2='LP - rho (kg/m3)'
    elif VAR2=='LP - New P (bar)':
        NAM2='LP - New P (bar)'
    elif VAR2=='LP - Mass initial (kg)':
        NAM2='LP - Mass initial (kg)'
    elif VAR2=='LP - Add mass (kg)':
        NAM2='LP - Add mass (kg)'
    elif VAR2=='LP - Mass Final (kg)':
        NAM2='LP - Mass Final (kg)'
    elif VAR2=='HP (bar)':
        NAM2='HP (bar)'
    elif VAR2=='HP - rho (kg/m3)':
        NAM1='HP - rho (kg/m3)'
    elif VAR2=='HP - New P (bar)':
        NAM2='HP - New P (bar)'
    elif VAR2=='HP - Mass initial (kg)':
        NAM2='HP - Mass initial (kg)'
    elif VAR2=='HP - Add mass (kg)':
        NAM2='HP - Add mass (kg)'       
    else:
        NAM2='HP - Mass Final (kg)'  
    
        
fig = make_subplots(specs=[[{"secondary_y": True}]])
fig.add_trace(go.Scatter(y=df_n[VAR1],x=df_n.index,name=NAM1),secondary_y=False)
fig.add_trace(go.Scatter(y=df_n[VAR2],x=df_n.index,name=NAM2),secondary_y=True)
fig.update_layout(height=600, width=800, title_text="Flow Rate - High Pressure - Slurry INEOS")
fig.update_xaxes(title_text="Time(s)")

st.plotly_chart(fig, use_container_width=True)


col11, col12 = st.columns(2)

with col11:
    VAR_01=st.selectbox('Select the variable 1 - Fig 2:', ['Flow Rate (kg/h)','LP (bar)','LP - rho (kg/m3)','LP - New P (bar)'
        'LP - Mass initial (kg)', 'LP - Add mass (kg)','LP - Mass Final (kg)','HP (bar)','HP - rho (kg/m3)',
        'HP - New P (bar)','HP - Mass initial (kg)','HP - Add mass (kg)','HP - Mass Final (kg)'])
with col12:
    VAR_02=st.selectbox('Select the variable 2 - Fig 2:', ['LP (bar)','Flow Rate (kg/h)','LP - rho (kg/m3)','LP - New P (bar)'
        'LP - Mass initial (kg)', 'LP - Add mass (kg)','LP - Mass Final (kg)','HP (bar)','HP - rho (kg/m3)',
        'HP - New P (bar)','HP - Mass initial (kg)','HP - Add mass (kg)','HP - Mass Final (kg)'])

    if VAR_01=='LP (bar)':
        NAM_01='LP (bar)'
    elif VAR_01=='Flow Rate (kg/h)':
        NAM_01='Flow Rate (kg/h)'
    elif VAR_01=='LP - rho (kg/m3)':
        NAM_01='LP - rho (kg/m3)'
    elif VAR_01=='LP - New P (bar)':
        NAM_01='LP - New P (bar)'
    elif VAR_01=='LP - Mass initial (kg)':
        NAM_01='LP - Mass initial (kg)'
    elif VAR_01=='LP - Add mass (kg)':
        NAM_01='LP - Add mass (kg)'
    elif VAR_01=='LP - Mass Final (kg)':
        NAM_01='LP - Mass Final (kg)'
    elif VAR_01=='HP (bar)':
        NAM_01='HP (bar)'
    elif VAR_01=='HP - rho (kg/m3)':
        NAM_01='HP - rho (kg/m3)'
    elif VAR_01=='HP - New P (bar)':
        NAM_01='HP - New P (bar)'
    elif VAR_01=='HP - Mass initial (kg)':
        NAM_01='HP - Mass initial (kg)'
    elif VAR_01=='HP - Add mass (kg)':
        NAM_01='HP - Add mass (kg)'       
    else:
        NAM_01='HP - Mass Final (kg)'


    if VAR_02=='LP (bar)':
        NAM_02='LP (bar)'
    elif VAR_02=='Flow Rate (kg/h)':
        NAM_02='Flow Rate (kg/h)'
    elif VAR_02=='LP - rho (kg/m3)':
        NAM_02='LP - rho (kg/m3)'
    elif VAR_02=='LP - New P (bar)':
        NAM_02='LP - New P (bar)'
    elif VAR_02=='LP - Mass initial (kg)':
        NAM_02='LP - Mass initial (kg)'
    elif VAR_02=='LP - Add mass (kg)':
        NAM_02='LP - Add mass (kg)'
    elif VAR_02=='LP - Mass Final (kg)':
        NAM_02='LP - Mass Final (kg)'
    elif VAR_02=='HP (bar)':
        NAM_02='HP (bar)'
    elif VAR_02=='HP - rho (kg/m3)':
        NAM1='HP - rho (kg/m3)'
    elif VAR_02=='HP - New P (bar)':
        NAM_02='HP - New P (bar)'
    elif VAR_02=='HP - Mass initial (kg)':
        NAM_02='HP - Mass initial (kg)'
    elif VAR_02=='HP - Add mass (kg)':
        NAM_02='HP - Add mass (kg)'       
    else:
        NAM_02='HP - Mass Final (kg)'
   
    
        
fig_p = make_subplots(specs=[[{"secondary_y": True}]])
fig_p.add_trace(go.Scatter(y=df_n[VAR_01],x=df_n.index,name=NAM_01),secondary_y=False)
fig_p.add_trace(go.Scatter(y=df_n[VAR_02],x=df_n.index,name=NAM_02),secondary_y=True)
fig_p.update_layout(height=600, width=800, title_text="Flow Rate - Low Pressure - Slurry INEOS")
fig_p.update_xaxes(title_text="Time (s)")

st.plotly_chart(fig_p, use_container_width=True)
