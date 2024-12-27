##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------LIBRERIAS---------------------------------------------------
##-----------------------------------------------------------------------------------------------------##
import handcalcs.render
from handcalcs.decorator import handcalc
import forallpeople as fap
fap.environment('alejandrov6',top_level="False")
import math
from math import pi,sqrt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import func_struct_abr

##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------FUNCIONES IMPRIMIR PROPIEDAES----------------------------

#Función que muestra propiedaes mecánicas de un acero
@handcalc(override="params",precision=2,jupyter_display=True)
def display_propiedades_material(prop_mat):
    nombre_material=prop_mat[0]
    Fy=prop_mat[1]*1000000*fap.kg/(fap.s**2*fap.m)
    Fu=prop_mat[2]*1000000*fap.kg/(fap.s**2*fap.m)
    Ry=prop_mat[3]
    Rt=prop_mat[4]
    Es=prop_mat[5]*1000000*fap.kg/(fap.s**2*fap.m)
    return nombre_material,Fy,Fu,Ry,Rt,Es

#Función que muestra propiedaes geométricas de perfiles tubulares cuadrados/rectangulares usando handcalcs y forallpeople
@handcalc(override="params",precision=2,jupyter_display=True)
def display_geomet_perfil_tub_rect(prop_perfil_tub):   
    name_perfil=prop_perfil_tub[0]
    d_perfil=((prop_perfil_tub[1]/1000 * fap.m).prefix("m"))
    tw_perfil=((prop_perfil_tub[2]/1000 * fap.m).prefix("m"))
    b_perfil=((prop_perfil_tub[3]/1000 * fap.m).prefix("m"))
    tf_perfil=((prop_perfil_tub[4]/1000 * fap.m).prefix("m"))
    Ag_perfil=((prop_perfil_tub[5]/1000000 * fap.m**2).prefix("c"))
    Ix_perfil=((prop_perfil_tub[6]/1000000000000 * fap.m**4).prefix("c"))  
    Iy_perfil=((prop_perfil_tub[7]/1000000000000 * fap.m**4).prefix("c"))  
    Sx_perfil=((prop_perfil_tub[8]/1000000000 * fap.m**3).prefix("c"))  
    Sy_perfil=((prop_perfil_tub[9]/1000000000 * fap.m**3).prefix("c"))  
    Zx_perfil=((prop_perfil_tub[10]/1000000000 * fap.m**3).prefix("c"))  
    Zy_perfil=((prop_perfil_tub[11]/1000000000 * fap.m**3).prefix("c"))  
    rx_perfil=((prop_perfil_tub[12]/1000 * fap.m).prefix("c"))
    ry_perfil=((prop_perfil_tub[13]/1000 * fap.m).prefix("c"))
    peso_perfil=((prop_perfil_tub[14]*9.807* fap.kg/fap.s**2).to("kgf_m"))

    return name_perfil,d_perfil,tw_perfil,b_perfil,tf_perfil,Ag_perfil,Ix_perfil,Iy_perfil,Sx_perfil,Sy_perfil,Zx_perfil,Zy_perfil,rx_perfil,ry_perfil,peso_perfil

#Función que muestra propiedaes geométricas de perfiles i usando handcalcs y forallpeople
@handcalc(override="params",precision=2,jupyter_display=True)
def display_geomet_perfil_i(prop_perfil_i):
    name_perfil=prop_perfil_i[0]
    d_perfil=((prop_perfil_i[1]/1000 * fap.m).prefix("m"))
    tw_perfil=((prop_perfil_i[2]/1000 * fap.m).prefix("m"))
    bf_perfil=((prop_perfil_i[3]/1000 * fap.m).prefix("m"))
    tf_perfil=((prop_perfil_i[4]/1000 * fap.m).prefix("m"))
    h_perfil=((prop_perfil_i[5]/1000 * fap.m).prefix("m"))    
    Ag_perfil=((prop_perfil_i[6]/1000000 * fap.m**2).prefix("c"))  
    Ix_perfil=((prop_perfil_i[7]/1000000000000 * fap.m**4).prefix("c"))  
    Iy_perfil=((prop_perfil_i[8]/1000000000000 * fap.m**4).prefix("c"))  
    Sx_perfil=((prop_perfil_i[9]/1000000000 * fap.m**3).prefix("c"))  
    Sy_perfil=((prop_perfil_i[10]/1000000000 * fap.m**3).prefix("c"))  
    Zx_perfil=((prop_perfil_i[11]/1000000000 * fap.m**3).prefix("c"))  
    Zy_perfil=((prop_perfil_i[12]/1000000000 * fap.m**3).prefix("c"))  
    rx_perfil=((prop_perfil_i[13]/1000 * fap.m).prefix("c"))
    ry_perfil=((prop_perfil_i[14]/1000 * fap.m).prefix("c"))
    J_perfil=((prop_perfil_i[15]/1000000000000 * fap.m**4).prefix("c"))
    Cw_perfil=((prop_perfil_i[16]/1000000000000000000 * fap.m**6).prefix("c"))

    peso_perfil=((prop_perfil_i[17]*9.807* fap.kg/fap.s**2).to("kgf_m"))

    return name_perfil,d_perfil,tw_perfil,bf_perfil,tf_perfil,h_perfil,Ag_perfil,Ix_perfil,Iy_perfil,Sx_perfil,Sy_perfil,Zx_perfil,Zy_perfil,rx_perfil,ry_perfil,J_perfil,Cw_perfil,peso_perfil

##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------FUNCIONES CÁLCULO RESISTENCIA----------------------------
##-----------------------------------------------------------------------------------------------------##
#Función resistencia a tensión por fluencia del area bruta 
@handcalc(override="long",precision=2,jupyter_display=True)
def resist_tens_fluencia(Ag,Fy):
    phi_tens=0.90 #Coeficiente reducción resistencia para tensión
    Tr = phi_tens*Ag*Fy #Resistencia fluencia tensión
    return Tr
##-----------------------------------------------------------------------------------------------------##
#Función resistencia a tensión esperada
@handcalc(override="long",precision=2,jupyter_display=True)
def resist_tens_esperada(Ag,Fy,Ry):
    Tre = Ry*Fy*Ag
    return Tre
##-----------------------------------------------------------------------------------------------------##
#Función resistencia a compresión 
@handcalc(override="long", precision=2,jupyter_display=True)
def resist_comp(K,L,r,Fy,Es,Ag):
    phi_comp=0.90 #Coeficiente reducción resistencia para compresión
    esbelt_glob = K*L/r #Esbeltez global del elemento a compresión
    Fe=(pi**2*Es)/(K*L/r)**2 #Esfuerzo crítico de pandeo elástico
    if Fe>=0.44*Fy:
        Fcr=(0.658**(Fy/Fe))*Fy
    elif Fe<0.44*Fy:
        Fcr = 0.877*Fe
    Cr = phi_comp*Fcr*Ag #Resistencia a compresión
    return Cr
##-----------------------------------------------------------------------------------------------------##
#Función resistencia a compresión esperada
@handcalc(override="long", precision=2,jupyter_display=True)
def resist_comp_esperada(K,L,r,Fy,Ry,Es,Ag):

    Cre_1 = Ry*Fy*Ag
    
    esbelt_glob = K*L/r #Esbeltez global del elemento a compresión
    Fe=(pi**2*Es)/(K*L/r)**2 #Esfuerzo crítico de pandeo elástico
    if Fe>=0.44*Fy*Ry:
        Fcre=(0.658**(Fy*Ry/Fe))*Fy*Ry
    elif Fe<0.44*Fy*Ry:
        Fcre = 0.877*Fe
    Cre_2 = 1.14*Fcre*Ag 
    Cre = min(Cre_1,Cre_2) #Resistencia esperada a compresión
    return Cre

#Función resistencia a flexión tubulares
@handcalc(override="long", precision=2,jupyter_display=True)
def resist_flex_tub_rect(prop_mat,prop_geom_tub,esb_perfil):
    Es=prop_mat[5]
    Fy=prop_mat[1]
    b_perfil=prop_geom_tub[3]
    tf_perfil=prop_geom_tub[4]
    d_perfil=prop_geom_tub[1]
    tw_perfil=prop_geom_tub[2]
    Sx_perfil=prop_geom_tub[8]
    Zx_perfil=prop_geom_tub[10]
    esb_patin=esb_perfil[0]
    esb_alma=esb_perfil[1]    
    bef_perfil = min(1.92*tf_perfil*sqrt(Es/Fy)*(1-0.38/(b_perfil/tf_perfil)*sqrt(Es/Fy)),b_perfil)
    Se_perfil = min(((bef_perfil*d_perfil**3)-((bef_perfil-(2*tw_perfil))*(d_perfil-(2*tf_perfil))**3))/(6*d_perfil)*mm**3,Sx_perfil)
    Mp_perfil = Fy*Zx_perfil #Momento plastifiación
    if esb_patin=="Compacto": Mn1 = Mp_perfil
    elif esb_patin=="No compacto": Mn1 = min(Mp_perfil-(Mp_perfil-Fy*Sx_perfil)*(3.57*b_perfil/tf_perfil*sqrt(Fy/Es)-4),Mp_perfil)
    elif esb_patin=="Esbelto": Mn1 = Fy*Se_perfil
    if esb_alma=="Compacto": Mn2 = Mp_perfil
    elif esb_alma=="No compacto": Mn2 = min(Mp_perfil-(Mp_perfil-Fy*Sx_perfil)*(0.305*d_perfil/tw_perfil*sqrt(Fy/Es)-0.738),Mp_perfil)
    Mn = min(Mp_perfil,Mn1,Mn2)
    phi_flex=0.90
    Mres = phi_flex*Mn
    return Mres




##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------FUNCIONES COMPACIDAD PERFILES----------------------------
##-----------------------------------------------------------------------------------------------------##
#Función para calcular la compacidad de perfiles tubulares rectangulares (Tabla F.2.2.4-1b NSR10)/(Table B4.1b AISC 360-16)
@handcalc(override="long", precision=2,jupyter_display=True)
def compacidad_tub_rect_nsr(prop_mat,prop_geom_tub):
    Es=prop_mat[5]
    Fy=prop_mat[1]
    b_perfil=prop_geom_tub[3]
    tf_dis=prop_geom_tub[4]
    d_perfil=prop_geom_tub[1]
    tw_dis=prop_geom_tub[2]
    b_int = b_perfil-3*tw_dis
    h_int = d_perfil-3*tf_dis
    lamb_p_patin = 1.12*sqrt(Es/Fy) #Esbeltez límite patines-compacto
    lamb_r_patin = 1.40*sqrt(Es/Fy) #Esbeltez límite patines-no compacto
    lamb_patin = b_int/tf_dis #Esbeltez patin
    if lamb_patin <= lamb_p_patin: esb_patin = "Compacto" #Clasificación esbeltez patín
    elif lamb_patin <= lamb_r_patin: esb_patin = "No compacto" #Clasificación esbeltez patín
    elif lamb_patin > lamb_r_patin: esb_patin = "Esbelto" #Clasificación esbeltez patín
    lamb_p_alma = 2.42*sqrt(Es/Fy) #Esbeltez límite alma-compacto
    lamb_r_alma = 5.70*sqrt(Es/Fy) #Esbeltez límite alma-no compacto
    lamb_alma = h_int/tw_dis #Esbeltez alma
    if lamb_alma <= lamb_p_alma: esb_alma = "Compacto" #Clasificación esbeltez patín
    elif lamb_alma <= lamb_r_alma: esb_alma = "No compacto" #Clasificación esbeltez patín
    elif lamb_alma > lamb_r_alma: esb_alma = "Esbelto" #Clasificación esbeltez patín
    esb_perfil=[esb_patin,esb_alma]
    return esb_perfil

