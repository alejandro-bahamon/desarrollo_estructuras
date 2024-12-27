##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------LIBRERIAS---------------------------------------------------
##-----------------------------------------------------------------------------------------------------##
import handcalcs
#import handcalcs.render
#from handcalcs.decorator import handcalc
import forallpeople
forallpeople.environment('alejandrov6',top_level="True")
import math
from math import pi,sqrt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------FUNCIONES PERFILES Y MATERIALES----------------------------
##-----------------------------------------------------------------------------------------------------##

#Función para obtener las propiedades mecánicas del material
#Se retornan los valores numéricos en unidades Mpa
#Argumentos: 
#    nombre_material (str):
#                            A36
#                            A572-G50
#                            A500-GB-CIRC
#                            A500-GB-RECT
#                            A500-GC-CIRC
#                            A500-GC-RECT
#                            A325
#                            A490
#                            A193-GB7-und2.5"
#                            A193-GB7-und4"
#                            SAE1020
#                            SAE1045
#Retorna:
#   prop_mat=[nombre_material,Fy,Fu,Ry,Rt,Es] (list[str])
def importar_propiedades_material(nombre_material):
    bd_material = pd.read_csv('bdmaterial.csv', sep=';') #Importa la base de datos
    material_select=bd_material[bd_material["Acero"]==nombre_material] #Filtra los datos del material
    nombre_material=str(material_select[("Acero")].iloc[0]) #Seleciona el nombre del perfil
    Fy=float(material_select[("Fy (Mpa)")].iloc[0]) #Seleciona el nombre del perfil
    Fu=float(material_select[("Fu (Mpa)")].iloc[0]) #Seleciona el nombre del perfil
    Ry=float(material_select[("Ry")].iloc[0]) #Seleciona el nombre del perfil
    Rt=float(material_select[("Rt")].iloc[0]) #Seleciona el nombre del perfil
    Es=float(material_select[("Es(Mpa)")].iloc[0]) #Modulo elasticidad
    prop_mat=[nombre_material,Fy,Fu,Ry,Rt,Es]
    
    return prop_mat
##-----------------------------------------------------------------------------------------------------##

#Función para obtener las propiedades geométricas de un perfil laminado tipo I
#Se retornan las unidades en mm,mmm2,mm3 y mm4 y el peso en kgf/m
#Argumentos: 
#    nombre_perfil (str):
#Retorna:
#   prop_perfil_i=[name_perfil,d_perfil,tw_perfil,bf_perfil,tf_perfil,h_perfil,Ag_perfil,Ix_perfil,Iy_perfil,Sx_perfil,Sy_perfil,Zx_perfil,Zy_perfil,rx_perfil,ry_perfil,J_perfil,Cw_perfil,peso_perfil]
def importar_propiedades_perfil_i(nombre_perfil):
    #Función para obtener las propiedades geométricas del perfil
    #Se retornan las unidades en mm,mmm2,mm3 y mm4 y el peso en kgf/m

    bd_material = pd.read_csv('bdperfilesi.csv', sep=';') #Importa la base de datos
    perfil_select=bd_material[bd_material["Referencia"]==nombre_perfil] #Filtra los datos del material
    name_perfil=str(perfil_select[("Referencia")].iloc[0]) #Seleciona el nombre del perfil
    d_perfil=float(perfil_select[("d(mm)")].iloc[0])
    tw_perfil=float(perfil_select[("tw(mm)")].iloc[0])
    bf_perfil=float(perfil_select[("bf(mm)")].iloc[0])
    tf_perfil=float(perfil_select[("tf(mm)")].iloc[0])
    h_perfil =float(perfil_select[("T(mm)")].iloc[0])

    Ag_perfil=float(perfil_select[("A(cm2)")].iloc[0])*100
    Ix_perfil=float(perfil_select[("Ixx(cm4)")].iloc[0])*10000
    Iy_perfil=float(perfil_select[("Iyy(cm4)")].iloc[0])*10000
    Sx_perfil=float(perfil_select[("Sx(cm3)")].iloc[0])*1000
    Sy_perfil=float(perfil_select[("Sy(cm3)")].iloc[0])*1000
    Zx_perfil=float(perfil_select[("Zx(cm3)")].iloc[0])*1000
    Zy_perfil=float(perfil_select[("Zy(cm3)")].iloc[0])*1000
    rx_perfil=float(perfil_select[("rx(cm)")].iloc[0])*10
    ry_perfil=float(perfil_select[("ry(cm)")].iloc[0])*10
    J_perfil=float(perfil_select[("J(cm4)")].iloc[0])*10000
    Cw_perfil=float(perfil_select[("Cw(cm6)")].iloc[0])*1000000

    peso_perfil =float(perfil_select[("Peso(Kg/m)")].iloc[0])
    
    prop_perfil_i=[name_perfil,d_perfil,tw_perfil,bf_perfil,tf_perfil,h_perfil,Ag_perfil,Ix_perfil,Iy_perfil,Sx_perfil,Sy_perfil,Zx_perfil,Zy_perfil,rx_perfil,ry_perfil,J_perfil,Cw_perfil,peso_perfil]

    return prop_perfil_i

#Función para obtener las propiedades geométricas de un perfil tubular
#Se retornan las unidades en mm,mmm2,mm3 y mm4 y el peso en kgf/m
#Argumentos: 
#   nombre_perfil (str):
#       Ejemplo:TUB.300X150X8
#Retorna:
#   prop_perfil_tub=[name_perfil,d_perfil,tw_perfil,b_perfil,tf_perfil,Ag_perfil,Ix_perfil,Iy_perfil,Sx_perfil,Sy_perfil,Zx_perfil,Zy_perfil,rx_perfil,ry_perfil,peso_perfil]
def importar_propiedades_perfil_tub_rect(nombre_perfil):
    #Función para obtener las propiedades geométricas del perfil
    #Se retornan las unidades en mm,mmm2,mm3 y mm4 y el peso en kgf/m

    bd_material = pd.read_csv('bd_tub_rect_cuad.csv', sep=';') #Importa la base de datos
    perfil_select=bd_material[bd_material["Seccion_tub"]==nombre_perfil] #Filtra los datos del material
    name_perfil=str(perfil_select[("Seccion_tub")].iloc[0]) #Seleciona el nombre del perfil

    d_perfil=float(perfil_select[("H (mm)")].iloc[0])
    tw_perfil=float(perfil_select[("tdis (mm)")].iloc[0])
    b_perfil=float(perfil_select[("B (mm)")].iloc[0])
    tf_perfil=float(perfil_select[("tdis (mm)")].iloc[0])

    Ag_perfil=float(perfil_select[("Ag (mm2)")].iloc[0])
    Ix_perfil=float(perfil_select[("Ix (mm4)")].iloc[0])
    Iy_perfil=float(perfil_select[("Iy (mm4)")].iloc[0])
    Sx_perfil=float(perfil_select[("Sx (mm3)")].iloc[0])
    Sy_perfil=float(perfil_select[("Sy (mm3)")].iloc[0])
    Zx_perfil=float(perfil_select[("Zx (mm3)")].iloc[0])
    Zy_perfil=float(perfil_select[("Zy (mm3)")].iloc[0])
    rx_perfil=float(perfil_select[("rx (mm)")].iloc[0])
    ry_perfil=float(perfil_select[("ry (mm)")].iloc[0])
    
    peso_perfil =float(perfil_select[("Peso (kg/m)")].iloc[0])

    prop_perfil_tub=[name_perfil,d_perfil,tw_perfil,b_perfil,tf_perfil,Ag_perfil,Ix_perfil,Iy_perfil,Sx_perfil,Sy_perfil,Zx_perfil,Zy_perfil,rx_perfil,ry_perfil,peso_perfil]
    
    return prop_perfil_tub


#Función para calcular las propiedades geométricas de un perfil "I" armado
#Argumentos: 
#   d_perfil (float)
#   tw_perfil (float)
#   bf_perfil (float)
#   tf_perfil (float)
#Retorna:
#   prop_perfil_i=[name_perfil,d_perfil,tw_perfil,bf_perfil,tf_perfil,h_perfil,Ag_perfil,Ix_perfil,Iy_perfil,Sx_perfil,Sy_perfil,Zx_perfil,Zy_perfil,rx_perfil,ry_perfil,J_perfil,Cw_perfil,peso_perfil]
def propiedades_geom_perfil_i_pref(d_perfil,tw_perfil,bf_perfil,tf_perfil):
    name_perfil="HS"+str(d_perfil)
    h_perfil=d_perfil-2*tf_perfil
    Ag_perfil=(h_perfil*tw_perfil)+2*(bf_perfil*tf_perfil)
    Ix_perfil=(1/12*bf_perfil*d_perfil**3)-2*(1/12*(bf_perfil-tw_perfil)/2*h_perfil**3)
    Iy_perfil=(1/12*h_perfil*tw_perfil**3)+2*(1/12*tf_perfil*bf_perfil**3)
    Sx_perfil=Ix_perfil/(d_perfil/2)
    Sy_perfil=Iy_perfil/(bf_perfil/2)
    Zx_perfil=(bf_perfil*tf_perfil*(d_perfil-tf_perfil))+(1/4*tw_perfil*(d_perfil-2*tf_perfil)**2*mm*2)
    Zy_perfil=((bf_perfil**2*tf_perfil)/2)+(1/4*tw_perfil**2*(d_perfil-2*tf_perfil)*mm**2)
    rx_perfil=math.sqrt(Ix_perfil/Ag_perfil)
    ry_perfil=math.sqrt(Iy_perfil/Ag_perfil)
    peso_perfil=Ag_perfil/1000000*7850
    
    prop_perfil_i=[name_perfil,d_perfil,tw_perfil,bf_perfil,tf_perfil,h_perfil,Ag_perfil,Ix_perfil,Iy_perfil,Sx_perfil,Sy_perfil,Zx_perfil,Zy_perfil,rx_perfil,ry_perfil,J_perfil,Cw_perfil,peso_perfil]

    return prop_perfil_i
    

##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------FUNCIONES CÁLCULO RESISTENCIA----------------------------
##-----------------------------------------------------------------------------------------------------##
#Función resistencia a tensión por fluencia del area bruta 

def resist_tens_fluencia(Ag,Fy):
    phi_tens=0.90 #Coeficiente reducción resistencia para tensión
    Tr = phi_tens*Ag*Fy #Resistencia fluencia tensión
    return Tr
##-----------------------------------------------------------------------------------------------------##
#Función resistencia a tensión esperada
def resist_tens_esperada(Ag,Fy,Ry):
    Tre = Ry*Fy*Ag
    return Tre
##-----------------------------------------------------------------------------------------------------##
#Función resistencia a compresión 
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
##-----------------------------------------------------------------------------------------------------##
#Función resistencia a flexión tubulares
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
#Función pandeo lateral torsional
def Pandeo_lat_tors_seccion_i(Es,Fy,d_perfil,tw_perfil,bf_perfil,tf_perfil,h_perfil,Iy_perfil,ry_perfil,Sx_perfil,Zx_perfil,J_perfil,Cw_perfil,Cb,Lb,esb_alma,esb_patin,resist_comp_esperada):
    
    Mp_perfil=resist_comp_esperada(Fy,Zx_perfil)

    if esb_alma == "Compacto":

        #Alma Compacta- Aletas Compactas-F.2.6.2
        #Alma Compacta- Aletas No Compactas o esbeltas-F.2.6.3
        rts = (sqrt(sqrt(Iy_perfil.value*Cw_perfil.value)/Sx_perfil.value)*m) #Radio efectivo de giro
        Lp = (1.76*ry_perfil*sqrt(Es/Fy)) #Longitud arriostramiento plastificación
        ho_perfil = d_perfil-tf_perfil  #Distancia entre centroides de aletas
        c_perfil=1 #Factor para perfil "I" doble simetría
        Lr = 1.95*rts*Es/(0.7*Fy)*sqrt((J_perfil*c_perfil)/(Sx_perfil*ho_perfil)+sqrt(((J_perfil*c_perfil)/(Sx_perfil*ho_perfil))**2+6.76*((0.7*Fy)/Es)**2)) #Máx long arriostram pandeo inelástico
        Fcr= (Cb*pi**2*Es)/(Lb/rts)**2*sqrt(1+0.078*(J_perfil*c_perfil)/(Sx_perfil*ho_perfil)*(Lb/rts)**2)
        
        if Lb<=Lp : Mn=Mp_perfil
        elif Lb<=Lr : Mn = min(Cb*(Mp_perfil-(Mp_perfil-0.70*Fy*Sx_perfil)*((Lb-Lp)/(Lr-Lp))),Mp_perfil)
        elif Lb>Lr : Mn = min(Fcr*Sx_perfil,Mp_perfil)

    elif esb_alma == "No compacto":

        b_fc = bf_perfil #Ancho aleta compresión
        t_fc = tf_perfil #Espesor aleta compresión
        hc = h_perfil
        aw = hc*tw_perfil/(b_fc*t_fc)
        rt = b_fc/(sqrt(12*(1+(1/6)*aw)))

        Mp_2 = min(Zx_perfil*Fy,1.6*Sx_perfil*Fy)

        #1
        Myc = Fy*Sx
        #2
        Iyc=1/12*tf_perfil*bf_perfil**3
        if Iyc/Iy_perfil<=0.23:
            J_perfil=0
        Fcr=(Cb*pi**2*Es)/(Lb/rt)**2*sqrt(1+0.078*(J_perfil)/(Sx_perfil*ho_perfil)*(Lb/rts)**2)
        #3        
        Fl = 0.7*Fy
        



        if Lb<=Lp : Mn=Mp_perfil
        elif Lb<=Lr : Mn = min(Cb*(Mp_perfil-(Mp_perfil-0.70*Fy*Sx_perfil)*((Lb-Lp)/(Lr-Lp))),Mp_perfil)
        elif Lb>Lr : Mn = min(Fcr*Sx_perfil,Mp_perfil)



    #Alma No Compacta-F.2.6.4


    #Alma Esbelta-F.2.6.5


    esb_alma = "Compacto"
    esb_alma = "No compacto"
    esb_alma = "Esbelto"

    esb_patin = "Compacto"
    esb_patin = "No compacto"
    esb_patin = "Esbelto"




    return

##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------FUNCIONES COMPACIDAD PERFILES----------------------------
##-----------------------------------------------------------------------------------------------------##

##NSR-10 TABLA F.2.2.4-1b

def compacidad_i_nsr(Es,Fy,bf_perfil,tf_perfil,h_perfil,tw_perfil):

    #Esbeltez patines
    lamb_p_patin = 0.35*sqrt(Es/Fy) #Esbeltez límite patines-compacto
    lamb_r_patin = 1.0*sqrt(Es/Fy) #Esbeltez límite patines-no compacto

    b_patin = bf_perfil/2
    lamb_patin = b_patin/tf_perfil #Esbeltez patin

    if lamb_patin <= lamb_p_patin: esb_patin = "Compacto" #Clasificación esbeltez patín
    elif lamb_patin <= lamb_r_patin: esb_patin = "No compacto" #Clasificación esbeltez patín
    elif lamb_patin > lamb_r_patin: esb_patin = "Esbelto" #Clasificación esbeltez patín

    #Esbeltez alma
    lamb_p_alma = 3.76*sqrt(Es/Fy) #Esbeltez límite alma-compacto
    lamb_r_alma = 5.70*sqrt(Es/Fy) #Esbeltez límite alma-no compacto

    lamb_alma = h_perfil/tw_perfil #Esbeltez alma

    if lamb_alma <= lamb_p_alma: esb_alma = "Compacto" #Clasificación esbeltez patín
    elif lamb_alma <= lamb_r_alma: esb_alma = "No compacto" #Clasificación esbeltez patín
    elif lamb_alma > lamb_r_alma: esb_alma = "Esbelto" #Clasificación esbeltez patín

    return esb_patin,esb_alma

#Función para calcular la compacidad de perfiles tubulares rectangulares (Tabla F.2.2.4-1b NSR10)/(Table B4.1b AISC 360-16)
#Argumentos: 
#   prop_mat(list[float])
#       prop_mat=[nombre_material,Fy,Fu,Ry,Rt,Es]
#   prop_geom_tub(list[float])
#       prop_perfil_tub=[name_perfil,d_perfil,tw_perfil,b_perfil,tf_perfil,Ag_perfil,Ix_perfil,Iy_perfil,Sx_perfil,Sy_perfil,Zx_perfil,Zy_perfil,rx_perfil,ry_perfil,peso_perfil]
#Retorna:
#   prop_mat=esb_patin,esb_alma (str)
def compacidad_tub_rect_nsr(prop_mat,prop_geom_tub):
    #Extrae las propiedades necesarias de las listas de argumentos
    Es=prop_mat[5]
    Fy=prop_mat[1]
    b_perfil=prop_geom_tub[3]
    tf_dis=prop_geom_tub[4]
    d_perfil=prop_geom_tub[1]
    tw_dis=prop_geom_tub[2]

    #Calcula el espesor de diseño y dimensiones internas según AISC360-16 B4.1b-d
    b_int=b_perfil-3*tw_dis
    h_int=d_perfil-3*tf_dis

    #Esbeltez patines
    lamb_p_patin = 1.12*sqrt(Es/Fy) #Esbeltez límite patines-compacto
    lamb_r_patin = 1.40*sqrt(Es/Fy) #Esbeltez límite patines-no compacto

    lamb_patin = b_int/tf_dis #Esbeltez patin

    if lamb_patin <= lamb_p_patin: esb_patin = "Compacto" #Clasificación esbeltez patín
    elif lamb_patin <= lamb_r_patin: esb_patin = "No compacto" #Clasificación esbeltez patín
    elif lamb_patin > lamb_r_patin: esb_patin = "Esbelto" #Clasificación esbeltez patín

    #Esbeltez alma
    lamb_p_alma = 2.42*sqrt(Es/Fy) #Esbeltez límite alma-compacto
    lamb_r_alma = 5.70*sqrt(Es/Fy) #Esbeltez límite alma-no compacto

    lamb_alma = h_int/tw_dis #Esbeltez alma

    if lamb_alma <= lamb_p_alma: esb_alma = "Compacto" #Clasificación esbeltez patín
    elif lamb_alma <= lamb_r_alma: esb_alma = "No compacto" #Clasificación esbeltez patín
    elif lamb_alma > lamb_r_alma: esb_alma = "Esbelto" #Clasificación esbeltez patín

    return esb_patin,esb_alma

##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------FUNCIONES ESBELTEZ PERFILES----------------------------
##-----------------------------------------------------------------------------------------------------##
#Función para calcular las esbeltez local de riostras tipo I según NSR-10
##NSR-10 TABLA F.3.4-1
def riostras_i_esbeltez_nsr(Es,Fy,h_perfil,tw_perfil,bf_perfil,tf_perfil):
    
    lamb_da_patin = 0.30*sqrt(Es/Fy)#Límite esbeltez patín ductilidad alta
    lamb_dm_patin = 0.38*sqrt(Es/Fy)#Límite esbeltez patín ductilidad moderada
    
    lamb_patin = bf_perfil/(2*tf_perfil) #Esbeltez del patín

    if lamb_patin <= lamb_da_patin: duct_patin = "Ductilidad alta"
    elif lamb_patin <= lamb_dm_patin: duct_patin = "Ductilidad moderada"
    elif lamb_patin > lamb_dm_patin: duct_patin = "Ductilidad baja"

    lamb_da_alma = 1.49*sqrt(Es/Fy)#Límite esbeltez alma ductilidad alta
    lamb_dm_alma = 1.49*sqrt(Es/Fy)#Límite esbeltez alma ductilidad moderada

    lamb_alma = h_perfil/tw_perfil #Esbeltez del alma

    if lamb_alma <= lamb_da_alma: duct_alma = "Ductilidad alta"
    elif lamb_alma <= lamb_dm_alma: duct_alma = "Ductilidad moderada"
    elif lamb_alma > lamb_dm_alma: duct_alma = "Ductilidad baja"

    return duct_patin,duct_alma


##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------FUNCIONES PÓRTICOS ARRIOSTRADOS PAC-DES--------------------
##-----------------------------------------------------------------------------------------------------##
def pac_des_x(Lport,n_pisos,alturas_piso,resistencias_esperadas,cargas_axiales_muertas,cargas_axiales_vivas):
    #DATOS DE ENTRADA
        #1.Geometría_pórtico
            #Lport=7 #Luz del pórtico
            #n_pisos=4 #Número de niveles en los que inician o finalizan las riostras, se inicia la cuenta desde 0
        #2.Diccionario de alturas de piso
            #alturas_piso={
            #"Hp1":4.00, #Altura del piso 1 medida desde el nivel 0 al nivel 1
            #"Hp2":3.85, #Altura del piso 2 medida desde el nivel 1 al nivel 2
            #"Hp3":3.175, #Altura del piso 2 medida desde el nivel 2 al nivel 3
            #"Hp4":3.175, #Altura del piso 2 medida desde el nivel 3 al nivel 4
            #}
        #3.Diccionario con resistencia esperada a tensión y resistencia esperada a compresión de las riostras en cada nivel
            #resistencias_esperadas={
            #"Tre_p1":2396.25, #Resistencia esperada a tensión de riostras del piso 1
            #"Cre_p1":1858.9, #Resistencia esperada a compresión de riostras del piso 1
            #"Tre_p2":2396.25, #Resistencia esperada a tensión de riostras del piso 2
            #"Cre_p2":1858.9, #Resistencia esperada a compresión de riostras del piso 2
            #"Tre_p3":2396.25,#Resistencia esperada a tensión de riostras del piso 3
            #"Cre_p3":1858.9, #Resistencia esperada a compresión de riostras del piso 3
            #"Tre_p4":2396.25,#Resistencia esperada a tensión de riostras del piso 4
            #"Cre_p4":1858.9 #Resistencia esperada a compresión de riostras del piso 4
        #4.Diccionario con cargas axiales por carga muerta sin mayorar en cada piso de la columna
            #cargas_axiales_muertas={
            #"Pcm_p1":2396.25, #Resistencia esperada a tensión de riostras del piso 1
            #"Pcm_p2":2396.25, #Resistencia esperada a tensión de riostras del piso 2
            #"Pcm_p3":2396.25,#Resistencia esperada a tensión de riostras del piso 3
            #"Pcm_p4":2396.25,#Resistencia esperada a tensión de riostras del piso 4
            #}
        #5.Diccionario con cargas axiales vivas muertas sin mayorar en cada piso de la columna
            #cargas_axiales_vivas={
            #"Pcv_p1":2396.25, #Resistencia esperada a tensión de riostras del piso 1
            #"Pcv_p2":2396.25, #Resistencia esperada a tensión de riostras del piso 2
            #"Pcv_p3":2396.25,#Resistencia esperada a tensión de riostras del piso 3
            #"Pcv_p4":2396.25,#Resistencia esperada a tensión de riostras del piso 4
            #}

    #CREA DATAFRAME CON DATOS DE ENTRADA
    vector_nombre_piso=[] #vector para guardar los nombres de los pisos
    vector_h_piso=[0] #vector para guardar las alturas de cada piso
    vector_Tre=[] #vector para guardar las resistencias esperadas a tensión de todos los pisos
    vector_Cre=[] #vector para guardar las resistencias esperadas a compresión de todos los pisos
    vector_angulo_grados=[] #Vector para guerdar el ángulo en grados que que forman las riostras con la horizontal
    vector_angulo_rad=[] #Vector para guerdar el ángulo en radianes que que forman las riostras con la horizontal
    vector_cm_piso=[0] #Vector para guardar las cargas axiales por carga muerta en las columnas
    vector_cv_piso=[0] #Vector para guardar las cargas axiales por carga viva en las columnas
    
    for i in range (0,n_pisos+1): #Ciclo que alimenta el vector de nombres de piso
        nombre_piso="piso"+str(i)
        vector_nombre_piso.append(nombre_piso)
    
    for nombre, valor in alturas_piso.items(): #Ciclo que alimenta el vector de alturas de piso
        vector_h_piso.append(valor)
    for nombre, valor in cargas_axiales_muertas.items(): #Ciclo que alimenta el vector de cargfas mnuertas
         vector_cm_piso.append(-1*valor)
     #vector_cm_piso.append(0)
    for nombre, valor in cargas_axiales_vivas.items(): #Ciclo que alimenta el vector de cargfas mnuertas
         vector_cv_piso.append(-1*valor)
    
    for nombre, valor in resistencias_esperadas.items():#Ciclo que alimenta los vectores de resistencias esperadas a compresión y tensión
        if nombre.startswith('Tre'):
            vector_Tre.append(valor)
        elif nombre.startswith('Cre'):
            vector_Cre.append(valor)
    vector_Tre.append(0)
    vector_Cre.append(0)
    
    for i in range (0,n_pisos): #Vector que calcula calcula los ángulos de las riostras y alimenta los vectores de ángulos
        vector_angulo_rad.append(math.atan((vector_h_piso[i+1]) / (Lport/2)))
        vector_angulo_grados.append(math.atan((vector_h_piso[i+1]) / (Lport/2)) * 180 / math.pi)
    vector_angulo_rad.append(0)
    vector_angulo_grados.append(0)
    
    #Crea el dataframe a partir de los vectores definidos
    df_pac_capacidad = pd.DataFrame({
        "Piso": vector_nombre_piso,
        "Altura (m)":vector_h_piso,
        "Ango_grados":vector_angulo_grados,
        "Ang_rad":vector_angulo_rad,
        "Tre":vector_Tre,
        "Cre":vector_Cre})
    df_pac_capacidad["0.3Cre"]=df_pac_capacidad["Cre"]*0.3
    
    #Calcula el vector de cargas máximas de compresión en la columna
    vector_comp_col=[]
    for num, valor in df_pac_capacidad["Piso"].items():
        if num == 0:
            vector_comp_col.append(-1*df_pac_capacidad["Cre"][0]*math.sin(df_pac_capacidad["Ang_rad"][0]))
        elif num % 2 == 0:
            vector_comp_col.append(-1*df_pac_capacidad["Tre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1])-
                                        df_pac_capacidad["Cre"][num]*math.sin(df_pac_capacidad["Ang_rad"][num]))
        else:
            vector_comp_col.append((-1*df_pac_capacidad["Tre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1])+
                                         df_pac_capacidad["Cre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1])-
                                         df_pac_capacidad["Cre"][num]*math.sin(df_pac_capacidad["Ang_rad"][num])+
                                         df_pac_capacidad["Tre"][num]*math.sin(df_pac_capacidad["Ang_rad"][num]))/2)
    
    df_pac_capacidad["comp_col_sismo"]=vector_comp_col
    
    #Calcula el vector de cargas máximas de tensión en la columna
    vector_ten_col=[]
    for num, valor in df_pac_capacidad["Piso"].items():
        if num == 0:
            vector_ten_col.append(df_pac_capacidad["Tre"][0]*math.sin(df_pac_capacidad["Ang_rad"][0]))
        elif num % 2 == 0:
            vector_ten_col.append(df_pac_capacidad["Cre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1])+
                                        df_pac_capacidad["Tre"][num]*math.sin(df_pac_capacidad["Ang_rad"][num]))
        else:
            vector_ten_col.append((-1*df_pac_capacidad["Tre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1])+
                                         df_pac_capacidad["Cre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1])-
                                         df_pac_capacidad["Cre"][num]*math.sin(df_pac_capacidad["Ang_rad"][num])+
                                         df_pac_capacidad["Tre"][num]*math.sin(df_pac_capacidad["Ang_rad"][num]))/2)
    df_pac_capacidad["tens_col_sismo"]=vector_ten_col
    
    #Crea la suma acumulada de la carga a compresión
    df_pac_capacidad["comp_col_acum_sismo"]=df_pac_capacidad.loc[::-1, 'comp_col_sismo'].cumsum()[::-1]
    
    #crea la suma acumulada de la carga en tensión
    df_pac_capacidad["tens_col_acum_sismo"]=df_pac_capacidad.loc[::-1, 'tens_col_sismo'].cumsum()[::-1]

    df_pac_capacidad["comp_cm"]=vector_cm_piso
    df_pac_capacidad["comp_cv"]=vector_cv_piso
    df_pac_capacidad["1.2*cm+0.5cv+1.0Eq_comp"]=\
 df_pac_capacidad["comp_cm"]*1.2+df_pac_capacidad["comp_cv"]*0.5+df_pac_capacidad["comp_col_acum_sismo"]
    df_pac_capacidad["0.9*cm+1.0Eq_tens"]=df_pac_capacidad["comp_cm"]*0.9+df_pac_capacidad["tens_col_acum_sismo"]
    
        
    df_pac_capacidad=df_pac_capacidad.round(2)

    return df_pac_capacidad
##-----------------------------------------------------------------------------------------------------##
def pac_des_chevron(Lport,n_pisos,alturas_piso,resistencias_esperadas,cargas_axiales_muertas,cargas_axiales_vivas):
 
    #DATOS DE ENTRADA
        #1.Geometría_pórtico
            #Lport=7 #Luz del pórtico
            #n_pisos=4 #Número de niveles en los que inician o finalizan las riostras, se inicia la cuenta desde 0
        #2.Diccionario de alturas de piso
            #alturas_piso={
            #"Hp1":4.00, #Altura del piso 1 medida desde el nivel 0 al nivel 1
            #"Hp2":3.85, #Altura del piso 2 medida desde el nivel 1 al nivel 2
            #"Hp3":3.175, #Altura del piso 2 medida desde el nivel 2 al nivel 3
            #"Hp4":3.175, #Altura del piso 2 medida desde el nivel 3 al nivel 4
            #}
        #3.Diccionario con resistencia esperada a tensión y resistencia esperada a compresión de las riostras en cada nivel
            #resistencias_esperadas={
            #"Tre_p1":2396.25, #Resistencia esperada a tensión de riostras del piso 1
            #"Cre_p1":1858.9, #Resistencia esperada a compresión de riostras del piso 1
            #"Tre_p2":2396.25, #Resistencia esperada a tensión de riostras del piso 2
            #"Cre_p2":1858.9, #Resistencia esperada a compresión de riostras del piso 2
            #"Tre_p3":2396.25,#Resistencia esperada a tensión de riostras del piso 3
            #"Cre_p3":1858.9, #Resistencia esperada a compresión de riostras del piso 3
            #"Tre_p4":2396.25,#Resistencia esperada a tensión de riostras del piso 4
            #"Cre_p4":1858.9 #Resistencia esperada a compresión de riostras del piso 4
            #}
        #4.Diccionario con cargas axiales por carga muerta sin mayorar en cada piso de la columna
            #cargas_axiales_muertas={
            #"Pcm_p1":2396.25, #Resistencia esperada a tensión de riostras del piso 1
            #"Pcm_p2":2396.25, #Resistencia esperada a tensión de riostras del piso 2
            #"Pcm_p3":2396.25,#Resistencia esperada a tensión de riostras del piso 3
            #"Pcm_p4":2396.25,#Resistencia esperada a tensión de riostras del piso 4
            #}
        #5.Diccionario con cargas axiales vivas muertas sin mayorar en cada piso de la columna
            #cargas_axiales_vivas={
            #"Pcv_p1":2396.25, #Resistencia esperada a tensión de riostras del piso 1
            #"Pcv_p2":2396.25, #Resistencia esperada a tensión de riostras del piso 2
            #"Pcv_p3":2396.25,#Resistencia esperada a tensión de riostras del piso 3
            #"Pcv_p4":2396.25,#Resistencia esperada a tensión de riostras del piso 4
            #}

    #CREA DATAFRAME CON DATOS DE ENTRADA
     vector_nombre_piso=[] #vector para guardar los nombres de los pisos
     vector_h_piso=[0] #vector para guardar las alturas de cada piso
     vector_Tre=[] #vector para guardar las resistencias esperadas a tensión de todos los pisos
     vector_Cre=[] #vector para guardar las resistencias esperadas a compresión de todos los pisos
     vector_angulo_grados=[] #Vector para guerdar el ángulo en grados que que forman las riostras con la horizontal
     vector_angulo_rad=[] #Vector para guerdar el ángulo en radianes que que forman las riostras con la horizontal
     vector_cm_piso=[0] #Vector para guardar las cargas axiales por carga muerta en las columnas
     vector_cv_piso=[0] #Vector para guardar las cargas axiales por carga viva en las columnas
     vector_cargas_vigas=[0] #Vector para guardar las cargas axiales por carga viva en las columnas
     vector_momentos_vigas=[0] #Vector para guardar las cargas axiales por carga viva en las columnas
     
     for i in range (0,n_pisos+1): #Ciclo que alimenta el vector de nombres de piso
         nombre_piso="piso"+str(i)
         vector_nombre_piso.append(nombre_piso)
     
     for nombre, valor in alturas_piso.items(): #Ciclo que alimenta el vector de alturas de piso
         vector_h_piso.append(valor)
     
     for nombre, valor in cargas_axiales_muertas.items(): #Ciclo que alimenta el vector de cargfas mnuertas
         vector_cm_piso.append(-1*valor)
     #vector_cm_piso.append(0)
     for nombre, valor in cargas_axiales_vivas.items(): #Ciclo que alimenta el vector de cargfas mnuertas
         vector_cv_piso.append(-1*valor)
     #vector_cv_piso.append(0)    
     for nombre, valor in resistencias_esperadas.items():#Ciclo que alimenta los vectores de resistencias esperadas a compresión y tensión
         if nombre.startswith('Tre'):
             vector_Tre.append(valor)
         elif nombre.startswith('Cre'):
             vector_Cre.append(valor)
     vector_Tre.append(0)
     vector_Cre.append(0)
     
     for i in range (0,n_pisos): #Vector que calcula calcula los ángulos de las riostras y alimenta los vectores de ángulos
         vector_angulo_rad.append(math.atan((vector_h_piso[i+1]) / (Lport/2)))
         vector_angulo_grados.append(math.atan((vector_h_piso[i+1]) / (Lport/2)) * 180 / math.pi)
     vector_angulo_rad.append(0)
     vector_angulo_grados.append(0)
     
     #Crea el dataframe a partir de los vectores definidos
     df_pac_capacidad = pd.DataFrame({
         "Piso": vector_nombre_piso,
         "Altura (m)":vector_h_piso,
         "Ango_grados":vector_angulo_grados,
         "Ang_rad":vector_angulo_rad,
         "Tre":vector_Tre,
         "Cre":vector_Cre})
     df_pac_capacidad["0.3Cre"]=df_pac_capacidad["Cre"]*0.3
     
     #----OPCIÓN 1 Tre y Cre----
     
     #Calcula el vector de cargas máximas de compresión en la columna
     vector_comp_col_1=[]
     vector_ten_col_1=[]
     for num, valor in df_pac_capacidad["Piso"].items():
         if num == 0:
             vector_comp_col_1.append(-1*df_pac_capacidad["Cre"][0]*math.sin(df_pac_capacidad["Ang_rad"][0]))
             vector_ten_col_1.append(df_pac_capacidad["Tre"][0]*math.sin(df_pac_capacidad["Ang_rad"][0]))
         
         else:
             vector_comp_col_1.append(-1*df_pac_capacidad["Cre"][num]*math.sin(df_pac_capacidad["Ang_rad"][num])+
                                       (df_pac_capacidad["Cre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1])-
                                       df_pac_capacidad["Tre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1]))/2)
     
             vector_ten_col_1.append(df_pac_capacidad["Tre"][num]*math.sin(df_pac_capacidad["Ang_rad"][num])+
                                       (df_pac_capacidad["Cre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1])-
                                       df_pac_capacidad["Tre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1]))/2)
     
     df_pac_capacidad["comp_col_1_sismo"]=vector_comp_col_1
     df_pac_capacidad["tens_col_1_sismo"]=vector_ten_col_1
     
     
     #----OPCIÓN 2 Tre y 0.3*Cre----
     
     #Calcula el vector de cargas máximas de compresión en la columna
     vector_comp_col_2=[]
     vector_ten_col_2=[]
     for num, valor in df_pac_capacidad["Piso"].items():
         if num == 0:
             vector_comp_col_2.append(-1*df_pac_capacidad["0.3Cre"][0]*math.sin(df_pac_capacidad["Ang_rad"][0]))
             vector_ten_col_2.append(df_pac_capacidad["Tre"][0]*math.sin(df_pac_capacidad["Ang_rad"][0]))
         
         else:
             vector_comp_col_2.append(-1*df_pac_capacidad["0.3Cre"][num]*math.sin(df_pac_capacidad["Ang_rad"][num])+
                                       (df_pac_capacidad["0.3Cre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1])-
                                       df_pac_capacidad["Tre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1]))/2)
     
             vector_ten_col_2.append(df_pac_capacidad["Tre"][num]*math.sin(df_pac_capacidad["Ang_rad"][num])+
                                       (df_pac_capacidad["0.3Cre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1])-
                                       df_pac_capacidad["Tre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1]))/2)
     
     df_pac_capacidad["comp_col_2_sismo"]=vector_comp_col_2
     df_pac_capacidad["tens_col_2_sismo"]=vector_ten_col_2
         
     if df_pac_capacidad["comp_col_1_sismo"].sum()<df_pac_capacidad["comp_col_2_sismo"].sum():
         #Crea la suma acumulada de la carga a compresión
         #df_pac_capacidad["comp_col_acum_sismo"]=df_pac_capacidad["comp_col_1_sismo"].cumsum()
         df_pac_capacidad["comp_col_acum_sismo"]=df_pac_capacidad.loc[::-1, 'comp_col_1_sismo'].cumsum()[::-1]
     else:
         #Crea la suma acumulada de la carga a compresión
         #df_pac_capacidad["comp_col_acum_sismo"]=df_pac_capacidad["comp_col_2_sismo"].cumsum()
         df_pac_capacidad["comp_col_acum_sismo"]=df_pac_capacidad.loc[::-1, 'comp_col_2_sismo'].cumsum()[::-1]
         
     if df_pac_capacidad["tens_col_1_sismo"].sum()>df_pac_capacidad["tens_col_2_sismo"].sum():
         #crea la suma acumulada de la carga en tensión
         #df_pac_capacidad["tens_col_acum_sismo"]=df_pac_capacidad["tens_col_1_sismo"].cumsum()
         df_pac_capacidad["tens_col_acum_sismo"]=df_pac_capacidad.loc[::-1, 'tens_col_1_sismo'].cumsum()[::-1]      
     else:
         #crea la suma acumulada de la carga en tensión
         #df_pac_capacidad["tens_col_acum_sismo"]=df_pac_capacidad["tens_col_2_sismo"].cumsum()
         df_pac_capacidad["tens_col_acum_sismo"]=df_pac_capacidad.loc[::-1, 'tens_col_2_sismo'].cumsum()[::-1]       
     
     
     df_pac_capacidad["comp_cm"]=vector_cm_piso
     df_pac_capacidad["comp_cv"]=vector_cv_piso
     df_pac_capacidad["1.2*cm+0.5cv+1.0Eq_comp"]=\
df_pac_capacidad["comp_cm"]*1.2+df_pac_capacidad["comp_cv"]*0.5+df_pac_capacidad["comp_col_acum_sismo"]
     df_pac_capacidad["0.9*cm+1.0Eq_tens"]=df_pac_capacidad["comp_cm"]*0.9+df_pac_capacidad["tens_col_acum_sismo"]


     for num, valor in df_pac_capacidad["Piso"].items():
         if num != 0:
             vector_cargas_vigas.append((-1*df_pac_capacidad["Tre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1]))+\
                                        (df_pac_capacidad["0.3Cre"][num-1]*math.sin(df_pac_capacidad["Ang_rad"][num-1])))

             vector_momentos_vigas.append(vector_cargas_vigas[num]*Lport/4)
          

     df_pac_capacidad["carga_puntual_viga"]=vector_cargas_vigas    
     df_pac_capacidad["momento_max_viga"]=vector_momentos_vigas


 
     df_pac_capacidad=df_pac_capacidad.round(2)
     return df_pac_capacidad