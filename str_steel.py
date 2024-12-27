"""
str_steel: Structural Steel Calculations Module


This module provides functions to perform structural calculations for steel elements.

Functions available:
- prop_mat:
    Retrieve the mechanical properties of a specified material.
- prop_i:
    Retrieve geometric properties of an I-shaped rolled section.
- prop_hss_rect:
    Retrieve geometric properties of a rectangular tubular section.
- prop_i_bu:
    Calculate the geometric properties of a custom I-section (assembled).
-stren_flex_yield
-stren_flex_c_c_ltb
-stren_flex_ncs_c_flb
-stren_flex_cncs_nc_cfy
-slend_i_comp_ns
-slend_i_flex_ns
-slend_i_s
"""

##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------LIBRARIES---------------------------------------------------
##-----------------------------------------------------------------------------------------------------##
# Standard library imports
import os  # For file validation

# Third-party library imports
import math
from math import pi, sqrt  # Specific math functions
import numpy as np  # For numerical operations
import pandas as pd  # For handling CSV data
import matplotlib.pyplot as plt  # For plotting
import handcalcs  # For hand calculations
import forallpeople  # For physical units and engineering calculations

# Initialize forallpeople environment
forallpeople.environment('alejandrov6', top_level="True")

##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------MATERIAL & SECTION PROPERTIES----------------------------
##-----------------------------------------------------------------------------------------------------##

def prop_mat(mat_name):
    """
    Retrieve the mechanical properties of a specified material.

    This function retrieves the mechanical properties of a specified material 
    from a database stored in a CSV file. All properties are returned as numerical
    values in MPa units.

    Arguments:
        mat_name (str): The name of the material. Accepted values are:
            - A36
            - A572-G50
            - A500-GB-CIRC
            - A500-GB-RECT
            - A500-GC-CIRC
            - A500-GC-RECT
            - A325
            - A490
            - A193-GB7-und2.5"
            - A193-GB7-und4"
            - SAE1020
            - SAE1045

    Returns:
        dict: A dictionary containing the following mechanical properties:
            - mat_name (str)
            - Fy (float)
            - Fu (float)
            - Ry (float)
            - Rt (float)
            - Es (float)

    Raises:
        ValueError: If the material is not found in the database.
        FileNotFoundError: If the database file is missing.
    """
 # Define the path to the database file
    current_dir = os.path.dirname(__file__)
    db_path = os.path.join(current_dir, 'bdmaterial.csv')

    # Check if the database file exists
    if not os.path.exists(db_path):
        # Raise an error if the file is missing
        raise FileNotFoundError(f"The database file '{db_path}' was not found.")

    # Read the database into a pandas DataFrame
    bd_material = pd.read_csv(db_path, sep=';')

    # Filter the DataFrame to select rows matching the specified material name
    material_select = bd_material[bd_material["Acero"] == mat_name]

    # Check if the filtered DataFrame is empty (material not found)
    if material_select.empty:
        # Raise an error if the material is not found in the database
        raise ValueError(f"Material '{mat_name}' not found in the database.")

    # Extract material properties
    mat_name = str(material_select["Acero"].iloc[0])
    Fy = float(material_select["Fy (Mpa)"].iloc[0])
    Fu = float(material_select["Fu (Mpa)"].iloc[0])
    Ry = float(material_select["Ry"].iloc[0])
    Rt = float(material_select["Rt"].iloc[0])
    Es = float(material_select["Es(Mpa)"].iloc[0])

    # Group all properties into a dictionary
    prop_mat = {
        "mat_name": mat_name,
        "Fy": Fy,
        "Fu": Fu,
        "Ry": Ry,
        "Rt": Rt,
        "Es": Es}   
    # Return the dictionary of properties
    return prop_mat

def prop_i(sect_name):
    """
    Retrieve geometric properties of an I-shaped rolled section.

    This function fetches geometric properties from a database stored in a CSV file. 
    Units are returned in mm, mm², mm³, mm⁴, and weight in kgf/m.

    Arguments:
        sect_name (str): The reference name of the I-shaped section. Some accepted values are:
            - IPE450
            - HE550A
            - W30x108

    Returns:
        dict: A dictionary containing the following geometric properties:
            - sect_name (str)
            - d_sect (float)
            - tw_sect (float)
            - bf_sect (float)
            - tf_sect (float)
            - h_sect (float)
            - Ag_sect (float)
            - Ix_sect (float)
            - Iy_sect (float)
            - Sx_sect (float)
            - Sy_sect (float)
            - Zx_sect (float)
            - Zy_sect (float)
            - rx_sect (float)
            - ry_sect (float)
            - J_sect (float)
            - Cw_sect (float)
            - weight_sect (float)

    Raises:
        ValueError: If the section is not found in the database.
        FileNotFoundError: If the database file is missing.
    """

    # Define the path to the database file
    current_dir = os.path.dirname(__file__)
    db_path = os.path.join(current_dir, 'bdperfilesi.csv')

    # Check if the database file exists
    if not os.path.exists(db_path):
        # Raise an error if the file is missing
        raise FileNotFoundError(f"The database file '{db_path}' was not found.")

    # Read the database into a pandas DataFrame
    bd_sections = pd.read_csv(db_path, sep=';')

    # Filter the DataFrame to select rows matching the specified section name
    sect_select = bd_sections[bd_sections["Referencia"] == sect_name]

    # Check if the filtered DataFrame is empty (section not found)
    if sect_select.empty:
        # Raise an error if the section is not found in the database
        raise ValueError(f"section '{sect_name}' not found in the database.")

    # Extract geometric properties
    sect_name = str(sect_select["Referencia"].iloc[0])
    d_sect = float(sect_select["d(mm)"].iloc[0])
    tw_sect = float(sect_select["tw(mm)"].iloc[0])
    bf_sect = float(sect_select["bf(mm)"].iloc[0])
    tf_sect = float(sect_select["tf(mm)"].iloc[0])
    h_sect = float(sect_select["T(mm)"].iloc[0])

    Ag_sect = float(sect_select["A(cm2)"].iloc[0]) * 100
    Ix_sect = float(sect_select["Ixx(cm4)"].iloc[0]) * 10000
    Iy_sect = float(sect_select["Iyy(cm4)"].iloc[0]) * 10000
    Sx_sect = float(sect_select["Sx(cm3)"].iloc[0]) * 1000
    Sy_sect = float(sect_select["Sy(cm3)"].iloc[0]) * 1000
    Zx_sect = float(sect_select["Zx(cm3)"].iloc[0]) * 1000
    Zy_sect = float(sect_select["Zy(cm3)"].iloc[0]) * 1000
    rx_sect = float(sect_select["rx(cm)"].iloc[0]) * 10
    ry_sect = float(sect_select["ry(cm)"].iloc[0]) * 10
    J_sect = float(sect_select["J(cm4)"].iloc[0]) * 10000
    Cw_sect = float(sect_select["Cw(cm6)"].iloc[0]) * 1000000
    weight_sect = float(sect_select["Peso(Kg/m)"].iloc[0])
    sect_type = "Rolled"

    # Group all properties into a dictionary
    prop_sect_i = {
        "sect_name": sect_name, 
        "d_sect": d_sect,
        "tw_sect": tw_sect,
        "bf_sect": bf_sect,
        "tf_sect": tf_sect, 
        "h_sect": h_sect,
        "Ag_sect": Ag_sect, 
        "Ix_sect": Ix_sect, 
        "Iy_sect": Iy_sect, 
        "Sx_sect": Sx_sect, 
        "Sy_sect": Sy_sect, 
        "Zx_sect": Zx_sect, 
        "Zy_sect": Zy_sect, 
        "rx_sect": rx_sect, 
        "ry_sect": ry_sect, 
        "J_sect": J_sect, 
        "Cw_sect": Cw_sect, 
        "weight_sect": weight_sect,
        "sect_type": sect_type
    }

    # Return the dictionary of properties
    return prop_sect_i

def prop_hss_rect(sect_name):
    """
    Retrieve geometric properties of a rectangular tubular section.

    This function fetches geometric properties from a database stored in a CSV file. 
    Units are returned in mm, mm², mm³, mm⁴, and weight in kgf/m.

    Arguments:
        sect_name (str): The reference name of the rectangular tubular section. 
            Example values:
                - TUB.300X150X8
                - TUB.400X200X10

    Returns:
        list: A dictionary containing the following geometric properties:
            - sect_name (str)
            - d_sect (float)
            - tw_sect (float)
            - b_sect (float)
            - tf_sect (float)
            - Ag_sect (float)
            - Ix_sect (float)
            - Iy_sect (float)
            - Sx_sect (float)
            - Sy_sect (float)
            - Zx_sect (float)
            - Zy_sect (float)
            - rx_sect (float)
            - ry_sect (float)
            - weight_sect (float)
        ```

    Raises:
        ValueError: If the section is not found in the database.
        FileNotFoundError: If the database file is missing.
    """

    # Define the path to the database file
    current_dir = os.path.dirname(__file__)
    db_path = os.path.join(current_dir, 'bd_tub_rect_cuad.csv')

    # Check if the database file exists
    if not os.path.exists(db_path):
        # Raise an error if the file is missing
        raise FileNotFoundError(f"The database file '{db_path}' was not found.")

    # Read the database into a pandas DataFrame
    bd_material = pd.read_csv(db_path, sep=';')

    # Filter the DataFrame to select rows matching the specified section name
    sect_select = bd_material[bd_material["Seccion_tub"] == sect_name]

    # Check if the filtered DataFrame is empty (section not found)
    if sect_select.empty:
        # Raise an error if the section is not found in the database
        raise ValueError(f"Section '{sect_name}' not found in the database.")

    # Extract geometric properties
    sect_name = str(sect_select["Seccion_tub"].iloc[0])
    d_sect = float(sect_select["H (mm)"].iloc[0])
    tw_sect = float(sect_select["tdis (mm)"].iloc[0])
    b_sect = float(sect_select["B (mm)"].iloc[0])
    tf_sect = float(sect_select["tdis (mm)"].iloc[0])

    Ag_sect = float(sect_select["Ag (mm2)"].iloc[0])
    Ix_sect = float(sect_select["Ix (mm4)"].iloc[0])
    Iy_sect = float(sect_select["Iy (mm4)"].iloc[0])
    Sx_sect = float(sect_select["Sx (mm3)"].iloc[0])
    Sy_sect = float(sect_select["Sy (mm3)"].iloc[0])
    Zx_sect = float(sect_select["Zx (mm3)"].iloc[0])
    Zy_sect = float(sect_select["Zy (mm3)"].iloc[0])
    rx_sect = float(sect_select["rx (mm)"].iloc[0])
    ry_sect = float(sect_select["ry (mm)"].iloc[0])

    weight_sect = float(sect_select["Peso (kg/m)"].iloc[0])

    # Group all properties into a dictionary
    prop_sect_tub = {
        "sect_name": sect_name,
        "d_sect": d_sect,
        "tw_sect": tw_sect, 
        "b_sect": b_sect, 
        "tf_sect": tf_sect, 
        "Ag_sect": Ag_sect, 
        "Ix_sect": Ix_sect, 
        "Iy_sect": Iy_sect, 
        "Sx_sect": Sx_sect, 
        "Sy_sect": Sy_sect, 
        "Zx_sect": Zx_sect, 
        "Zy_sect": Zy_sect, 
        "rx_sect": rx_sect, 
        "ry_sect": ry_sect, 
        "weight_sect": weight_sect
    }

    # Return the dictionary of properties
    return prop_sect_tub

def prop_i_bu(section_dims):
    """
    Calculate the geometric properties of a custom I-section (assembled).

    This function calculates the geometric properties of a custom I-section based on given dimensions 
    provided in a dictionary. It returns the properties in various units such as mm, mm², mm³, mm⁴, and 
    weight in kgf/m.

    Arguments:
        section_dims (dict): A dictionary containing the dimensions of the I-section with the following keys:
            - "d_sect" (float): Height of the I-section (mm).
            - "tw_sect" (float): Thickness of the web (mm).
            - "bf_sect" (float): Width of the flange (mm).
            - "tf_sect" (float): Thickness of the flange (mm).

    Returns:
        list: A list containing the following geometric properties:
            - sect_name (str)
            - d_sect (float)
            - tw_sect (float)
            - bf_sect (float)
            - tf_sect (float)
            - h_sect (float)
            - Ag_sect (float)
            - Ix_sect (float)
            - Iy_sect (float)
            - Sx_sect (float)
            - Sy_sect (float)
            - Zx_sect (float)
            - Zy_sect (float)
            - rx_sect (float)
            - ry_sect (float)
            - J_sect (float)
            - Cw_sect (float)
            - weight_sect (float)

    Raises:
        ValueError: If any of the input dimensions are non-positive or if required keys are missing.
    """
    # Validate dictionary keys
    required_keys = {"d_sect", "tw_sect", "bf_sect", "tf_sect"}
    if not required_keys.issubset(section_dims):
        raise ValueError(f"The dictionary must contain the keys: {required_keys}")

    # Extract dimensions
    d_sect = section_dims["d_sect"]
    tw_sect = section_dims["tw_sect"]
    bf_sect = section_dims["bf_sect"]
    tf_sect = section_dims["tf_sect"]

    # Validate the input dimensions
    if d_sect <= 0 or tw_sect <= 0 or bf_sect <= 0 or tf_sect <= 0:
        raise ValueError("All input dimensions must be positive.")

    # Calculate intermediate properties
    sect_name = f"HS{d_sect}x{tw_sect}x{bf_sect}x{tf_sect}"
    h_sect = d_sect - 2 * tf_sect
    Ag_sect = (h_sect * tw_sect) + 2 * (bf_sect * tf_sect)
    
    Ix_sect = (1/12 * bf_sect * d_sect**3) - 2 * (1/12 * (bf_sect - tw_sect) / 2 * h_sect**3)
    Iy_sect = (1/12 * h_sect * tw_sect**3) + 2 * (1/12 * tf_sect * bf_sect**3)
    
    Sx_sect = Ix_sect / (d_sect / 2)
    Sy_sect = Iy_sect / (bf_sect / 2)
    
    Zx_sect = (bf_sect * tf_sect * (d_sect - tf_sect)) + (1/4 * tw_sect * (d_sect - 2 * tf_sect)**2)
    Zy_sect = ((bf_sect**2 * tf_sect) / 2) + (1/4 * tw_sect**2 * (d_sect - 2 * tf_sect))
    
    rx_sect = math.sqrt(Ix_sect / Ag_sect)
    ry_sect = math.sqrt(Iy_sect / Ag_sect)
    
    J_sect = 1/3*(2*bf_sect*tf_sect**3+(d_sect-tf_sect)*tw_sect**3)
    Cw_sect = tf_sect*bf_sect**3*(d_sect-tf_sect)**2/24
    
    weight_sect = Ag_sect / 1000000 * 7850  # Weight in kgf/m (assuming steel with a density of 7850 kg/m³)

    sect_type = "Built_up"


    # Group all properties into a dictionary
    prop_sect_i = {
        "sect_name": sect_name, 
        "d_sect": d_sect,
        "tw_sect": tw_sect,
        "bf_sect": bf_sect,
        "tf_sect": tf_sect, 
        "h_sect": h_sect,
        "Ag_sect": Ag_sect, 
        "Ix_sect": Ix_sect, 
        "Iy_sect": Iy_sect, 
        "Sx_sect": Sx_sect, 
        "Sy_sect": Sy_sect, 
        "Zx_sect": Zx_sect, 
        "Zy_sect": Zy_sect, 
        "rx_sect": rx_sect, 
        "ry_sect": ry_sect, 
        "J_sect": J_sect, 
        "Cw_sect": Cw_sect, 
        "weight_sect": weight_sect,
        "sect_type": sect_type
    }

    # Return the dictionary of properties
    return prop_sect_i

##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------STRENGTH----------------------------------------------------
##-----------------------------------------------------------------------------------------------------##

def stren_flex_yield(prop_mat, prop_sect):
    """
    AISC 360-22 F2.1
    Calculate the flexural strength of a section based on yield criteria.

    This function computes the nominal flexural strength of a section using the plastic modulus
    and material yield stress, reduced by a resistance reduction coefficient.

    Arguments:
        prop_mat (dict): Material properties containing:
            - "Fy" (float): Yield stress of the material (MPa).

        prop_sect (dict): Section properties containing:
            - "Zx_sect" (float): Plastic modulus of the section (mm³).

    Returns:
        float: The reduced plastic moment capacity (Mp_sect), in the same unit system as the inputs (N*mm).

    Example:
        ```python
        prop_mat = {"Fy": 250}  # Yield stress in MPa
        prop_sect = {"Zx_sect": 500000}  # Plastic modulus in mm³
        Mp_sect = stren_flex_yield(prop_mat, prop_sect)
        print(Mp_sect)  # Output: 112500000.0 (N·mm)
        ```
    """
   
    # Extract material and section properties
    Fy = prop_mat["Fy"]  # Yield stress
    Zx_sect = prop_sect["Zx_sect"]  # Plastic modulus

    # Calculate the plastic moment capacity
    Mp_sect = Fy * Zx_sect

    return Mp_sect

def stren_flex_c_c_ltb(prop_mat, prop_sect, Lb, Cb):
    
    #AISC 360-22 F2
      
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    d_sect = prop_sect["d_sect"]
    tw_sect = prop_sect["tw_sect"]
    bf_sect = prop_sect["bf_sect"]
    tf_sect = prop_sect["tf_sect"]
    h_sect = prop_sect["h_sect"]
    Ag_sect = prop_sect["Ag_sect"]
    Ix_sect = prop_sect["Ix_sect"]
    Iy_sect = prop_sect["Iy_sect"]
    Sx_sect = prop_sect["Sx_sect"]
    Sy_sect = prop_sect["Sy_sect"]
    Zx_sect = prop_sect["Zx_sect"]
    Zy_sect = prop_sect["Zy_sect"]
    rx_sect = prop_sect["rx_sect"]
    ry_sect = prop_sect["ry_sect"]
    J_sect = prop_sect["J_sect"]
    Cw_sect = prop_sect["Cw_sect"]

    
    Mp = stren_flex_yield(prop_mat, prop_sect) #Calls function that calculates plastic moment
    
    
    c_sect = 1
    ho_sect = d_sect-tf_sect
    rts = sqrt(sqrt(Iy_sect*Cw_sect)/Sx_sect) #F2-7


    Lp = 1.76*ry_sect*sqrt(Es/Fy) #F2-5
    Lr = 1.95*rts*Es/(0.7*Fy)*sqrt(J_sect*c_sect/(Sx_sect*ho_sect)+sqrt((J_sect*c_sect/(Sx_sect*ho_sect))**2+6.76*(0.7*Fy/Es)**2)) #F2-6

    Fcr=Cb*math.pi**2*Es/((Lb/rts)**2)*sqrt(1+0.078*J_sect*c_sect/(Sx_sect*ho_sect)*(Lb/rts)**2) #F2-4

    if Lb<=Lp:
        Mn = Mp
    elif Lb<=Lr:
        Mn = min(Cb*(Mp-(Mp-0.7*Fy*Sx_sect)*((Lb-Lp)/(Lr-Lp))),Mp) #F2-2
    else:
        Mn = min(Fcr*Sx_sect,Mp) #F2-3

    output_flex ={
        "rts": rts,
        "Lp": Lp,
        "Lr": Lr,
        "Fcr": Fcr,
        "Mn": Mn
    }
    return output_flex

def stren_flex_ncs_c_flb(prop_mat, prop_sect, Lb, Cb,):
    
    #AISC 360-22 F3
      
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    d_sect = prop_sect["d_sect"]
    tw_sect = prop_sect["tw_sect"]
    bf_sect = prop_sect["bf_sect"]
    tf_sect = prop_sect["tf_sect"]
    h_sect = prop_sect["h_sect"]
    Ag_sect = prop_sect["Ag_sect"]
    Ix_sect = prop_sect["Ix_sect"]
    Iy_sect = prop_sect["Iy_sect"]
    Sx_sect = prop_sect["Sx_sect"]
    Sy_sect = prop_sect["Sy_sect"]
    Zx_sect = prop_sect["Zx_sect"]
    Zy_sect = prop_sect["Zy_sect"]
    rx_sect = prop_sect["rx_sect"]
    ry_sect = prop_sect["ry_sect"]
    J_sect = prop_sect["J_sect"]
    Cw_sect = prop_sect["Cw_sect"]
    
    Mp = stren_flex_yield(prop_mat, prop_sect) #Calls function that calculates plastic moment
    
    slend_flex = slend_i_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    
    kc = min(max(4/(sqrt(h_sect/tw_sect)),0.35),0.76)

    if slend_flange == "Noncompact":
        Mn = Mp-(Mp-0.7*Fy*Sx_sect)*((lamb_flange-lamb_p_flange)/(lamb_r_flange-lamb_p_flange))
    elif slend_flange == "Slender":
        Mn = 0.9*Es*kc*Sx_sect/(lamb_flange**2)

    Mn1 = Mp-(Mp-0.7*Fy*Sx_sect)*((lamb_flange-lamb_p_flange)/(lamb_r_flange-lamb_p_flange))

    output_flex ={"Mn": Mn,
                  "Mn1": Mn1}
    
    return output_flex

def stren_flex_cncs_nc_cfy(prop_mat, prop_sect, Lb, Cb,):
    
    #AISC 360-22 F4.1
      
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    d_sect = prop_sect["d_sect"]
    tw_sect = prop_sect["tw_sect"]
    bf_sect = prop_sect["bf_sect"]
    tf_sect = prop_sect["tf_sect"]
    h_sect = prop_sect["h_sect"]
    Ag_sect = prop_sect["Ag_sect"]
    Ix_sect = prop_sect["Ix_sect"]
    Iy_sect = prop_sect["Iy_sect"]
    Sx_sect = prop_sect["Sx_sect"]
    Sy_sect = prop_sect["Sy_sect"]
    Zx_sect = prop_sect["Zx_sect"]
    Zy_sect = prop_sect["Zy_sect"]
    rx_sect = prop_sect["rx_sect"]
    ry_sect = prop_sect["ry_sect"]
    J_sect = prop_sect["J_sect"]
    Cw_sect = prop_sect["Cw_sect"]
    
    Myc = Fy*Sx_sect
    Rpc = Mp/Myc
    Mn = Rpc*Myc
    
    output_flex ={"Mn": Mn}
    
    return output_flex

def stren_flex_cncs_nc_cfy(prop_mat, prop_sect, Lb, Cb,):
    
    #AISC 360-22 F4.1
      
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    d_sect = prop_sect["d_sect"]
    tw_sect = prop_sect["tw_sect"]
    bf_sect = prop_sect["bf_sect"]
    tf_sect = prop_sect["tf_sect"]
    h_sect = prop_sect["h_sect"]
    Ag_sect = prop_sect["Ag_sect"]
    Ix_sect = prop_sect["Ix_sect"]
    Iy_sect = prop_sect["Iy_sect"]
    Sx_sect = prop_sect["Sx_sect"]
    Sy_sect = prop_sect["Sy_sect"]
    Zx_sect = prop_sect["Zx_sect"]
    Zy_sect = prop_sect["Zy_sect"]
    rx_sect = prop_sect["rx_sect"]
    ry_sect = prop_sect["ry_sect"]
    J_sect = prop_sect["J_sect"]
    Cw_sect = prop_sect["Cw_sect"]
    
    Myc = Fy*Sx_sect
    Rpc = Mp/Myc
    Mn = Rpc*Myc
    
    output_flex ={"Mn": Mn}
    
    return output_flex

def stren_flex_cncs_nc_ltb(prop_mat, prop_sect, Lb, Cb,):
    
    #AISC 360-22 F4.1
      
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    d_sect = prop_sect["d_sect"]
    tw_sect = prop_sect["tw_sect"]
    bf_sect = prop_sect["bf_sect"]
    tf_sect = prop_sect["tf_sect"]
    h_sect = prop_sect["h_sect"]
    Ag_sect = prop_sect["Ag_sect"]
    Ix_sect = prop_sect["Ix_sect"]
    Iy_sect = prop_sect["Iy_sect"]
    Sx_sect = prop_sect["Sx_sect"]
    Sy_sect = prop_sect["Sy_sect"]
    Zx_sect = prop_sect["Zx_sect"]
    Zy_sect = prop_sect["Zy_sect"]
    rx_sect = prop_sect["rx_sect"]
    ry_sect = prop_sect["ry_sect"]
    J_sect = prop_sect["J_sect"]
    Cw_sect = prop_sect["Cw_sect"]
    
    Myc = Fy*Sx_sect
    Rpc = Mp/Myc
    Mn = Rpc*Myc
    FL = 0.7*Fy
    
    aw = h_sect*tw_sect/(bf_sect*tf_sect)
    rt = bf_sect/sqrt(12*(1+1/6*aw))

    Lp = 1.1*rt*sqrt(Es/Fy)
    Lr = 1.95*rt*Es/FL*sqrt(J_sect/(Sx_sect*h_sect)+sqrt((J_sect/(Sx_sect*h_sect))**2+6.76*(FL/Es)**2))

    Mn = Cb[Rpc*Myc-(Rpc*Myc-FL*Sx_sect)*(()/())]


    output_flex ={"Myc": Myc,
                  "aw": aw,
                  "Lp": Lp,
                  "Lr": Lr}
    
    return output_flex



##-----------------------------------------------------------------------------------------------------##
##--------------------------------------------SLENDERNESS FUNCTIONS----------------------------
##-----------------------------------------------------------------------------------------------------##

def slend_i_comp_ns(prop_mat, prop_sect):
    """
    Local buckling for I sections for members in compression - non seismic applications AISC360-22

    This function calculates the slenderness ratios of the flange and web of an I-section
    based on the material properties and section geometry. It classifies the section as
    "Nonslender" or "Slender" according to the criteria for compactness.

    Arguments:
        prop_mat (dict): Material properties containing:
            - "Fy" (float): Yield stress of the material (MPa or ksi).
            - "Es" (float): Modulus of elasticity of the material (MPa or ksi).

        prop_sect (dict): Section properties containing:
            - "bf_sect" (float): Flange width (mm or in).
            - "tf_sect" (float): Flange thickness (mm or in).
            - "h_sect" (float): Clear height of the web (mm or in).
            - "tw_sect" (float): Web thickness (mm or in).

    Returns:
        tuple: A tuple containing:
            - slend_flange (str): Classification of flange slenderness ("Nonslender" or "Slender").
            - slend_web (str): Classification of web slenderness ("Nonslender" or "Slender").

    Example:
        ```python
        prop_mat = {"Fy": 250, "Es": 200000}  # Yield stress in MPa, Modulus of elasticity in MPa
        prop_sect = {"bf_sect": 300, "tf_sect": 15, "h_sect": 500, "tw_sect": 10}  # Dimensions in mm
        slend_flange, slend_web = slend_i_comp_ns(prop_mat, prop_sect)
        print(slend_flange, slend_web)  # Output: "Nonslender", "Slender"
        ```
    """
    # Extract material properties
    Fy = prop_mat["Fy"]  # Yield stress
    Es = prop_mat["Es"]  # Modulus of elasticity

    # Extract section properties
    bf_sect = prop_sect["bf_sect"]  # Flange width
    tf_sect = prop_sect["tf_sect"]  # Flange thickness
    h_sect = prop_sect["h_sect"]  # Clear height of the web
    tw_sect = prop_sect["tw_sect"]  # Web thickness

    # Slenderness of flange
    lamb_r_flange = 0.56 * math.sqrt(Es / Fy)  # Limiting slenderness for non-compact flanges
    b_flange = bf_sect / 2  # Half flange width
    lamb_flange = b_flange / tf_sect  # Flange slenderness ratio

    slend_flange = "Nonslender" if lamb_flange <= lamb_r_flange else "Slender"

    # Slenderness of web
    lamb_r_web = 1.49 * math.sqrt(Es / Fy)  # Limiting slenderness for compact webs
    lamb_web = h_sect / tw_sect  # Web slenderness ratio

    slend_web = "Nonslender" if lamb_web <= lamb_r_web else "Slender"
    
    slend_comp_ns = {
        "lamb_flange": round(lamb_flange,2),
        "lamb_r_flange": round(lamb_r_flange,2),
        "slend_flange": slend_flange,
        "lamb_web": round(lamb_web,2),
        "lamb_r_web": round(lamb_r_web,2),
        "slend_web": slend_web
    } 

    return slend_comp_ns

def slend_i_flex_ns(prop_mat, prop_sect):
    """
    Local buckling for I sections in flexure - non seismic applications AISC360-22

    This function calculates the slenderness ratios of the flange and web of an I-section
    based on the material properties and section geometry. It classifies the section as
    "Nonslender" or "Slender" according to the criteria for compactness.

    Arguments:
        prop_mat (dict): Material properties containing:
            - "Fy" (float): Yield stress of the material (MPa or ksi).
            - "Es" (float): Modulus of elasticity of the material (MPa or ksi).

        prop_sect (dict): Section properties containing:
            - "bf_sect" (float): Flange width (mm or in).
            - "tf_sect" (float): Flange thickness (mm or in).
            - "h_sect" (float): Clear height of the web (mm or in).
            - "tw_sect" (float): Web thickness (mm or in).

    Returns:
        tuple: A tuple containing:
            - slend_flange (str): Classification of flange slenderness ("Nonslender" or "Slender").
            - slend_web (str): Classification of web slenderness ("Nonslender" or "Slender").

    Example:
        ```python
        prop_mat = {"Fy": 250, "Es": 200000}  # Yield stress in MPa, Modulus of elasticity in MPa
        prop_sect = {"bf_sect": 300, "tf_sect": 15, "h_sect": 500, "tw_sect": 10}  # Dimensions in mm
        slend_flange, slend_web = slend_i_comp_ns(prop_mat, prop_sect)
        print(slend_flange, slend_web)  # Output: "Nonslender", "Slender"
        ```
    """
    # Extract material properties
    Fy = prop_mat["Fy"]  # Yield stress
    Es = prop_mat["Es"]  # Modulus of elasticity

    # Extract section properties
    bf_sect = prop_sect["bf_sect"]  # Flange width
    tf_sect = prop_sect["tf_sect"]  # Flange thickness
    h_sect = prop_sect["h_sect"]  # Clear height of the web
    tw_sect = prop_sect["tw_sect"]  # Web thickness
    sect_type = prop_sect["sect_type"]

    # Slenderness of flange
    if sect_type == "Rolled":
        lamb_p_flange = 0.38 * math.sqrt(Es / Fy)  # Limiting slenderness for non-compact flanges
        lamb_r_flange = 1.0 * math.sqrt(Es / Fy)  # Limiting slenderness for non-compact flanges
    elif sect_type == "Built_up":
        kc = min(max(4/(sqrt(h_sect/tw_sect)),0.35),0.76)
        FL = 0.7*Fy
        lamb_p_flange = 0.38 * math.sqrt(Es / Fy)  # Limiting slenderness for non-compact flanges
        lamb_r_flange = 0.95 * math.sqrt(kc*Es / FL)  # Limiting slenderness for non-compact flang


    b_flange = bf_sect / 2  # Half flange width
    lamb_flange = b_flange / tf_sect  # Flange slenderness ratio

    if lamb_flange <= lamb_p_flange: slend_flange = "Compact" 
    elif lamb_flange <= lamb_r_flange: slend_flange = "Noncompact" 
    elif lamb_flange > lamb_r_flange: slend_flange = "Slender" 


    # Slenderness of web
    lamb_p_web = 3.76 * math.sqrt(Es / Fy)  # Limiting slenderness for compact webs
    lamb_r_web = 5.70 * math.sqrt(Es / Fy)  # Limiting slenderness for compact webs
    
    lamb_web = h_sect / tw_sect  # Web slenderness ratio

    if lamb_web <= lamb_p_web: slend_web = "Compact" 
    elif lamb_web <= lamb_r_web: slend_web = "Noncompact" 
    elif lamb_web > lamb_r_web: slend_web = "Slender"
  
    slend_flex_ns = {
        "lamb_flange": round(lamb_flange,2),
        "lamb_r_flange": round(lamb_r_flange,2),
        "lamb_p_flange": round(lamb_p_flange,2),
        "slend_flange": slend_flange,
        "lamb_web": round(lamb_web,2),
        "lamb_r_web": round(lamb_r_web,2),
        "lamb_p_web": round(lamb_p_web,2),
        "slend_web": slend_web
    } 

    return slend_flex_ns

def slend_i_s(prop_mat, prop_sect,element_type,struct_syst,axial_load=0):
    """
    Local buckling for I sections in flexure - non seismic applications AISC360-22

    This function calculates the slenderness ratios of the flange and web of an I-section
    based on the material properties and section geometry. It classifies the section as
    "Nonslender" or "Slender" according to the criteria for compactness.

    Arguments:
        prop_mat (dict): Material properties containing:
            - "Fy" (float): Yield stress of the material (MPa or ksi).
            - "Es" (float): Modulus of elasticity of the material (MPa or ksi).

        prop_sect (dict): Section properties containing:
            - "bf_sect" (float): Flange width (mm or in).
            - "tf_sect" (float): Flange thickness (mm or in).
            - "h_sect" (float): Clear height of the web (mm or in).
            - "tw_sect" (float): Web thickness (mm or in).

        element_type (str): Element Type between this options:
            - "Diagonal_brace"
            - "Colum/beam"
        
        struct_syst (str): Structural system between this options:
            - "Moment_frame"
            - "Other"

        axial_load (float): value of axial load if available in [N]

    Returns:
        tuple: A tuple containing:
            - slend_flange (str): Classification of flange slenderness ("Nonslender" or "Slender").
            - slend_web (str): Classification of web slenderness ("Nonslender" or "Slender").

    Example:
        ```python
        prop_mat = {"Fy": 250, "Es": 200000}  # Yield stress in MPa, Modulus of elasticity in MPa
        prop_sect = {"bf_sect": 300, "tf_sect": 15, "h_sect": 500, "tw_sect": 10}  # Dimensions in mm
        slend_flange, slend_web = slend_i_comp_ns(prop_mat, prop_sect)
        print(slend_flange, slend_web)  # Output: "Nonslender", "Slender"
        ```
    """
    # Extract material properties
    Fy = prop_mat["Fy"]  # Yield stress
    Es = prop_mat["Es"]  # Modulus of elasticity
    Ry = prop_mat["Ry"]  

    # Extract section properties
    bf_sect = prop_sect["bf_sect"]  # Flange width
    tf_sect = prop_sect["tf_sect"]  # Flange thickness
    h_sect = prop_sect["h_sect"]  # Clear height of the web
    tw_sect = prop_sect["tw_sect"]  # Web thickness
    Ag_sect = prop_sect["Ag_sect"] #Gross area


    # Slenderness of flange
    lamb_hd_flange = 0.30 * math.sqrt(Es / (Ry*Fy))  # Limiting slenderness for highly ductile flanges
    lamb_md_flange = 0.38 * math.sqrt(Es / (Ry*Fy))  # Limiting slenderness for moderately ductile flanges

    b_flange = bf_sect / 2  # Half flange width
    lamb_flange = b_flange / tf_sect  # Flange slenderness ratio

    if lamb_flange <= lamb_hd_flange: slend_flange = "Highly_Ductile" 
    elif lamb_flange <= lamb_md_flange: slend_flange = "Moderately_Ductile" 
    elif lamb_flange > lamb_md_flange: slend_flange = "Low_Ductility" 

    # Slenderness of web
    Pr = axial_load #Ultimate axial load 
    Ca = Pr/(Ry*Fy*Ag_sect)

    if element_type == "Diagonal_brace":
        lamb_hd_flange = 1.49 * math.sqrt(Es / (Ry*Fy))  # Limiting slenderness for highly ductile flanges
        lamb_md_flange = 1.49 * math.sqrt(Es / (Ry*Fy))  # Limiting slenderness for moderately ductile flanges
    elif struct_syst == "Moment_frame":
        Ca = 0.1
        lamb_hd_web_ax_low = 2.5*(1-Ca)**2.3*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for highly ductile flanges
        lamb_md_web_ax_low = 5.4*(1-Ca)**2.3*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for moderately ductile flanges
        Ca = 0.7
        lamb_hd_web_ax_high = 2.5*(1-Ca)**2.3*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for highly ductile flanges
        lamb_md_web_ax_high = 5.4*(1-Ca)**2.3*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for moderately ductile flanges
        Ca = Pr/(Ry*Fy*Ag_sect)
        lamb_hd_web_ax_calc = 2.5*(1-Ca)**2.3*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for highly ductile flanges
        lamb_md_web_ax_calc = 5.4*(1-Ca)**2.3*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for moderately ductile flanges
    elif struct_syst == "Other":
        Ca = 0.1
        if  Ca<=0.113:
            lamb_hd_web_ax_low = 2.45*(1-1.04*Ca)*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for highly ductile flanges
            lamb_md_web_ax_low = 3.76*(1-3.05*Ca)*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for moderately ductile flanges
        else:
            lamb_hd_web_ax_low = max(2.26*(1-0.38*Ca)*math.sqrt(Es / (Ry*Fy)),1.56*math.sqrt(Es / (Ry*Fy))) # Limiting slenderness for highly ductile flanges
            lamb_md_web_ax_low = max(2.61*(1-0.49*Ca)*math.sqrt(Es / (Ry*Fy)),1.56*math.sqrt(Es / (Ry*Fy))) # Limiting slenderness for moderately ductile flanges
        Ca = 0.7
        if  Ca<=0.113:
            lamb_hd_web_ax_high = 2.45*(1-1.04*Ca)*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for highly ductile flanges
            lamb_md_web_ax_high = 3.76*(1-3.05*Ca)*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for moderately ductile flanges
        else:
            lamb_hd_web_ax_high = max(2.26*(1-0.38*Ca)*math.sqrt(Es / (Ry*Fy)),1.56*math.sqrt(Es / (Ry*Fy))) # Limiting slenderness for highly ductile flanges
            lamb_md_web_ax_high = max(2.61*(1-0.49*Ca)*math.sqrt(Es / (Ry*Fy)),1.56*math.sqrt(Es / (Ry*Fy))) # Limiting slenderness for moderately ductile flanges        
        Ca = Pr/(Ry*Fy*Ag_sect)
        if  Ca<=0.113:
            lamb_hd_web_ax_calc = 2.45*(1-1.04*Ca)*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for highly ductile flanges
            lamb_md_web_ax_calc = 3.76*(1-3.05*Ca)*math.sqrt(Es / (Ry*Fy)) # Limiting slenderness for moderately ductile flanges
        else:
            lamb_hd_web_ax_calc = max(2.26*(1-0.38*Ca)*math.sqrt(Es / (Ry*Fy)),1.56*math.sqrt(Es / (Ry*Fy))) # Limiting slenderness for highly ductile flanges
            lamb_md_web_ax_calc = max(2.61*(1-0.49*Ca)*math.sqrt(Es / (Ry*Fy)),1.56*math.sqrt(Es / (Ry*Fy))) # Limiting slenderness for moderately ductile flanges 
    
    lamb_web = h_sect / tw_sect  # Web slenderness ratio

    #Axial load available
    if lamb_web <= lamb_hd_web_ax_calc: slend_web_axial_calc = "Highly_Ductile" 
    elif lamb_web <= lamb_md_web_ax_calc: slend_web_axial_calc = "Moderately_Ductile" 
    elif lamb_web > lamb_md_web_ax_calc: slend_web_axial_calc = "Low_Ductility" 

    #Low Axial load
    if lamb_web <= lamb_hd_web_ax_low: slend_web_axial_low = "Highly_Ductile" 
    elif lamb_web <= lamb_md_web_ax_low: slend_web_axial_low = "Moderately_Ductile" 
    elif lamb_web > lamb_md_web_ax_low: slend_web_axial_low = "Low_Ductility" 

    #High Axial load
    if lamb_web <= lamb_hd_web_ax_high: slend_web_axial_high = "Highly_Ductile" 
    elif lamb_web <= lamb_md_web_ax_high: slend_web_axial_high = "Moderately_Ductile" 
    elif lamb_web > lamb_md_web_ax_high: slend_web_axial_high = "Low_Ductility" 

    slend_i_s = {
        "lamb_flange": round(lamb_flange,2),
        "lamb_hd_flange": round(lamb_hd_flange,2),
        "lamb_md_flange": round(lamb_md_flange,2),
        "slend_flange": slend_flange,

        "lamb_web": round(lamb_web,2),
        "lamb_hd_web_ax_calc": round(lamb_hd_web_ax_calc,2),
        "lamb_md_web_ax_calc": round(lamb_md_web_ax_calc,2),
        "slend_web_axial_calc": slend_web_axial_calc,

        "lamb_hd_web_ax_low": round(lamb_hd_web_ax_low,2),
        "lamb_md_web_ax_low": round(lamb_md_web_ax_calc,2),
        "slend_web_axial_low": slend_web_axial_calc,
        
        "lamb_hd_web_ax_high": round(lamb_hd_web_ax_low,2),
        "lamb_md_web_ax_high": round(lamb_md_web_ax_high,2),
        "slend_web_axial_high": slend_web_axial_calc
    } 

    return slend_i_s