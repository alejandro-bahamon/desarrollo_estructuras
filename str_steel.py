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
    Calculate the geometric properties of a built-up I-section.
- stren_flex_i
    Calculate the flexural strength of a I-section.
- stren_flex_hss_rect
    Calculate the flexural strength of a Rectangular-HSS section.

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
            - A1011-G50

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
    bf_sect = float(sect_select["B (mm)"].iloc[0])
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
    J_sect = float(sect_select["Js (mm4)"].iloc[0])
    weight_sect = float(sect_select["Peso (kg/m)"].iloc[0])
    sect_type = "Rolled"

    # Group all properties into a dictionary
    prop_sect_tub = {
        "sect_name": sect_name,
        "d_sect": d_sect,
        "tw_sect": tw_sect, 
        "bf_sect": bf_sect, 
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
        "J_sect":J_sect,
        "weight_sect": weight_sect,
        "sect_type": sect_type
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

def stren_flex_i_y_F2_1(prop_mat, prop_sect):
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
        Mp_sect = stren_flex_i_y_F2_1(prop_mat, prop_sect)
        print(Mp_sect)  # Output: 112500000.0 (N·mm)
        ```
    """
   
    # Extract material and section properties
    Fy = prop_mat["Fy"]  # Yield stress
    Zx_sect = prop_sect["Zx_sect"]  # Plastic modulus

    # Calculate the plastic moment capacity
    Mn = Fy * Zx_sect

    return Mn

def stren_flex_i_ltb_F2_2(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the flexural strength of an I-section considering lateral-torsional buckling (LTB)
    per AISC 360-22, Section F2.2.

    Arguments:
        prop_mat (dict): Material properties containing:
            - Fy (float): Yield stress (MPa)
            - Es (float): Elastic modulus (MPa)
        prop_sect (dict): Section properties containing:
            - d_sect (float): Section depth (mm)
            - tw_sect (float): Web thickness (mm)
            - bf_sect (float): Flange width (mm)
            - tf_sect (float): Flange thickness (mm)
            - h_sect (float): Clear web height (mm)
            - Ag_sect (float): Gross area (mm²)
            - Ix_sect (float): Moment of inertia about x-axis (mm⁴)
            - Iy_sect (float): Moment of inertia about y-axis (mm⁴)
            - Sx_sect (float): Section modulus about x-axis (mm³)
            - Sy_sect (float): Section modulus about y-axis (mm³)
            - Zx_sect (float): Plastic section modulus about x-axis (mm³)
            - Zy_sect (float): Plastic section modulus about y-axis (mm³)
            - rx_sect (float): Radius of gyration about x-axis (mm)
            - ry_sect (float): Radius of gyration about y-axis (mm)
            - J_sect (float): Torsional constant (mm⁴)
            - Cw_sect (float): Warping constant (mm⁶)
        Lb (float): Unbraced length of the beam (mm)
        Cb (float): Lateral-torsional buckling modification factor

    Returns:
        float: Nominal flexural strength Mn (kN·m)
    """
     # Material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    # Section properties
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

    # Calculate plastic moment (Mp) using a separate function
    Mp = stren_flex_i_y_F2_1(prop_mat, prop_sect) #Calls function that calculates plastic moment
    
    # Calculate lateral-torsional buckling parameters
    c_sect = 1
    ho_sect = d_sect-tf_sect # Distance between flange centroids 
    
    # Radius of gyration for torsional buckling (Eq. F2-7)
    rts = sqrt(sqrt(Iy_sect*Cw_sect)/Sx_sect) #F2-7 

    # Limiting unbraced lengths
    Lp = 1.76*ry_sect*sqrt(Es/Fy) #F2-5
    Lr = 1.95*rts*Es/(0.7*Fy)*sqrt(J_sect*c_sect/(Sx_sect*ho_sect)+sqrt((J_sect*c_sect/(Sx_sect*ho_sect))**2+6.76*(0.7*Fy/Es)**2)) #F2-6

    # Elastic critical moment for LTB (Eq. F2-4)
    Fcr=Cb*math.pi**2*Es/((Lb/rts)**2)*sqrt(1+0.078*J_sect*c_sect/(Sx_sect*ho_sect)*(Lb/rts)**2) #F2-4

    # Determine nominal flexural strength Mn
    if Lb<=Lp:
        Mn = Mp # No lateral-torsional buckling
    elif Lb<=Lr:
        Mn = min(Cb*(Mp-(Mp-0.7*Fy*Sx_sect)*((Lb-Lp)/(Lr-Lp))),Mp) #F2-2 Inelastic buckling
    else:
        Mn = min(Fcr*Sx_sect,Mp) #F2-3 Elastic buckling

    return Mn

def stren_flex_i_flb_F3_2(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the flexural strength of an I-section considering flange local buckling (FLB)
    per AISC 360-22, Section F3.2.

    Arguments:
        prop_mat (dict): Material properties containing:
            - Fy (float): Yield stress (MPa)
            - Es (float): Elastic modulus (MPa)
        prop_sect (dict): Section properties containing:
            - d_sect (float): Section depth (mm)
            - tw_sect (float): Web thickness (mm)
            - bf_sect (float): Flange width (mm)
            - tf_sect (float): Flange thickness (mm)
            - h_sect (float): Clear web height (mm)
            - Ag_sect (float): Gross area (mm²)
            - Ix_sect (float): Moment of inertia about x-axis (mm⁴)
            - Iy_sect (float): Moment of inertia about y-axis (mm⁴)
            - Sx_sect (float): Section modulus about x-axis (mm³)
            - Sy_sect (float): Section modulus about y-axis (mm³)
            - Zx_sect (float): Plastic section modulus about x-axis (mm³)
            - Zy_sect (float): Plastic section modulus about y-axis (mm³)
            - rx_sect (float): Radius of gyration about x-axis (mm)
            - ry_sect (float): Radius of gyration about y-axis (mm)
            - J_sect (float): Torsional constant (mm⁴)
            - Cw_sect (float): Warping constant (mm⁶)
        Lb (float): Unbraced length of the beam (mm) (Not used in this function)
        Cb (float): Lateral-torsional buckling modification factor (Not used in this function)

    Returns:
        float: Nominal flexural strength Mn (N·mm)
    """
    # Material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    # Section properties
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
    
    # Calculate plastic moment (Mp) using a separate function
    Mp = stren_flex_i_y_F2_1(prop_mat, prop_sect) #Calls function that calculates plastic moment
    
    # Compute slenderness parameters
    slend_flex = slend_i_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    
    # Compute the flange buckling coefficient kc
    kc = min(max(4/(sqrt(h_sect/tw_sect)),0.35),0.76)

    # Determine nominal flexural strength Mn
    if slend_flange == "Noncompact":
        Mn = Mp-(Mp-0.7*Fy*Sx_sect)*((lamb_flange-lamb_p_flange)/(lamb_r_flange-lamb_p_flange)) #F3-1
    elif slend_flange == "Slender":
        Mn = 0.9*Es*kc*Sx_sect/(lamb_flange**2) #F3-2
           
    return Mn

def stren_flex_i_cfy_F4_1(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the flexural strength of an I-section considering compression flange yielding (CFY)
    per AISC 360-22, Section F4.1.

    Arguments:
        prop_mat (dict): Material properties containing:
            - Fy (float): Yield stress (MPa)
            - Es (float): Elastic modulus (MPa)
        prop_sect (dict): Section properties containing:
            - d_sect (float): Section depth (mm)
            - tw_sect (float): Web thickness (mm)
            - bf_sect (float): Flange width (mm)
            - tf_sect (float): Flange thickness (mm)
            - h_sect (float): Clear web height (mm)
            - Ag_sect (float): Gross area (mm²)
            - Ix_sect (float): Moment of inertia about x-axis (mm⁴)
            - Iy_sect (float): Moment of inertia about y-axis (mm⁴)
            - Sx_sect (float): Section modulus about x-axis (mm³)
            - Sy_sect (float): Section modulus about y-axis (mm³)
            - Zx_sect (float): Plastic section modulus about x-axis (mm³)
            - Zy_sect (float): Plastic section modulus about y-axis (mm³)
            - rx_sect (float): Radius of gyration about x-axis (mm)
            - ry_sect (float): Radius of gyration about y-axis (mm)
            - J_sect (float): Torsional constant (mm⁴)
            - Cw_sect (float): Warping constant (mm⁶)
        Lb (float): Unbraced length of the beam (mm) (Not used in this function)
        Cb (float): Lateral-torsional buckling modification factor (Not used in this function)

    Returns:
        float: Nominal flexural strength Mn (N·mm)
    """
    # Material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    # Section properties
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
    
    # Compute slenderness parameters
    slend_flex = slend_i_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]

    # Compute moment of inertia of the compression flange about its own axis
    Iyc_sect=1/12*tf_sect*bf_sect**3

    # Compute yield moment
    Myc=Fy*Sx_sect #F4-4

    # Compute plastic moment considering compactness
    Mp1 = min(Zx_sect*Fy,1.6*Sx_sect*Fy)

    # Compute Rpc factor for strength adjustment
    if Iyc_sect/Iy_sect>0.23:
        if h_sect/tw_sect<=lamb_p_web:
            Rpc=Mp1/Myc #F4-9a
        elif h_sect/tw_sect>lamb_p_web:
            Rpc=min(Mp1/Myc-((Mp1/Myc)-1)*((lamb_web-lamb_p_web)/(lamb_r_web-lamb_p_web)),Mp1/Myc) #F4-9b
    elif Iyc_sect/Iy_sect<=0.23:
        Rpc=1 #F4-10
    
    # Compute nominal flexural strength Mn
    Mn=Rpc*Myc #F4-1
    
    return Mn

def stren_flex_i_ltb_F4_2(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the flexural strength of an I-section considering lateral-torsional buckling (LTB)
    per AISC 360-22, Section F4.2.

    Arguments:
        prop_mat (dict): Material properties containing:
            - Fy (float): Yield stress (MPa)
            - Es (float): Elastic modulus (MPa)
        prop_sect (dict): Section properties containing:
            - d_sect (float): Section depth (mm)
            - tw_sect (float): Web thickness (mm)
            - bf_sect (float): Flange width (mm)
            - tf_sect (float): Flange thickness (mm)
            - h_sect (float): Clear web height (mm)
            - Ag_sect (float): Gross area (mm²)
            - Ix_sect (float): Moment of inertia about x-axis (mm⁴)
            - Iy_sect (float): Moment of inertia about y-axis (mm⁴)
            - Sx_sect (float): Section modulus about x-axis (mm³)
            - Sy_sect (float): Section modulus about y-axis (mm³)
            - Zx_sect (float): Plastic section modulus about x-axis (mm³)
            - Zy_sect (float): Plastic section modulus about y-axis (mm³)
            - rx_sect (float): Radius of gyration about x-axis (mm)
            - ry_sect (float): Radius of gyration about y-axis (mm)
            - J_sect (float): Torsional constant (mm⁴)
            - Cw_sect (float): Warping constant (mm⁶)
        Lb (float): Unbraced length of the beam (mm)
        Cb (float): Lateral-torsional buckling modification factor

    Returns:
        float: Nominal flexural strength Mn (N·mm)
    """
    # Material properties  
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]
    
    # Section properties
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

    # Compute slenderness parameters
    slend_flex = slend_i_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]
    
    # Compute lateral-torsional buckling parameters
    aw = h_sect*tw_sect/(bf_sect*tf_sect) #F4-12
    rt = bf_sect/sqrt(12*(1+1/6*aw)) #F4-11
    FL = 0.7*Fy #F4-6a

    Lp = 1.1*rt*sqrt(Es/Fy) #F4-7
    Lr = 1.95*rt*Es/FL*sqrt(J_sect/(Sx_sect*h_sect)+sqrt((J_sect/(Sx_sect*h_sect))**2+6.76*(FL/Es)**2)) #F4-8
    Fcr = Cb*pi**2*Es/(Lb/rt)**2*sqrt(1+0.078*(J_sect/(Sx_sect*h_sect*(Lb/rt)**2))) #F4-5

    # Compute moment of inertia of the compression flange about its own axis (Eq. F4-8)
    Iyc_sect=1/12*tf_sect*bf_sect**3

    # Compute yield moment
    Myc=Fy*Sx_sect #F4-4
    
    # Compute plastic moment considering compactness
    Mp1 = min(Zx_sect*Fy,1.6*Sx_sect*Fy)

    # Compute Rpc factor for strength adjustment
    if Iyc_sect/Iy_sect>0.23:
        if h_sect/tw_sect<=lamb_p_web:
            Rpc=Mp1/Myc #F4-9a
        elif h_sect/tw_sect>lamb_p_web:
            Rpc=min(Mp1/Myc-((Mp1/Myc)-1)*((lamb_web-lamb_p_web)/(lamb_r_web-lamb_p_web)),Mp1/Myc) #F4-9b
    elif Iyc_sect/Iy_sect<=0.23:
        Rpc=1 #F4-10

    # Compute nominal flexural strength Mn
    if Lb<=Lp:
        Mn=Mp1
    elif Lb<=Lr:
        Mn=min(Cb*(Rpc*Myc-(Rpc*Myc-FL*Sx_sect)*((Lb-Lp)/(Lr-Lp))),Rpc*Myc) #F4-2
    elif Lb>Lr:
        Mn=min(Fcr*Sx_sect,Rpc*Myc) #F4-3
  
    return Mn

def stren_flex_i_ltb_F4_3(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the flexural strength of an I-section per AISC 360-22, Section F4.3.
    
    Parameters:
        prop_mat (dict): Material properties including:
            - Fy: Yield strength (MPa)
            - Es: Elastic modulus (MPa)
        
        prop_sect (dict): Section properties including:
            - bf_sect: Flange width (mm)
            - tf_sect: Flange thickness (mm)
            - h_sect: Web depth (mm)
            - tw_sect: Web thickness (mm)
            - Sx_sect: Elastic section modulus about the x-axis (mm³)
            - Zx_sect: Plastic section modulus about the x-axis (mm³)
        
        Lb (float): Unbraced length of the beam (mm).
        Cb (float): Bending coefficient.

    Returns:
        float: Nominal flexural strength Mn (N·mm).
    """
    # Material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]
    
    # Section properties
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

    # Compute slenderness parameters
    slend_flex = slend_i_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]

    # Compute additional parameters
    kc = min(max(4/(sqrt(h_sect/tw_sect)),0.35),0.76)
    FL = 0.7*Fy
    
    # Compute plastic moment considering compactness
    Mp1 = min(Zx_sect*Fy,1.6*Sx_sect*Fy)

    # Compute moment of inertia of the compression flange about its own axis
    Iyc_sect=1/12*tf_sect*bf_sect**3
    
    # Compute yield moment
    Myc=Fy*Sx_sect

    # Compute Rpc factor for strength adjustment
    if Iyc_sect/Iy_sect>0.23:
        if h_sect/tw_sect<=lamb_p_web:
            Rpc=Mp1/Myc
        elif h_sect/tw_sect>lamb_p_web:
            Rpc=min(Mp1/Myc-((Mp1/Myc)-1)*((lamb_web-lamb_p_web)/(lamb_r_web-lamb_p_web)),Mp1/Myc)
    elif Iyc_sect/Iy_sect<=0.23:
        Rpc=1

    # Compute nominal flexural strength Mn
    if slend_flange == "Compact":
        Mn=Mp1
    elif slend_flange == "Noncompact":
        Mn=Rpc*Myc-(Rpc*Myc-FL*Sx_sect)*((lamb_flange-lamb_p_flange)/(lamb_r_flange-lamb_p_flange)) #F4-13
    elif slend_flange == "Slender":
        Mn=0.90*Es*kc*Sx_sect/(lamb_flange**2) #F4-14

    return Mn

def stren_flex_i_cfy_F5_1(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the flexural strength of an I-section per AISC 360-22, Section F5.1.
    
    Parameters:
        prop_mat (dict): Material properties including:
            - Fy: Yield strength (MPa)
            - Es: Elastic modulus (MPa)
        
        prop_sect (dict): Section properties including:
            - d_sect: Section depth (mm)
            - tw_sect: Web thickness (mm)
            - bf_sect: Flange width (mm)
            - tf_sect: Flange thickness (mm)
            - h_sect: Clear web height (mm)
            - Ag_sect: Gross area (mm²)
            - Ix_sect: Moment of inertia about x-axis (mm⁴)
            - Iy_sect: Moment of inertia about y-axis (mm⁴)
            - Sx_sect: Section modulus about x-axis (mm³)
            - Sy_sect: Section modulus about y-axis (mm³)
            - Zx_sect: Plastic section modulus about x-axis (mm³)
            - Zy_sect: Plastic section modulus about y-axis (mm³)
            - rx_sect: Radius of gyration about x-axis (mm)
            - ry_sect: Radius of gyration about y-axis (mm)
            - J_sect: Torsional constant (mm⁴)
            - Cw_sect: Warping constant (mm⁶)
        
        Lb (float): Unbraced length of the beam (mm).
        Cb (float): Lateral-torsional buckling modification factor.

    Returns:
        float: Nominal flexural strength Mn (N·mm).
    """
    # Material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    # Section properties
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

    # Compute slenderness parameters
    slend_flex = slend_i_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]

    # Compute adjustment factors
    aw = min(h_sect*tw_sect/(bf_sect*tf_sect),10) #F4-12
    Rpg = min(1-aw/(1200+300*aw)*((h_sect/tw_sect)-5.7*sqrt(Es/Fy)),1) #F5-6

    # Compute nominal flexural strength
    Mn = Rpg*Fy*Sx_sect
    
    return Mn

def stren_flex_i_ltb_F5_2(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the flexural strength of an I-section per AISC 360-22, Section F5.2.
    
    Parameters:
        prop_mat (dict): Material properties including:
            - Fy: Yield strength (MPa)
            - Es: Elastic modulus (MPa)
        
        prop_sect (dict): Section properties including:
            - d_sect: Section depth (mm)
            - tw_sect: Web thickness (mm)
            - bf_sect: Flange width (mm)
            - tf_sect: Flange thickness (mm)
            - h_sect: Clear web height (mm)
            - Ag_sect: Gross area (mm²)
            - Ix_sect: Moment of inertia about x-axis (mm⁴)
            - Iy_sect: Moment of inertia about y-axis (mm⁴)
            - Sx_sect: Section modulus about x-axis (mm³)
            - Sy_sect: Section modulus about y-axis (mm³)
            - Zx_sect: Plastic section modulus about x-axis (mm³)
            - Zy_sect: Plastic section modulus about y-axis (mm³)
            - rx_sect: Radius of gyration about x-axis (mm)
            - ry_sect: Radius of gyration about y-axis (mm)
            - J_sect: Torsional constant (mm⁴)
            - Cw_sect: Warping constant (mm⁶)
        
        Lb (float): Unbraced length of the beam (mm).
        Cb (float): Lateral-torsional buckling modification factor.

    Returns:
        float: Nominal flexural strength Mn (N·mm).
    """
    # Material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    # Section properties
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

    # Compute slenderness parameters
    slend_flex = slend_i_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]

    # Limiting laterally unbraced length for full plastic strength
    Lp = 1.1*rt*sqrt(Es/Fy) #F4-7
    aw = h_sect*tw_sect/(bf_sect*tf_sect) #F4-12
    rt = bf_sect/sqrt(12*(1+1/6*aw)) #F4-11
    Rpg = min(1-aw/(1200+300*aw)*((h_sect/tw_sect)-5.7*sqrt(Es/Fy)),1) #F5-6
    # Limiting laterally unbraced length for inelastic buckling
    Lr = pi*rt*sqrt(Es/(0.7*Fy)) #F5-5
    # Plastic moment capacity considering compactness
    Mp1 = min(Zx_sect*Fy,1.6*Sx_sect*Fy)

    # Compute nominal flexural strength
    if Lb<= Lp:
        Mn =  Mp1
    elif Lb<=Lr:
        Fcr = min(Cb*(Fy-(0.3*Fy)*((Lb-Lp)/(Lr-Lp))),Fy) #F5-3
        Mn = Rpg*Fcr*Sx_sect
    elif Lb>Lr:
        Fcr = min(Cb*pi**2*Es/((Lb/rt)**2),Fy)
        Mn = Rpg*Fcr*Sx_sect
    
    return Mn

def stren_flex_i_flb_F5_3(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the flexural strength of an I-section per AISC 360-22, Section F5.3.
    
    Parameters:
        prop_mat (dict): Material properties including:
            - Fy: Yield strength (MPa)
            - Es: Elastic modulus (MPa)
        
        prop_sect (dict): Section properties including:
            - d_sect: Section depth (mm)
            - tw_sect: Web thickness (mm)
            - bf_sect: Flange width (mm)
            - tf_sect: Flange thickness (mm)
            - h_sect: Clear web height (mm)
            - Ag_sect: Gross area (mm²)
            - Ix_sect: Moment of inertia about x-axis (mm⁴)
            - Iy_sect: Moment of inertia about y-axis (mm⁴)
            - Sx_sect: Section modulus about x-axis (mm³)
            - Sy_sect: Section modulus about y-axis (mm³)
            - Zx_sect: Plastic section modulus about x-axis (mm³)
            - Zy_sect: Plastic section modulus about y-axis (mm³)
            - rx_sect: Radius of gyration about x-axis (mm)
            - ry_sect: Radius of gyration about y-axis (mm)
            - J_sect: Torsional constant (mm⁴)
            - Cw_sect: Warping constant (mm⁶)
        
        Lb (float): Unbraced length of the beam (mm).
        Cb (float): Lateral-torsional buckling modification factor.

    Returns:
        float: Nominal flexural strength Mn (N·mm).
    """
    # Extract material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    # Extract section properties
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

    # Compute slenderness parameters
    slend_flex = slend_i_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]

    
    kc = min(max(4/(sqrt(h_sect/tw_sect)),0.35),0.76)
    # Plastic moment capacity considering compactness
    Mp1 = min(Zx_sect*Fy,1.6*Sx_sect*Fy)

    # Compute nominal flexural strength based on flange slenderness category
    if slend_flange == "Compact":
        Mn = Mp1
    elif slend_flange == "Noncompact":
        Fcr = Fy-(0.3*Fy)*((lamb_flange-lamb_p_flange)/(lamb_r_flange-lamb_p_flange)) #F5-8
        Mn = Rpg*Fcr*Sx_sect
    elif slend_flange == "Slender":
        Fcr = 0.9*Es*kc/(lamb_flange**2) #F5-9
        Mn = Rpg*Fcr*Sx_sect
    
    return Mn

def stren_flex_i(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the nominal flexural strength of an I-section per AISC 360-22, Section F.
    
    Parameters:
        prop_mat (dict): Material properties including:
            - Fy: Yield strength (MPa)
            - Es: Elastic modulus (MPa)
        
        prop_sect (dict): Section properties including:
            - d_sect, tw_sect, bf_sect, tf_sect, h_sect, Ag_sect
            - Ix_sect, Iy_sect, Sx_sect, Sy_sect, Zx_sect, Zy_sect
            - rx_sect, ry_sect, J_sect, Cw_sect
        
        Lb (float): Unbraced length of the beam (mm).
        Cb (float): Lateral-torsional buckling modification factor.

    Returns:
        float: Factored nominal flexural strength (N·mm).
    """
    # Extract material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    # Extract section properties
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

    # Compute slenderness parameters
    slend_flex = slend_i_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    slend_web= slend_flex["slend_web"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]

    # Determine flexural strength based on slenderness classification
    if slend_web == "Compact":
        if slend_flange == "Compact":
            #2.1
            Mn_F2_1 = stren_flex_i_y_F2_1(prop_mat, prop_sect)
            #2.2
            Mn_F2_2 = stren_flex_i_ltb_F2_2(prop_mat, prop_sect, Lb, Cb)

            Mn_fin=min(Mn_F2_1,Mn_F2_2)
        else:
            #2.2
            Mn_F2_2 = stren_flex_i_ltb_F2_2(prop_mat, prop_sect, Lb, Cb)
            #3.2
            Mn_F3_2 = stren_flex_i_flb_F3_2(prop_mat, prop_sect, Lb, Cb)

            Mn_fin=min(Mn_F2_2,Mn_F3_2)
        
    elif slend_web == "Noncompact":
        #4.1
        Mn_F4_1 = stren_flex_i_cfy_F4_1(prop_mat, prop_sect, Lb, Cb)
        #4.2
        Mn_F4_2 = stren_flex_i_ltb_F4_2(prop_mat, prop_sect, Lb, Cb)
        #4.3
        Mn_F4_3 = stren_flex_i_ltb_F4_3(prop_mat, prop_sect, Lb, Cb)
        
        Mn_fin=min(Mn_F4_1,Mn_F4_2,Mn_F4_3)

    elif slend_web == "Slender":
        #5.1
        Mn_F5_1 = stren_flex_i_cfy_F5_1(prop_mat, prop_sect, Lb, Cb)
        #5.2
        Mn_F5_2 = stren_flex_i_ltb_F5_2(prop_mat, prop_sect, Lb, Cb)
        #5.3
        Mn_F5_3 = stren_flex_i_flb_F5_3(prop_mat, prop_sect, Lb, Cb)

        Mn_fin=min(Mn_F5_1,Mn_F5_2,Mn_F5_3)
    
    # Apply resistance factor
    phi_b = 0.90 #
    Mres = phi_b*Mn_fin
    
    return Mres

def stren_flex_hss_rect_y_F7_1(prop_mat, prop_sect, Lb, Cb):   
    """
    Compute the nominal flexural strength of a rectangular HSS section per AISC 360-22, Section F7.1.
    
    Parameters:
        prop_mat (dict): Material properties including:
            - Fy: Yield strength (MPa)
            - Es: Elastic modulus (MPa)
        
        prop_sect (dict): Section properties including:
            - d_sect, tw_sect, bf_sect, tf_sect, Ag_sect
            - Ix_sect, Iy_sect, Sx_sect, Sy_sect, Zx_sect, Zy_sect
            - rx_sect, ry_sect, J_sect, sect_type
        
        Lb (float): Unbraced length of the beam (mm).
        Cb (float): Lateral-torsional buckling modification factor.

    Returns:
        float: Nominal flexural strength (N·mm).
    """
    # Extract material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    # Extract section properties
    d_sect = prop_sect["d_sect"]
    tw_sect = prop_sect["tw_sect"]
    bf_sect = prop_sect["bf_sect"]
    tf_sect = prop_sect["tf_sect"]
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
    sect_type = prop_sect["sect_type"]

    # Compute slenderness parameters
    slend_flex = slend_hss_rect_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation   
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    slend_web = slend_flex["slend_web"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]

    # Compute plastic moment capacity
    Mp = Fy*Zx_sect #F7-1
    Mn = Mp

    return Mn

def stren_flex_hss_rect_flb_F7_2(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the nominal flexural strength of a rectangular HSS per AISC 360-22, Section F7.2.
    
    Parameters:
        prop_mat (dict): Material properties including:
            - Fy: Yield strength (MPa)
            - Es: Elastic modulus (MPa)
        
        prop_sect (dict): Section properties including:
            - d_sect, tw_sect, bf_sect, tf_sect, Ag_sect
            - Ix_sect, Iy_sect, Sx_sect, Sy_sect, Zx_sect, Zy_sect
            - rx_sect, ry_sect, J_sect, sect_type
        
        Lb (float): Unbraced length of the beam (mm).
        Cb (float): Lateral-torsional buckling modification factor.

    Returns:
        float: Nominal flexural strength (N·mm).
    """
    # Extract material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    # Extract section properties
    d_sect = prop_sect["d_sect"]
    tw_sect = prop_sect["tw_sect"]
    bf_sect = prop_sect["bf_sect"]
    tf_sect = prop_sect["tf_sect"]
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
    sect_type = prop_sect["sect_type"]

    # Compute slenderness parameters
    slend_flex = slend_hss_rect_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation   
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    slend_web = slend_flex["slend_web"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]

    # Plastic moment
    Mp = Fy*Zx_sect #F7-1

    if slend_flange == "Compact":
        Mn = Mp
    elif slend_flange == "Noncompact":
        Mn = min(Mp-(Mp-Fy*Sx_sect)*(3.57*lamb_flange*sqrt(Fy/Es)-4),Mp)
    elif slend_flange == "Slender":
        if sect_type == "Rolled":
            be = min(1.95*tf_sect*sqrt(Es/Fy)*(1-0.38/lamb_flange*sqrt(Es/Fy)),bf_sect)
        elif sect_type == "Built_up":
            be = min(1.92*tf_sect*sqrt(Es/Fy)*(1-0.34/lamb_flange*sqrt(Es/Fy)),bf_sect)

        if sect_type == "Rolled":
            be_int = be-3*tdis #Clear distance between webs minus radius
            d_int = d_sect-3*tdis #Clear distance between flanges minus radius
        elif sect_type == "Built_up":
            be_int = be-2*tw_sect #Clear distance between webs
            d_int = d_sect-2*tf_sect #Clear distance between flanges
        Se = (be*d_sect**2/6)-(be_int*d_int**3/(6*d_sect))
        Mn = Fy*Se #F7-3

    return Mn

def stren_flex_hss_rect_wlb_F7_3(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the nominal flexural strength of a rectangular HSS per AISC 360-22, Section F7.3.
    
    Parameters:
        prop_mat (dict): Material properties including:
            - Fy: Yield strength (MPa)
            - Es: Elastic modulus (MPa)
        
        prop_sect (dict): Section properties including:
            - d_sect, tw_sect, bf_sect, tf_sect, Ag_sect
            - Ix_sect, Iy_sect, Sx_sect, Sy_sect, Zx_sect, Zy_sect
            - rx_sect, ry_sect, J_sect, sect_type
        
        Lb (float): Unbraced length of the beam (mm).
        Cb (float): Lateral-torsional buckling modification factor.

    Returns:
        float: Nominal flexural strength (N·mm).
    """
    # Extract material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    # Extract section properties
    d_sect = prop_sect["d_sect"]
    tw_sect = prop_sect["tw_sect"]
    bf_sect = prop_sect["bf_sect"]
    tf_sect = prop_sect["tf_sect"]
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
    sect_type = prop_sect["sect_type"]

    slend_flex = slend_hss_rect_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation   
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    slend_web = slend_flex["slend_web"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]

    Mp = Fy*Zx_sect #F7-1

    if slend_web == "Compact":
        Mn = Mp
    elif slend_flange == "Noncompact":
        Mn = min(Mp-(Mp-Fy*Sx_sect)*(0.305*lamb_web*sqrt(Fy/Es)-0.738),Mp) #F7-6

    elif slend_flange == "Slender":
        if sect_type == "Rolled":
            b_int = bf_sect-3*tdis #Clear distance between webs minus radius
            d_int = d_sect-3*tdis #Clear distance between flanges minus radius
        elif sect_type == "Built_up":
            b_int = bf_sect-2*tw_sect #Clear distance between webs
            d_int = d_sect-2*tf_sect #Clear distance between flanges
        aw = 2*d_int*tw_sect/(b_int*tf_sect)
        kc = 4.0
        Rpg = min(1-aw/(1200+300*aw)*(lamb_web-5.7*sqrt(Es/Fy)),1)
        Mn1 = Rpg*Fy*Sx_sect #F7-7
        
        Fcr = 0.9*Es*kc/(lamb_flange**2)
        Mn2 = Rpg*Fcr*Sx_sect

        Mn = min(Mn1,Mn2)
    return Mn

def stren_flex_hss_rect_ltb_F7_4(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the nominal flexural strength of a rectangular HSS considering lateral-torsional buckling per AISC 360-22, Section F7.4.
    
    Parameters:
        prop_mat (dict): Material properties including:
            - Fy: Yield strength (MPa)
            - Es: Elastic modulus (MPa)
        
        prop_sect (dict): Section properties including:
            - d_sect, tw_sect, bf_sect, tf_sect, Ag_sect
            - Ix_sect, Iy_sect, Sx_sect, Sy_sect, Zx_sect, Zy_sect
            - rx_sect, ry_sect, J_sect, sect_type
        
        Lb (float): Unbraced length of the beam (mm).
        Cb (float): Lateral-torsional buckling modification factor.

    Returns:
        float: Nominal flexural strength (N·mm).
    """
    # Extract material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    # Extract section properties
    d_sect = prop_sect["d_sect"]
    tw_sect = prop_sect["tw_sect"]
    bf_sect = prop_sect["bf_sect"]
    tf_sect = prop_sect["tf_sect"]
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
    sect_type = prop_sect["sect_type"]

   
    slend_flex = slend_hss_rect_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation   
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    slend_web = slend_flex["slend_web"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]

    # Compute plastic moment
    Mp = Fy*Zx_sect #F7-1

    # Compute lateral-torsional buckling limits
    Lp = 0.13*Es*ry_sect*sqrt(J_sect*Ag_sect)/Mp #F7-12
    Lr = 2*Es*ry_sect*sqrt(J_sect*Ag_sect)/(0.7*Fy*Sx_sect) #F7-13

    # Compute nominal flexural strength based on Lb range
    if Lb<=Lp:
        Mn = Mp
    elif Lb<=Lr:
        Mn = min(Cb*(Mp-(Mp-0.7*Fy*Sx_sect)*((Lb-Lp)/(Lr-Lp))),Mp)
    elif Lb>Lr:
        Mn = min(2*E*Cb*sqrt(J_sect*Ag_sect)/(Lb/ry_sect),Mp)

    return Mn

def stren_flex_hss_rect(prop_mat, prop_sect, Lb, Cb):
    """
    Compute the nominal flexural strength of a rectangular HSS per AISC 360-22, Section F.
    
    Parameters:
        prop_mat (dict): Material properties including:
            - Fy: Yield strength (MPa)
            - Es: Elastic modulus (MPa)
        
        prop_sect (dict): Section properties including:
            - d_sect, tw_sect, bf_sect, tf_sect, h_sect, Ag_sect
            - Ix_sect, Iy_sect, Sx_sect, Sy_sect, Zx_sect, Zy_sect
            - rx_sect, ry_sect, J_sect, Cw_sect
        
        Lb (float): Unbraced length of the beam (mm).
        Cb (float): Lateral-torsional buckling modification factor.

    Returns:
        float: Factored flexural strength (N·mm).
    """
    # Extract material properties
    Fy = prop_mat["Fy"]
    Es = prop_mat["Es"]

    d_sect = prop_sect["d_sect"]
    tw_sect = prop_sect["tw_sect"]
    bf_sect = prop_sect["bf_sect"]
    tf_sect = prop_sect["tf_sect"]
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
    sect_type = prop_sect["sect_type"]

    # Compute slenderness parameters
    slend_flex = slend_hss_rect_flex_ns(prop_mat, prop_sect) #Calls function slenderness calculation 
    slend_flange = slend_flex["slend_flange"]
    lamb_flange= slend_flex["lamb_flange"]
    lamb_r_flange= slend_flex["lamb_r_flange"]
    lamb_p_flange= slend_flex["lamb_p_flange"]
    slend_web= slend_flex["slend_web"]
    lamb_web= slend_flex["lamb_web"]
    lamb_r_web= slend_flex["lamb_r_web"]
    lamb_p_web= slend_flex["lamb_p_web"]
    
    # Compute nominal moment capacities based on different limit states
    Mn_F7_1 = stren_flex_hss_rect_y_F7_1(prop_mat, prop_sect, Lb, Cb)
    Mn_F7_2 = stren_flex_hss_rect_flb_F7_2(prop_mat, prop_sect, Lb, Cb)
    Mn_F7_3 = stren_flex_hss_rect_wlb_F7_3(prop_mat, prop_sect, Lb, Cb)
    Mn_F7_4 = stren_flex_hss_rect_ltb_F7_4(prop_mat, prop_sect, Lb, Cb)

    # Determine controlling nominal moment
    Mn_fin=min(Mn_F7_1,Mn_F7_2,Mn_F7_3,Mn_F7_4)
    
    # Apply resistance factor
    phi_b = 0.90
    Mres = phi_b*Mn_fin

    return Mres

def stren_comp_flex_buck_ns(prop_mat,prop_sect,Kef,Long):   
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

    rmin=min(rx_sect,ry_sect)

    glob_slend = Kef*Long/rmin #Esbeltez global del elemento a compresión
    Fe=(pi**2*Es)/(Kef*Long/rmin)**2 #Esfuerzo crítico de pandeo elástico
    if Fe>=0.44*Fy:
        Fcr=(0.658**(Fy/Fe))*Fy
    elif Fe<0.44*Fy:
        Fcr = 0.877*Fe
    Cr = Fcr*Ag_sect #Resistencia a compresión

    return Cr


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
        "lamb_flange": lamb_flange,
        "lamb_r_flange": lamb_r_flange,
        "lamb_p_flange": lamb_p_flange,
        "slend_flange": slend_flange,
        "lamb_web": lamb_web,
        "lamb_r_web": lamb_r_web,
        "lamb_p_web": lamb_p_web,
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

def slend_hss_rect_flex_ns(prop_mat, prop_sect):

    # Extract material properties
    Fy = prop_mat["Fy"]  # Yield stress
    Es = prop_mat["Es"]  # Modulus of elasticity

    # Extract section properties
    bf_sect = prop_sect["bf_sect"]  # Flange width
    tf_sect = prop_sect["tf_sect"]  # Flange thickness
    d_sect = prop_sect["d_sect"]  # Clear height of the web
    tw_sect = prop_sect["tw_sect"]  # Web thickness
    sect_type = prop_sect["sect_type"]

    tdis = 0.93*tf_sect #Design wall thickness B.4.2

    if sect_type == "Rolled":
        bf_int = bf_sect-3*tdis #Clear distance between webs minus radius
        d_int = d_sect-3*tdis #Clear distance between flanges minus radius
    elif sect_type == "Built_up":
        bf_int = bf_sect-2*tw_sect #Clear distance between webs
        d_int = d_sect-2*tf_sect #Clear distance between flanges

    # Slenderness of flange
    lamb_p_flange = 1.12 * math.sqrt(Es / Fy)  # Limiting slenderness for compact flanges
    lamb_r_flange = 1.40 * math.sqrt(Es / Fy)  # Limiting slenderness for non compact flanges

    lamb_flange = bf_int/tf_sect # Flange slenderness ratio

    if lamb_flange <= lamb_p_flange: slend_flange = "Compact" 
    elif lamb_flange <= lamb_r_flange: slend_flange = "Noncompact" 
    elif lamb_flange > lamb_r_flange: slend_flange = "Slender" 

    # Slenderness of web
    lamb_p_web = 2.42 * math.sqrt(Es / Fy)  # Limiting slenderness for compact webs
    lamb_r_web = 5.70 * math.sqrt(Es / Fy)  # Limiting slenderness for compact webs
    
    lamb_web = d_int/tw_sect # Web slenderness ratio

    if lamb_web <= lamb_p_web: slend_web = "Compact" 
    elif lamb_web <= lamb_r_web: slend_web = "Noncompact" 
    elif lamb_web > lamb_r_web: slend_web = "Slender"
  
    slend_flex_ns = {
        "lamb_flange": lamb_flange,
        "lamb_r_flange": lamb_r_flange,
        "lamb_p_flange": lamb_p_flange,
        "slend_flange": slend_flange,
        "lamb_web": lamb_web,
        "lamb_r_web": lamb_r_web,
        "lamb_p_web": lamb_p_web,
        "slend_web": slend_web
    } 

    return slend_flex_ns