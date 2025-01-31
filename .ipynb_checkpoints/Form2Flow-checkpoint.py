# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# This program creates an app using Python tkinter. The purpose of this app is to automate the process of file creation inputs for the MultiFlow3D digital twin. The app displays 4 tabs - Control, Infodom, Mdmap, and LPT. Each tab displays a form for the user to fill out and the backend of the app uses the inputs to create files in the correct format. The user can save the files to their desktop.

import pandas as pd
import numpy as np
import customtkinter as ctk
from customtkinter import CTkEntry
import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
import os
import math
import tempfile

#create a global list to hold error messages which will be displayed to the user on submit
error_messages = []
#create a variable to check if there are user inputs to restore after a refresh
saved_user_inputs = None
# Global variable to track whether the splash screen has been shown
splash_shown = False

# Step 1: Dictionaries are created to keep track of variables for each question in each section in each tab. These dictionaries are used to create the forms for each tab.

#Control Sections with keys, questions, input_types, and options in the case of menu or radio button input types
control_sections = [
    {
        "section_name": "Numeric Parameters",
        "questions": [
            {"key":"keyword","question": "What kind of scenario do you want to simulate?", 
             "input_type": "dropdown","options": ["channel", "cavity", "column"]},
            {"key":"ubulk","question":"What is the ubulk value?","input_type": "float"},
            {"key":"dx","question":"What is the mesh size for x?","input_type": "float"},
            {"key":"dy","question":"What is the mesh size for y?","input_type": "float"},
            {"key":"dz","question":"What is the mesh size for z?","input_type": "float"},
            {"key":"fluid","question":"What will be the main fluid, water or air?",
             "input_type": "radio","options": ["water", "air"]},
            {"key":"differencing","question":"What discretisation scheme do you choose?",
             "input_type": "dropdown","options":["2nd order CDS", "4th order CDS","WENO"]},
            {"key":"dt","question":"Which timestep do you want?","input_type":"float"},
            {"key":"itime_end","question":"How many timesteps do you want to run?","input_type":"int"},
            {"key":"restart","question": "Is this a restart?","input_type":"radio","options":["yes","no"]},
            {"key":"reinitmean","question":"Do you want to calculate statistics from scratch?",
             "input_type":"radio","options":["YES","NO"]},
            {"key":"n_out","question":"How often do you want to produce outputs (tecbins)?","input_type":"int"},
            {"key":"l_transient","question":"Do you want to produce transient outputs?",
             "input_type":"radio","options":["yes","no"]},
            {"key":"results","question":"How often do you want to produce transient outputs? (If none, enter 0)",
             "input_type":"int"}
        ]
    },
    {
        "section_name": "Flow Boundary Conditions",
        "questions": [
            {"key":"west_boundary_condition","question":"Which west boundary condition option do you want?","input_type":
                                 "dropdown","options": ["inflow - uniform","inflow - 1/7th power law","inflow - prescribed",
                                                        "inflow - SEM", "wall - no-slip", "wall - smooth log law",
                                                        "wall - rough log law",
                                 "slip","periodic condition","1/6th law", "1/7th law","1/8th law"]},
            {"key":"east_boundary_condition","question":"Which east boundary condition option do you want?", 
                                 "input_type":"dropdown","options": ["outflow - Neumann", "outflow - convective", 
                                                                     "wall - no-slip","wall - smooth log law", 
                                                                     "wall - rough log law","slip","periodic condition",
                                                                    "1/6th law", "1/7th law", "1/8th law"]},
            {"key":"south_boundary_condition","question":"Which south boundary condition option do you want?",
                                 "input_type": "dropdown","options": ["wall - no-slip", "wall - smooth log law",
                                                                      "wall - rough log law", "slip","periodic condition", 
                                "1/6th law", "1/7th law", "1/8th law"]},
            {"key":"north_boundary_condition","question":"Which north boundary condition option do you want?", 
                                 "input_type":"dropdown","options": ["wall - no-slip", "wall - smooth log law",
                                                                      "wall - rough log law", "slip","periodic condition", 
                                "1/6th law", "1/7th law", "1/8th law"]},
            {"key":"bottom_boundary_condition","question":"Which bottom boundary condition option do you want?",
                                 "input_type": "dropdown","options": ["wall - no-slip", "wall - smooth log law",
                                                                      "wall - rough log law", "slip","periodic condition", 
                                "1/6th law", "1/7th law", "1/8th law"]},
            {"key":"top_boundary_condition","question":"Which top boundary condition option do you want?", 
                                 "input_type":"dropdown","options": ["wall - no-slip", "wall - smooth log law",
                                                                      "wall - rough log law", "slip","periodic condition", 
                                "1/6th law", "1/7th law", "1/8th law"]}
        ]
    },
    {
        "section_name": "Modelling Options",
        "questions": [
            {"key":"ctime_averaging","question":"At which ctime do you want to start collecting statistics?",
             "input_type": "float"}, 
            {"key":"SGS_model_choice","question":"Which SGS model do you want?", "input_type":"dropdown",
             "options": ["None","Smagorinsky","WALE","OneEqnModel","k-eps model (RANS)"]},
            {"key":"limb","question":"Do you want immersed boundaries to be on or off?","input_type": "radio",
             "options": ["on", "off"]},
            {"key":"lenergy","question":"Do you want energy equation?","input_type":"radio","options": ["yes", "no"]},
            {"key":"lpt","question":"Do you want to create lagrangian particles?","input_type":"radio",
             "options": ["yes", "no"]},
            {"key":"LSM","question":"Do you want to use LSM?","input_type":"radio","options":["yes","no"]},
            {"key":"L_LSMbase","question":"Do you want to use LSM base?","input_type":"radio","options":["yes","no"]},
            {"key":"lscalar","question":"Passive scalar?","input_type":"radio","options": ["yes", "no"]},
            {"key":"l_active_scalar","question":"Active scalar?","input_type":"radio","options":["yes","no"]},
            {"key":"l_non_newt","question":"L Non Newtonian?","input_type":"radio","options":["yes","no"]}
        ]
    }
]

#Advanced control sections additionally have defaults, which will store a value if the user doesn't choose to edit advanced
advanced_control_sections = [
    {
        "section_name": "Numeric Parameters",
        "questions": [
            {"key":"pr","question":"pr:","input_type":"float","default":0.72},
            {"key":"turb_schmidt","question":"turb_schmidt:","input_type":"float","default":0.6},
            {"key":"beta","question":"beta:","input_type": "float","default":0.0034},
            {"key":"gx","question":"gx:","input_type":"float","default":'0.0'},
            {"key":"gy","question":"gy:", "input_type":"float","default":'0.0'},
            {"key":"gz","question":"gz:", "input_type":"float","default":-9.81},
            {"key":"convection_scheme","question":"convection scheme:","input_type": "dropdown",
             "options":["Exp.Euler","AdamsBashfort","RK2","RK3"],"default":"RK2"},
            {"key":"diffusion_scheme","question":"diffusion scheme:", "input_type":"dropdown",
             "options":["Imp.Euler","CrankNicholson","ExplicitRungeKuttaDiff"],"default":"ExplicitRungeKuttaDiff"},
            {"key":"solver","question":"solver:", "input_type":"dropdown","options":["sip","mg"],"default":"mg"},
            {"key":"multigrid_step","question":"multigrid step (ngrid):","input_type": "int","default":1},
            {"key":"multigrid_iteration_scheme","question":"multigrid iteration scheme:", "input_type":"dropdown",
             "options":["GSS","TDMA","RBGSS"],"default":"TDMA"},
            {"key":"multigrid_max_iter","question":"multigrid maximum iteration per timestep:","input_type": "int",
             "default":30},
            {"key":"restriction_iter","question":"restriction iter:","input_type":"int","default":2},
            {"key":"prolongation_iter","question":"prolongation iter:", "input_type":"int","default":1},
            {"key":"variable_fixed","question":"variable or fixed dt:","input_type": "radio","options":["variable","fixed"],
            "default":"fixed"},
            {"key":"sweeps","question":"sweeps:","input_type": "int","default":25},
            {"key":"safety_factor","question":"safety factor:","input_type": "float","default":0.2},
            {"key":"niter","question":"niter:", "input_type":"int","default":20},
            {"key":"eps","question":"eps:","input_type": "float","default":1.0E-5},
            {"key":"nswp_1","question":"nswp(1):","input_type": "int","default":5},
            {"key":"nswp_2","question":"nswp(2):", "input_type":"int","default":5},
            {"key":"nswp_3","question":"nswp(3):", "input_type":"int","default":5},
            {"key":"nswp_4","question":"nswp(4):","input_type": "int","default":20}
        ]
    },
    {
        "section_name": "Flow Boundary Conditions",
                "questions": [
            {"key":"friction_coeff","question":"Friction coefficient:","input_type":"radio",
             "options": ["manning", "equivalent sand"],"default":"manning"},
            {"key":"friction_coeff_value","question":"What is the value of the friction coefficient?", "input_type":"float",
            "default":0.03},
            {"key":"save_inflow_data","question":"Save inflow data (precursor sim.):","input_type":"radio",
             "options":["True","False"],"default":"False"},
            {"key":"inlet_num","question":"Number of inlets:","input_type":"int","default":5000}
        ]
    },
    {
        "section_name": "Synthetic Eddy Method",
        "questions": [
            {"key":"velocity_profile","question":"Velocity Profile:","input_type":"radio",
             "options": ["Uniform", "1/7th PL"],"default":"1/7th PL"},
            {"key":"turbulence_intensity","question":"Turbulence intensity:","input_type":"float","default":0.1},
            {"key":"inlet_profile_num","question":"Number inlet profiles:","input_type":"int","default":1000}
        ]
    },
    {
        "section_name": "Modelling Options",
        "questions": [
            {"key":"time_averaging","question":"time_averaging:","input_type":"radio","options": ["True", "False"],
             "default":"True"},
            {"key":"noise","question":"noise:","input_type":"float","default":'0.0'},
            {"key":"LMR","question":"LMR:","input_type":"radio",
             "options":["old ghost cell approach","new ghost cell approach"],"default":"new ghost cell approach"},
            {"key":"normal_ghost_velocity_interpolation","question":"normal_ghost_velocity_interpolation",
             "input_type":"radio","options": ["2nd order", "4th order"],"default":"4th order"},
            {"key":"lrough","question":"lrough:","input_type":"radio","options":["True","False"],"default":"False"},
            {"key":"pl_ex","question":"pl_ex (# of extra ghost planes. pl_ex=0 -> blocks have only 1 ghost layer):", 
             "input_type":"int","default":2},
            {"key":"th","question":"Th:","input_type":"float","default":293},
            {"key":"tc","question":"Tc:","input_type":"float","default":293},
            {"key":"tinit","question":"Tinit:","input_type":"float","default":273}
        ]
    },
    {
        "section_name": "Energy Boundary Conditions",
        "questions": [
            {"key":"west_temp_boundary","question":"West Temp Boundary Conditions:","input_type":"dropdown",
             "options": ["Periodic", "Adiabatic","Cold Surface","Hot Surface"],"default":"Adiabatic"},
            {"key":"east_temp_boundary","question":"East Temp Boundary Conditions:","input_type":"dropdown",
             "options": ["Periodic", "Adiabatic","Cold Surface","Hot Surface"],"default":"Adiabatic"},
            {"key":"south_temp_boundary","question":"South Temp Boundary Conditions:","input_type":"dropdown",
             "options": ["Periodic", "Adiabatic","Cold Surface","Hot Surface"],"default":"Adiabatic"},
            {"key":"north_temp_boundary","question":"North Temp Boundary Conditions:","input_type":"dropdown",
             "options": ["Periodic", "Adiabatic","Cold Surface","Hot Surface"],"default":"Adiabatic"},
            {"key":"bottom_temp_boundary","question":"Bottom Temp Boundary Conditions:","input_type":"dropdown",
             "options": ["Periodic", "Adiabatic","Cold Surface","Hot Surface"],"default":"Cold Surface"},
            {"key":"top_temp_boundary","question":"Top Temp Boundary Conditions:","input_type":"dropdown",
             "options": ["Periodic", "Adiabatic","Cold Surface","Hot Surface"],"default":"Hot Surface"}
        ]
    },
    {
        "section_name": "Time Series",
        "questions": [
            {"key":"num_time_series_points","question":"num of time series points:","input_type": "int","default":1},
            {"key":"block_pos_1","question":"Block number position #1:","input_type": "int","default":1},
            {"key":"block_pos_2","question":"Block number position #2:","input_type": "int","default":10},
            {"key":"block_pos_3","question":"Block number position #3:","input_type": "int","default":10},
            {"key":"block_pos_4","question":"Block number position #4:","input_type": "int","default":10}
        ]
    }
]

#General LPT sections control the first few parameters in LPT file
lpt_sections = [
    {
        "section_name": "General Parameters",
        "questions": [
            {"key":"lagrangian_num",
             "question":"How many lagrangian fractions do you want? (Press tab or click on the next box to see fraction options)",
             "input_type":"int"},
            {"key":"particle_to_particle","question":"Do you want particle to particle collisions on or off?",
             "input_type":"radio","options":['on','off']},
            {"key":"particle_to_wall","question":"Do you want particle to wall collisions on or off?",
             "input_type":"radio","options":['on','off']},
            {"key":"PSI_cell_ball","question":"PSI cell or PSI ball?","input_type":"radio","options":["PSIcell","PSIball"],
             "default":"PSIball"},
            {"key":"delta_function","question":"Delta Function:","input_type":"int","default":5},
            {"key":"k_n","question":"Stiffness Constant (k_n):","input_type":"float","default":0.00001}
        ]
    }
]

# Fraction Sections in LPT are dynamic, where depending on the user choice of lagrangian fractions, the form will be created 
# multiple times

fraction_lpt_sections = [
    {
        "section_name": "Fraction",
        "questions": [
            {"key":"release_frequency","question":"What release frequency do you want?","input_type":"int"},
            {"key":"particles_per_release","question":"How many particles per release do you want?","input_type":"int"},
            {"key":"particle_diameter","question":"Particle Diameter:","input_type":"float","default":0.005},
            {"key":"std_dev_diameter","question":"Diameter Standard Deviation:","input_type":"float","default":'0.0'},
            {"key":"particle_density","question":"Particle Density:","input_type":"int","default":1320},
            {"key":"std_dev_density","question":"Density Standard Deviation:","input_type":"int","default":100},
            {"key":"release_volume_lx","question":"Release Volume Lx:","input_type":"float","default":0.2},
            {"key":"release_volume_ly","question":"Release Volume Ly:","input_type":"float","default":0.2},
            {"key":"release_volume_lz","question":"Release Volume Lz:","input_type":"float","default":0.05},
            {"key":"spherical_volume_release","question":"Spherical volume of release?","input_type":"radio","options":
             ["True","False"],"default":"False"},
            {"key":"generate_on_surface","question":"Generate on Surface?","input_type":"radio","options":["True","False"],
             "default":"False"},
            {"key":"radius","question":"Radius in m:","input_type":"float","default":0.1},
            {"key":"spherical_release_options","question":"Which spherical release option do you choose?",
             "input_type":"dropdown","options":["3D","2D plane-XY","2D plane-ZY","2D plane-ZX"],"default":"2D plane-XY"},
            {"key":"random_release","question":"Random release? (If choosing True you must manually input coordinates)",
             "input_type":"radio","options":["True","False"],"default":"True"},
            {"key":"xp","question":"xp:","input_type":"float","default":0.1},
            {"key":"yp","question":"yp:","input_type":"float","default":0.125},
            {"key":"zp","question":"zp:","input_type":"float","default":0.125},
            {"key":"up","question":"up:","input_type":"float","default":'0.0'},
            {"key":"vp","question":"vp:","input_type":"float","default":'0.0'},
            {"key":"wp","question":"wp:","input_type":"float","default":-0.05},
        ]
    }
]

#Infodom sections contain 9 questions
infodom_sections = [
    {
        "section_name": "Infodom Parameters",
        "questions": [
            {"key": "grid_size_x","question": "What is the grid size for x?","input_type":"float"},
            {"key": "grid_size_y","question": "What is the grid size for y?","input_type":"float"},
            {"key": "grid_size_z","question": "What is the grid size for z?","input_type":"float"},
            {"key": "domain_length","question": "How long is your domain in meters?","input_type":"float"},
            {"key": "domain_width","question": "How wide is your domain in meters?","input_type":"float"},
            {"key": "domain_height","question": "How high is your domain in meters?","input_type":"float"},
            {"key": "sub_domain_i","question": "How many sub-domains do you want in the x direction?","input_type":"int"},
            {"key": "sub_domain_j","question": "How many sub-domains do you want in the y direction?","input_type":"int"},
            {"key": "sub_domain_k","question": "How many sub-domains do you want in the z direction?","input_type":"int"},
        ]
    }
]

#Mdmap only requires one question, but is based on the dependency total_block_num which is from infodom. If the user tries to 
#create MDmap before infodom, an error message will appear
mdmap_sections = [
    {
        "section_name":"Processors Assignment",
        "questions": [
            {"key":"num_processors","question":"How many processors do you want to assign to each domain? (Default of 1).",
             "input_type":"int","default":1}
        ]
    }
]

# Some values are chosen from dropdown menus but need to be entered in the files as integers. This dictionary provides a mapping to convert the user choices to the appropriate integers

#value_mapping so user input is entered in files correctly
mapping_dict = {"Exp.Euler":1,
                        "AdamsBashfort":2,
                        "RK2":3,
                        "RK3":4,
                        "Imp.Euler":1,
                        "CrankNicholson":2,
                        "ExplicitRungeKuttaDiff":3,
                        "2nd order CDS":1,
                        "4th order CDS":2,
                        "WENO":3,
                        "sip":1,
                        "mg":2,
                        "GSS":1,
                        "TDMA":2,
                        "RBGSS":3,
                        "variable":"T",
                        "fixed":"F",
                        "inflow - uniform":1,
                        "inflow - 1/7th power law":12,
                        "inflow - prescribed":7,
                        "inflow - SEM":8,
                        "outflow - Neumann":2,
                        "outflow - convective":21,
                        "slip":3,
                        "wall - no-slip":4,
                        "periodic condition":5,
                        "wall - smooth log law":61,
                        "wall - rough log law":62,
                        "1/6th law":63,
                        "1/7th law":64,
                        "1/8th law":65,
                        "manning":"n",
                        "equivalent sand":"k",
                        "Uniform":1,
                        "1/7th PL":12,
                        "Smagorinsky":1,
                        "WALE":2,
                        "OneEqnModel":3,
                        "k-eps model (RANS)":4,
                        "old ghost cell approach":1,
                        "new ghost cell approach":2,
                        "2nd order":1,
                        "4th order":2,
                        "Periodic":5,
                        "Adiabatic":7,
                        "Cold Surface":8,
                        "Hot Surface":9,
                        "HJ_ENO":11,
                        "HJ_WENO":22,
                        "YES":"F",
                        "NO":"T",
                        "yes":"T",
                        "no":"F",
                        "on":"T",
                        "off":"F",
                        "True":"T",
                        "False":"F",
                        "PSIcell":"T",
                        "PSIball":"F",
                        "3D":0,
                        "2D plane-XY":1,
                        "2D plane-ZY":2,
                        "2D plane-ZX":3
                       }


# Step 2: Dataframes are created by collecting user input, adding it to a dictionary, passing it through the mapping dictionary, and following appropriate formats for each tab. The dataframes are returned as lists to be input into a file creation function.

#Create control dataframes
def create_control_dfs():
    
    #logic for density and kinematic viscosity
    density = 0
    kinematic_visc = 0
    if (control_entries['fluid'] == "water"):
        density = 1000
        kinematic_visc = 10e-6
    elif (control_entries['fluid'] == "air"):
        density = 1.2
        kinematic_visc = 1.6e-5

    #All dataframes for control using values from advanced and regular sections
    
    #NUMERIC PARAMETERS
    numeric_params = pd.DataFrame({
        "Value": [
            f"{control_entries['keyword']} {control_entries['ubulk']}",
            f"{control_entries['dx']} {control_entries['dy']} {control_entries['dz']}",
            f"{density} {kinematic_visc} {control_entries['pr']} {control_entries['turb_schmidt']} {control_entries['beta']}",
            f"{control_entries['gx']} {control_entries['gy']} {control_entries['gz']}",
            f"{control_entries['convection_scheme']}",
            f"{control_entries['diffusion_scheme']}",
            f"{control_entries['differencing']}",
            f"{control_entries['solver']}",
            f"{control_entries['multigrid_step']} {control_entries['multigrid_iteration_scheme']}",
            f"{control_entries['multigrid_max_iter']} {control_entries['restriction_iter']} {control_entries['prolongation_iter']}",
            f"{control_entries['dt']} {control_entries['variable_fixed']} {control_entries['sweeps']} {control_entries['safety_factor']}",
            f"{control_entries['itime_end']} {control_entries['restart']} {control_entries['reinitmean']} {control_entries['n_out']}",
            f"{control_entries['l_transient']} {control_entries['results']}",
            f"{control_entries['niter']} {control_entries['eps']} {control_entries['nswp_1']} {control_entries['nswp_2']} {control_entries['nswp_3']} {control_entries['nswp_4']}"
        ],
        "Parameter": [
            f"{'Keyword'}, {'Ubulk'}",
            f"{'dx'},{'dy'},{'dz'}",
            f"{'dens'}, {'kinematic visc'}, {'Pr'}, {'turb Schmidt'}, {'beta'}",
            f"{'gx'},{'gy'},{'gz'}",
            f"{'convection_scheme(1=Exp.Euler,2=AdamsBashfort,3=RK2,4=RK3)'}",
            f"{'diffusion_scheme(1=Imp.Euler,2=CrankNicholson,3=ExplicitRungeKuttaDiff)'}",
            f"{'differencing(1=2ndCDS,2=4thCDS,3=WENO)'}",
            f"{'solver(1=sip,2=mg)'}",
            f"{'multigrid step (ngrid)'}, {'multigrid iteration scheme(1=GSS,2=TDMA,3=RBGSS)'}",
            f"{'multigrid maximum iteration per time step'}, {'restriction iter'}, {'prolongation iter'}",
            f"{'dt'}, {'variable(T)/fixed(F) dt'}, {'sweeps'}, {'safety_factor'}",
            f"{'itime_end'}, {'restart'}, {'reinitmean'}, {'n_out'}",
            f"{'LTRANSIENT'}, {'results output'}",
            f"{'niter'}, {'eps'}, {'nswp(1)'}, {'nswp(2)'}, {'nswp(3)'}, {'nswp(4)'}",
        ]
    })
    
    # FLOW BOUNDS
    flow_bounds = pd.DataFrame({
        "Value": [
            f"{control_entries['west_boundary_condition']}",
            f"{control_entries['east_boundary_condition']}",
            f"{control_entries['south_boundary_condition']}",
            f"{control_entries['north_boundary_condition']}",
            f"{control_entries['bottom_boundary_condition']}",
            f"{control_entries['top_boundary_condition']}",
            f"{control_entries['friction_coeff']} {control_entries['friction_coeff_value']}",
            f"{control_entries['save_inflow_data']} {control_entries['inlet_num']}"
        ],
        "Parameter": [
            f"{'West Boundary Condition (1= Inflow; 12=1/7th power law; 7=prescribed inflow; 8=SEM)'}",
            f"{'East Boundary Condition (2= Outflow, 2=NeumannBC(NBC) 21= ConvectiveBC(CBC))'}",
            f"{'South Boundary Condition (3= Slip Condition)'}",
            f"{'North Boundary Condition (4= No-Slip Condition)'}",
            f"{'Bottom Boundary Condition (5= Periodic Condition)'}",
            f"{'Top Boundary Condition (61=smooth log law; 62=rough log law; 63=1/6th law; 64=1/7th law; 65=1/8th law)'}",
            f"{'Friction coefficient (n:manning;k:equivalent sand; only if bc=62)'}",
            f"{'save inflow data (precursor sim.)'} ; {'number of inlets'}"
        ]
    })
    
    # SYNTHETIC EDDY
    synthetic_eddy = pd.DataFrame({
        "Value": [
            f"{control_entries['velocity_profile']}",
            f"{control_entries['turbulence_intensity']}",
            f"{control_entries['inlet_profile_num']}"
        ],
        "Parameter": [
            f"{'Velocity profile:1=Uniform; 12=1/7th PL'}",
            f"{'Turbulence intensity'}",
            f"{'Number inlet profiles'}"
        ]
    })
    
    # MODEL OPTIONS
    t_start_averaging_2 = float(control_entries['ctime_averaging'])-1
    SGS_model = "T"
    if (control_entries['SGS_model_choice']=="None"):
        SGS_model = "F"
    model_opts = pd.DataFrame({
        "Value": [
            f"{control_entries['time_averaging']} {control_entries['ctime_averaging']} {t_start_averaging_2} {control_entries['noise']}",
            f"{SGS_model} {control_entries['SGS_model_choice']}",
            f"{control_entries['LMR']} {control_entries['normal_ghost_velocity_interpolation']}",
            f"{control_entries['limb']} {control_entries['lenergy']} {control_entries['lrough']}",
            f"{control_entries['lpt']} {control_entries['LSM']} {control_entries['L_LSMbase']}",
            f"{control_entries['lscalar']} {control_entries['l_active_scalar']} {control_entries['l_non_newt']}",
            f"{control_entries['pl_ex']}",
            f"{control_entries['th']} {control_entries['tc']} {control_entries['tinit']}"
        ],
        "Parameter": [
            f"{'time_averaging'}, {'t_start_averaging1'}, {'t_start_averaging2'}, {'noise'}",
            f"{'SGS-model,1=Smagorinsky,2=WALE,3=OneEqnModel,4=k-eps model (RANS)'}",
            f"{'LMR (1=old ghost cell approach, 2=new ghost cell approach)'}, {'normal ghost velocity interpolation(1=2nd-order,2=4th-order)'}",
            f"{'LIMB'},{'LENERGY'},{'LROUGH'}",
            f"{'LPT'}, {'LSM'}, {'L_LSMbase'}",
            f"{'LSCALAR'}, {'LActiveScalar'} {'LNonNewt'}",
            f"{'pl_ex (# of extra ghost planes. pl_ex=0 -> blocks have only 1 ghost layer)'}",
            f"{'Th'}, {'Tc'}, {'Tinit'}"
        ]
    })
    
    # ENERGY BOUNDS
    energy_bounds = pd.DataFrame({
        "Value": [
            f"{control_entries['west_temp_boundary']}",
            f"{control_entries['east_temp_boundary']}",
            f"{control_entries['south_temp_boundary']}",
            f"{control_entries['north_temp_boundary']}",
            f"{control_entries['bottom_temp_boundary']}",
            f"{control_entries['top_temp_boundary']}"
        ],
        "Parameter": [
            f"{'West Temp Boundary Condition (5= Periodic)'}",
            f"{'East Temp Boundary Condition (7= Adiabatic)'}",
            f"{'South Temp Boundary Condition (8= Cold Surface)'}",
            f"{'North Temp Boundary Condition (9= Hot Surface)'}",
            f"{'Bottom Temp Boundary Condition'}",
            f"{'Top Temp Boundary Condition'}"
        ]
    })
    
    # TIME SERIES
    time_series = pd.DataFrame({
        "Value": [
            f"{control_entries['num_time_series_points']}",
            f"{control_entries['block_pos_1']} {control_entries['block_pos_2']} {control_entries['block_pos_3']} {control_entries['block_pos_4']}"
        ],
        "Parameter": [
            f"{'num of time series points'}",
            f"{'point #1'}"
        ]
    })

    return [numeric_params, flow_bounds, synthetic_eddy, model_opts, energy_bounds,time_series]


# LPT is a special case for dataframe creation, since it requires a dynamic number of dataframes based on the user's choice of lagrangian fractions. When the user presses submit, dataframes are created depending on the number of fractions in the form.

#create LPT dataframe
def create_lpt_dfs():

    #General Parameters
    lpt_general = pd.DataFrame({"Values":[f"{lpt_entries['PSI_cell_ball']} {lpt_entries['delta_function']}",
                                          f"{lpt_entries['lagrangian_num']}",f"{lpt_entries['particle_to_particle']}",
                                         f"{lpt_entries['particle_to_wall']}",f"{lpt_entries['k_n']}"],
                           "Parameters":[f"{'PSIcell(T)/PSIball(F)'}, {'delta function'}",
                           f"{'number of Lagrangian fractions'}",f"{'Inter-particle collisions'}",
                           f"{'Particle-wall collisions'}", f"{'k_n'}"]})



    lagrangian_num = int(lpt_entries['lagrangian_num'])
    duplicated_dfs = [lpt_general]  # List to store the duplicated dataframes
    
    # Create new DataFrames based on lagrangian_num
    for i in range(1, lagrangian_num + 1):
        lpt_fraction = pd.DataFrame({"Values":[f"{lpt_entries[f'release_frequency_{i}']}",f"{lpt_entries[f'particles_per_release_{i}']}",
                                      f"{lpt_entries[f'particle_diameter_{i}']} {lpt_entries[f'std_dev_diameter_{i}']}",
                                        f"{lpt_entries[f'particle_density_{i}']} {lpt_entries[f'std_dev_density_{i}']}",
                                      f"{lpt_entries[f'release_volume_lx_{i}']} {lpt_entries[f'release_volume_ly_{i}']} {lpt_entries[f'release_volume_lz_{i}']}",
                                      f"{lpt_entries[f'spherical_volume_release_{i}']} {lpt_entries[f'generate_on_surface_{i}']} {lpt_entries[f'radius_{i}']}",
                                      f"{lpt_entries[f'spherical_release_options_{i}']}",f"{lpt_entries[f'random_release_{i}']}",
                                      f"{lpt_entries[f'xp_{i}']} {lpt_entries[f'yp_{i}']} {lpt_entries[f'zp_{i}']} {lpt_entries[f'up_{i}']} {lpt_entries[f'vp_{i}']} {lpt_entries[f'wp_{i}']}"],
                             "Parameters":[f"{'release frequency (i.e., how many tsteps between releases)'}",
                                          f"{'number of parts/release'}",
                                          f"{'Diameter of particle in meters;'}, {'standard deviation (0 if constant)'}",
                                          f"{'Density of particle in kg/m3;'}, {'standard deviation (0 if constant)'}",
                                          f"{'Volume of release (Lx,Ly,Lz)'}",
                                          f"{'Spherical volume of release'}, {'generate on surface'}, {'radius in m'}",
                                          f"{'Spherical release options (0=3D, 1-3=2D plane --> 1=XY, 2=ZY, 3=ZX plane)'}",
                                          f"{'Random release (if T) within the above volume'}",
                                          f"{'xp'},{'yp'},{'zp'},{'up'},{'vp'},{'wp'}"]})
        
        duplicated_dfs.append(lpt_fraction)

    return duplicated_dfs


# LPT also requires dynamic dictionaries. One section for fractions is not enough, since it depends on user's choice. the create_fractions function is used to dynamically create dictionaries and key names to store user input. When the form is created, it uses create_fractions to dynamically insert questions.

#Create dynamic dictionaries
def create_fractions(sections,lagrangian_num):
    duplicated_sections = []

     #Iterate through the base sections
    for section in sections:
        #Create separate sections for each fraction (based on user input)
        for i in range(1, lagrangian_num + 1):
            #Create a clean list for the questions of the current fraction
            modified_questions = []

            #Modify the keys of each question for the current fraction
            for question in section.get("questions", []):
                new_key = f"{question['key']}_{i}"
                #Create a copy of the original question and update its key
                duplicated_question = question.copy()
                duplicated_question["key"] = new_key
                modified_questions.append(duplicated_question)

            #Create a new section with an updated name and modified questions
            new_section = {
                "section_name": f"{section['section_name']} {i}",
                "questions": modified_questions,
            }

            #Append the new section to the duplicated_sections list
            duplicated_sections.append(new_section)

    return duplicated_sections


# The infodom tab requires some special logic. Since cells per block and cell num must both be integers, this logic iteratively adjusts user choices until it converges on an integer. This process involves rounding decimals in domain and mesh sizes, increasing and decreasing sub-domains, and increasing and decreasing mesh sizes. The process will only change the sub domain block number by up to 5 before quitting. If it converges on an integer for cells per block and cell num, it will ask the user to accept the changes. If not, it will ask the user to re-enter values.

#whole process for finding an integer value
def get_integer_values(var_name,isinteger,og_sub_domain,og_grid_size,grid_size,sub_domain,cell_num,domain_size,cells_per_block):

    #steps to be taken with different parameters each time
    def converge_on_integer(adjust_grid_size=False,increase_sub_domain=False,decrease_sub_domain=False,use_iterator=False):

        global iterator,error_message
        nonlocal var_name,isinteger,grid_size,sub_domain,cell_num,domain_size,cells_per_block

        #adjust grid size up to 1
        while not isinteger:
            if (adjust_grid_size):
                grid_size += 0.001
                cell_num = domain_size/grid_size

            #adjust grid size and increase/decrease sub_domain
            if (increase_sub_domain):
                sub_domain+=1
            if (decrease_sub_domain):
                sub_domain-=1
                if (sub_domain <= 0):
                    error_message = print("Sub-domain cannot be 0.")
                    return error_message
            cells_per_block = cell_num/sub_domain
            
            if (iterator):
                iterator+=1

            #Break loop for errors
            if (sub_domain <= 0):
                error_message = print("Sub-domains cannot be 0 or less.")
                return error_message
            if (grid_size <= 0):
                error_message = print("Grid size cannot be 0 or less.")
                return error_message
            if (cells_per_block >= 100):
                error_message = print("Cells per block cannot be more than 100.")
                return error_message

            #Stop if integer is found, ask user to confirm choices
            if cells_per_block.is_integer():
                if cells_per_block > 0:
                    print("Converged. Found integer.")
                    print("New cells per block:",cells_per_block,"New grid size:",grid_size,"New cell num:",cell_num,
                          "New sub domain:",sub_domain,"New domain length:",domain_size)
                    # Display confirmation dialog box
                    response = messagebox.askyesno(
                        "Confirm Changes",
                        f"The values you entered could not be used for {var_name}. Values have been updated:\n\n"
                        f"Old Grid Size: {og_grid_size}\nNew Grid Size: {grid_size}\n"
                        f"Old Sub-domain: {og_sub_domain}\nNew Sub-domain: {sub_domain}\n\n"
                        f"Do you want to save these changes?",
                    )

                    #Allow user to choose yes or no in confirmation (if no, it will try the process until failure)
                    if response:
                        print("Changes confirmed by user. Saving new values.")
                        isinteger = True
                        return cells_per_block
                    else:
                        print("Changes rejected by user. Reverting to original values.")
                        return None
                        isinteger = True
                        
                elif cells_per_block <= 0:
                    error_message = print("Cells per block cannot be 0. Try again.")
                    return error_message
            
            if (decrease_sub_domain):
                if (sub_domain == 1):
                    eror_message = print("Number of blocks cannot be 0.")
                    return error_message
                if (iterator == 6):
                    error_message = "Tried up to 5 blocks less."
                    return error_message
                    
            #must choose between 5 blocks less or 5 blocks more
            if (increase_sub_domain):
                if (iterator == 6):
                    error_message = print("Tried up to 5 blocks more.")
                    return error_message

            if (adjust_grid_size):
                if (grid_size == 1):
                    error_message = print("Mesh size too large.")
                    return error_message
            if (cells_per_block < 1):
                error_message = print("Cell num too low.")
                return error_message
    
        return grid_size,cell_num,cells_per_block,sub_domain

    #If adjusting mesh size and block num doesn't work, the program will round domain sizes and mesh sizes
    def rounding(domain=False,grid_size_10ths=False,grid_size_100ths=False):
    
        global iterator
        nonlocal isinteger,grid_size,sub_domain,cell_num,domain_size,cells_per_block
        
        if (not (cells_per_block.is_integer())):
    
            if domain:
                domain_size = round(domain_size)
                print("Rounded domain size:",domain_size)
            if grid_size_100ths:
                grid_size = math.ceil(og_grid_size * 100) / 100
                print("Rounded grid size:",grid_size)
            if grid_size_10ths:
                grid_size = math.ceil(og_grid_size * 10) / 10
                print("Rounded grid size:",grid_size)
                
            try:
            # Calculate number of cells for the whole domain
                cell_num_x = domain_size / grid_size
                print("Cell_num after rounding", cell_num)
            
                # Check if the result is less than or equal to zero
                if grid_size <= 0:
                    raise ValueError("Cell number cannot be zero or negative.")
            
                # Calculate cells per block
                cells_per_block = cell_num / sub_domain
                print("Cells per block after rounding:", cells_per_block)
            
                # Check if cells per block is less than or equal to zero
                if sub_domain <= 0:
                    raise ValueError("Cells per block cannot be zero or negative.")

            #control for errors
            except ZeroDivisionError:
                print("Error: Division by zero. Ensure grid size and subdomain values are valid.")
            except ValueError as ve:
                print(f"Error: {ve}")
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
    
        return domain_size,grid_size
    
    #print statements for debugging. These will not be visible to the user.
    process_statements = ["Trying to increase number of sub-domains.","Trying to number of sub-domains.",
                          "Trying to increase number of sub-domains and update mesh size.",
                          "Trying to decrease number of sub-domains and update mesh size.","Rounding Domain size.",
                          "Trying to increase number of sub-domains and change domain size.",
                          "Trying to decrease number of sub-domains and change domain size.","Rounding mesh size.",
                          "Trying to increase number of sub-domains and round grid size to 100ths place.",
                          "Trying to decrease number of sub-domains and round grid size to 100ths place.",
                          "Rounding mesh size","Trying to increase number of sub-domains and round grid size to 10ths place.",
                          "Trying to decrease number of sub-domains and round grid size to 10ths place."]

    #functions for doing different steps in the process
    function_calls = [lambda: converge_on_integer(increase_sub_domain=True,use_iterator=True),
                      lambda: converge_on_integer(decrease_sub_domain=True,use_iterator=True),
                      lambda: converge_on_integer(adjust_grid_size=True,increase_sub_domain=True),
                      lambda: converge_on_integer(adjust_grid_size=True,decrease_sub_domain=True),
                      lambda: rounding(domain=True), lambda: converge_on_integer(increase_sub_domain=True,use_iterator=True),
                      lambda: converge_on_integer(decrease_sub_domain=True,use_iterator=True),
                      lambda: rounding(grid_size_100ths=True), 
                      lambda: converge_on_integer(increase_sub_domain=True,use_iterator=True),
                      lambda: converge_on_integer(decrease_sub_domain=True,use_iterator=True),
                      lambda: rounding(grid_size_10ths=True),
                      lambda: converge_on_integer(increase_sub_domain=True,use_iterator=True),
                      lambda: converge_on_integer(decrease_sub_domain=True,use_iterator=True)]

    #function to try all the steps in order
    def do_steps(print_statement,function_call):
        global iterator
        nonlocal isinteger, grid_size, sub_domain
        iterator = 0
        grid_size = og_grid_size
        sub_domain = og_sub_domain
        print(print_statement)
        result = function_call()
        
        # If an integer is found, return success to terminate the loop
        if isinteger:
            return True  # Convergence achieved
        return False  # Continue loop

    #do steps
    for i, (statement, func) in enumerate(zip(process_statements, function_calls)):
        if do_steps(statement, func):  # Pass each statement and function
            print(f"Process completed successfully after step {i + 1}.")
            break  # Exit loop if successful
            
    #display error message if convergence is not successful       
    else:
        error_msg = f"The values you entered could not be used for {var_name}: cells per block must be a whole number. Please try again."
        error_messages.append(error_msg)
        print("Failed to converge on an integer after all steps. Please try different values.")


    return grid_size,cell_num,cells_per_block,sub_domain,domain_size,isinteger


# Step 3: File creation. Since infodom and mdmap are small, their dataframe creation and file creation are handled in one function. The infodom entries are used to make calculations, check for integer as above, and finally processed into coordinates for the final file. Mdmap uses two values to create lists. LPT and Control use the same file creation functions to create headers and data is imported from dataframes. All files are saved to .cin based on user's choice.

#create the infodom file using user input to make calculations for blocks, find coordinates for each block
def create_infodom(filename,overwrite=True):
    
    #access user inputs
    domain_length = float(infodom_entries['domain_length'])
    domain_width = float(infodom_entries['domain_width'])
    domain_height = float(infodom_entries['domain_height'])
    grid_size_x = float(infodom_entries['grid_size_x'])
    grid_size_y = float(infodom_entries['grid_size_y'])
    grid_size_z = float(infodom_entries['grid_size_z'])
    sub_domain_i = int(infodom_entries['sub_domain_i'])
    sub_domain_j = int(infodom_entries['sub_domain_j'])
    sub_domain_k = int(infodom_entries['sub_domain_k'])
    
    try:
    #Calculate number of cells for the whole domain
        cell_num_x = domain_length / grid_size_x
        cell_num_y = domain_width / grid_size_y
        cell_num_z = domain_height / grid_size_z
    
        #Check if any of the values are less than or equal to zero
        if cell_num_x <= 0 or cell_num_y <= 0 or cell_num_z <= 0:
            raise ValueError("Cell count cannot be zero or negative.")
    
        #Calculate cells per block
        cells_per_block_x = cell_num_x / sub_domain_i
        cells_per_block_y = cell_num_y / sub_domain_j
        cells_per_block_z = cell_num_z / sub_domain_k
    
        #Check if the cells per block are less than or equal to zero
        if cells_per_block_x <= 0 or cells_per_block_y <= 0 or cells_per_block_z <= 0:
            raise ValueError("Cells per block cannot be zero or negative.")
            
    except ZeroDivisionError:
        messagebox.showerror(
        "Division by Zero Error",
        "Error: Division by zero occurred. Ensure that grid size and subdomain values are valid File not created."
    )
        print("Error: Division by zero. Ensure grid size and subdomain values are valid.")
        
    except ValueError as ve:
        messagebox.showerror(
        "Value Error",
        f"Error: {ve}"
    )
        print(f"Error: {ve}")

    except Exception as e:
        messagebox.showerror(
            "Unexpected Error",
            f"An unexpected error occurred: {e}"
        )
        print(f"An unexpected error occurred: {e}")

    #check if user input results in integers
    def check_integer(cells_per_block):
        if cells_per_block.is_integer():
            isinteger = True
        else:
            isinteger = False
            
        return isinteger

    #if integers are not achieved, call function to try to adjust user input until integer convergence
    isinteger_x = check_integer(cells_per_block_x)
    isinteger_y = check_integer(cells_per_block_y)
    isinteger_z = check_integer(cells_per_block_z)

    if (isinteger_x == False):
        grid_size_x,cell_num_x,cells_per_block_x,sub_domain_i,domain_length,isinteger_x = get_integer_values('x',isinteger_x,
                                                                                             sub_domain_i,grid_size_x,
                                                                                             grid_size_x,sub_domain_i,
                                                                                             cell_num_x,domain_length,
                                                                                             cells_per_block_x)
    
    if (isinteger_y == False):
        grid_size_y,cell_num_y,cells_per_block_y,sub_domain_j,domain_width,isinteger_y = get_integer_values('y',isinteger_y,
                                                                                                sub_domain_j,
                                                                                                grid_size_y,grid_size_y,
                                                                                                sub_domain_j,cell_num_y,
                                                                                                domain_width,cells_per_block_y)
    if (isinteger_z == False):
        grid_size_z,cell_num_z,cells_per_block_z,sub_domain_k,domain_height,isinteger_z = get_integer_values('z',isinteger_z,sub_domain_k,
                                                                                                 grid_size_z,grid_size_z,
                                                                                                 sub_domain_k,cell_num_z,
                                                                                                 domain_height,cells_per_block_z)
    
    if not (isinteger_x==True and isinteger_y==True and isinteger_z==True):
        print("Integer validation failed. File will not be written.")
        return

    print("Creating infodom df")
    
    #create df with 6 columns and 4 rows from infodom.cin example
    infodom_initial_df = pd.DataFrame({'Grid Size': [grid_size_x,grid_size_y,grid_size_z,np.nan],
                                       'Overall Domain Size': [domain_length,domain_width,domain_height,np.nan],
                                       'Number of Cells for Whole Domain': [cell_num_x,cell_num_y,cell_num_z,cell_num_x*cell_num_y*cell_num_z],
                                       'Sub-domains': [sub_domain_i,sub_domain_j,sub_domain_k,sub_domain_i*sub_domain_j*sub_domain_k],
                                       'Cells per Block': [cells_per_block_x,cells_per_block_y,cells_per_block_z,cells_per_block_x*cells_per_block_y*cells_per_block_z]})

    print(infodom_initial_df)
    
    #individual block dimension calculation
    block_dimension_x = domain_length/sub_domain_i
    block_dimension_y = domain_width/sub_domain_j
    block_dimension_z = domain_height/sub_domain_k
    
    #get total block num
    global total_block_num 
    total_block_num = sub_domain_i*sub_domain_j*sub_domain_k
    
    #column 1 is a list iterating total_block_num times
    domain = list(range(total_block_num))
    rdiv = [1]*total_block_num
    print(block_dimension_x)

    for floats in domain:
        int(floats)
    print("Domain after casting:",domain)
    for floats in rdiv:
        int(floats)
    print("rdiv after casting:",rdiv)
    
    #get initial coordinates for infodom
    def get_coords(sub_domain,block_dimension):
        counter = 0
        v1 = []
        v2 = []
        for i in range(sub_domain):
            v1.append(round(counter,2))
            counter = counter+block_dimension
            v2.append(round(counter,2))
    
        return v1,v2
    
    #repeat coordinates for infodom
    x1,x2 = get_coords(sub_domain_i,block_dimension_x)
    y1,y2 = get_coords(sub_domain_j,block_dimension_y)
    z1,z2 = get_coords(sub_domain_k,block_dimension_z)
    
    #Repeat coordinates for all blocks so they are sitting in different spots in the 3D space
    x1_repeated = x1 * (len(y1) * len(z1))
    x2_repeated = x2 * (len(y1) * len(z1))
    
    y1_repeated = [y for _ in range(len(z1)) for y in y1 for _ in range(len(x1))]
    y2_repeated = [y for _ in range(len(z1)) for y in y2 for _ in range(len(x1))]
    
    z1_repeated = [z for z in z1 for _ in range(len(y1) * len(x1))]
    z2_repeated = [z for z in z2 for _ in range(len(y1) * len(x1))]
    
    infodom_coord_df = pd.DataFrame({'Domain':domain,'rdiv':rdiv,'x1':x1_repeated,'x2':x2_repeated,
                                    'y1':y1_repeated,'y2':y2_repeated,'z1':z1_repeated,'z2':z2_repeated})
    print(infodom_coord_df)
    
    #Hard-code the length of the equal signs line to 49 to create headers
    equals_line = "=" * 49
    header = f"{' '*8}{total_block_num} number of domains\n"
    footer = f"{' '*8}{sub_domain_i} number of divisions in i\n{' '*8}{sub_domain_j} number of divisions in j\n{' '*8}{sub_domain_k} number of divisions in k\n"

    #Ask user for file name
    if not overwrite and os.path.exists(filename):
        print(f"File '{filename}' already exists. Overwrite is disabled.")
        return
        
    #Writing to the .cin file
    try:
        with open(filename, "w") as f:
            
            #Write the header
            f.write(header)
            f.write(equals_line + "\n")
            #Write the DataFrame as a space-separated format with no column or row names, with centered rows
            for index, row in infodom_coord_df.iterrows():
                f.write(f"{int(row['Domain'])}  "
                        f"{int(row['rdiv'])}  "
                        f"{row['x1']}  "
                        f"{row['x2']}  "
                        f"{row['y1']}  "
                        f"{row['y2']}  "
                        f"{row['z1']}  "
                        f"{row['z2']}\n")
                print(row)
            #Add the line of equals signs at the end of the file
            f.write(equals_line + "\n")
            f.write(footer)
    
        #If the file is written successfully, print a success message
        print(f"File '{'infodom_test.cin'}' created successfully!")
    
    except Exception as e:
        #Handle any errors that occur during file creation
        print(f"An error occurred while creating the file: {e}")
        
    return 


#Create MDmap file
def create_mdmap_file_with_columns(filename, overwrite=True):
    
    global total_block_num

    #Display error if the user has not submitted the infodom file first
    if 'total_block_num' not in globals() or total_block_num is None:
        error_msg = f"Number of domains is not defined. Please ensure you have created the infodom file before submitting this tab."
        error_messages.append(error_msg)
        return

    try:
        #Access the processor count
        num_processors = int(mdmap_entries['num_processors'])
        print("Processors assigned.")
    
        #Create processor_num column for file
        num_processor_list = []
        current_processor_num = 0
        while current_processor_num < total_block_num:
            num_processor_list.append(num_processors)
            current_processor_num += 1

        #calculate processors 
        total_processors_num = total_block_num * num_processors
        domains_processors_ID_list = list(range(total_processors_num))
        ID_list = [
            " ".join(map(str, domains_processors_ID_list[i:i + num_processors]))
            for i in range(0, len(domains_processors_ID_list), num_processors)
        ]

        #create dataframe
        domain_processor_df = pd.DataFrame({'Sub-Domain ID': ID_list,
                                            'Num Processors': num_processor_list, 
                                            'Processor ID': ID_list})

        #Hard-code the length of the equal signs line to 49 for header
        equals_line = "=" * 49
        header = f"{' '*8}{total_block_num} number of domains\n{' '*8}{total_processors_num} number of processors\n{equals_line}\n"

        #Ask user for file name
        if not overwrite and os.path.exists(filename):
            print(f"File '{filename}' already exists. Overwrite is disabled.")
            return
        
        #write file
        with open(filename, "w") as f:
            f.write(header)
            for index, row in domain_processor_df.iterrows():
                f.write(f"{str(row['Sub-Domain ID']).rjust(4)}  "
                        f"{str(row['Num Processors']).rjust(4)}  "
                        f"{str(row['Processor ID']).rjust(4)}\n")
            f.write(equals_line + "\n")
        
        print(f"File '{filename}' created successfully!")
    
    except Exception as e:
        # Catch all exceptions and print error message
        errror_msg = f"Error encountered: {str(e)}"
        error_messages.append(error_msg)
        print("Error message shown. File not created.")


#concat dfs into single .cin file - call for control.cin and lpt.cin
def create_file(dfs,headers,filename,overwrite=True):

    #Ask user for file name
    if not overwrite and os.path.exists(filename):
        print(f"File '{filename}' already exists. Overwrite is disabled.")
        return
        
    #Function to format rows with alignment
    def format_aligned_row(values, names, spacing=30):
        value_str = ' '.join(f"{val:<{spacing}}" for val in values)  # Align values
        name_str = ' '.join(f"{name:<{spacing}}" for name in names)  # Align names
        return f"{value_str} {name_str}"

    #Writing to the .cin file
    with open(filename, "w") as f:
        for header, df in zip(headers, dfs):
            # Write the header
            f.write(header + "\n")
        
            # Write each row of the DataFrame
            for index, row in df.iterrows():
                values = row[:-1].values  
                names = row[-1:].values 
                formatted_row = format_aligned_row(values, names)
                f.write(formatted_row + "\n")


# Headers saved in this list are used for the create file for Control and LPT tabs only. LPT also requires dynamic headers for fractions, which are handled later.

# +
control_headers = [
        "=" * 30 + "Numeric Parameters" + "=" * 30,
        "=" * 27 + "Flow Boundary Conditions" + "=" * 26,
        "=" * 28 + "Synthetic Eddy Method" + "=" * 28,
        "=" * 30 + "Modelling Options" + "=" * 30,
        "=" * 25 + "Energy Boundary Conditions" + "=" * 26,
        "=" * 33 + "Time Series" + "=" * 33
    ]

lpt_headers = [
        "#" * 30 + "General Parameters" + "#" * 30,
    ]


# -

# Step 4: Front end. Functions are defined to create sections and questions. Sections are first created by iterating through each dictionary. Create_question handles form appearance depending on question types and adds user input to the appropriate entries dictionary. Create_section creates headers for each section within the dictionary and adds the questions to each section.

# +
#Function to handle creating a question of any type
def create_question(key,question_text, input_type, options, parent_frame, default_value,entries):
    label = ctk.CTkLabel(parent_frame, text=question_text, anchor="w").pack(pady=5, anchor="w")
    
    if input_type == "dropdown":
        var = ctk.StringVar(value=default_value if default_value else options[0])
        dropdown = ctk.CTkOptionMenu(parent_frame, variable=var, values=options)
        dropdown.pack(pady=5, anchor="w")
        entries[key] = var

    elif input_type == "radio":
        var = ctk.StringVar(value=default_value if default_value else options[0])
        radio_button_1 = ctk.CTkRadioButton(parent_frame, text=options[0], variable=var, value=options[0])
        radio_button_2 = ctk.CTkRadioButton(parent_frame, text=options[1], variable=var, value=options[1])
        radio_button_1.pack(pady=5, anchor="w")
        radio_button_2.pack(pady=5, anchor="w")
        entries[key] = var

    elif input_type in ["int","float"]:
        entry = ctk.CTkEntry(parent_frame,placeholder_text=f"Enter {input_type}")
        if default_value:
            entry.insert(0, default_value)
        entry.pack(pady=5, anchor="w")
        entries[key] = entry
        
        #Attach a listener for lagrangian_num, updates form on LPT to add fraction inputs for number of fractions user chooses
        if key == "lagrangian_num":
            entry.bind("<FocusOut>", lambda event: update_lpt_sections(entries, parent_frame))


#Function to handle creating a section with its header and questions
def create_section(section, parent_frame,entries):
    
    # Create a frame for the section
    section_frame = ctk.CTkFrame(parent_frame)
    section_frame.pack(fill="x", pady=10, padx=10)

    # Add a label for the section name
    section_name = section.get("section_name", "Unnamed Section")
    section_label = ctk.CTkLabel(section_frame, text=section_name, font=("Arial", 18, "bold"))
    section_label.pack(anchor="w", pady=5)

    # Create questions within this section
    for question in section.get("questions", []):
        key = question["key"]  # Access 'key' field
        question_text = question["question"]  # Access 'question' field
        input_type = question["input_type"]  # Access 'input_type' field
        options = question.get("options", [])  # Use .get() in case there are no options
        default_value = question.get("default", None)
        create_question(key,question_text, input_type, options, section_frame, default_value,entries)


# -

# In the control tab, an advanced popup is used to decrease the burden on the user. Some inputs have defaults that the user can choose to modify or not. These two functions handle the creation of the advanced popup button (for control tab only) and the opening of the advanced popup. On submit, the window is closed.

# +
def open_advanced_popup(sections, entries):
    
    popup_window = ctk.CTkToplevel()
    popup_window.title("Advanced Options")
    popup_window.geometry("500x400")

    popup_scrollable_frame = ctk.CTkScrollableFrame(popup_window, width=480, height=360)
    popup_scrollable_frame.pack(fill="both", expand=True)

    for section in sections:
        create_section(section, popup_scrollable_frame, entries)
        
    def handle_submit(entries):
        # Hide the popup window after submission
        popup_window.withdraw()
    
    # Submit button to confirm inputs
    submit_button = ctk.CTkButton(
        popup_scrollable_frame, 
        text="Submit", 
        command=lambda: handle_submit(entries) 
    )
    submit_button.pack(pady=10)

def add_advanced_button_to_tab(tab_name, advanced_sections):
    if not advanced_sections:  # Skip adding the button if no advanced sections
        return
    tab_frame = tabview.tab(tab_name)
    advanced_button = ctk.CTkButton(
        tab_frame,
        text="Show Advanced Options",
        command=lambda: open_advanced_popup(advanced_sections, entries_dict[tab_name])
    )
    advanced_button.pack(pady=10)


# -

# Since the LPT tab requires dynamic creation, these two functions handle both the dynamic headers (Fraction 1,2..etc. for each Fraction section) and the dynamic sections.

# +
#function to generate more headers based on lagrangian_num
def lpt_fraction_headers(entries):
    fraction_headers = []
    counter = 1
    for fraction in range(1,int(entries["lagrangian_num"].get())+1):
        fraction_string = "#" * 30 + "Fraction" + str(counter) + "#" * 38
        counter+=1
        fraction_headers.append(fraction_string)
        
    return fraction_headers

def update_lpt_sections(lpt_entries, parent_frame):

    try:
        # Get the value of lagrangian_num
        lagrangian_num = int(lpt_entries["lagrangian_num"].get())
        print("How many fractions the user entered:", lagrangian_num)
        if lagrangian_num < 1:
            raise ValueError("lagrangian_num must be greater than 0.")
    except ValueError as e:
        error_messages.append(f"Invalid lagrangian_num: {e}")
        error_label.configure(text="\n".join(error_messages))
        return

    # Create new dynamic sections for fractions
    lpt_dict = create_fractions(fraction_lpt_sections, lagrangian_num)
    for section in lpt_dict:
        create_section(section, parent_frame, lpt_entries)
    


# -

# If the user submits and decides they want to keep editing and submit again, these functions allow the page to be refreshed (and the dictionaries to be refreshed) while saving the previous inputs. This allows for small edits on files after submitting. 
#

# +
#Function to save current user input values
def save_user_inputs(entries_dict):
    saved_values = {}
    for tab_name, entries in entries_dict.items():
        saved_values[tab_name] = {}
        for key, field in entries.items():
            if hasattr(field, "get"):  #If it's a widget with a 'get' method
                saved_values[tab_name][key] = field.get()
            else:
                saved_values[tab_name][key] = field  #For plain string values
    return saved_values

#Function to restore user inputs into fields
def restore_user_inputs(saved_values, entries_dict):
    for tab_name, entries in entries_dict.items():
        if tab_name in saved_values:
            for key, value in saved_values[tab_name].items():
                if key in entries:
                    field = entries[key]
                    if hasattr(field, "set"):  # If it's a widget with a 'set' method
                        field.set(value)
                    elif hasattr(field, "insert"):  # If it's an entry widget
                        field.delete(0, "end")  # Clear existing value
                        field.insert(0, value)

def refresh_page():
    global saved_user_inputs,error_messages
    saved_user_inputs = save_user_inputs(entries_dict)  # Save current user inputs
    saved_error_messages = list(error_messages)
    root.destroy()  # Destroy the current window
    main()  # Reinitialise the main function
    restore_user_inputs(saved_user_inputs, entries_dict)  # Restore saved inputs
    #After the page is refreshed, restore the error messages
    error_messages = saved_error_messages


# -

# A submit button is added to each tab. The logic for submit is typically calling 1-2 functions, handled for each tab in dictionaries within main

def add_submit_button(parent_frame, entries, dfs_function=None, headers=None, file_creation_function=None):

    # Define a wrapper function to call submit_inputs with the correct arguments
    def submit_for_tab():
        # Call the headers if it's a callable (for dynamic generation)
        actual_headers = headers() if callable(headers) else headers
        submit_inputs(
            entries=entries,
            file_creation_function=file_creation_function,
            dfs_function=dfs_function,
            headers=actual_headers,
        )

    # Create the Submit button and attach the wrapper function
    submit_button = ctk.CTkButton(parent_frame, text="Submit", command=submit_for_tab)
    submit_button.pack(pady=10)


# Submit inputs requires complex logic to handle errors, create file names, call file creation functions, and refresh the page. First, user inputs are collected and passed through the mapping dictionary. If the user has missed or incorrectly entered a value, it is displayed as an error. A temporary filename is created and the file creation function is passed. More errors may be displayed (particularly from integer convergence on Infodom tab) if errors arise during file creation. If there are no errors, a file dialog pops up asking the user to save the file. Upon save, the page is refreshed.

def submit_inputs(entries,file_creation_function,dfs_function=None, headers=None):
            
    error_messages.clear()
    
     #Step 1: Validate inputs and populate the dictionary
    for key, var in entries.items():
        #Check if the variable is a widget (has a .get method), otherwise it's a plain string
        if hasattr(var, "get"):
            value = var.get().strip()  #Retrieve value from the widget
        else:
            value = var  #If it's not a widget, it's already a string
        print(f"Key: {key} var {value}")
        
        if isinstance(var, CTkEntry):
            try:
                float(value)
            except ValueError:
                error_messages.append(f"Input Error: The field '{key}' must be a valid number")
             
        #Step 2: Pass inputs through mapping dictionary
        if not value:
            error_messages.append(f"The field '{key}' is required.")
        else:
            if value in mapping_dict:
                mapped_value = mapping_dict[value]
                entries[key] = mapped_value  # Replace with mapped value
            else:
                entries[key] = value  # Retain original value if no mapping

    #Step 3: Create a temporary placeholder filename to create the file and check for errors before user is asked to save
    with tempfile.NamedTemporaryFile(delete=False, suffix=".cin") as temp_file:
        placeholder_path = temp_file.name  # Temporary filename

    #Step 4: Try creating the file using the placeholder filename
    try:
        if file_creation_function:
            if dfs_function and headers:  # If dfs_function and headers are provided
                dfs = dfs_function()  # Generate the dataframe
                file_creation_function(dfs, headers, placeholder_path)  # Attempt to create the file
            else:  # No dfs_function or headers provided
                file_creation_function(placeholder_path)  # Call file creation with the placeholder
        else:
            print("No valid function provided for file creation.")

    except Exception as e:
        # Step 5: Handle any errors during file creation
        os.remove(placeholder_path)  #Clean up the temporary file
        error_label_text = f"Error during file creation: {e}"
        error_label.configure(text=error_label_text)
        print("Error during file creation:", e)

        #Step 5a: If there are validation errors, display them and exit
    if error_messages:
        error_label_text = "\n".join(error_messages)
        error_label.configure(text=error_label_text)  # Display errors
        error_label.pack(anchor="w")
        print("Errors:", error_messages)
        refresh_page()
        return

    #Event binding to clear errors dynamically
    def clear_error(event, key):
        if key in error_messages:
            error_messages.remove(key)
            error_label.configure(text="\n".join(error_messages))

    # Step 5b: If no errors, prompt the user to save the file
    file_path = filedialog.asksaveasfilename(
        defaultextension=".cin",
        filetypes=[("Input File", "*.cin"), ("All Files", "*.*")],
        title="Save File",
        parent=root
    )

    if not file_path:
        #If the user cancels the file dialog, show a warning and exit the function
        os.remove(placeholder_path)  #Clean up the temporary file
        messagebox.showwarning("File Save Cancelled", "File save was cancelled. Your inputs were not processed.")
        refresh_page()
        return

    #Save the temporary file to the user-selected path
    os.rename(placeholder_path, file_path)

    #Refresh the page, saving user input
    refresh_page()


# These three functions are used to initialize parts of the code. Initialize entries adds keys and default values to the user_entries dictionary, which are later updated from user inputs. Load_sections_to_tab creates each tab with the proper sections and questions, using functions created above. Initialize_tabs calls load_sections_to_tabs, and adds advanced and submit buttons

# +
#Add keys and defaults to empty user_input dictionaries (different for each tab)
def initialize_entries(sections,entries,advanced_sections=None):

    all_sections = sections + (advanced_sections or [])
    
    for section in all_sections:
        for question in section.get("questions", []):
            key = question["key"]
            default_value = question.get("default", "")
            entries[key] = ctk.StringVar(value=str(default_value))

#Function to load sections into a tab
def load_sections_to_tab(sections, tab_name):
    tab_frame = tabview.tab(tab_name)  #Get the frame for the tab
    #Create a scrollable frame for the current tab (Control, LPT, or Infodom)
    scrollable_frame = ctk.CTkScrollableFrame(tab_frame, width=480, height=360)
    scrollable_frame.pack(pady=10, anchor="w",fill="both",expand=True)
    entries = entries_dict[tab_name]
    for section in sections:
        create_section(section, scrollable_frame,entries)

#Load sections and add buttons on each tab
def initialize_tabs(sections_dict, advanced_sections_dict, entries_dict, file_creation_functions,dfs_functions, headers_dict):
    for tab_name, sections in sections_dict.items():
        load_sections_to_tab(sections, tab_name)
        
        tab_frame = tabview.tab(tab_name)  #Get the frame for the current tab

        #Add the Submit button with the corresponding entries and parameters
        add_submit_button(
            parent_frame=tab_frame,
            entries=entries_dict[tab_name],         #Pass the correct entries dictionary
            dfs_function=dfs_functions.get(tab_name),  #Pass the DataFrame creation function (if applicable)
            headers=headers_dict.get(tab_name),        #Pass the headers (if applicable)
            file_creation_function=file_creation_functions.get(tab_name),  #Pass the specialized file creation function
        )

        #Add the Advanced Options button for this tab if it exists (for control only)
        add_advanced_button_to_tab(tab_name, advanced_sections_dict.get(tab_name, []))



# -

# The main function keeps track of proper inputs for each tab and the app layout. It then calls the intialize_entries function for each tab, and calls initialize_tabs to begin the program. On refresh, it also handles the displaying of error messages if they exist and the repopulation of old user inputs if they have entered them previously.

# +
def show_splash_screen():
    global splash_shown
    if splash_shown:
        return  # Skip if the splash screen has already been displayed

    splash_shown = True
    splash = tk.Tk()
    splash.overrideredirect(True)

    # Centre the splash screen
    screen_width = splash.winfo_screenwidth()
    screen_height = splash.winfo_screenheight()
    x = (screen_width // 2) - (300 // 2)
    y = (screen_height // 2) - (200 // 2)
    splash.geometry(f"300x200+{x}+{y}")
    splash.configure(bg="black")  # Set background to black

    # Loading text (no background, matches the dark theme)
    label = tk.Label(
        splash, text="Loading...", font=("Arial", 14), fg="white", bg="black"
    )
    label.pack(pady=20)

    # Rounded progress bar using CustomTkinter
    progress = ctk.CTkProgressBar(splash, width=250, height=20, corner_radius=10)
    progress.pack(pady=20)
    progress.set(0)  # Set initial value to 0

    # Simulate loading
    def update_progress():
        for i in range(101):
            progress.set(i / 100)  # Update progress (value between 0 and 1)
            splash.update_idletasks()
            splash.after(30)
        splash.destroy()  # Destroy the splash screen after loading

    update_progress()


splash = show_splash_screen()

def main():
    
    global saved_user_inputs,error_messages,splash
    global root, tabview, sections_dict,advanced_sections_dict,entries_dict,dfs_functions,headers_dict,file_creation_functions,error_label
    global control_entries,lpt_entries,infodom_entries,mdmap_entries# Use global if needed elsewhere
    
    control_entries = {}  #For Control tab entries
    lpt_entries = {}      #For LPT tab entries
    infodom_entries = {}  #For Infodom tab entries
    mdmap_entries = {}

    #Data configuration for each tab (Control, LPT, Infodom)
    sections_dict = {
        "Control": control_sections,
        "LPT (optional)": lpt_sections,
        "Infodom": infodom_sections,
        "Mdmap": mdmap_sections
    }
    
    #Advanced options configuration for each tab
    advanced_sections_dict = {
        "Control": advanced_control_sections
    }

    #Entries dictionaries for each tab
    entries_dict = {
        "Control": control_entries,
        "LPT (optional)": lpt_entries,
        "Infodom": infodom_entries,
        "Mdmap": mdmap_entries
    }
    
    #DataFrame creation functions for each tab
    dfs_functions = {
        "Control": create_control_dfs,
        "LPT (optional)": create_lpt_dfs,
        "Infodom": None,
        "Mdmap": None
    }
    
    #Headers for each tab
    headers_dict = {
        "Control": control_headers,
        "LPT (optional)": lambda: lpt_headers + lpt_fraction_headers(lpt_entries),  #Dynamic headers
        "Infodom": None,
        "Mdmap": None
    }
    
    #File creation functions for each tab
    file_creation_functions = {
        "Control": create_file, 
        "LPT (optional)": create_file,
        "Infodom": create_infodom,  #Specialized function
        "Mdmap": create_mdmap_file_with_columns   #Specialized function
    }
    
    #Initialize the CustomTkinter app Main Page
    ctk.set_appearance_mode("dark") 
    ctk.set_default_color_theme("green") 
    
    root = ctk.CTk()
    root.title("Automated File Creation")
    
    #Create a CTkTabview widget for tabs
    tabview = ctk.CTkTabview(root)
    tabview.pack(pady=10, padx=10, fill="both", expand=True)
    
    #Add tabs: Control, LPT, and Infodom
    tabview.add("Control")
    tabview.add("Infodom")
    tabview.add("Mdmap")
    tabview.add("LPT (optional)")

    #Initialise error_label
    error_label = ctk.CTkLabel(
        root,
        text_color="#FF6666",
        fg_color="transparent",
        anchor="w",# Red background for errors
        wraplength=400  # Adjust as needed
    )
    
    initialize_entries(control_sections, control_entries,advanced_control_sections)
    initialize_entries(lpt_sections, lpt_entries)
    initialize_entries(infodom_sections,infodom_entries)
    initialize_entries(mdmap_sections,mdmap_entries)

    
    #Initialize the tabs with their sections, advanced options, entries, and file creation functions
    initialize_tabs(
        sections_dict=sections_dict,
        advanced_sections_dict=advanced_sections_dict,
        entries_dict=entries_dict,
        file_creation_functions=file_creation_functions,
        dfs_functions=dfs_functions,
        headers_dict=headers_dict
    )

    #Check for saved inputs. 
    #On refresh, when the main function is called again, this will evaluate to True if the user has entered values
    if saved_user_inputs:
        print("Inputs were saved.")
        print(saved_user_inputs)
        restore_user_inputs(saved_user_inputs, entries_dict) 

    if error_messages:  #Show error message if there are any
        error_label_text = "\n".join(error_messages)
        error_label.configure(text=error_label_text)  #Display errors
        error_label.pack()
    

    #Start the main loop
    root.mainloop()

# Show the splash screen and schedule the main application if not shown already
if not splash_shown:
    show_splash_screen()

main()
# -




