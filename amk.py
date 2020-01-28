#!/usr/bin/python3 
# Created by Rodrigo García-Muelas 
# Rev.2020Jan28
# 
"""This package generates the Input of microkinetic models for Maple.  
   
Possible expansions:  
    * Non-isothermal reactors, T dependent on time. 
    * Non-isobaric   reactors, P dependent on time.  
    * Cycle the model making the energies depend on two or more parameters (PCA). 
    * Consider coverage effects.    
    * Put area of active site. 
    * Unidimensionah diffusions, taking stationary state conditions in Fick's law. 
Security checks to implement: 
    * Check if float(conf['Catalyst']['areaactivesite']) gives no error. 
    * Idem float(conf['Catalyst']['secondlayerthickness']) if gas-phase species. 
"""

# Load libraries 
import amklib       

# Initialize variables
ltp={}                    # List-to-print dictionary of lists

# Read configuration file 
conf=amklib.readconf("./parameters.txt") 

# Read the input files int&rxn as dictionary of dictionaries. 
itm=amklib.read('./itm.csv')
rxn=amklib.read('./rxn.csv')
#print('\n \n', int, '\n \n' , rxn, '\n \n')

# Prepare site balance equation, solver for SODE, and initial conditions. 
# Also initialize the list of differential equations.
itm,sbalance,sodesolv,initialc,rhsparse=amklib.process_intermediates(conf,itm,ltp)

# Prepare kinetic constants and rates of all chemical steps. 
# Also expand list of differential equations in "itm" to include chemical steps. 
itm,rxn=amklib.process_rxn(conf,itm,rxn,ltp)

# Print Maple input. 
amklib.printtxt(conf,itm,rxn,sbalance,initialc,sodesolv,rhsparse,ltp)


