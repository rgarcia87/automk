#!/usr/bin/python3 
# Created by Rodrigo García-Muelas 
# Rev.2019Oct25
# 
"""This package generates the Input of a microkinetic model to be used in Maple.  
   
Future modules: 
    * Reactor behaviour as a function of time. 
    * Inherit initial conditions (preconverge).  
    * Non-isothermal reactors, T dependent on time. 
    * Non-isobaric   reactors, P dependent on time.  
    * Deactivation of key intermediates: Sensibility analysis.  
    * Cycle the model making the energies depend on two or more parameters (PCA). 
    * Consider coverage effects.    
"""

# Load libraries 
import amklib       

# Initialize variables
ltp={}                    # List to print dictionary

# Read configuration file 
conf=amklib.readconf("./parameters.txt") 

# Read the input files: gas, int, and rxn. 
# Format: dictionary of dictionaries. 
#cat=amklib.amklib.read('./cat.csv') ##### CHANGE THIS ASAP!!!
cat='iO' ################################# Correct when splitting into functions
gas=amklib.read('./gas.csv')
int=amklib.read('./int.csv')
rxn=amklib.read('./rxn.csv')
#print(cat, '\n \n', gas,'\n \n', int, '\n \n' , rxn)

# Prepare site balance equation, solver for SODE, and initial conditions. 
# Also initialize the list of differential equations.
int,sbalance,sodesolv,initialc,rhsparse=amklib.fint(conf,int,cat,ltp)

# Prepare kinetic constants and rates of adsorption/desorption.
# Also expand list of differential equations in "int" to include adsorption/desorptions. 
gas,int=amklib.fgas(conf,gas,int,cat,ltp) 

# Prepare kinetic constants and rates of all chemical steps. 
# Also expand list of differential equations in "int" to include chemical steps. 
rxn,int=amklib.frxn(conf,rxn,int,cat,ltp)

# Print Maple input. 
amklib.printtxt(conf,gas,int,rxn,cat,sbalance,initialc,sodesolv,rhsparse,ltp)

#amklib.rxntime(conf)

