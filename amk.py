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

import amklib       
import os

#Constants 
kbh="20836612225.1252"    # Boltzmann constant divided by Planck constant, s^-1, string.  
kbev="8.617333262145E−5"  # Boltzmann constant in eV·K−1, string. 

#try: 
#    from var import time1  
#except: 
#    time1=[10800.0]

time1="10800" 
time1=[1.0,2.0,4.0,8.0] 
print(type(time1))   



# Read the input files: gas, int, and rxn. 
# Format: dictionary of dictionaries. 
#cat=amklib.amklib.read('./cat.csv') ##### CHANGE THIS ASAP!!!
cat='iO' ################################# Correct when splitting into functions
gas=amklib.read('./gas.csv')
int=amklib.read('./int.csv')
rxn=amklib.read('./rxn.csv')
#print(cat, '\n \n', gas,'\n \n', int, '\n \n' , rxn)
#print('\n \n', int['iP']['diff'], '\n \n')    

# Prepare site balance equation, solver for SODE, and initial conditions. 
# Also initialize the list of differential equations. 
int,sbalance,sodesolv,initialc,rhsparse = amklib.fint(int,cat)

# Prepare kinetic constants and rates of adsorption/desorption.
# Also expand list of differential equations in "int" to include adsorption/desorptions. 
gas,int=amklib.fgas(gas,int,cat) 

# Prepare kinetic constants and rates of all chemical steps. 
# Also expand list of differential equations in "int" to include chemical steps. 
rxn,int=amklib.frxn(rxn,int,cat)

#print("\n",int,"\n\n",sbalance,"\n\n",sodesolv,"\n\n",initialc,"\n")

# Print Maple input. 
printtxt(gas,int,rxn,cat,time1,sbalance,initialc,sodesolv,rhsparse) 


