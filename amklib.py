# -*- coding: utf-8 -*-
import pandas as pd
import os
import configparser, ast
import copy  
     
#Constants 
kbh="20836612225.1252"    # Boltzmann constant divided by Planck constant, s^-1, string.  
kbev="8.617333262145E-5"  # Boltzmann constant in eV·K−1, string. 
avogadro=6.02214199E23    # Avogadro's constant. 
     
def readconf(filename='./parameters.txt'):  
    """This function reads the input parameters from a file
     
    Args: 
        filename: Input file in Windows .ini format. Comments should be provided as "#" 
        
    Returns: 
        conf: Configuration data. 
        
    """
    
    conf=configparser.ConfigParser(inline_comment_prefixes=('#'))
    conf.read(filename)
    return conf   
    
def rxntime(conf) : 
    """Subroutine that interpretes the time.   
        Requires the ast package 
    
    Args: 
        conf: Configuration data. 
    """
    
    time1raw=conf['Reactor']['time1']
    if time1raw.find(',') >0 : # If it is a list 
        timel=True 
        time1=ast.literal_eval(time1raw)
    else :                     # If it is an unique value  
        timel=False 
        time1=time1raw
    return time1, timel 
    
def read(filename='./itm.csv') :  
    """This function reads a file containing information of species in 
    gas, aqu(eous), or adsorbed on cat(alyst). 
    It can also read the reactions file. . 
    It requires pandas to be loaded.  
      
    Args: 
        filename: Input file. The columns are separated by one or more spaces.
            Energies must be provided in eV and frequencies in cm-1 Format:
      
    Returns: 
        dicint: a dictionary containing at least the tags, energies, and frequencies of all species.     
      
    """
     
    dic=pd.read_csv(filename, delim_whitespace=True, index_col='label').T.to_dict()
    return dic  
     
     
def get_damptime(conf) :  
    """Parse pressure damp from configuration file
     
    Args: 
        conf: Configuration data. 
     
    Returns:      
        dampt1: Pressure damp in processing. 
        dampt2: Pressure damp in post-processing. 
    """
         
    try :          
        damptime=float(conf['Reactor']['damptime'])   
    except :   
        damptime=1.0   
    
    if  damptime>1E-13 : 
        dampt1="*(1-exp(-"+"{:.6E}".format(damptime)+"*t))^2"
        dampt2="*(1-exp(-"+"{:.6E}".format(damptime)+"*timei))^2"
    else : 
        dampt1=""
        dampt2=""
    return dampt1, dampt2 
     
      
def get_elecpot(conf) : 
    """This function extracts the electric potential vs RHE from the configuration file. 
    returns the electric potential vs SHE. 
    """
    try :  
        elecpot=(float(conf['Electrochemistry']['electricpotentialrhe'])-
                 float(conf['Electrochemistry']['pH'])*float(kbeV)*  
                 float(conf['Reactor']['reactortemp'])*ln(10.0) ) 
    except :  
        elecpot=0.0 
    return elecpot 
        
    
def get_nelect_for_itm(itm,item,label) : 
    if item==None : 
        nelect=0.0  
    else : 
        try : 
            nelect=float(itm[item][label]) 
        except : 
            nelect=0.0 
            print("Nasty problem found in ",item) 
    return nelect
       
       
def get_nelect_for_rxn(conf,itm,rxn) :  
    """Get the number of electrons for a particular transition state
    from alpha values
    """ 
    try :   
        label=conf['Electrochemistry']['nelectronslabel']  
    except :   
        label="ne"  
    for item in sorted(rxn) : 
        rxn[item][label]=((1-float(rxn[item][alpha]))*(
                          get_nelect_for_itm(itm,rxn[item]['is1'],label)+
                          get_nelect_for_itm(itm,rxn[item]['is2'],label))+
                          float(rxn[item][alpha])*(
                          get_nelect_for_itm(itm,rxn[item]['fs1'],label)+
                          get_nelect_for_itm(itm,rxn[item]['fs2'],label)))  
          
      
def adjust_energy_with_potential(conf,itm,elecpot) : 
    """Adds the electric potential component to the Gibbs energy. 
    """ 
    try :  
        label=conf['Electrochemistry']['nelectronslabel']  
    except :  
        label="ne"  
    for item in sorted(itm) :  
        try :  
            itm[item]['G']=float(itm[item]['G'])+float(itm[item][label])*elecpot
        except :   
            print("Error found adjusting the potential of ", item) 
            print(item)
            print(label) 
            print(itm[item]['G'])
            print(float(itm[item]['G'])) 
            print(float(itm[item][label]))
            exit()  
        
      
def process_intermediates(conf,itm,ltp) :
    """This function process the "intermediates" dataframe to generate 
    the site-balance equation, the SODE-solver, and the initial conditions as clean surface. 
    It also initializes the list of differential equations.  
    
    Args: 
        conf: Configuration data.
        itm: Dict of dicts containing at least a list of intermediates as index. 
    
    Returns:
        itm:      Expanded dict of dicts containing also the list of differential equations. (Mutable)
        sbalance: Site-balance equation. (Unmutable)
        sodesolv: Input of the SODE-solver. (Unmutable)   
        initialc: Initial conditions as clean surface. (Unmutable) 
    """
    
    # Initialize variables related to intermediates. 
    sbalance="c"+conf['Catalyst']['sitebalancespecies']+":=(t)-> 1.0"  
    sodesolv="Solution:=dsolve({"   
    initialc="IC0:="
    rhsparse="" 
    index=1  
    # Initialize list-to-print for postprocessing
    ltp['prs']=[] # ltp of pressures and concentrations-in-second-layer.     
    #ltp['itm']=["sc"+conf['Catalyst']['sitebalancespecies']] # ltp of interm.: init w/ s-b species     
    ltp['itm']=[conf['Catalyst']['sitebalancespecies']] # ltp of interm.: init w/ s-b species
      
    # Process intermediates, starting by adsorbed (cat), then gas. 
    for item in sorted(itm) :  
    # SERGIO: 
    # for key,value in sorted(itm).items() : #key~item ; value~itM[item] (all line)             
    # so the input of the sub-function will be the key and value
        if  itm[item]['phase']=='cat' and item!=conf['Catalyst']['sitebalancespecies'] :      
            # A surface species 
               
            # Initialize diff equations to count in which reactions each species participate     
            itm[item]['diff']="eqd"+item+":=diff(c"+item+"(t),t)="         
            #value['diff']="eqd"+key+":=diff(c"+key+"(t),t)="
             
            # Prepare site balance  
            sbalance+=" -c"+item+"(t)"
             
            # Prepare list of differential equations for the SODE solver  
            sodesolv+="eqd"+item+", "
             
            # Prepare list of default initial conditions as clean surface
            initialc+=" c"+item+"(0.0)=0.0,"
              
            # Prepare parser of concentrations after SODE is solved  
            index+=1 # First element should be 1+1=2. Do not touch.  
            rhsparse+="sc"+item+":=rhs(S["+str(index)+"]) : "
             
            # List of reactions for fprintf function in Maple 
            ltp['itm'].append(item)
              
        elif itm[item]['phase']=='gas' : 
            # Get partial pressures 
            try : 
                itm[item]['pressure']=conf['Pressures'][item] 
            except : 
                itm[item]['pressure']=0 
            # Generate list of pressures
            ltp['prs'].append("P"+item)
             
        elif itm[item]['phase']=='aqu' : 
            # Get concentrations and convert to molecules/activesite. 
            try :  
                itm[item]['concentration']=(float(conf['Concentrations'][item])*
                                            float(conf['Catalyst']['areaactivesite'])*
                                            float(conf['Catalyst']['secondlayerthickness'])*
                                            avogadro*1E-27) 
            except :   
                itm[item]['concentration']=0.0  
            # Generate list-to-print of concentrations-in-the-second-layer; put along pressures. 
            ltp['prs'].append("CSL"+item)   
                  
        elif item!=conf['Catalyst']['sitebalancespecies'] : 
            print("Unknown phase for ",item,itm[item]['phase'],
                  "\n I only recognize 'aqu', 'cat', and 'gas'") 
            exit()
                                            
    # Close the site-balance equation     
    sbalance=sbalance+" : " 
     
    # Close the sodesolv
    sodesolv=sodesolv+"IC0}, numeric, method=rosenbrock, maxfun=0, abserr=1E-16, interr=false);"
     
    # In the initial conditions, replace the last comma by a colon 
    initialc=initialc[:-1]+" : "
      
    return itm, sbalance, sodesolv, initialc, rhsparse  
     
     
def is_gas(itm,rxn,item,state) :  
    """ Returns 1 if a given (initial/final) state of rxn #item is gas.
    Returns 0 otherwise.  """
    #print(item,state,rxn[item][state],itm['gP']['phase']) 
    if   rxn[item][state]=='None' or rxn[item][state]==None : 
        gas=0 
    else : 
        if itm[rxn[item][state]]['phase']=='gas' : 
            gas=1  
        elif itm[rxn[item][state]]['phase']=='cat' or itm[rxn[item][state]]['phase']=='aqu' : 
            gas=0  
        else : 
            print("Phase of rxn#",item," intermediate ",rxn[item][state],":", 
                  itm[rxn[item][state]]['phase'],"Not recognized" ) 
    return gas   
        
        
def mw_gas(itm,rxn,item,state) :  
    """ Returns the mass weight of a gas-phase intermediate. Zero if adsorbed """
    if   rxn[item][state]=='None' or rxn[item][state]==None :
        mw=0
    else :
        if itm[rxn[item][state]]['phase']=='gas' :
            mw=float(itm[rxn[item][state]]['mw'])
        elif itm[rxn[item][state]]['phase']=='cat' :
            mw=0              
        elif itm[rxn[item][state]]['phase']=='aqu' :
            mw=0 
        else :
            print("Phase of rxn#",item," intermediate ",rxn[item][state],":",
                  itm[rxn[item][state]]['phase'],"Not recognized" )
    return mw 
     
     
def kinetic_constants(conf,itm,rxn,item) : 
    """ Prepares the kinetic constants for direct and (i)reverse semireactions 
    depending on the number of gas-phase intermediates. 
    Returns error if there are more than two species in gas for a given semirxn. """
    rxn[item]['kd']="" 
    rxn[item]['ki']="" 
    howmanygasd=is_gas(itm,rxn,item,'is1')+is_gas(itm,rxn,item,'is2')
    howmanygasi=is_gas(itm,rxn,item,'fs1')+is_gas(itm,rxn,item,'fs2')
    area="{:.6f}".format( float(conf['Catalyst']['areaactivesite']) ) # Site area in Å²
    # Direct semireaction:      
    if   howmanygasd==0 :
        # If semireaction on surface: use Arrhenius kb*T/h*exp(-Ga/kB*T)
        rxn[item]['kd']="k"+item+"d:=evalf("+kbh+"*T*exp(-max(0.0,"+\
                        "{:.6f}".format( rxn[item]['aGd'])+","+\
                        "{:.6f}".format( rxn[item]['dGd'])+\
                        ")/("+kbev+"*T)) ) : "
    elif howmanygasd==1 :
        mw="{:.6f}".format(mw_gas(itm,rxn,item,'is1')+mw_gas(itm,rxn,item,'is2'))
                                           # (atm=>Pa)*Area*(Å²=>m²)
        rxn[item]['kd']="k"+item+"d:=evalf((101325*"+area+"*1E-20"+\
                        "*exp(-max(0.0,"+\
                        "{:.6f}".format( rxn[item]['aGd'])+","+\
                        "{:.6f}".format( rxn[item]['dGd'])+\
                        ")/("+kbev+"*T)))"+\
                        "/sqrt(2*Pi*1.6605390400E-27*"+mw+"*1.3806485200E-23*T )) : "
                        # Denominator: sqrt(2Pi(elemmass@kg)*massweight*kB(SI)*T
    else :
        print("WARNING! direct reaction #",item,"has",howmanygasd,"gas/aq reactants.")
        print("Abnormal termination")
        exit()
    # Reverse (i) semireaction:      
    if   howmanygasi==0 : 
        rxn[item]['ki']="k"+item+"i:=evalf("+kbh+"*T*exp(-max(0.0,"+\
                        "{:.6f}".format( rxn[item]['aGi'])+","+\
                        "{:.6f}".format(-rxn[item]['dGd'])+\
                        ")/("+kbev+"*T)) ) : "
    elif howmanygasi==1 : 
        mw="{:.6f}".format(mw_gas(itm,rxn,item,'fs1')+mw_gas(itm,rxn,item,'fs2')) 
        rxn[item]['ki']="k"+item+"i:=evalf((101325*"+area+"*1E-20"+\
                        "*exp(-max(0.0,"+\
                        "{:.6f}".format( rxn[item]['aGi'])+","+\
                        "{:.6f}".format(-rxn[item]['dGd'])+\
                        ")/("+kbev+"*T)))"+\
                        "/sqrt(2*Pi*1.6605390400E-27*"+mw+"*1.3806485200E-23*T )) : "
    else :
        print("WARNING! reverse reaction #",item,"has",howmanygasd,"gas/aq reactants.")
        print("Abnormal termination")
        exit()
        
      
def process_itm_on_rxn(conf,itm,rxn,item,state,dampt1,dampt2) :  
    """ Use the is/fs states for each reaction to get their activation energies. 
    Then write the formula for reaction rate according to their intermediates. 
        This formula is split between rtd (direct part) and rti (inverse part). 
    Then update the differential equations in which each adsorbed species participates. 
    If a rectant is in gas phase, include a Hertz-Knudsen term in the constant 
        and its pressure as variable in the reaction rate.
    """
     
    if   state=='is1' or state=='is2' : 
        semirxn='d'
        sign='-' # Consume reactants 
    elif state=='fs1' or state=='fs2' : 
        semirxn='i'
        sign='+' # Increase products 
    else : 
        print("Wrong state for reaction", item, "\nOnly 'is1', 'is2', 'fs1', and 'fs2' supported") 
        exit() 
         
    # Get energy of the (initial/final) state "i" 
    if rxn[item][state]=='None' or rxn[item][state]==None :
        G=0.0
    else:
        try:
            G=itm[rxn[item][state]]['G'] 
        except: 
            print("\n Error!, reaction ",item, " comes from ",state, rxn[item][state],
                  " whose energy was not found.")
            exit()
        # If (initial/final) state "i" is on catalyst, include concentration in rxn equation
        # and add rxn to respective differential equation. 
        if  itm[rxn[item][state]]['phase']=='cat':
            rxn[item]['rt'+semirxn]+="*c"+rxn[item][state]+"(t)"
            rxn[item]['srt'+semirxn]+="*sc"+rxn[item][state]
            # Do not generate differential equation for site-balance species (empty site?). 
            if rxn[item][state]!=conf['Catalyst']['sitebalancespecies'] :
                itm[rxn[item][state]]['diff']+=sign+"r"+item+"(t)"
        # If (initial/final) state "i" is "gas" (or aqueous) use P instead of c(t) 
        # and do not generate any differential equation.  
        elif itm[rxn[item][state]]['phase']=='gas':
            rxn[item]['rt'+semirxn]+= dampt1+"*P"+rxn[item][state]
            rxn[item]['srt'+semirxn]+=dampt2+"*P"+rxn[item][state] 
        # If (initial/final) state "i" is "aqu" (or aqueous) use CSL instead of c(t) 
        # and do not generate any differential equation.  
        elif itm[rxn[item][state]]['phase']=='aqu': 
            rxn[item]['rt'+semirxn]+= dampt1+"*CSL"+rxn[item][state]    
            rxn[item]['srt'+semirxn]+=dampt2+"*CSL"+rxn[item][state]   
    return G  
       
        
def process_rxn(conf,itm,rxn,ltp) : 
    """Subroutine that expands the "rxn" dictionary of dictionaries to include 
    the kinetic constants and rates of all chemical reactions. 
    It also expands the list of differential equations in "itm"
    and the list of reactions in which each intermediate participates. 
    
    Args: 
        conf: Configuration data.
        rxn: Dict of dicts containing the reactions. (Mutable)
        itm: Dict of dicts containing at least a list of intermediates as index. (Mutable)
    
    Returns: 
        rxn: Expanded dict of dicts containing adsorption/desorption constants and rates. (Mutable)
        itm: Expanded dict of dicts with list of differential equations updated with chemical reactions. (Mutable)
    """ 
          
    # Get pressure damp     
    dampt1,dampt2=get_damptime(conf)     
         
    # Initialize list-to-print: reactions, for postprocessing. 
    ltp['rxn']=[]  
        
    for item in sorted(rxn) : 
        # Initialize variables 
        rxn[item]['dGd']=""
        rxn[item]['aGd']=""
        rxn[item]['aGi']=""
        rxn[item]['rtd']="r"+item+":=(t)-> k"+item+"d"
        rxn[item]['rti']="-k"+item+"i"
        rxn[item]['srtd']="sr"+item+":= k"+item+"d" 
        rxn[item]['srti']="-k"+item+"i"
        #howmanygasd=0
        #howmanygasi=0        
           
        # Use the is/fs states for each reaction to get their activation energies. 
        # Then write the formula for reaction rate according to their intermediates. 
        #     This formula is split between rtd (direct part) and rti (inverse part). 
        # Then update the differential equations in which each adsorbed species participates. 
        # If a rectant is in gas phase, include a Hertz-Knudsen term in the constant 
        #     and its pressure as variable in the reaction rate.
        # For gaseous and aqueous species, include a damping term for numerical stability.       
        Gi1=process_itm_on_rxn(conf,itm,rxn,item,'is1',dampt1,dampt2)  
        Gi2=process_itm_on_rxn(conf,itm,rxn,item,'is2',dampt1,dampt2)  
        Gf1=process_itm_on_rxn(conf,itm,rxn,item,'fs1',dampt1,dampt2)  
        Gf2=process_itm_on_rxn(conf,itm,rxn,item,'fs2',dampt1,dampt2)  
         
        # Get reaction (dG) and (aG) activation energies for each reaction, both direct and inverse.  
        #print("\n",Gdi1,Gdi2,Gdf1,Gdf2)    
        rxn[item]['dGd']=Gf1+Gf2       -Gi1-Gi2 
        rxn[item]['aGd']=rxn[item]['G']-Gi1-Gi2 
        rxn[item]['aGi']=rxn[item]['G']-Gf1-Gf2 
             
        # Close the reaction formula with ":"
        rxn[item]['rti']=rxn[item]['rti']+" : "
          
        # Get kinetic constants 
        kinetic_constants(conf,itm,rxn,item)        
           
        # List of reactions for fprintf function in Maple 
        ltp['rxn'].append(item)
           
    return itm, rxn 
        
          
def printtxt(conf,itm,rxn,sbalance,initialc,sodesolv,rhsparse,ltp) :  
    # Before called printtxtsr
    """Subroutine that prints a given calculation for Maple, just a 's'ingle 'r'un 
      
    Args: 
        conf: Configuration data. 
            time1: Tiime or times for which the concentrations and reactions shall be printed (str/int/float, or list).
        itm: Dict of dicts listing the intermediates by an index. (Mutable)
        rxn: Dict of dicts for reactions. (Mutable)
        sbalance: Site-balance equetion, string. 
        initialc: Initial conditions, string. 
        sodesolv: Calls SODE solver in Maple, string.  
        rhsparse: Parser of surface concentrations, string. 
    """
    
    print("# Heading " )
    print("restart : \n " )  
        
    # Open file and print labels
    print('filename1:=FileTools[Text][Open]("',
          conf['General']['mapleoutput'].replace('"','').replace("'","").replace(" ",""),
          '",create,overwrite) : ') # Remove " ' and spaces from name of files. 
    print('fprintf(filename1,"%q %q\\n",','catalyst, "timei", "T",', 
          ', '.join(['"'+item+'"' for item in ltp['prs']]) ,",", 
          ', '.join(['"'+item+'"' for item in ltp['itm']]) ,",",
          ', '.join(['"'+item+'"' for item in ltp['rxn']]) ,   
          " ): " )   
    print('FileTools[Flush](filename1) : \n ')  
      
    # Temperature, pressures, and concentration.  
    print("T:=", conf.get("Reactor","reactortemp"), " : " )
    for item in sorted(itm) : 
        if itm[item]['phase']=='gas' :  
            print('P'+item+":=",itm[item]['pressure']," : ") 
    for item in sorted(itm) : 
        if itm[item]['phase']=='aqu' :   
            print('CSL'+item+":=",itm[item]['concentration']," : ") 
      
    print("\n# Kinetic constants")
    for item in sorted(rxn) :
        print(rxn[item]['kd']+"\n"+rxn[item]['ki']  )
      
    print("\n# Reaction rates:")
    for item in sorted(rxn) :
        print(rxn[item]['rtd'],rxn[item]['rti'])
      
    print("\n# Site-balance equation: ")
    print(sbalance)
     
    print("\n# Differential equations: ")
    for item in sorted(itm) :
        if  itm[item]['phase']=='cat' and item!=conf['Catalyst']['sitebalancespecies'] : 
            print(itm[item]['diff']," : ")
      
    print("\n# Initial conditions: ")
    print(initialc)
      
    print("\n# SODE Solver: ")
    print(sodesolv)
              
    # Time control: 
    time1,timel=rxntime(conf)
    if timel : 
        print("\n\nfor timei in " + str(time1) + " do ")
    else : 
        print("timei:= "+time1+" : ")
      
    print("S:=Solution(timei) : ")
    
    print("\n# Solution parser: ")
    print(rhsparse)
    
    print("\n# Site-balance equation after solver: ")
    tmp=os.popen(    "echo \""+sbalance+"\" | sed 's/(t)//g' "   ).read() # popen and read() are used to save in variable 
    tmp="s"+os.popen("echo \""+tmp[:-1]+"\" | sed 's/->//g' "    ).read() # rm last character (newline) 
    print(  os.popen("echo \""+tmp[:-1]+"\" | sed 's/-c/-sc/g' " ).read()[:-1] )
    
    print("\n# Reaction rates after solver: ")
    for item in sorted(rxn) :
        print(rxn[item]['srtd'],rxn[item]['srti']," : ")
                   
    # Print results 
    print("\nfprintf(filename1",',"%q %q\\n",',conf['Catalyst']['name'],', timei, T,',
          ', '.join([item for item in ltp['prs']]) ,",",
          ', '.join([item for item in ltp['itm']]) ,",", 
          ', '.join([item for item in ltp['rxn']]) , 
          " ): " )  
       
    print('\nFileTools[Flush](filename1) : ' )   
     
    if timel :     
        print("\nod: \n ") 
    
    # Print close file instruction  
    print('\nclose(filename1) : \n \n ')
     
    
