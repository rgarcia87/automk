import pandas as pd
import os
import configparser, ast
import copy  
     
#Constants 
kbh="20836612225.1252"    # Boltzmann constant divided by Planck constant, s^-1, string.  
kbev="8.617333262145E-5"  # Boltzmann constant in eV·K−1, string. 
     
def readconf(filename='./parameters.txt'):  
    """This function reads the input parameters from a file
     
    Args: 
        filename: Input file in Windows .ini format. Comments should be provided as "#" 
        
    Returns: 
        conf: Configuration data. 
        
    """
    
    conf=configparser.ConfigParser(inline_comment_prefixes=('#'))
    conf.read(filename)
    return(conf)  
    
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
    return(time1,timel)
  
def read(filename='./int.csv'): 
    """This function reads a file containing information of catalysts, gas, intermediates, or reactions. 
    It requires pandas to be loaded.  
    
    Args: 
        filename: Input file. The columns are separated by one or more spaces.
            Energies must be provided in eV and frequencies in cm-1 Format:
             label  phase  formula     G      e   mw    frq                     
              gR     gas   CH3CHO      0.000  1  52.0   [99,500,2000]    
              gP     gas   CH2OCH2     0.100  1  52.0   [100,200,1000]  
              gU     gas   CH2CHOH     0.100  1  52.0   [120,350,1500] 
              iO     cat   EmptySurf   0.000  1   0.0   []              
              iR     cat   CH3CHO     -1.000  1  52.0   [99,500,2000]            
              iI1    cat   CH2OCH2    -1.050  1  52.0   [100,200,1000] 
              iI2    cat   CH2CHOH    -0.950  1  52.0   [100,200,1000] 
              iP     cat   Unknown    -1.000  1  52.0   [100,200,1000]            
              iU     cat   Unknown    -2.000  1  52.0   [120,350,1500]            
    Returns: 
        dicint: a dictionary containing at least the tags, energies, and frequencies of all species. 
    
    """
    
    dic=pd.read_csv(filename, delim_whitespace=True, index_col='label').T.to_dict()
    return(dic) 
     
     
def process_intermediates(conf,itm,ltp):
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
    sbalance="c"+conf["Reactor"]["sitebalancespecies"]+":=(t)-> 1.0"  
    sodesolv="Solution:=dsolve({"   
    initialc="IC0:="
    rhsparse=""
    index=1
    
    # List to print from postprocessing: intermediates: Initialize with site-balance species. 
    ltp['itm']=["sc"+conf["Reactor"]["sitebalancespecies"]]     
     
    # Get pressure damp 
    try :      
        pressuredamptime=float(conf['Reactor']['pressuredamptime'])   
    except :   
        pressuredamptime=1.0   
    
    if pressuredamptime>1E-13 : 
        pdamp1="(1-exp(-"+"{:.6E}".format(pressuredamptime)+"*t))^2*"
        pdamp2="(1-exp(-"+"{:.6E}".format(pressuredamptime)+"*timei))^2*"
    else : 
        pdamp1=""
        pdamp2=""
        
    ltp['prs']=[]
      
    # Process intermediates, starting by adsorbed (cat), then gas. 
    for item in sorted(itm) :  
    # SERGIO: 
    # for key,value in sorted(itm).items() : #key~item ; value~itM[item] (all line)             
    # so the input of the sub-function will be the key and value
        if  itm[item]['phase']=='cat' and item!=conf["Reactor"]["sitebalancespecies"] :      
            # A surface species 
               
            # Initialize diff equations to count in which reactions each species participate.     
            itm[item]['diff']="eqd"+item+":=diff(c"+item+"(t),t)="         
            #value['diff']="eqd"+key+":=diff(c"+key+"(t),t)="
             
            # Prepare site balance  
            sbalance+=" -c"+item+"(t)"
             
            # Prepare list of differential equations for the SODE solver 
            sodesolv+="eqd"+item+", "
             
            # Prepare list of default initial conditions as clean surface
            initialc+=" c"+item+"(0.0)=0.0,"
              
            # Prepare parser of concentrations after SODE is solved  
            index+=1 
            #print("index: ",index, "rhsparse: ", rhsparse , "sodesolv: ", sodesolv)
            rhsparse+="sc"+item+":=rhs(S["+str(index)+"]) : "
            #print(rhsparse)  
             
            # List of reactions for fprintf function in Maple 
            #ltp['itm']+="sc"+item+", " 
            ltp['itm'].append("sc"+item)
              
            # Initialize list of reactions in which each intermediate participate. 
            # Deprecated: No longer needed.              
            #itm[item]['rxnlst']=[] 
                 
        elif itm[item]['phase']=='gas' : 
            # Get partial pressures 
            try : 
                itm[item]['pressure']=conf['Pressures'][item] 
            except : 
                itm[item]['pressure']=0 
            # Generate list of pressures
            #ltp['prs']+="P"+item+", "
            ltp['prs'].append("P"+item)
              
        elif item!=conf["Reactor"]["sitebalancespecies"] : 
            print("Unknown phase for ",item," \n I only recognize 'cat' and 'gas'") 
            exit()
                                            
    # Close the site-balance equation     
    sbalance=sbalance+":" 
     
    # Close the sodesolv
    sodesolv=sodesolv+"IC0}, numeric, method=rosenbrock, maxfun=0, abserr=1E-16, interr=false);"
     
    # In the initial conditions, replace the last comma by a colon 
    initialc=initialc[:-1]+" : "
      
    #print("\n",sbalance,"\n",initialc,"\n",sodesolv) 
    return(itm,sbalance,sodesolv,initialc,rhsparse) 
     
     
def is_gas(rxn,itm,item,state): 
    """ Returns 1 if a given (initial/final) state of rxn #item is gas.
    Returns 0 otherwise.  """
    #print(item,state,rxn[item][state],itm['gP']['phase']) 
    if   rxn[item][state]=='None' or rxn[item][state]==None :
        gas=0 
    else : 
        if itm[rxn[item][state]]['phase']=='gas' :
            gas=1 
        elif itm[rxn[item][state]]['phase']=='cat' :
            gas=0
        else : 
            print("Phase of rxn#",item," intermediate ",rxn[item][state],":",
                  itm[rxn[item][state]]['phase'],"Not recognized" ) 
    return(gas)  
     
     
def kinetic_constants(conf,itm,rxn,item) : 
    """ Prepares the kinetic constants for direct and (i)reverse semireactions 
    depending on the number of gas-phase intermediates. 
    Returns error if there are more than two species in gas for a given semirxn. 
    """ 
    # Kinetic constant for each reaction.         
    # If there are no species in gas-phase, multiply by kB*T/h. 
    # In any case, multiply by exp(-Ga/kB*T) (Arrhenius Eq, or sticking coefficient).  
    # Not yet implemented: reactions at interface and diffusions. 
    howmanygasd=is_gas(rxn,itm,item,'is1')+is_gas(rxn,itm,item,'is2')
    howmanygasi=is_gas(rxn,itm,item,'fs1')+is_gas(rxn,itm,item,'fs2')
    area=conf['Reactor']['areaactivesite']
    
    if   howmanygasd==0 :
        rxn[item]['kd']+=kbh+"*T*exp(-max(0.0,"+\
                        "{:.6f}".format( rxn[item]['aGd'])+","+\
                        "{:.6f}".format( rxn[item]['dGd'])+\
                        ")/("+kbev+"*T)) ): "
    elif howmanygasd==1 :
        rxn[item]['kd']+="exp(-max(0.0,"+\
                        "{:.6f}".format( rxn[item]['aGd'])+","+\
                        "{:.6f}".format( rxn[item]['dGd'])+\
                        ")/("+kbev+"*T)) ): "
    else :
        print("WARNING! direct reaction #",item,"has",howmanygasd,"gas/aq reactants.")
        print("Abnormal termination")
        exit()
     
    if   howmanygasi==0 :
        rxn[item]['ki']+=kbh+"*T*exp(-max(0.0,"+\
                        "{:.6f}".format( rxn[item]['aGi'])+","+\
                        "{:.6f}".format(-rxn[item]['dGd'])+\
                        ")/("+kbev+"*T)) ): "
    elif howmanygasi==1 :
        rxn[item]['ki']+="*exp(-max(0.0,"+\
                        "{:.6f}".format( rxn[item]['aGi'])+","+\
                        "{:.6f}".format(-rxn[item]['dGd'])+\
                        ")/("+kbev+"*T)) ): "
    else :
        print("WARNING! reverse reaction #",item,"has",howmanygasd,"gas/aq reactants.")
        print("Abnormal termination")
        exit()
     
     
def process_itm_on_rxn(conf,itm,rxn,item,state='is1'): 
    """ Use the is/fs states for each reaction to get their activation energies. 
    Then write the formula for reaction rate according to their intermediates. 
        This formula is split between rtd (direct part) and rti (inverse part). 
    Then update the differential equations in which each adsorbed species participates. 
    If a rectant is in gas phase, include a Hertz-Knudsen term in the constant 
        and its pressure as variable in the reaction rate.
    """ 
    if   state=='is1' or state=='is2' : 
        semirxn='d' 
    elif state=='fs1' or state=='fs2' : 
        semirxn='i'
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
            #rxn[item]['rt'+semirxn]=rxn[item]['rt'+semirxn]+"*c"+rxn[item][state]+"(t)"
            #rxn[item]['srt'+semirxn]=rxn[item]['srt'+semirxn]+"*sc"+rxn[item][state]
            rxn[item]['rt'+semirxn]+="*c"+rxn[item][state]+"(t)"
            rxn[item]['srt'+semirxn]+="*sc"+rxn[item][state]
            # Exclude the central species from site-balance equation from differential equations
            if rxn[item][state]!=conf["Reactor"]["sitebalancespecies"] :
                itm[rxn[item][state]]['diff']+="-r"+item+"(t)"
        # If (initial/final) state "i" is "gas" (or aqueous) use p instead of c(t) 
        # and do not generate any differential equation.  
        elif itm[rxn[item][state]]['phase']=='gas':
            rxn[item]['rt'+semirxn]+="*P"+rxn[item][state]
            rxn[item]['srt'+semirxn]+="*P"+rxn[item][state]
    return(G) 
     
     
def process_rxn(conf,itm,rxn,ltp):
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
        
    ltp['rxn']=[]  
     
    for item in sorted(rxn) : 
        # Initialize variables 
        rxn[item]['dGd']=""
        rxn[item]['aGd']=""
        rxn[item]['aGi']=""
        rxn[item]['kd']="k"+item+"d:=evalf("   
        rxn[item]['ki']="k"+item+"i:=evalf("    
        rxn[item]['rtd']="r"+item+":=(t)-> k"+item+"d"
        rxn[item]['rti']="-k"+item+"i"
        rxn[item]['srtd']="sr"+item+":= k"+item+"d" 
        rxn[item]['srti']="-k"+item+"i"
        howmanygasd=0
        howmanygasi=0        
        
        # Use the is/fs states for each reaction to get their activation energies. 
        # Then write the formula for reaction rate according to their intermediates. 
        #     This formula is splitt between rtd (direct part) and rti (inverse part). 
        # Then update the differential equations in which each adsorbed species participates. 
        # If a rectant is in gas phase, include a Hertz-Knudsen term in the constant 
        #     and its pressure as variable in the reaction rate.      
        Gi1=process_itm_on_rxn(conf,itm,rxn,item,'is1')  
        Gi2=process_itm_on_rxn(conf,itm,rxn,item,'is2')  
        Gf1=process_itm_on_rxn(conf,itm,rxn,item,'fs1')  
        Gf2=process_itm_on_rxn(conf,itm,rxn,item,'fs2')  
         
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
        ltp['rxn'].append("sr"+item)
           
    return(rxn,itm)
        
def printtxt(conf,itm,rxn,sbalance,initialc,sodesolv,rhsparse,ltp): 
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
    #print("restart: " )
     
    print("T:=", conf.get("Reactor","reactortemp"), " : " )
    for item in sorted(itm) : 
        if itm[item]['phase']=='gas' :  
            print('P'+item+":=",itm[item]['pressure']," : ") 
      
    print("\n# Kinetic constants")
    #for item in sorted(gas) :
    #    print(gas[item]['kads'+cat],gas[item]['kdes'+cat])
    for item in sorted(rxn) :
        print(rxn[item]['kd'],"\n",rxn[item]['ki']  )
      
    print("\n# Reaction rates:")
    #for item in sorted(gas) :
    #    print(gas[item]['rads'+cat])
    for item in sorted(rxn) :
        print(rxn[item]['rtd'],rxn[item]['rti'])
      
    print("\n# Site-balance equation: ")
    print(sbalance)
     
    print("\n# Differential equations: ")
    for item in sorted(itm) :
        if  itm[item]['phase']=='cat' and item!=conf["Reactor"]["sitebalancespecies"] : 
            print(itm[item]['diff']," : ")
      
    print("\n# Initial conditions: ")
    print(initialc)
      
    print("\n# SODE Solver: ")
    print(sodesolv)
       
    # Print labels, before the timeloop and flush. 
    print("\nfprintf(",conf['General']['mapleoutput'],',"%q %q\\n",','catalyst, "timei", "T",', 
          ', '.join(['"'+item+'"' for item in ltp['prs']]) ,",", 
          ', '.join(['"'+item+'"' for item in ltp['itm']]) ,",",
          ', '.join(['"'+item+'"' for item in ltp['rxn']]) ,   
          " ): \n " )  
    print('\nflush(',conf['General']['mapleoutput'],'): ')
         
    # Time control: 
    time1,timel=rxntime(conf)
    if timel : 
        print("for timei in " + str(time1) + " do ")
    else : 
        print("timei:= "+time1+" : ")
      
    print("S:=Solution(timei):")
    
    print("\n# Solution parser: ")
    print(rhsparse)
    
    print("\n# Site-balance equation after solver: ")
    tmp=os.popen(    "echo \""+sbalance+"\" | sed 's/(t)//g' "   ).read() # popen and read() are used to save in variable 
    tmp="s"+os.popen("echo \""+tmp[:-1]+"\" | sed 's/->//g' "    ).read() # rm last character (newline) 
    print(  os.popen("echo \""+tmp[:-1]+"\" | sed 's/-c/-sc/g' " ).read()[:-1] )
    
    print("\n# Reaction rates after solver: ")
    #for item in sorted(gas) :
    #    print(gas[item]['srads'+cat])
    for item in sorted(rxn) :
        print(rxn[item]['srtd'],rxn[item]['srti'],":")
                   
    # Print results 
    print("\nfprintf(",conf['General']['mapleoutput'],',"%q %q\\n",',conf['General']['catalyst'],', timei, T,',
          ', '.join([item for item in ltp['prs']]) ,",", 
          ', '.join([item for item in ltp['itm']]) ,",", 
          ', '.join([item for item in ltp['rxn']]) , 
          " ): " )  
       
    print('\nflush(',conf['General']['mapleoutput'],'): ')
     
    if timel :     
        print("\nod: \n ") 
    
    # Print close file instruction  
    print('\nclose(',conf['General']['mapleoutput'],'): \n\n')
     
    
