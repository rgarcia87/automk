import pandas as pd
import os
import configparser
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
    return conf 
      
def read(filename='./int.csv'): 
    """This function reads a file containing information of catalysts, gas, intermediates, or reactions. 
    It requires pandas to be loaded.  
    
    Args: 
        filename: Input file. The columns are separated by one or more spaces.
            Energies must be provided in eV and frequencies in cm-1 Format: 
            Label      E      frq     
            i0101   -1.205    [50,500,2000]   
            i0011   -0.701    [100,200,1000] 
    
    Returns: 
        dicint: a dictionary containing at least the tags, energies, and frequencies of all species. 
    
    """
    
    dic=pd.read_csv(filename, delim_whitespace=True, index_col='Label').T.to_dict()
    return(dic) 
     
def fint(conf,int,cat):
    """This function process the "intermediates" dataframe to generate 
    the site-balance equation, the SODE-solver, and the initial conditions as clean surface. 
    It also initializes the list of differential equations.  
    
    Args: 
        conf: Configuration data.
        int: Dict of dicts containing at least a list of intermediates as index. 
        cat: Name of the catalyst. Currently only string is supported. 
    
    Returns:
        int:      Expanded dict of dicts containing also the list of differential equations. (Mutable)
        sbalance: Site-balance equation. (Unmutable)
        sodesolv: Input of the SODE-solver. (Unmutable)   
        initialc: Initial conditions as clean surface. (Unmutable) 
    """
    
    # Initialize variables related to intermediates. 
    sbalance="c"+cat+":=(t)-> 1.0"  
    sodesolv="Solution:=dsolve({"   
    initialc="IC0:="
    rhsparse=""
    index=1
        
    # Process intermediates
    for item in sorted(int) : 
        
        # Initialize diff equations to count in which reactions each species participate.     
        int[item]['diff']="eqd"+item+":=diff(c"+item+"(t),t)="         
          
        # Initialize list of reactions in which each intermediate participate. 
        int[item]['gaslst']=[] 
        int[item]['rxnlst']=[] 
        
        # Prepare site balance  
        sbalance=sbalance+" -c"+item+"(t)"
        
        # Prepare list of differential equations for the SODE solver 
        sodesolv=sodesolv+"eqd"+item+", "
        
        # Prepare list of default initial conditions as clean surface
        initialc=initialc+" c"+item+"(0.0)=0.0,"
          
        # Prepare reader of concentrations from solver  
        index+=1 
        #print("index: ",index, "rhsparse: ", rhsparse , "sodesolv: ", sodesolv)
        rhsparse+="sc"+item+":=RHS["+str(index)+"] : "
        #print(rhsparse)  
        
    # Close the site-balance equation     
    sbalance=sbalance+":" 
    
    # Close the sodesolv
    sodesolv=sodesolv+"IC0}, numeric, method=rosenbrock, maxfun=0, abserr=1E-16, interr=false);"
    
    # In the initial conditions, replace the last comma by a colon 
    initialc=initialc[:-1]+" : "
      
    #print("\n",sbalance,"\n",initialc,"\n",sodesolv) 
    return(int,sbalance,sodesolv,initialc,rhsparse) 
     

 
def fgas(conf,gas,int,cat):
    """Subroutine that expands the "gas" dictionary of dictionaries to include 
    the kinetic constants and rates of adsorption/desorptions.  
    It also expands the list of differential equations in "int"
    and the list of adsorption/desorptions in which each intermediate participates (reaction-like). 
    
    Args: 
        conf: Configuration data. 
        gas: Dict of dicts containing the volatile species. (Mutable)
        int: Dict of dicts containing at least a list of intermediates as index. (Mutable)
        cat: Name of the catalyst. Currently only string is supported.
    
    Returns: 
        gas: Expanded dict of dicts containing adsorption/desorption constants and rates. (Mutable)
        int: Expanded dict of dicts with list of differential equations updated with adsorption/desorptions. (Mutable)
    """
    
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
        
    # Process adsorption and desorptions 
    for item in sorted(gas) : 
        
        # Get partial pressures 
        try : 
            gas[item]['pressure']=conf['Pressures'][item] 
        except : 
            gas[item]['pressure']=0  
          
        # Adsorption energy 
        DGads=int[item][cat]-gas[item][cat]
         
        # Activation energy of adsorption
        aGads=gas[item][cat]-gas[item]['gas'] 
         
        # Activation energy of desorption
        aGdes=gas[item][cat]-gas[item]['gas']
         
        # Kinetic constant: adsorption of gas to cat
        gas[item]['kads'+cat]="kads"+item+cat+":=evalf(101325*P"+item+"*exp(-max("+\
                              "0.00,"+"{:.6f}".format(aGads)+","+"{:.6f}".format(DGads)+\
                              ")/(0.861733E-4*T))/(1.7492150414*10^19*sqrt(1.6605390400*"+\
                              "(2*evalf(Pi)*1.3806485200)*10^(-23)*T*"+"{:.2f}".format(gas[item]['mw'])+\
                              "*10^(-27)))) : " 
         
        # Kinetic constant: desorption of volatile species from cat to gas. 
        gas[item]['kdes'+cat]="kdes"+item+cat+":=evalf("+kbh+"*T*exp(-max("+\
                              "0.00,"+"{:.6f}".format(aGdes)+","+"{:.6f}".format(-DGads)+\
                              ")/("+kbev+"*T)) ): "
        
        # Formula: adsorption/desorption rate of volatile adsorbates
        gas[item]['rads'+cat]="rads"+item+cat+":=(t)-> "+pdamp1+"kads"+item+cat+"*c"+cat+"(t)"+\
                              "-kdes"+item+cat+"*c"+item+"(t) : "  
         
        # Formula: adsorption/desorption rate of volatile adsorbates after solver
        gas[item]['srads'+cat]="srads"+item+cat+":= "   +pdamp2+"kads"+item+cat+"*sc"+cat+\
                               "-kdes"+item+cat+"*sc"+item+" : " 
             
        # Update differential equations for the species
        int[item]['diff']=int[item]['diff']+"+rads"+item+cat+"(t)"
         
        # Update list of adsorption/desorptions in which each intermediate participates. 
        int[item]['gaslst'].append(item)    
        # Example: int['P']['gaslst'] will append 'P' as item. 
        # This must be modified somewhat to differentiate the species in gas and adsorbed.  
        
        # Final formulaes
        #print("\n",gas[item]["kads"+cat],"\n\n",gas[item]["kdes"+cat],"\n\n",gas[item]['rads'+cat],"\n")
        
    return(gas,int)

def frxn(conf,rxn,int,cat):
    """Subroutine that expands the "rxn" dictionary of dictionaries to include 
    the kinetic constants and rates of all chemical reactions. 
    It also expands the list of differential equations in "int"
    and the list of reactions in which each intermediate participates. 
    
    Args: 
        conf: Configuration data.
        rxn: Dict of dicts containing the reactions. (Mutable)
        int: Dict of dicts containing at least a list of intermediates as index. (Mutable)
        cat: Name of the catalyst. Currently only string is supported.
    
    Returns: 
        rxn: Expanded dict of dicts containing adsorption/desorption constants and rates. (Mutable)
        int: Expanded dict of dicts with list of differential equations updated with chemical reactions. (Mutable)
    """ 
        
    for item in sorted(rxn) : 
        # Initialize variables 
        rxn[item]['aGd']=""
        rxn[item]['aGi']=""
        rxn[item]['kd']=""
        rxn[item]['ki']=""    
        rxn[item]['rtd']="r"+item+":=(t)-> k"+item+"d"
        rxn[item]['rti']="-k"+item+"i"
        rxn[item]['srtd']="sr"+item+":= k"+item+"d" 
        rxn[item]['srti']="-k"+item+"i"
        
        # Use the is/fs states for each reaction to get their activation energies. 
        # Then write the formula for reaction rate according to their intermediates. 
        #     This formula is splitt between rtd (direct part) and rti (inverse part). 
        # At the same time, do a list of reactions for each intermediate.     
        if rxn[item]['is1']=='None' :
            Gdi1=0.0
        else :   
            try: 
                Gdi1=int[rxn[item]['is1']][cat]
                rxn[item]['rtd']=rxn[item]['rtd']+"*c"+rxn[item]['is1']+"(t)"
                rxn[item]['srtd']=rxn[item]['srtd']+"*sc"+rxn[item]['is1'] 
                int[rxn[item]['is1']]['diff']=int[rxn[item]['is1']]['diff']+"-r"+item+"(t)"
            except: 
                Gdi1=0.0
                print("\n Error!, reaction ",item, " comes from IS1 ",rxn[item]['is1']," which was not found.")
         
        if rxn[item]['is2']=='None' :
            Gdi2=0.0
        else :        
            try: 
                Gdi2=int[rxn[item]['is2']][cat]
                rxn[item]['rtd']=rxn[item]['rtd']+"*c"+rxn[item]['is2']+"(t)"
                rxn[item]['srtd']=rxn[item]['srtd']+"*sc"+rxn[item]['is2']
                int[rxn[item]['is2']]['diff']=int[rxn[item]['is2']]['diff']+"-r"+item+"(t)"
            except: 
                Gdi2=0.0
                print("\n Error!, reaction ",item, " comes from IS2 ",rxn[item]['is2']," which was not found.")
         
        if rxn[item]['fs1']=='None' :
            Gdif=0.0
        else :
            try: 
                Gdf1=int[rxn[item]['fs1']][cat]
                rxn[item]['rti']=rxn[item]['rti']+"*c"+rxn[item]['fs1']+"(t)"
                rxn[item]['srti']=rxn[item]['srti']+"*sc"+rxn[item]['fs1'] 
                int[rxn[item]['fs1']]['diff']=int[rxn[item]['fs1']]['diff']+"+r"+item+"(t)"
            except: 
                Gdf1=0.0
                print("\n Error!, reaction ",item, " goes to FS1 ",rxn[item]['fs1']," which was not found.")
         
        if rxn[item]['fs2']=='None' :
            Gdf2=0.0
        else :        
            try: 
                Gdf2=int[rxn[item]['fs2']][cat]
                rxn[item]['rti']=rxn[item]['rti']+"*c"+rxn[item]['fs2']+"(t)"
                rxn[item]['srti']=rxn[item]['srti']+"*sc"+rxn[item]['fs2'] 
                int[rxn[item]['fs2']]['diff']=int[rxn[item]['fs2']]['diff']+"+r"+item+"(t)"
            except: 
                Gdf2=0.0        
                print("\n Error!, reaction ",item, " goes to FS2 ",rxn[item]['fs2']," which was not found.")
         
        # Close the formula with ":"
        rxn[item]['rti']=rxn[item]['rti']+" : "
          
        # Get reaction (dG) and (aG) activation energies for each reaction, both direct and inverse.  
        #print("\n",Gdi1,Gdi2,Gdf1,Gdf2)    
        rxn[item]['dGd']=Gdf1+Gdf2     -Gdi1-Gdi2    
        rxn[item]['aGd']=rxn[item][cat]-Gdi1-Gdi2
        rxn[item]['aGi']=rxn[item][cat]-Gdf1-Gdf2
        #print('\n \n' , rxn, '\n \n')    
         
        # Kinetic constant for each reaction. 
        rxn[item]['kd']="k"+item+"d:=evalf("+kbh+"*T*exp(-max(0.0,"+\
                        "{:.6f}".format( rxn[item]['aGd'])+","+\
                        "{:.6f}".format( rxn[item]['dGd'])+\
                        ")/("+kbev+"*T)) ): "
        rxn[item]['ki']="k"+item+"i:=evalf("+kbh+"*T*exp(-max(0.0,"+\
                        "{:.6f}".format( rxn[item]['aGi'])+","+\
                        "{:.6f}".format(-rxn[item]['dGd'])+\
                        ")/("+kbev+"*T)) ): "
        
        # Update list of reactions in which each intermediate participates. 
        try: 
            int[rxn[item]['is1']]['rxnlst'].append(item)    
        except: 
            pass 
        try: 
            int[rxn[item]['is2']]['rxnlst'].append(item)   
        except: 
            pass 
        try: 
            int[rxn[item]['fs1']]['rxnlst'].append(item)    
        except: 
            pass 
        try: 
            int[rxn[item]['fs2']]['rxnlst'].append(item)   
        except: 
            pass 
         
    return(rxn,int)
 
def printtxt(conf,gas,int,rxn,cat,sbalance,initialc,sodesolv,rhsparse): 
    """Subroutine that prints a given calculation for Maple 
      
    Args: 
        conf: Configuration data. 
            time1: Time or times for which the concentrations and reactions shall be printed (str/int/float, or list).
        gas: Dict of dicts containing molecules in gas phase. (Mutable) 
        int: Dict of dicts containing at least a list of intermediates as index. (Mutable)
        rxn: Dict of dicts containing the reactions. (Mutable)
        cat: Name of the catalyst. Currently only string is supported.
        sbalance: Site-balance equetion, string. 
        initialc: Initial conditions, string. 
        sodesolv: Calls SODE solver in Maple, string.  
        rhsparse: Parser of surface concentrations, string. 
    """
    
    print("# Heading " )
    print("restart: " )
    print("T:=", conf.get("Reactor","reactortemp"), " : " )
    for item in sorted(gas) : 
        print('P'+item+":=",gas[item]['pressure']," : ") 
    
    print("\n# Kinetic constants")
    for item in sorted(gas) :
        print(gas[item]['kads'+cat],gas[item]['kdes'+cat])
    for item in rxn :
        print(rxn[item]['kd'],  rxn[item]['ki']  )
    
    print("\n# Reaction rates:")
    for item in sorted(gas) :
        print(gas[item]['rads'+cat])
    for item in sorted(rxn) :
        print(rxn[item]['rtd'],rxn[item]['rti'])
    
    print("\n# Site-balance equation: ")
    print(sbalance)
    
    print("\n# Differential equations: ")
    for item in sorted(int) :
        print(int[item]['diff']," : ")
    
    print("\n# Initial conditions: ")
    print(initialc)
    
    print("\n# SODE Solver: ")
    print(sodesolv)
    
    time1=conf['Reactor']['time1'] 
    print("\n# Preparing postprocessing: ")
    if   type(time1) is str :
        print("timei:= "+time1+" : ")
    elif type(time1) is int :
        print("timei:= "+"g".format(time1)+" : ")  
    elif type(time1) is float :
        print("timei:= "+"g".format(time1)+" : ")
    elif type(time1) is list :
        print("for timei in " + str(time1) + " do ")
        # Lists are limited by [ ] and they should be printed as that. 
    else :
        print("\t Warning! time1 should be type string, float, integer, or list.")
        print("\t current type: "+type(time1))
        print("\t time1 contains: \n\t",time1) 
        print("\t automk abnormal termination. " ) 
        exit() 
    
    print("S:=Solution(timei):")
    
    print("\n# Solution parser: ")
    print(rhsparse)
    
    print("\n# Site-balance equation after solver: ")
    tmp=os.popen(    "echo \""+sbalance+"\" | sed 's/(t)//g' "   ).read() # popen and read() are used to save in variable 
    tmp="s"+os.popen("echo \""+tmp[:-1]+"\" | sed 's/->//g' "    ).read() # rm last character (newline) 
    print(  os.popen("echo \""+tmp[:-1]+"\" | sed 's/-c/-sc/g' " ).read()[:-1] )
    
    print("\n# Reaction rates after solver: ")
    for item in sorted(gas) :
        print(gas[item]['srads'+cat])
    for item in sorted(rxn) :
        print(rxn[item]['srtd'],rxn[item]['srti'])
    
    if type(time1) is list :  
        print("\nod: \n ") 
    
    print()
    
def printtxtpd(conf,gas,int,rxn,cat,sbalance,initialc,sodesolv,rhsparse): 
    """Subroutine that prints microkinetics in Maple removing an intermediate each time.  
        It cals subroutine printtxt. 
        
    Args: 
        conf: Configuration data. 
        gas: Dict of dicts containing molecules in gas phase. (Mutable) 
        int: Dict of dicts containing at least a list of intermediates as index. (Mutable)
        rxn: Dict of dicts containing the reactions. (Mutable)
        cat: Name of the catalyst. Currently only string is supported.
        sbalance: Site-balance equetion, string. 
        initialc: Initial conditions, string. 
        sodesolv: Calls SODE solver in Maple, string.  
        rhsparse: Parser of surface concentrations, string. 
    """ 
        
    for item in sorted(int) : 
        inttmp = copy.deepcopy(int)
        gastmp = copy.deepcopy(gas)  
        rxntmp = copy.deepcopy(rxn)  
         
        inttmp[item]['diff']="eqd"+item+":=diff(c"+item+"(t),t)=0.00 "
        #print(item) 
        #print(int[item]['gaslst']) 
        #print(int[item]['rxnlst']) 
        for jtem in int[item]['gaslst'] : 
            gastmp[jtem]['kads'+cat]="kads"+jtem+cat+":=evalf(0.0): "
            gastmp[jtem]['kdes'+cat]="kdes"+jtem+cat+":=evalf(0.0): "
            gastmp[jtem]['rads'+cat]="rads"+item+cat+":=(t)-> 0.00 :"
            gastmp[jtem]['srads'+cat]="srads"+jtem+cat+":=0.00 :" 
        for jtem in int[item]['rxnlst'] : 
            rxntmp[jtem]['kd']="k"+jtem+"d:=evalf(0.0): "   
            rxntmp[jtem]['ki']="k"+jtem+"i:=evalf(0.0): "   
            rxntmp[jtem]['rtd']="r"+jtem+":=(t)-> 0.00"
            rxntmp[jtem]['rti']=" : "
            rxntmp[jtem]['srtd']="sr"+jtem+":= 0.00"
            rxntmp[jtem]['srti']=" : "
        printtxt(conf,gastmp,inttmp,rxntmp,cat,sbalance,initialc,sodesolv,rhsparse)



