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
     
def fint(conf,int,ltp):
    """This function process the "intermediates" dataframe to generate 
    the site-balance equation, the SODE-solver, and the initial conditions as clean surface. 
    It also initializes the list of differential equations.  
    
    Args: 
        conf: Configuration data.
        int: Dict of dicts containing at least a list of intermediates as index. 
    
    Returns:
        int:      Expanded dict of dicts containing also the list of differential equations. (Mutable)
        sbalance: Site-balance equation. (Unmutable)
        sodesolv: Input of the SODE-solver. (Unmutable)   
        initialc: Initial conditions as clean surface. (Unmutable) 
    """
    
    # Initialize variables related to intermediates. 
    if conf.get("Reactor","sitebalancespecies")
    sbalance="c"+conf.get("Reactor","sitebalancespecies")+":=(t)-> 1.0"  
    sodesolv="Solution:=dsolve({"   
    initialc="IC0:="
    rhsparse=""
    index=1
    ltp['int']=""
        
    # Process intermediates
    for item in sorted(int) : 
        if  int[item]['phase']=='cat' and int[item]!=conf.get("Reactor","sitebalancespecies") :  
               
            # Initialize diff equations to count in which reactions each species participate.     
            int[item]['diff']="eqd"+item+":=diff(c"+item+"(t),t)="         
             
            # Prepare site balance  
            sbalance=sbalance+" -c"+item+"(t)"
             
            # Prepare list of differential equations for the SODE solver 
            sodesolv=sodesolv+"eqd"+item+", "
             
            # Prepare list of default initial conditions as clean surface
            initialc=initialc+" c"+item+"(0.0)=0.0,"
              
            # Prepare parser of concentrations after SODE is solved  
            index+=1 
            #print("index: ",index, "rhsparse: ", rhsparse , "sodesolv: ", sodesolv)
            rhsparse+="sc"+item+":=rhs(S["+str(index)+"]) : "
            #print(rhsparse)  
             
            # List of reactions for fprintf function in Maple 
            ltp['int']+="sc"+item+", " 
              
            # Initialize list of reactions in which each intermediate participate. 
            # Deprecated: No longer needed.              
            #int[item]['rxnlst']=[] 
                         
    # Close the site-balance equation     
    sbalance=sbalance+":" 
     
    # Close the sodesolv
    sodesolv=sodesolv+"IC0}, numeric, method=rosenbrock, maxfun=0, abserr=1E-16, interr=false);"
     
    # In the initial conditions, replace the last comma by a colon 
    initialc=initialc[:-1]+" : "
      
    #print("\n",sbalance,"\n",initialc,"\n",sodesolv) 
    return(int,sbalance,sodesolv,initialc,rhsparse) 
     
     
def frxn(conf,int,rxn,cat,ltp):
    """Subroutine that expands the "rxn" dictionary of dictionaries to include 
    the kinetic constants and rates of all chemical reactions. 
    It also expands the list of differential equations in "int"
    and the list of reactions in which each intermediate participates. 
    
    Args: 
        conf: Configuration data.
        rxn: Dict of dicts containing the reactions. (Mutable)
        int: Dict of dicts containing at least a list of intermediates as index. (Mutable)
    
    Returns: 
        rxn: Expanded dict of dicts containing adsorption/desorption constants and rates. (Mutable)
        int: Expanded dict of dicts with list of differential equations updated with chemical reactions. (Mutable)
    """ 
        
    ltp['rxn']=""  
     
    for item in sorted(rxn) : 
        # Initialize variables 
        rxn[item]['aGd']=""
        rxn[item]['aGi']=""
        rxn[item]['kd']=""
        rxn[item]['ki']=""
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
        # If a rectant is in gas phase, include a Hertz-Knudsen term in the constant.    
        if rxn[item]['is1']=='None' :
            Gdi1=0.0
        else:   
            try:
                Gdi1=int[rxn[item]['is1']]['G']
                if int[item]['phase']=='int':  
                    rxn[item]['rtd']=rxn[item]['rtd']+"*c"+rxn[item]['is1']+"(t)"
                    rxn[item]['srtd']=rxn[item]['srtd']+"*sc"+rxn[item]['is1'] 
                    int[rxn[item]['is1']]['diff']=int[rxn[item]['is1']]['diff']+"-r"+item+"(t)"
                elif int[item]['phase']=='gas':
                    howmanygasd+=1 
                    rxn[item]['rtd']=rxn[item]['rtd']+"*p"+rxn[item]['is1']
                    rxn[item]['srtd']=rxn[item]['srtd']+"*p"+rxn[item]['is1'] 
            except: 
                Gdi1=0.0
                print("\n Error!, reaction ",item, " comes from IS1 ",rxn[item]['is1']," which was not found.")
         
        if rxn[item]['is2']=='None' :
            Gdi2=0.0
        else:        
            try:
                Gdi2=int[rxn[item]['is2']]['G']
                if int[item]['phase']=='int':  
                    rxn[item]['rtd']=rxn[item]['rtd']+"*c"+rxn[item]['is2']+"(t)"
                    rxn[item]['srtd']=rxn[item]['srtd']+"*sc"+rxn[item]['is2']
                    int[rxn[item]['is2']]['diff']=int[rxn[item]['is2']]['diff']+"-r"+item+"(t)"
                elif int[item]['phase']=='gas':
                    howmanygasd+=1 
                    rxn[item]['rtd']=rxn[item]['rtd']+"*p"+rxn[item]['is2']
                    rxn[item]['srtd']=rxn[item]['srtd']+"*p"+rxn[item]['is2']
            except: 
                Gdi2=0.0
                print("\n Error!, reaction ",item, " comes from IS2 ",rxn[item]['is2']," which was not found.")
         
        if rxn[item]['fs1']=='None' :
            Gdif=0.0
        else:
            try: 
                Gdf1=int[rxn[item]['fs1']]['G']
                if int[item]['phase']=='int': 
                    rxn[item]['rti']=rxn[item]['rti']+"*c"+rxn[item]['fs1']+"(t)"
                    rxn[item]['srti']=rxn[item]['srti']+"*sc"+rxn[item]['fs1'] 
                    int[rxn[item]['fs1']]['diff']=int[rxn[item]['fs1']]['diff']+"+r"+item+"(t)"
                elif int[item]['phase']=='gas': 
                    howmanygasi+=1 
                    rxn[item]['rti']=rxn[item]['rti']+"*p"+rxn[item]['fs1']
                    rxn[item]['srti']=rxn[item]['srti']+"*p"+rxn[item]['fs1'] 
            except: 
                Gdf1=0.0
                print("\n Error!, reaction ",item, " goes to FS1 ",rxn[item]['fs1']," which was not found.")
         
        if rxn[item]['fs2']=='None' :
            Gdf2=0.0
        else:        
            try: 
                Gdf2=int[rxn[item]['fs2']]['G']
                if int[item]['phase']=='int':
                    rxn[item]['rti']=rxn[item]['rti']+"*c"+rxn[item]['fs2']+"(t)"
                    rxn[item]['srti']=rxn[item]['srti']+"*sc"+rxn[item]['fs2'] 
                    int[rxn[item]['fs2']]['diff']=int[rxn[item]['fs2']]['diff']+"+r"+item+"(t)"
                elif int[item]['phase']=='gas':
                    howmanygasi+=1 
                    rxn[item]['rti']=rxn[item]['rti']+"*p"+rxn[item]['fs2']
                    rxn[item]['srti']=rxn[item]['srti']+"*p"+rxn[item]['fs2'] 
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
        # If there are no species in gas-phase, multiply by kB*T/h. 
        # In any case, multiply by exp(-Ga/kB*T) (Arrhenius Eq, or sticking coefficient).  
        # Not yet implemented: reactions at interface and diffusions. 
        if   howmanygasd=0 :   
            rxn[item]['kd']+=kbh+"*T*exp(-max(0.0,"+\
                            "{:.6f}".format( rxn[item]['aGd'])+","+\
                            "{:.6f}".format( rxn[item]['dGd'])+\
                            ")/("+kbev+"*T)) ): "      
        elif howmanygasd=1 :
            rxn[item]['kd']+="*exp(-max(0.0,"+\
                            "{:.6f}".format( rxn[item]['aGd'])+","+\
                            "{:.6f}".format( rxn[item]['dGd'])+\
                            ")/("+kbev+"*T)) ): "   
        else : 
            print("WARNING! direct reaction #",item,"has",howmanygasd,"gas/aq reactants.") 
            print("Abnormal termination") 
            exit()    
          
        if   howmanygasi=0 :   
            rxn[item]['ki']+=kbh+"*T*exp(-max(0.0,"+\
                            "{:.6f}".format( rxn[item]['aGi'])+","+\
                            "{:.6f}".format(-rxn[item]['dGd'])+\
                            ")/("+kbev+"*T)) ): "             
        elif howmanygasi=1 :
            rxn[item]['ki']+="*exp(-max(0.0,"+\
                            "{:.6f}".format( rxn[item]['aGi'])+","+\
                            "{:.6f}".format(-rxn[item]['dGd'])+\
                            ")/("+kbev+"*T)) ): "   
        else : 
            print("WARNING! reverse reaction #",item,"has",howmanygasd,"gas/aq reactants.") 
            print("Abnormal termination") 
            exit()           
        
        # List of reactions for fprintf function in Maple 
        ltp['rxn']+="sr"+item+", "
         
    return(rxn,int)
      
def printtxt(conf,int,rxn,sbalance,initialc,sodesolv,rhsparse,ltp): 
    # Before called printtxtsr
    """Subroutine that prints a given calculation for Maple, just a 's'ingle 'r'un 
      
    Args: 
        conf: Configuration data. 
            time1: Time or times for which the concentrations and reactions shall be printed (str/int/float, or list).
        int: Dict of dicts listing the intermediates by an index. (Mutable)
        rxn: Dict of dicts for reactions. (Mutable)
        sbalance: Site-balance equetion, string. 
        initialc: Initial conditions, string. 
        sodesolv: Calls SODE solver in Maple, string.  
        rhsparse: Parser of surface concentrations, string. 
    """
    
    print("# Heading " )
    #print("restart: " )
    print("catal:=\""+cat+"\": ") 
     
    print("T:=", conf.get("Reactor","reactortemp"), " : " )
    for item in sorted(gas) : 
        print('P'+item+":=",gas[item]['pressure']," : ") 
    
    print("\n# Kinetic constants")
    for item in sorted(gas) :
        print(gas[item]['kads'+cat],gas[item]['kdes'+cat])
    for item in sorted(rxn) :
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
    for item in sorted(gas) :
        print(gas[item]['srads'+cat])
    for item in sorted(rxn) :
        print(rxn[item]['srtd'],rxn[item]['srti'],":")
    
    print("\nfprintf(",conf['General']['mapleoutput'],',"%q %q\\n", catal, T,',ltp['prs'],"timei,","sc"+cat+",",ltp['int'],ltp['gas'],ltp['rxn'][:-2]," ): ")
      
    if timel :     
        print("\nod: \n ") 
      
    print()
    
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
    
