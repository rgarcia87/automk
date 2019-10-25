import pandas as pd

#Constants 
kbh="20836612225.1252"    # Boltzmann constant divided by Planck constant, s^-1, string.  
kbev="8.617333262145E-5"  # Boltzmann constant in eV·K−1, string.

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

def fint(int,cat):
    """This function process the "intermediates" dataframe to generate 
    the site-balance equation, the SODE-solver, and the initial conditions as clean surface. 
    It also initializes the list of differential equations.  
    
    Args: 
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
        
        # Prepare site balance  
        sbalance=sbalance+" -c"+item+"(t)"
        
        # Prepare list of differential equations for the SODE solver 
        sodesolv=sodesolv+"eqd"+item+", "
        
        # Prepare list of default initial conditions as clean surface
        initialc=initialc+" c"+item+"(0.0)=0.0,"
          
        # Prepare reader of concentrations from solver  
        index+=1 
        print("index: ",index, "rhsparse: ", rhsparse , "sodesolv: ", sodesolv)
        rhsparse+=rhsparse+"sc"+item+":=RHS["+str(index)+"] : "
        print(rhsparse)  
        
    # Close the site-balance equation     
    sbalance=sbalance+":" 
    
    # Close the sodesolv
    sodesolv=sodesolv+"IC0}, numeric, method=rosenbrock, maxfun=0, abserr=1E-16, interr=false);"
    
    # In the initial conditions, replace the last comma by a colon 
    initialc=initialc[:-1]+" : "
      
    #print("\n",sbalance,"\n",initialc,"\n",sodesolv) 
    return(int,sbalance,sodesolv,initialc,rhsparse) 
     

def fgas(gas,int,cat):
    """Subroutine that expands the "gas" dictionary of dictionaries to include 
    the kinetic constants and rates of adsorption/desorptions, 
    as well as expanding the list of differential equations in "int"
    
    Args: 
        gas: Dict of dicts containing the volatile species. (Mutable)
        int: Dict of dicts containing at least a list of intermediates as index. (Mutable)
        cat: Name of the catalyst. Currently only string is supported.
    
    Returns: 
        gas: Expanded dict of dicts containing adsorption/desorption constants and rates. (Mutable)
        int: Expanded dict of dicts with list of differential equations updated with adsorption/desorptions. (Mutable)
    """
    
    # Process adsorption and desorptions 
    for item in sorted(gas) : 
          
        # Adsorption energy 
        DGads=int[item][cat]-gas[item][cat]
         
        # Activation energy of adsorption
        aGads=gas[item][cat]-gas[item]['gas'] 
         
        # Activation energy of desorption
        aGdes=gas[item][cat]-gas[item]['gas']
         
        # Kinetic constant: adsorption of gas to cat
        gas[item]["kads"+cat]="kads"+item+cat+":=evalf(101325*P"+item+"*exp(-max("+\
                              "0.00,"+"{:.6f}".format(aGads)+","+"{:.6f}".format(DGads)+\
                              ")/(0.861733E-4*T))/(1.7492150414*10^19*sqrt(1.6605390400*"+\
                              "(2*evalf(Pi)*1.3806485200)*10^(-23)*T*"+"{:.2f}".format(gas[item]['mw'])+\
                              "*10^(-27)))) : " 
         
        # Kinetic constant: desorption of volatile species from cat to gas. 
        gas[item]["kdes"+cat]="kdes"+item+cat+":=evalf("+kbh+"*T*exp(-max("+\
                              "0.00,"+"{:.6f}".format(aGdes)+","+"{:.6f}".format(-DGads)+\
                              ")/("+kbev+"*T)) ): "
        
        # Formula: adsorption/desorption rate of volatile adsorbates
        gas[item]['rads'+cat]="rads"+item+cat+":=(t)-> (1-exp(-1*t))*kads"+item+cat+"*c"+cat+"(t)"+\
                              "-kdes"+item+cat+"*c"+item+"(t) : "  
         
        # Formula: adsorption/desorption rate of volatile adsorbates after solver
        gas[item]['srads'+cat]="srads"+item+cat+":= (1-exp(-1*itime))*kads"+item+cat+"*sc"+cat+\
                               "-kdes"+item+cat+"*sc"+item+" : " 
             
        # Update differential equations for the species
        int[item]['diff']=int[item]['diff']+"+rads"+item+cat+"(t)"
        
        # Final formulaes
        #print("\n",gas[item]["kads"+cat],"\n\n",gas[item]["kdes"+cat],"\n\n",gas[item]['rads'+cat],"\n")
        
    return(gas,int)

def frxn(rxn,int,cat):
    """Subroutine that expands the "rxn" dictionary of dictionaries to include 
    the kinetic constants and rates of all chemical reactions,   
    as well as expanding the list of differential equations in "int"
    
    Args: 
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
         
        # Formula: Chemical reactions. The adsorption and desorptions are considered before. 
        #print("\n",rxn[item]['rtd']+rxn[item]['rti'],"\n",rxn[item]['kd'],rxn[item]['ki'])    
         
    return(rxn,int)
 
#def printtxt(

