def fgas(conf,gas,ltp):
    """Subroutine that expands the "ints" dictionary of dictionaries to include 
    the kinetic constants and rates of adsorption/desorptions of gas-phase species.  
    It also expands the list of differential equations in adsorbed species
    and the list of adsorption/desorptions in which each intermediate participates (reaction-like). 
    
    Args: 
        conf: Configuration data. 
        int: Dict of dicts containing at least a list of intermediates as index. (Mutable)
    
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
        
    ltp['prs']=""
    ltp['gas']=""
     
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
        
        # List of gas-phase species for the fprintf Maple function 
        ltp['prs']+="P"+item+", " 
        ltp['gas']+="srads"+item+cat+", "   
        
    return(gas,int)
