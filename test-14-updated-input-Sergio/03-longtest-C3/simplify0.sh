#/bin/bash 

sed -i "s/g000000/i000000/g" $1 

sed -i "/i4/d" $1 
#sed -i "/i33/d" $1 
#sed -i "/i32/d" $1 
#sed -i "/i31/d" $1 
#sed -i "/i30/d" $1
#sed -i "/i22/d" $1 
#sed -i "/i21/d" $1 
#sed -i "/i20/d" $1 

sed -i "/g4/d" $1  
#sed -i "/g33/d" $1  
#sed -i "/g32/d" $1  
#sed -i "/g31/d" $1  
#sed -i "/g30/d" $1
#sed -i "/i22/d" $1 
#sed -i "/i21/d" $1 
#sed -i "/i20/d" $1

