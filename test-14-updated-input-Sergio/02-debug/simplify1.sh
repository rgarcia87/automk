#/bin/bash 

sed -i "s/g000000/i000000/g" $1 

sed -i "/i44/d" $1 
sed -i "/i43/d" $1 
sed -i "/i42/d" $1 
sed -i "/i41/d" $1 
sed -i "/i40/d" $1 
sed -i "/i33/d" $1 
sed -i "/i32/d" $1 
sed -i "/i31/d" $1 
sed -i "/i30/d" $1
#sed -i "/i22/d" $1 
#sed -i "/i21/d" $1 
#sed -i "/i20/d" $1 

sed -i "/g44/d" $1  
sed -i "/g43/d" $1  
sed -i "/g42/d" $1  
sed -i "/g41/d" $1  
sed -i "/g40/d" $1  
sed -i "/g33/d" $1  
sed -i "/g32/d" $1  
sed -i "/g31/d" $1  
sed -i "/g30/d" $1
#sed -i "/i22/d" $1 
#sed -i "/i21/d" $1 
#sed -i "/i20/d" $1

