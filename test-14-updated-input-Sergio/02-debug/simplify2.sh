#/bin/bash 

sed -i "s/g000000/i000000/g" $1 

sed -i "/i4/d" $1 
sed -i "/i3/d" $1 
#sed -i "/i2/d" $1 

sed -i "/g4/d" $1  
sed -i "/g3/d" $1  
#sed -i "/g2/d" $1 

