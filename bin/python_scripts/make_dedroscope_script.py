x = ["255 0 0","255 153 0","255 222 0","0 255 0","0 0 255",
     "102 0 153","123 74 18","0 0 0","204 204 204","203 255 250","31 61 12","255 0 128"]

y = [0,1,2,3,4,5,6,7,8,9,10,11]


t = "find searchtext={0};select induced network ;set color={1};set labelcolor={1};deselect all;"

for ex,ey in zip(x,y): 
    print(t.format(ey,ex))





