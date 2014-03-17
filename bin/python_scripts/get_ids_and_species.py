


unip_file = "laccases_329.txt"


phylo_names = set( [ "Viridiplantae","Bacteria","Fungi","Metazoa" ] )


for line in open(unip_file,'r'):

    spl = line.split()
    
    if spl[0] == "ID":
        ac = spl[1]
    
    if spl[0] == "OC":

        phy_set = set(line.replace(";","").strip().split()[1:])
        is_name = phylo_names & phy_set        


        if is_name != set([]):        
            print( ac,list(is_name)[0]  )




