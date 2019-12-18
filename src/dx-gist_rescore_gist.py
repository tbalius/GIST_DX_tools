
import sys, os, math
import mol2
import dx_gist_lib as dxlib

# this is written by Trent Balius in the Shoichet Lab
# written in march, 2016
#

def intialize_vdw_parm(vdwfile):
    # this function reads in the vdw file. 
    fileh = open(vdwfile,'r')
    vdw_dict = {}
    for line in fileh:
        if (line[0] == '!'):
           continue
        splitline = line.split()
        t = int(splitline[0])
        a = float(splitline[1]) # A^(1/2)
        b = float(splitline[2]) # B^(1/2)
        #r = (((2A/B)^(1/2))^(1/3))/2
        if a == 0 or b == 0: 
           r = 0.5
        else: 
           r = (math.sqrt(2.0)*a/b)**(1.0/3.0)/2.0 # radius
        print t,a,b,r
        vdw_dict[t] = r
    return vdw_dict

def get_radius(type,vdw_dict):
    return vdw_dict[type]

def distance(v1,v2):
    if (len(v1)!=len(v2)):
       print "error" 
       exit()
    dist = 0.0
    for i in range(len(v1)):
        dist = dist + (v1[i]-v2[i])**2.0
    dist = math.sqrt(dist)
    return dist

def cal_grid(values,gridscale,xn,yn,zn,origin):
   # here the vector of values is transformed to a multidemitional array (grid)
    grid_old = []
    count = 0
    for i in range(xn):
        ydim = []
        for j in range(yn):
            zdem = []
            for k in range(zn):
                zdem.append(values[count])
                count = count + 1
            ydim.append(zdem)
        grid_old.append(ydim)
    return grid_old


######################################################################
# this function calculates gist score. 
# also calculates peratom controbutions distributing overlap 
# atoms equaly among all atoms contating the point
# also calculates volume displaced by molcule.
######################################################################

def calc_score(fileprefix,values,gridscale,xn,yn,zn,origin,mol,vdw_dict,fileh,return_vals,Hchoice):
# file - log file
# values gist values 
# gridscale - the grid spacing


    if Hchoice !=0 and Hchoice !=1 and Hchoice !=2: 
       print "Hchoice must be 0(use standard radius), 1 (use 0.8 for nonpolar H, instead of 1.540), or 2 (ignore H)"
       exit()

    dt   = mol2.convert_sybyl_to_dock(mol) # get dock atom types.

    # loop over each atom. 

    dict_gridpoint = {} # for each grid point have a list of atoms it is in

    
    sum_per_atom = []
    for atom_i,atom in enumerate(mol.atom_list):
        sum_per_atom.append(0.0)

        #if (Hchoice == 1 and dt[atom_i] == 6):
        if (Hchoice == 1 and dt[atom_i] == 7):
           radius = 0.8
        elif (Hchoice == 2 and (dt[atom_i] == 6 or dt[atom_i] == 7)):
           continue # go to nexed atom
        else:
           radius = vdw_dict[dt[atom_i]]
        #print atom.type, atom.X, atom.Y, atom.Z, atom.type, dt[atom_i], radius

        grid_i = round((atom.X - origin[0] ) / gridscale)
        grid_j = round((atom.Y - origin[1] ) / gridscale)
        grid_k = round((atom.Z - origin[2] ) / gridscale)

        if (grid_i<0.0 or grid_j<0.0 or grid_k<0.0):
           print "ERROR. . . "
           exit()

        radius_gridpoint = math.ceil(radius / gridscale) + 1
        #print radius, gridscale, radius_gridpoint
        #print  grid_i,grid_j, grid_k

        start_i = max([(grid_i - radius_gridpoint),0.0])
        stop_i  = min([(grid_i + radius_gridpoint),xn])
        start_j = max([(grid_j - radius_gridpoint),0.0])
        stop_j  = min([(grid_j + radius_gridpoint),yn])
        start_k = max([(grid_k - radius_gridpoint),0.0])
        stop_k  = min([(grid_k + radius_gridpoint),zn])

        '''
        print "centerpoint:",grid_i, grid_j, grid_k

        print "i", start_i, stop_i
        print "j", start_i, stop_i
        print "k", start_i, stop_i
        '''

        '''
        count = 0
        for i in range(xn):
            x = (i * gridscale) + origin[0]
            for j in range(yn):
                y = (j * gridscale) + origin[1]
                for k in range(zn):
                    z = (k * gridscale) + origin[2]
                    dist = distance([x,y,z],[atom.X,atom.Y,atom.Z])
                    if (dist <= radius):
                        if not (count in dict_gridpoint):
                           dict_gridpoint[count] = [atom_i]
                        else:
                           dict_gridpoint[count].append(atom_i)
                    count = count+1

        '''
        for i in range(int(start_i),int(stop_i)): # looping over the cube about the center of the atom 
            x = (i * gridscale) + origin[0] 
            for j in range(int(start_j),int(stop_j)):
                y = (j * gridscale) + origin[1] 
                for k in range(int(start_k),int(stop_k)):
                    z = (k * gridscale) + origin[2] 
                    dist = distance([x,y,z],[atom.X,atom.Y,atom.Z])
                    if (dist <= radius):
                        #count = i*((xn)**2)+j*(xn) + k
                        count = i*(yn*zn)+j*(zn) + k
                        if not (count in dict_gridpoint):
                           dict_gridpoint[count] = [atom_i]
                        else:
                           dict_gridpoint[count].append(atom_i)
    sum_val = 0
    sum_val_positive = 0
    sum_val_negative = 0
    voxel_vol =  gridscale**3.0
    # loop over all displaced grid points.
    for key in dict_gridpoint.keys():
        # sum up total value as well as positive, and negative contrabutions.   
        sum_val = sum_val + values[key]
        if values[key] > 0:
           sum_val_positive = sum_val_positive + values[key]
        else:
           sum_val_negative = sum_val_negative + values[key]

        # also calculate the per-atom gist decomposition. 
        # if a point is shared by multiple atoms, divied the value by the 
        # number of atoms that share that point and asign that fractional 
        # value to each atom
        N = len(dict_gridpoint[key])
        for atom_i in dict_gridpoint[key]:
            #print atom_i
            sum_per_atom[atom_i] = sum_per_atom[atom_i] + (values[key]/N)

    molN = len(dict_gridpoint.keys())
    boxN = len(values)
    boxV = xn*gridscale*yn*gridscale*zn*gridscale

    molV = float(molN)/float(boxN)*boxV

    print "gist_val:", sum_val*voxel_vol
    print "gist_val_positive:", sum_val_positive*voxel_vol
    print "gist_val_negative:", sum_val_negative*voxel_vol
    print "molN=",molN,"  boxN=",boxN,"  boxV=",boxV
    print "molV=",molV
    
    fileh.write('%s,%f\n'%("gist_val", sum_val*voxel_vol))
    fileh.write('%s,%f\n'%("gist_val_positive", sum_val_positive*voxel_vol))
    fileh.write('%s,%f\n'%("gist_val_negative", sum_val_negative*voxel_vol))
    fileh.write('%s,%f\n'%("molN",molN))
    fileh.write('%s,%f\n'%("boxN",boxN))
    fileh.write('%s,%f\n'%("boxV",boxV))
    fileh.write('%s,%f\n'%("molV",molV))

    for i,atom in enumerate(mol.atom_list):
        print i,":",sum_per_atom[i]*voxel_vol
        fileh.write('atom%d,%f\n'%(i+1,sum_per_atom[i]*voxel_vol))


    new_values = []
    # make a new grid with only the voxels in the ligand are non-zerro

    if (return_vals):
        new_values = []
        for i,val in enumerate(values):
            #if i in dict_gridpoint.keys():
            if i in dict_gridpoint:
               new_values.append(val)
            else:
               new_values.append(0.0)

    return new_values

def main():


   if len(sys.argv) != 4: # if no input
       print "ERORR:"
       print "syntex: dx-gist_rescore.py dx-file mol2"
       print "dx-input-file input file in dx formate produed by gist, may be disities or energies"
       print "mol2 containing docked poses. "
       print "Hchoice [0(use standard radius), 1 (use 0.8 for nonpolar H, instead of 1.540), or 2 (ignore H)]"
       return
 

   infiledx      = sys.argv[1]
   infilemol2    = sys.argv[2]
   Hchoice       = int(sys.argv[3])
   outfile       = "out"

   print infiledx
   print infilemol2

   DOCKpath = os.getenv('DOCKBASE')
   vdwdict = intialize_vdw_parm(DOCKpath+'/proteins/defaults/vdw.parms.amb.mindock') 
   xn,yn,zn,dx,dy,dz,origin,values = dxlib.read_in_dx_file(infiledx)

   gridscale = dx # assumes that they are all the same spaceing
   
   
   mols  = mol2.read_Mol2_file(infilemol2)

   file1 = open(outfile+'gist_values.txt','w')

   N = len(mols)
   count = 0
   for mol in mols:
      new_values = calc_score(outfile,values,gridscale,xn,yn,zn,origin, mol,vdwdict,file1,(N<3),Hchoice)
      if N < 3: # only print gist grids if there are a few poses in the mol2 file
         print "writting gist grid for overlap with ligand . . ."
         dxlib.write_out_dx_file(outfile +str(count) +"new_gist.dx",xn,yn,zn,dx,dy,dz,origin,new_values)
      count=count+1
      #exit()
   file1.close()
main()


