
import sys, os, math
import mol2

# this is written by Trent Balius & Reed Stein in the Shoichet Lab
# written in 2016 and 2017 
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
        print (t,a,b,r)
        #vdw_dict[t] = round(r,3)
        vdw_dict[t] = r
    return vdw_dict

def get_radius(type,vdw_dict):
    return vdw_dict[type]



def read_in_dx_file(file):

  fileh = open(file,'r')

  flag_read_dx = False

  count = 0

  values = []
  for line in fileh:
      splitline = line.split()

      #print splitline
      #if len(splitline) < 2:
      if len(splitline) == 0:
          #print line
          continue

      ## this should be line 1
      if (splitline[0] == "object" and splitline[1] == "1"):
         print ("count = ", count, " line = ", line)
         xn = int(splitline[5])
         yn = int(splitline[6])
         zn = int(splitline[7])
         
      ## this should be line 2       
      if (splitline[0] == "origin"):
         #print line
         print ("count = ", count, " line = ", line)
         origin = [float(splitline[1]), float(splitline[2]), float(splitline[3])] 

      ## this should be lines 3-5
      if (splitline[0] == "delta"):
         #print line
         print ("count = ", count, " line = ", line)
         if (float(splitline[2]) == 0 and  float(splitline[3]) ==0):
            dx = float(splitline[1]) 
         elif (float(splitline[1]) == 0 and  float(splitline[3]) ==0):
            dy = float(splitline[2]) 
         elif (float(splitline[1])== 0 and  float(splitline[2])==0):
            dz = float(splitline[3]) 
            print (dx, dy, dz )


      if (splitline[0] == "object" and splitline[1] == "2"):
         #print line
         print ("count = ", count, " line = ", line)
      if (splitline[0] == "object" and splitline[1] == "3"):
         #print line
         print ("count = ", count, " line = ", line)
         flag_read_dx = True
         continue # go to next line
      if (flag_read_dx):

         if (len(splitline) > 3): 
            print ("Error: dx formate problem. more than 3 colums")
            exit()

         for value in splitline:
             values.append(float(value))  

      count = count + 1


  print (len(values))
  fileh.close()
  return xn,yn,zn,dx,dy,dz,origin,values 

def write_out_dx_file(file,xn,yn,zn,dx,dy,dz,origin,values):

   
  fileh = open(file,'w')
#object 1 class gridpositions counts 40 40 40
#origin 35.31 27.576 18.265
#delta 0.5 0 0
#delta 0 0.5 0
#delta 0 0 0.5
#object 2 class gridconnections counts 40 40 40
#object 3 class array type float rank 0 items 64000 data follows

  fileh.write('object 1 class gridpositions counts %d %d %d\n' % (xn,yn,zn))
  fileh.write('origin %6.2f %6.2f %6.2f\n' % (origin[0],origin[1],origin[2]))
  fileh.write('delta %2.1f 0 0\n' % dx)
  fileh.write('delta 0 %2.1f 0\n' % dy)
  fileh.write('delta 0 0 %2.1f\n' % dz)
  fileh.write('object 2 class gridconnections counts %d %d %d\n' % (xn,yn,zn))
  fileh.write('object 3 class array type float rank 0 items %d data follows\n' % len(values))

  count = 1
  for value in values:
       if (value == 0.0): 
          fileh.write('%d' % 0)
       else:
          fileh.write('%f' % value)
       # print newline after 3rd number.
       if (count == 3): 
            fileh.write('\n')
            count = 0
       # print space after number but not at the end of the line.
       else:
            fileh.write(' ')
       count = count + 1

  # if the last line has less than 3 numbers then print the a newline.
  if (count < 3):
       fileh.write('\n')
  fileh.close()


def distance(v1,v2):
    if (len(v1)!=len(v2)):
       print ("error" )
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

def cal_gausian(dist,sigma):
    #print dist, sigma
    var = sigma**2.0
    val = 1.0/math.sqrt(2.0*math.pi*var) * math.exp(-1.0*(dist**2.0)/(2.0*var))
    return val

######################################################################
# this function calculates gist score. 
# also calculates peratom controbutions distributing overlap 
# atoms equaly among all atoms contating the point
# also calculates volume displaced by molcule.
######################################################################

def calc_score(fileprefix,values,gridscale,xn,yn,zn,origin,mol,vdw_dict,fileh,return_vals,Hchoice,div):
# file - log file
# values gist values 
# gridscale - the grid spacing


    mol_name = mol.name
    if Hchoice !=0 and Hchoice !=1 and Hchoice !=2: 
       print ("Hchoice must be 0(use standard radius), 1 (use 0.8 for nonpolar H, instead of 1.540), or 2 (ignore H)")
       exit()

    dt   = mol2.convert_sybyl_to_dock(mol) # get dock atom types.

    # loop over each atom. 

    dict_gridpoint = {} # for each grid point have a list of atoms it is in

    gaus_val = 0.0 
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
        print( atom.type, dt[atom_i], radius)

        grid_i = round((atom.X - origin[0] ) / gridscale)
        grid_j = round((atom.Y - origin[1] ) / gridscale)
        grid_k = round((atom.Z - origin[2] ) / gridscale)

        if (grid_i<0.0 or grid_j<0.0 or grid_k<0.0):
           print ("ERROR. . . ")
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
                        count2 = i*((yn)*(zn))+j*(zn) + k
                        print (count,count2,count-count2)
                        #if not (count in dict_gridpoint):
                        #   dict_gridpoint[count] = [atom_i]
                        #else:
                        #   dict_gridpoint[count].append(atom_i)
                    count = count+1

        '''
        atom_val = 0.0
        for i in range(int(start_i),int(stop_i)): # looping over the cube about the center of the atom 
            x = (i * gridscale) + origin[0] 
            for j in range(int(start_j),int(stop_j)):
                y = (j * gridscale) + origin[1] 
                for k in range(int(start_k),int(stop_k)):
                    z = (k * gridscale) + origin[2] 
                    dist = distance([x,y,z],[atom.X,atom.Y,atom.Z])
                    if (dist <= radius):
                    #if (dist < radius):
                        #count = i*((xn)**2)+j*(xn) + k
			#scale = cal_gausian(dist, radius/2.0)
			scale = cal_gausian(dist, radius/div)
                        count = i*((yn)*(zn))+j*(zn) + k
                        gaus_val = gaus_val + scale*values[count] 
                        atom_val = atom_val + scale*values[count] 
                        if not (count in dict_gridpoint):
                           dict_gridpoint[count] = [[atom_i, scale]]
                        else:
                           #print("append")
                           dict_gridpoint[count].append([atom_i, scale])
        print ("atom_val = %f"%atom_val)
    print ("gaus_val = %f"%gaus_val)
    voxel_vol =  gridscale**3.0
    print ("vol*gaus_val = %f"%(voxel_vol*gaus_val))

    sum_val = 0
    sum_val_positive = 0
    sum_val_negative = 0
    #voxel_vol =  gridscale**3.0
    # loop over all displaced grid points.
    for key in dict_gridpoint.keys():
        #print(key)
	key_sum = 0
        for atom_i, scale  in dict_gridpoint[key]:
		key_sum += -voxel_vol*values[key]*scale
		sum_val = sum_val + -1.*values[key]*scale
        #print ("key %d (%d) : val = %f, grid_val = %f"%(key,len(dict_gridpoint[key]),key_sum,values[key])) 

    print ("gist_val:", sum_val*voxel_vol)
    
    fileh.write('%f\n'%(sum_val*voxel_vol))

    new_values = []
    # make a new grid with only the voxels in the ligand are non-zerro

    if (return_vals):
        #new_values = []
        for i,val in enumerate(values):
            #if i in dict_gridpoint.keys():
            if i in dict_gridpoint:
               new_values.append(val)
            else:
               new_values.append(0.0)

    # write out the scale for visulization of blur
    new_weights = []
    if (return_vals):
        for i in range(len(values)):
            #if i in dict_gridpoint.keys():
            if i in dict_gridpoint:
               sum_scale = 0.0
               #new_values.append(1.0)
               for atom_i, scale  in dict_gridpoint[i]:
                   #print(i,atom_i, scale)
                   sum_scale = sum_scale + scale
               #print("sum",sum_scale)
               new_weights.append(sum_scale)
            else:
               new_weights.append(0.0)



    return new_values, new_weights

def main():


   if len(sys.argv) != 6: # if no input
       print ("ERORR:")
       print ("syntex: dx-gist_rescore.py dx-file mol2 hchoice div output")
       print ("dx-input-file input file in dx formate produed by gist, may be disities or energies")
       print ("mol2 containing docked poses. ")
       print ("Hchoice [0(use standard radius), 1 (use 0.8 for nonpolar H, instead of 1.540), or 2 (ignore H)]")
       print ("div is the value to divided the radius by to determine the standard deviation. ")
       print ("output prefix")
       return
 

   infiledx      = sys.argv[1]
   infilemol2    = sys.argv[2]
   Hchoice       = int(sys.argv[3])
   div           = float(sys.argv[4])
   outfile       = sys.argv[5]

   print (infiledx)
   print (infilemol2)

   #vdwdict = intialize_vdw_parm('/nfs/home/tbalius/zzz.github/DOCK/proteins/defaults/vdw.parms.amb.mindock') 
   #vdwdict = intialize_vdw_parm('/home/baliuste/zzz.github/DOCK/ucsfdock/proteins/defaults/vdw.parms.amb.mindock') 
   DOCKpath = os.getenv('DOCKBASE')
   vdwdict = intialize_vdw_parm(DOCKpath+'/proteins/defaults/vdw.parms.amb.mindock') 
   xn,yn,zn,dx,dy,dz,origin,values = read_in_dx_file(infiledx)

   gridscale = dx # assumes that they are all the same spaceing
   
   mols  = mol2.read_Mol2_file(infilemol2)

   file1 = open(outfile+'_gist_values.txt','w')

   N = len(mols)
   count = 0
   for mol in mols:
      file1.write(mol.name+"\t")
      new_values,weights = calc_score(outfile,values,gridscale,xn,yn,zn,origin, mol,vdwdict,file1,(N<3),Hchoice,div)
      if N < 3: # only print gist grids if there are a few poses in the mol2 file
       #  print "writting gist grid for overlap with ligand . . ."
         write_out_dx_file(outfile +str(count) +"new_gist.dx",xn,yn,zn,dx,dy,dz,origin,new_values)
         write_out_dx_file(outfile +str(count) +"new_weights.dx",xn,yn,zn,dx,dy,dz,origin,weights)
      count=count+1
      #exit()
   file1.close()
main()


