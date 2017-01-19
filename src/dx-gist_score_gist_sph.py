
import sys, os, math
import mol2
import sph_lib

# this is written by Trent Balius in the Shoichet Lab
# written in march, 2016
#


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
         print "count = ", count, " line = ", line
         xn = int(splitline[5])
         yn = int(splitline[6])
         zn = int(splitline[7])
         
      ## this should be line 2       
      if (splitline[0] == "origin"):
         #print line
         print "count = ", count, " line = ", line
         origin = [float(splitline[1]), float(splitline[2]), float(splitline[3])] 

      ## this should be lines 3-5
      if (splitline[0] == "delta"):
         #print line
         print "count = ", count, " line = ", line
         if (float(splitline[2]) == 0 and  float(splitline[3]) ==0):
            dx = float(splitline[1]) 
         elif (float(splitline[1]) == 0 and  float(splitline[3]) ==0):
            dy = float(splitline[2]) 
         elif (float(splitline[1])== 0 and  float(splitline[2])==0):
            dz = float(splitline[3]) 
            print dx, dy, dz 


      if (splitline[0] == "object" and splitline[1] == "2"):
         #print line
         print "count = ", count, " line = ", line
      if (splitline[0] == "object" and splitline[1] == "3"):
         #print line
         print "count = ", count, " line = ", line
         flag_read_dx = True
         continue # go to next line
      if (flag_read_dx):

         if (len(splitline) > 3): 
            print "Error: dx formate problem. more than 3 colums"
            exit()

         for value in splitline:
             values.append(float(value))  

      count = count + 1


  print len(values)
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
  fileh.write('origin %6.3f %6.3f %6.3f\n' % (origin[0],origin[1],origin[2]))
  fileh.write('delta %6.3f 0 0\n' % dx)
  fileh.write('delta 0 %6.3f 0\n' % dy)
  fileh.write('delta 0 0 %6.3f\n' % dz)
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

def calc_score(fileprefix,values,gridscale,xn,yn,zn,origin,sphs,fileh,return_vals):
# file - log file
# values gist values 
# gridscale - the grid spacing


    # loop over each sphere. 

    dict_gridpoint = {} # for each grid point have a list of atoms it is in

    
    sum_per_atom = []
    for atom_i,atom in enumerate(sphs):
        sum_per_atom.append(0.0)

        radius = atom.radius
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
        stop_j  = min([(grid_j + radius_gridpoint),xn])
        start_k = max([(grid_k - radius_gridpoint),0.0])
        stop_k  = min([(grid_k + radius_gridpoint),xn])

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

    for i,atom in enumerate(sphs):
        print i,":",sum_per_atom[i]*voxel_vol
        fileh.write('atom%d,%f\n'%(i+1,sum_per_atom[i]*voxel_vol))


    new_values = []
    # make a new grid with only the voxels in the ligand are non-zerro

    if (return_vals):
        new_values = []
        for i,val in enumerate(values):
            #print i 
            if i in dict_gridpoint:
            #if i in dict_gridpoint.keys():
               new_values.append(val)
            else:
               new_values.append(0.0)

    return new_values

def main():


   if len(sys.argv) != 3: # if no input
       print "ERORR:"
       print "syntex: dx-gist_rescore.py dx-file sphere_file"
       print "dx-input-file input file in dx formate produed by gist, may be disities or energies"
       print "sphere file that defines the binding site (low dielectric sphere recomended) "
       return
 

   infiledx      = sys.argv[1]
   infilesph     = sys.argv[2]
   outfile       = "out"

   print infiledx
   print infilesph

   xn,yn,zn,dx,dy,dz,origin,values = read_in_dx_file(infiledx)

   gridscale = dx # assumes that they are all the same spaceing
   
   print gridscale ,xn,yn,zn,dx,dy,dz,origin   
   sphs  = sph_lib.read_sph(infilesph,'A','A')

   file1 = open(outfile+'gist_values.txt','w')

   count = 0
   new_values = calc_score(outfile,values,gridscale,xn,yn,zn,origin, sphs, file1, True)
   write_out_dx_file(outfile +str(count) +"new_gist.dx",xn,yn,zn,dx,dy,dz,origin,new_values)
   file1.close()
main()


