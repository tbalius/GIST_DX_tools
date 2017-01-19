
import sys, os, math

# this is written by Trent Balius in the Shoichet Lab
# written in Oct, 2014
#
# here we read in, manpulate and output dx file.
# the dx file are produed by GIST calculations.

# We place an sphere at each grid point. 
# rows and columns for which the sphere is not completely 
# acomicated in the grid are removed. 
# all points which are containd in the sphere are sumed and stored that the grid point
# this will per-compute the atomic displacement. 
# the space is no longer non-overlaping partiationing of the space.
# the spheres overlap.  

## THIS IS NOT COMPLETE.  

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
       
######################################################################
# this function compute the atomic displacement grids.  
######################################################################

def pre_compute(file,xn,yn,zn,dx,dy,dz,origin,values,radius):
# file - log file
# values gist values 
# gridscale - the grid spacing
# radius - sphere radius

    fileh = open(file,'w')


    padx = math.ceil(radius/dx)
    pady = math.ceil(radius/dy)
    padz = math.ceil(radius/dz)

    if not ( padx == pady and pady == padz and padx == padz):
       print "dx, dy and dz are not the same: ", dx, dy, dz
       exit()

    pad = int(padx)
    gspace = dx

    # convert array to grid

    #   inialize grid
    grid = []
    for i in range(xn):
        matrix = []
        for j in range(yn):
            row = []
            for k in range(zn):
                row.append(0)
            matrix.append(row)
        grid.append(matrix)

    index = 0
    for i in range(xn):
        for j in range(yn):
            for k in range(zn):
                grid[i][j][k] = values[index]
                index = index + 1
    

    # transform grid to precomput grid
    ngrid = []
    for i in range(xn-2*pad):
        matrix = []
        for j in range(yn-2*pad):
            row = []
            for k in range(zn-2*pad):
                row.append(0)
            matrix.append(row)
        ngrid.append(matrix)

    for i in range(xn-2*pad):
        for j in range(yn-2*pad):
            for k in range(zn-2*pad):
                oldgi = i + pad 
                oldgj = j + pad 
                oldgk = k + pad
                centerind = [oldgi, oldgj, oldgk]
                points = get_points_in_Sph(origin,xn,yn,zn,gspace,centerind,radius) 
                for p in points:
                    ngrid[i][j][k] = ngrid[i][j][k]  + grid[p[0]][p[1]][p[2]]

    # convert new grid to array
    nvalues = []

    index = 0
    for i in range(xn-2*pad):
        for j in range(yn-2*pad):
            for k in range(zn-2*pad):
                  nvalues.append(ngrid[i][j][k])

    new_xn = xn-2*pad
    new_yn = yn-2*pad
    new_zn = zn-2*pad
    
    new_origin = [ 0.0, 0.0, 0.0]
    new_origin[0] = origin[0] + gspace*pad
    new_origin[1] = origin[1] + gspace*pad
    new_origin[2] = origin[2] + gspace*pad

    fileh.close()
    return new_xn,new_yn,new_zn,new_origin,nvalues

######################################
def InSPH(center,radius,point):

    dist2 = 0.0
    for i in range(3):
        dist2 = dist2 + (center[i] - point[i])**2

    if dist2 == 0.0:
        #print center,  point 
        return False

    return (radius**2 > dist2)
    #return (radius**2 >= dist2)
    #if radius**2 > dist2:
    #   #print "%6.3f > %6.3f" % (radius**2,  dist2)
    #   return True
    #else:
    #   return False
    
######################################
def get_neighbors(cgpind):

    # cgpind - grid point index 

    n1 = [ cgpind[0], cgpind[1], cgpind[2] + 1 ]
    n2 = [ cgpind[0], cgpind[1], cgpind[2] - 1 ]
    n3 = [ cgpind[0], cgpind[1] + 1, cgpind[2] ]
    n4 = [ cgpind[0], cgpind[1] - 1, cgpind[2] ]
    n5 = [ cgpind[0] + 1, cgpind[1], cgpind[2] ]
    n6 = [ cgpind[0] - 1, cgpind[1], cgpind[2] ]

    return [n1, n2, n3, n4, n5, n6]
    

######################################
def get_points_in_Sph(origin,nx,ny,nz,gspace,centerind,radius):
    print " IN get_points_in_Sph"
    #cgpind = [0.0, 0.0, 0.0]
    # center is the center of the sphere, in this it is a grid point already
    #cgpind[0] = round((center[0]-origin[0])/nx) # closesed gridpoint index
    #cgpind[1] = round((center[1]-origin[1])/ny) # closesed gridpoint index
    #cgpind[2] = round((center[2]-origin[2])/nz) # closesed gridpoint index

    print origin,nx,ny,nz,gspace,centerind,radius

    #center = [0.0,0.0,0.0]
    #center[0] = centerind[0] * gspace + origin[0]
    #center[1] = centerind[1] * gspace + origin[1]
    #center[2] = centerind[2] * gspace + origin[2]

    #taged  = # taged # this is to be used 
    gpoints = [] #
    gpoints.append(centerind) 
    nieghbors = []
    nieghbors = nieghbors + get_neighbors(centerind) # combine list
    size = len(nieghbors)
    print "size =", size
    # while we still have neighbors to get.
    point = [0.0, 0.0, 0.0] 
    
    while (size > 0 ): 
        niegh_new = []
        for index in nieghbors: 
            if index in  gpoints: 
               continue
            # if point in taged
            #    continue
            #point[0] = index[0] * gspace + origin[0]# transform to coordenates 
            #point[1] = index[1] * gspace + origin[1] 
            #point[2] = index[2] * gspace + origin[2] 
            #if InSPH(center,radius,point):
            if InSPH(centerind,(radius/gspace),index):
               gpoints.append(index)
               #if index in niegh_new:
               #   print "index in niegh_new"
               #   continue
               nblist = get_neighbors(index) # here we get the neighbors
               #print len(nblist)
               for nb in nblist: 
                   flag = True
                   for val in nb: # donot append if negetibe value.
                       if val < 0.0: 
                          flag = False
                   if not flag:
                      print "**", nb
                      continue
                   if nb in nieghbors:
                      continue
                   if nb in niegh_new:
                      continue
                   if nb in gpoints:
                      continue
                   niegh_new.append(nb)
        nieghbors = niegh_new
        size = len(nieghbors)
        print "size =", size
        print "size of points = ", len(gpoints)
        
        #for nb in niegh_new:
        #for nb in gpoints:
        #    print nb
        if len(gpoints)>400:
           print "len(gpoints)>400"
           exit()

    return gpoints

##########################################

def write_pdb_threshold(file,xn,yn,zn,values,origin,dx,dy,dz):

    fileh = open(file,'w')
    #   inialize grid
    #grid = []
    #for i in range(xn):
    #    matrix = []
    #    for j in range(yn):
    #        row = []
    #        for k in range(zn):
    #            row.append(0)
    #        matrix.append(row)
    #    grid.append(matrix)

    index = 0
    for i in range(xn):
        for j in range(yn):
            for k in range(zn):
                #grid[i][j][k] = values[index]

                if (math.fabs(values[index]) > 0.5): # print pdb formate it greater than a threshold, for debuging
                    #convert grid index into xyz positions:
                    ic = i*dx + origin[0]
                    jc = j*dy + origin[1]
                    kc = k*dz + origin[2]
                    fileh.write("ATOM      1   O  WAT A   1     %7.3f %7.3f %7.3f %5.2f  0.00           O \n" % (ic,jc,kc,values[index]))
                index = index + 1
    fileh.close()

##########################################

def main():

   if len(sys.argv) != 4: # if no input
       print "ERORR:"
       print "syntex: dx-gist_precalculate_sphere.py infile sphradius outfileprefix"
       return
 

   infile1 = sys.argv[1]
   sphradius = float(sys.argv[2])
   outfile = sys.argv[3]

   print infile1
   print sphradius
   print outfile

   #exit()

   xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1 = read_in_dx_file(infile1)
   write_pdb_threshold(outfile+'old.pdb',xn1,yn1,zn1,values1,origin1,dx1,dy1,dz1)

   new_xn,new_yn,new_zn,new_origin,nvalues = pre_compute(outfile+'.dat',xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,sphradius)

   write_out_dx_file(outfile+'.dx',new_xn,new_yn,new_zn,dx1,dy1,dz1,new_origin,nvalues)

   write_pdb_threshold(outfile+'new.pdb',new_xn,new_yn,new_zn,nvalues,new_origin,dx1,dy1,dz1)

main()


