
import sys, os, math
import dx_gist_lib as dxlib

# this is written by Trent Balius in the Shoichet Lab
# written in Oct, 2014
# modifed by Reed Stein, 2016
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

######################################################################
# this function compute the atomic displacement grids.  
######################################################################

def pre_compute(file,xn,yn,zn,dx,dy,dz,origin,values,radius,divider,pad_radius):
# file - log file
# values gist values 
# gridscale - the grid spacing
# radius - sphere radius

    fileh = open(file,'w')

    padx = math.ceil(pad_radius/dx)
    pady = math.ceil(pad_radius/dy)
    padz = math.ceil(pad_radius/dz)

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
                points,dist = get_points_in_Sph(origin,xn,yn,zn,gspace,centerind,radius) 
                for indp,p in enumerate(points):
                    scale = cal_gausian(dist[indp],(radius/divider)) 
                    #print "I AM HERE ",indp,dist[indp],scale
                    ngrid[i][j][k] = ngrid[i][j][k]  + scale*grid[p[0]][p[1]][p[2]]
		#    if scale*grid[p[0]][p[1]][p[2]] > 0.0:
		#	print "indices are", (p[0], p[1], p[2])
		#	print "scale is",scale
		#	print "distance is", dist[indp]
		#	print "new_val is", scale*grid[p[0]][p[1]][p[2]]
		#	print "orig_val is", grid[p[0]][p[1]][p[2]]

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
def calc_distance(center,point):
    dist2 = 0.0
    for i in range(3):
        dist2 = dist2 + (center[i] - point[i])**2
    return math.sqrt(dist2)

######################################
def cal_gausian(dist,sigma):
    #print dist, sigma
    var = sigma**2.0
    val = 1.0/math.sqrt(2.0*math.pi*var) * math.exp(-1.0*(dist**2.0)/(2.0*var))
    return val
    
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
    #print " IN get_points_in_Sph"
    #cgpind = [0.0, 0.0, 0.0]
    # center is the center of the sphere, in this it is a grid point already
    #cgpind[0] = round((center[0]-origin[0])/nx) # closesed gridpoint index
    #cgpind[1] = round((center[1]-origin[1])/ny) # closesed gridpoint index
    #cgpind[2] = round((center[2]-origin[2])/nz) # closesed gridpoint index

    #print origin,nx,ny,nz,gspace,centerind,radius

    #center = [0.0,0.0,0.0]
    #center[0] = centerind[0] * gspace + origin[0]
    #center[1] = centerind[1] * gspace + origin[1]
    #center[2] = centerind[2] * gspace + origin[2]

    #taged  = # taged # this is to be used 
    gpoints = [] #
    distances = [0.0] #
    gpoints.append(centerind) 
    nieghbors = []
    nieghbors = nieghbors + get_neighbors(centerind) # combine list
    size = len(nieghbors)
    #print "size =", size
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
               dist = calc_distance(centerind,index)
               dist = gspace * dist  # converts from grid point to angstroms. 
               distances.append(dist)
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
                      #print "**", nb
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
        #print "size =", size
        #print "size of points = ", len(gpoints)
        
        if len(gpoints)>400:
           #print "len(gpoints)>400"
           exit()

    return gpoints, distances

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

   if len(sys.argv) != 6: # if no input
       print "ERORR:"
       print "syntex: dx-gist_precalculate_sphere.py infile sphradius divider pad_radius outfileprefix"
       return
 

   infile1 = sys.argv[1]
   sphradius = float(sys.argv[2])
   divider = float(sys.argv[3])
   pad_radius = float(sys.argv[4])
   outfile = sys.argv[5]

   if (pad_radius < sphradius): 
       print "warning pap_radius (%f) is less than sphradius (%f).\n" % (pad_radius,sphradius)
       print "This might cause boundary problems.\n"

   print infile1
   print sphradius
   print outfile

   #exit()

   xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1 = dxlib.read_in_dx_file(infile1)
   #write_pdb_threshold(outfile+'old.pdb',xn1,yn1,zn1,values1,origin1,dx1,dy1,dz1)

   new_xn,new_yn,new_zn,new_origin,nvalues = pre_compute(outfile+'.dat',xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,sphradius,divider,pad_radius)

   dxlib.write_out_dx_file(outfile+'.dx',new_xn,new_yn,new_zn,dx1,dy1,dz1,new_origin,nvalues)

   #write_pdb_threshold(outfile+'new.pdb',new_xn,new_yn,new_zn,nvalues,new_origin,dx1,dy1,dz1)

main()


