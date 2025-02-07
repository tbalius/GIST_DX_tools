
import sys, os, math
import dx_gist_lib as dxlib

# this is written by Trent Balius in the Shoichet Lab
# written in Oct, 2014
# modifed by Reed Stein, 2016
# modified in Jan, 2025
#
# here we read in sites with a value, manpulate and output dx file.
# the sites file are produed by  GIST calculations.

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

# pre compute displacement from sites
def pre_compute_site(file1,sx,sy,sz,svalue,xn,yn,zn,dx,dy,dz,origin,radius,divider):
#   print " in pre_compute "
# file - log file
# values gist values 
# gridscale - the grid spacing
# radius - sphere radius

    fileh = open(file1,'w')

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

    
    if (len(sx)!=len(sy) or len(sx)!=len(sz) or len(sx)!=len(svalue)):
        print("Error. ")
        exit()

    sN = len(sx)

    corner = [0.0, 0.0, 0.0]

    #corner[0] = origin[0] - (xn * dx / 2.0)
    #corner[1] = origin[1] - (yn * dy / 2.0)
    #corner[2] = origin[2] - (zn * dz / 2.0)
    corner[0] = origin[0] 
    corner[1] = origin[1] 
    corner[2] = origin[2] 
    fileh.write('origin = %f, %f, %f\n'%(origin[0],origin[1],origin[2]))
    fileh.write('corner = %f, %f, %f\n'%(corner[0],corner[1],corner[2]))

    # compute grid from sites
    for si in range(sN):
      # find closed grid point
      fileh.write('site x,y,z coordenate = %f, %f, %f\n'%(sx[si],sy[si],sz[si]))
      gi = round((sx[si] - corner[0])/dx)
      gj = round((sy[si] - corner[1])/dy)
      gk = round((sz[si] - corner[2])/dz)
      fileh.write('site grid point position i,j,k = = %d, %d, %d\n'%(gi,gj,gk))
      centerind = [gi, gj, gk]
      centercord = [sx[si],sy[si],sz[si]]
      points,dist = get_points_in_Sph(origin,xn,yn,zn,gspace,centerind,centercord,radius) 
      #points,dist = get_points_in_Sph(origin,xn,yn,zn,gspace,centerind,2*radius) # get more points 
      for indp,p in enumerate(points):
         scale = cal_gausian(dist[indp],(radius/divider)) 
         print(p)
         #grid[i][j][k] = grid[i][j][k]  + scale*grid[p[0]][p[1]][p[2]]
         grid[p[0]][p[1]][p[2]] = grid[p[0]][p[1]][p[2]] + scale*svalue[si]

    # convert new grid to array
    nvalues = []

    index = 0
    for i in range(xn):
        for j in range(yn):
            for k in range(zn):
                  nvalues.append(grid[i][j][k])

    new_xn = xn
    new_yn = yn
    new_zn = zn
    
    new_origin = [ 0.0, 0.0, 0.0]
    new_origin[0] = origin[0] 
    new_origin[1] = origin[1]
    new_origin[2] = origin[2] 

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
def get_points_in_Sph(origin,nx,ny,nz,gspace,centerind,centercord,radius):
#    print " IN get_points_in_Sph"
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
            point[0] = index[0] * gspace + origin[0]# transform to coordenates 
            point[1] = index[1] * gspace + origin[1] 
            point[2] = index[2] * gspace + origin[2] 
            #if InSPH(center,radius,point):
            #if InSPH(centerind,(radius/gspace),index):
            if InSPH(centercord,radius,point):
               gpoints.append(index)
               #dist = calc_distance(centerind,index)
               #dist = gspace * dist  # converts from grid point to angstroms. 
               dist = calc_distance(centercord,point)
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
        print ("size =", size)
        print ("size of points = ", len(gpoints))
        print (gpoints[0])
        
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
def read_in_site_file(infile1,margin):

    fh = open(infile1,'r')
    x = [];y=[];z=[];v=[]
    gcenter1 = [0.0,0.0,0.0]
    maxcorner = [0.0,0.0,0.0]
    mincorner = [1000.0,1000.0,1000.0]
    for line in fh:
        sline = line.split()
        tx=float(sline[0])
        ty=float(sline[1])
        tz=float(sline[2])
        tv=float(sline[3])
        x.append(tx)
        y.append(ty)
        z.append(tz)
        v.append(tv)
        gcenter1[0] = gcenter1[0] + tx
        gcenter1[1] = gcenter1[1] + ty
        gcenter1[2] = gcenter1[2] + tz
        if (maxcorner[0] < tx):
            maxcorner[0] = tx
        if (maxcorner[1] < ty):
            maxcorner[1] = ty
        if (maxcorner[2] < tz):
            maxcorner[2] = tz

        if (mincorner[0] > tx):
            mincorner[0] = tx
        if (mincorner[1] > ty):
            mincorner[1] = ty
        if (mincorner[2] > tz):
            mincorner[2] = tz

    N = len(x)
    gcenter1[0] = gcenter1[0]/float(N) 
    gcenter1[1] = gcenter1[1]/float(N) 
    gcenter1[2] = gcenter1[2]/float(N) 
    print("center=",gcenter1[0],gcenter1[1],gcenter1[2])
       
 
    gdx1 = 0.5; gdy1 = 0.5; gdz1 = 0.5

    maxcorner[0] = maxcorner[0]+margin
    maxcorner[1] = maxcorner[1]+margin
    maxcorner[2] = maxcorner[2]+margin

    mincorner[0] = mincorner[0]-margin
    mincorner[1] = mincorner[1]-margin
    mincorner[2] = mincorner[2]-margin

    #gorigin1 = mincorner
    gorigin1 = [0.0, 0.0, 0.0]
    gorigin1[0] = round(mincorner[0],4)
    gorigin1[1] = round(mincorner[1],4)
    gorigin1[2] = round(mincorner[2],4)

    print("origin=",gorigin1[0],gorigin1[1],gorigin1[2])

    diffx = maxcorner[0] - mincorner[0]
    diffy = maxcorner[1] - mincorner[1]
    diffz = maxcorner[2] - mincorner[2]

    gxn1 = math.ceil(diffx/gdx1)
    gyn1 = math.ceil(diffy/gdy1)
    gzn1 = math.ceil(diffz/gdz1)
    print("xn,yn,zn=",gxn1,gyn1,gzn1)
 

    return x,y,z,v,gxn1,gyn1,gzn1,gdx1,gdy1,gdz1,gorigin1
##########################################
##########################################

def main():

   if len(sys.argv) != 6: # if no input
       print ("ERORR:")
       print ("syntex: site_to_dx-gist_precalculate_sphere.py infile margin sphradius divider outfileprefix")
       print ("divider is used to calculate sigma:  sigma = radius/divider ")
       return
 

   infile1 = sys.argv[1]
   margin = float(sys.argv[2])
   sphradius = float(sys.argv[3])
   divider = float(sys.argv[4])
   outfile = sys.argv[5]

   if ( margin < sphradius): 
       print ("warning margin (%f) is less than sphradius (%f).\n" % (margin,sphradius))
       print ("This might cause boundary problems.\n")

   print (infile1)
   print (sphradius)
   print (divider)
   print (outfile)

   #exit()

   x,y,z,value,xn1,yn1,zn1,dx1,dy1,dz1,origin1 = read_in_site_file(infile1,margin)
   #write_pdb_threshold(outfile+'old.pdb',xn1,yn1,zn1,values1,origin1,dx1,dy1,dz1)

   new_xn,new_yn,new_zn,new_origin,nvalues = pre_compute_site(outfile+'.dat',x,y,z,value,xn1,yn1,zn1,dx1,dy1,dz1,origin1,sphradius,divider)

   dxlib.write_out_dx_file(outfile+'.dx',new_xn,new_yn,new_zn,dx1,dy1,dz1,new_origin,nvalues)

   #write_pdb_threshold(outfile+'new.pdb',new_xn,new_yn,new_zn,nvalues,new_origin,dx1,dy1,dz1)

main()


