
import sys, os, math

# this is written by Trent Balius in the Shoichet Lab
# written in june, 2015
#
# here we read in dx grid of some values (energies, denisties) and write out the center of intenities.
# frist apply a cut of.
# then clusters are determined of conneted grid points 
# (a point is in a cluster if that point is next to at least one other point in that cluster)
# the weighted center of each cluster is then calculated. 

class Point:
    def __init__(self,X,Y,Z,I,J,K,V,C):
        self.X = float(X) # coordenates
        self.Y = float(Y)
        self.Z = float(Z)
        self.I = int(I) # grid index
        self.J = int(J)
        self.K = int(K)
        self.V = float(V) # value
        self.C = int(C) # cluster



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
# this function takes as input a list of points and
# returns a dictionary of clusters.
# the lookup value is the point number, the definition is the cluster number
######################################################################
def cluster_points(points,dist):

    bonds = []
    bonded_p = {} # rember which are bonded. 
    N = len(points)
    # loop over list of points look at which are bonded. 
    # if 
    for i in range(N):
        ind_i1 = points[i].I 
        ind_j1 = points[i].J 
        ind_k1 = points[i].K 
        for j in range(i+1,N):
            ind_i2 = points[j].I 
            ind_j2 = points[j].J 
            ind_k2 = points[j].K 
            if (math.fabs(ind_i2-ind_i1) + math.fabs(ind_j2-ind_j1) + math.fabs(ind_k2-ind_k1)) == 0:
                print "Error. same point. duplicate found."
                exit()
            if (math.fabs(ind_i2-ind_i1) + math.fabs(ind_j2-ind_j1) + math.fabs(ind_k2-ind_k1)) <= dist: # if 0, they are the same, if grater than dist (one or two, must be an integer) they are not connected.  
               bond = [i, j]
               bonds.append(bond)
               bonded_p[i] = 1
               bonded_p[j] = 1
               print bond
    
    clusters = {}
    count = 0
    # 
    for bond in bonds:
        if bond[0] in clusters: # if the starting point is in a cluster then put the ending point in that same cluster
           clusters[bond[1]] = clusters[bond[0]]
        elif bond[1] in clusters: # if the ending point is in a cluster then put the starting point in that same cluster
           clusters[bond[0]] = clusters[bond[1]]
        else: # nether point is in a cluster
           clusters[bond[0]] = count 
           clusters[bond[1]] = count 
           count = count +1
    for i in range(N):
        if not i in bonded_p:
           # then singlton.
           clusters[i] = count
           count = count+1


    print clusters           
    for key in clusters.keys():
        print key, clusters[key]

    return clusters, count
        
######################################################################
######################################################################

def write_pdb(points,fileh):

    count = 0
    for point in points:
        fileh.write("ATOM  %5d  %-3s %3s %1s %3d %11.3f %7.3f %7.3f %5.2f\n" % (count, "CA", "CLU", "A", point.C, point.X, point.Y, point.Z, point.V/100.0 ))
        count = count + 1
       
######################################################################
# this function makes the cluster centers. 
# first find all all connected grid point after applying a threshold
# then calculate the weighted mean. 
# then print the  
######################################################################

def make_centers(fileprefix,values,gridscale,xn,yn,zn,origin,sign,threshold,dist):
# file - log file
# values gist values 
# gridscale - the grid spacing
# number -- number of grid spacings to combine
# sign -- direction to cluster


    fileh = open(fileprefix+'_cluster_centers.pdb','w')
    fileh2 = open(fileprefix+'_clusters.pdb','w')

    new_values = []
    countx = 0
    county = 0
    countz = 0
    #yflag = False
    #zflag = False
    #gridscale_new = gridscale*number

    ## here is how to calculate the new orgin: 
    ## x_new = x + gridscale / number * (1 + 2 + ... + number) 
    #sum = 0
    #for i in range (1,number): 
    #     sum = sum + i 
    #print "num=",number, 'sum=',sum
    #
    #ox_n =  origin[0] + gridscale / number * sum
    #oy_n =  origin[1] + gridscale / number * sum
    #oz_n =  origin[2] + gridscale / number * sum
    #origin_new = [ox_n,oy_n,oz_n]
#
    #print origin_new

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

    points = []

    # here we idenify grid points to be clustered
    for i in range(xn):
        x = (i * gridscale) + origin[0] 
        for j in range(yn):
            y = (j * gridscale) + origin[1] 
            for k in range(zn):
                z = (k * gridscale) + origin[2] 
                if (sign*grid_old[i][j][k] > threshold): # we only use the sign in the comparison with Threshold. 
                    print x,y,z, grid_old[i][j][k]
                    point = Point(x,y,z,i,j,k,grid_old[i][j][k],-1)
                    points.append(point)

    # here we cluster the points.  
    clusters, count = cluster_points(points,dist) 

    # put which cluster the points belong to in date structure.  
    for key in clusters.keys():
         i = key
         c = clusters[key]
         points[i].C = c
    write_pdb(points,fileh2)

    print "there are "+ str(count)+" clusters"

    # initalize a weight mean cluster point
    weight_mean_points = []
    
    num_cluster_points = []
    for c in range(count): # loop over the number of clusters
          weight_mean_points.append(Point(0.0,0.0,0.0,0,0,0,0.0,0))
          num_cluster_points.append(0)

    # calculate weighted mean of the cluters.
    for key in clusters.keys():
         i = key
         c = clusters[key]
         #print key, clusters[key]
         #weight_mean_points[c].X = weight_mean_points[c].X + points[i].X
         #weight_mean_points[c].Y = weight_mean_points[c].Y + points[i].Y
         #weight_mean_points[c].Z = weight_mean_points[c].Z + points[i].Z
         weight_mean_points[c].X = weight_mean_points[c].X + points[i].V * points[i].X
         weight_mean_points[c].Y = weight_mean_points[c].Y + points[i].V * points[i].Y
         weight_mean_points[c].Z = weight_mean_points[c].Z + points[i].V * points[i].Z
         weight_mean_points[c].V = weight_mean_points[c].V + points[i].V
         weight_mean_points[c].C = c
         num_cluster_points[c] = num_cluster_points[c] + 1  
    for c in range(count): # loop over number of clusters
         #weight_mean_points[c].X = weight_mean_points[c].X / num_cluster_points[c]
         #weight_mean_points[c].Y = weight_mean_points[c].Y / num_cluster_points[c]
         #weight_mean_points[c].Z = weight_mean_points[c].Z / num_cluster_points[c]
         weight_mean_points[c].X = weight_mean_points[c].X / weight_mean_points[c].V 
         weight_mean_points[c].Y = weight_mean_points[c].Y / weight_mean_points[c].V 
         weight_mean_points[c].Z = weight_mean_points[c].Z / weight_mean_points[c].V 
         print c, weight_mean_points[c].X, weight_mean_points[c].Y, weight_mean_points[c].Z, weight_mean_points[c].V


    write_pdb(weight_mean_points,fileh)

    fileh.close()
    fileh2.close()
    return 





def main():

   if len(sys.argv) != 6: # if no input
       print "ERORR:"
       print "syntex: dx-gist_make_centers.py dx-input-file sign cutoff distint pdb-output-file"
       print "sign can be 1 or -1: if one should always be used for densities,"
       print "                     for energies, weight of one means pdb contains cordenates where "
       print "                     waters do not want to be(anti-waters)."
       print "                     for energies, weight of negative one means pdb contains cordenates "
       print "                     where waters do want to be (negative energies)."
       print "cutoff can only be a positive floating point number. "
       print "                     this point is used to deterint clusters of points" 
       print "distint howmeny grid points away considered conected. should be 1 (0.5 angstroms) or 2(1.0 A), 3 (1.5 A)"
       print "dx-input-file input file in dx formate produed by gist, may be disities or energies"
       print "pdb-output-file is an output file in pdb formate contains the weighted mean poition of clusters"
       return
 

   infile1    = sys.argv[1]
   sign       = float(sys.argv[2])
   threshold  = float(sys.argv[3])
   dist       = int(sys.argv[4])
   outfile    = sys.argv[5]

   xn,yn,zn,dx,dy,dz,origin,values = read_in_dx_file(infile1)

   gridscale = dx # assumes that they are all the same spaceing

   make_centers(outfile,values,gridscale,xn,yn,zn,origin, sign, threshold, dist)


main()


