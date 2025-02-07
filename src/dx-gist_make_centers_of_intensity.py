
import sys, os, math
import dx_gist_lib as dxlib

# this is written by Trent Balius in the Shoichet Lab
# written in june, 2015
#
# changed names, MF -- March 2017
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
                print ("Error. same point. duplicate found.")
                exit()
            if (math.fabs(ind_i2-ind_i1) + math.fabs(ind_j2-ind_j1) + math.fabs(ind_k2-ind_k1)) <= dist: # if 0, they are the same, if grater than dist (one or two, must be an integer) they are not connected.  
               bond = [i, j]
               bonds.append(bond)
               bonded_p[i] = 1
               bonded_p[j] = 1
               print (bond)
    
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


    print (clusters)           
    for key in clusters.keys():
        print (key, clusters[key])

    return clusters, count
        
######################################################################
######################################################################

def write_pdb(points,fileh):

    count = 0
    for point in points:
        #fileh.write("ATOM  %5d  %-3s %3s %1s %3d %11.3f %7.3f %7.3f %5.2f\n" % (count, "CA", "CLU", "A", point.C, point.X, point.Y, point.Z, point.V/100.0 ))
        fileh.write("ATOM  %5d  %-3s %3s %1s%4d %11.3f %7.3f %7.3f %5.2f\n" % (count, "CA", "CLU", "A", point.C, point.X, point.Y, point.Z, point.V/100.0 ))
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


    fileh = open(fileprefix+'_clusters.pdb','w')
    fileh2 = open(fileprefix+'_points.pdb','w')

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
                    print (x,y,z, grid_old[i][j][k])
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

    print ("there are "+ str(count)+" clusters")

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
         print (c, weight_mean_points[c].X, weight_mean_points[c].Y, weight_mean_points[c].Z, weight_mean_points[c].V)


    write_pdb(weight_mean_points,fileh)

    fileh.close()
    fileh2.close()
    return 





def main():

   if len(sys.argv) != 6: # if no input
       print ("ERORR:"                                                                                                    )
       print ("syntex: dx-gist_make_centers.py dx-input-file sign cutoff distint pdb-output-file"                         )
       print ("sign can be 1 or -1: if one should always be used for densities,"                                          )
       print ("                     for energies, weight of one means pdb contains cordenates where "                     )
       print ("                     waters do not want to be(anti-waters)."                                               )
       print ("                     for energies, weight of negative one means pdb contains cordenates "                  )
       print ("                     where waters do want to be (negative energies)."                                      )
       print ("cutoff can only be a positive floating point number. "                                                     )
       print ("                     this point is used to deterint clusters of points"                                    )
       print ("distint howmeny grid points away considered conected. should be 1 (0.5 angstroms) or 2(1.0 A), 3 (1.5 A)"  )
       print ("dx-input-file input file in dx formate produed by gist, may be disities or energies"                       )
       print ("pdb-output-file is an output file in pdb formate contains the weighted mean poition of clusters"           )
       return
 

   infile1    = sys.argv[1]
   sign       = float(sys.argv[2])
   threshold  = float(sys.argv[3])
   dist       = int(sys.argv[4])
   outfile    = sys.argv[5]

   xn,yn,zn,dx,dy,dz,origin,values = dxlib.read_in_dx_file(infile1)

   gridscale = dx # assumes that they are all the same spaceing

   make_centers(outfile,values,gridscale,xn,yn,zn,origin, sign, threshold, dist)


main()


