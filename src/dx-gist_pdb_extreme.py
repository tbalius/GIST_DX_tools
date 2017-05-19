
import sys, os, math
import dx_gist_lib as dxlib

# this is written by Trent Balius in the Shoichet Lab
# written in june, 2015
#
# changed names, MF -- March 2017
# based on centers_of_intensity  code May, 2017
######################################################################
######################################################################

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



def write_pdb(points,fileh):

    count = 0
    for point in points:
        fileh.write("ATOM  %5d  %-3s %3s %1s %3d %11.3f %7.3f %7.3f %5.2f\n" % (count, "CA", "CLU", "A", point.C, point.X, point.Y, point.Z, point.V/100.0 ))
        count = count + 1
       
######################################################################
# this function writes extreme points to file. 
######################################################################

def get_extreme(fileprefix,values,gridscale,xn,yn,zn,origin,threshold):
# file - log file
# values gist values 
# gridscale - the grid spacing
# number -- number of grid spacings to combine
# sign -- direction to cluster


    fileh = open(fileprefix+'_extreme_points.pdb','w')

    new_values = []
    countx = 0
    county = 0
    countz = 0

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

    maxval = -1000.0
    minval = 1000.0

    # here we idenify grid points to be clustered
    for i in range(xn):
        x = (i * gridscale) + origin[0] 
        for j in range(yn):
            y = (j * gridscale) + origin[1] 
            for k in range(zn):
                z = (k * gridscale) + origin[2] 
                if (grid_old[i][j][k] > threshold or grid_old[i][j][k] < (-1.0 * threshold)): # sign in the comparison with Threshold. 
                    print i,j,k, x,y,z, grid_old[i][j][k]
                    point = Point(x,y,z,i,j,k,grid_old[i][j][k],1)
                    points.append(point)
                if (maxval < grid_old[i][j][k]):
                    maxval =  grid_old[i][j][k]
                if (minval > grid_old[i][j][k]):
                    minval =  grid_old[i][j][k]
                
    print "max value: %f"%(maxval)
    print "min value: %f"%(minval)
    write_pdb(points,fileh)

    fileh.close()
    return 





def main():

   if len(sys.argv) != 4: # if no input
       print "ERORR:"
       print "syntex: dx-gist_extreme.py dx-input-file sign cutoff pdb-output-file"
       print "cutoff can only be a positive floating point number. "
       print "                     this point is used to deterint clusters of points" 
       print "dx-input-file input file in dx formate produed by gist, may be disities or energies"
       print "pdb-output-file is an output file in pdb formate contains the weighted mean poition of clusters"
       return
 

   infile1    = sys.argv[1]
   threshold  = float(sys.argv[2])
   outfile    = sys.argv[3]

   xn,yn,zn,dx,dy,dz,origin,values = dxlib.read_in_dx_file(infile1)

   gridscale = dx # assumes that they are all the same spaceing

   get_extreme(outfile,values,gridscale,xn,yn,zn,origin, threshold)


main()


