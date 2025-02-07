
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

    # calculate mean, standared deviation, skew.  

    sum1 = 0
    sum2 = 0 # sum2 = sum over i of x_i^2
    sum3 = 0 # sum3 = sum over i of x_i^3
    sum1no0 = 0
    sum2no0 = 0 # sum2 = sum over i of x_i^2
    sum3no0 = 0 # sum3 = sum over i of x_i^3
    countno0     = 0
    lenv = len(values)
    for i in range(lenv):
        sum1 = sum1 + values[i]
        sum2 = sum2 + values[i]**2
        sum3 = sum3 + values[i]**3
        if values[i] != 0.0:
           sum1no0 = sum1no0 + values[i]
           sum2no0 = sum2no0 + values[i]**2
           sum3no0 = sum3no0 + values[i]**3
           countno0 = countno0 + 1
    mean = sum1/lenv
    std  = (sum2/lenv - mean**2.0)**(1.0/2.0)
    print ("mean1=%f,mean2=%f,mean3=%f"%(sum1/lenv,sum2/lenv,sum3/lenv))
    print ("mean=%f,standdev=%f"%(mean,std))

    meanno0 = sum1no0/countno0
    stdno0  = (sum2no0/countno0 - meanno0**2.0)**(1.0/2.0)
    print ("non-zero points = %d"%countno0)
    print ("no0:mean1=%f,mean2=%f,mean3=%f"%(sum1no0/countno0,sum2no0/countno0,sum3no0/countno0))
    print ("no0:mean=%f,standdev=%f"%(meanno0,stdno0))

    val    = 5.0
    countp = 0 # above mu + val* sigma
    countm = 0 # below mu - val* sigma
    countb = 0 # between mu +- val* sigma
    for i in range(lenv):
        if ((mean+val*std) < values[i]):
           countp = countp+1
        elif ((mean-val*std) > values[i]):
           countm = countm+1
        else:
           countb = countb+1
    print ("mu-%dsigma = %f"%(val,(mean-val*std))       )
    print ("mu+%dsigma = %f"%(val,(mean+val*std))       )
    print ("below mu-%dsigma = %d"%(int(val),countm)    )
    print ("above mu+%dsigma = %d"%(int(val),countp)    )
    print ("between mu+-%dsigma = %d"%(int(val),countb) )
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
                    print (i,j,k, x,y,z, grid_old[i][j][k])
                    point = Point(x,y,z,i,j,k,grid_old[i][j][k],1)
                    points.append(point)
                if (maxval < grid_old[i][j][k]):
                    maxval =  grid_old[i][j][k]
                if (minval > grid_old[i][j][k]):
                    minval =  grid_old[i][j][k]
                
    print ("max value: %f"%(maxval))
    print ("min value: %f"%(minval))
    write_pdb(points,fileh)

    fileh.close()
    return 





def main():

   if len(sys.argv) != 4: # if no input
       print ("ERORR:"                                                                                          )
       print ("syntex: dx-gist_extreme.py dx-input-file cutoff pdb-output-file"                                 )
       print ("cutoff can only be a positive floating point number. "                                           )
       print ("                     this point is used to deterint clusters of points"                          )
       print ("dx-input-file input file in dx formate produed by gist, may be disities or energies"             )
       print ("pdb-output-file is an output file in pdb formate contains the weighted mean poition of clusters" )
       return
 

   infile1    = sys.argv[1]
   threshold  = float(sys.argv[2])
   outfile    = sys.argv[3]

   xn,yn,zn,dx,dy,dz,origin,values = dxlib.read_in_dx_file(infile1)

   gridscale = dx # assumes that they are all the same spaceing

   get_extreme(outfile,values,gridscale,xn,yn,zn,origin, threshold)


main()


