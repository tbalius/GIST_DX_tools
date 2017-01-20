
import sys, os, math
import dx_gist_lib as dxlib

# this is written by Trent Balius in the Shoichet Lab
# written in 2014
#
# here we read in, manpulate and output dx file.
# the dx file are produed by GIST calculations.


# calculate the square root of the sum of the squares. 
def square_root_sum_of_squares(values):
    sum2 = 0 
    for ele in values:
        sum2 = sum2 + ele**2
    return math.sqrt(sum2)


def main():

   if len(sys.argv) != 4: # if no input
       print "ERORR:"
       print "this script caluculates the difference grid and the second norm over the grid. "
       print "syntex: dx-compare_2norm.py gist-Eww-dens.dx gist-gO.dx diff_grid"
       return
 

   infile1 = sys.argv[1]
   infile2 = sys.argv[2]
   outfile = sys.argv[3]
   weight1 = 1
   weight2 = -1
   const3  = 0 

   xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1 = dxlib.read_in_dx_file(infile1)
   xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2 = dxlib.read_in_dx_file(infile2)

   # caluclate the differents grid
   new_values = dxlib.combine_values(outfile+'.dat',xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,weight1,xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2,weight2,const3)
   dxlib.write_out_dx_file(outfile+'.dx',xn1,yn1,zn1,dx1,dy1,dz1,origin1,new_values)
   norm = square_root_sum_of_squares(new_values)
   print "norm = ", norm

main()


