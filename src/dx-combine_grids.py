
import sys, os
import dx_gist_lib as dxlib 

# this is written by Trent Balius in the Shoichet Lab
# written in 2014
#
# here we read in, manpulate and output dx file.
# the dx file are produed by GIST calculations.


def main():

   if len(sys.argv) != 7: # if no input
       print "ERORR:"
       print "syntex: dx-combine_grids.py gist-Eww-dens.dx 1 gist-gO.dx 1 0.0 temp"
       print "v[i] = w1*v1[i] + w2*v2[i] + c3"
       return
 

   infile1 = sys.argv[1]
   weight1 = float(sys.argv[2])
   infile2 = sys.argv[3]
   weight2 = float(sys.argv[4])
   const3 = float(sys.argv[5])
   outfile = sys.argv[6]

   xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1 = dxlib.read_in_dx_file(infile1)
   xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2 = dxlib.read_in_dx_file(infile2)

   new_values = dxlib.combine_values(outfile+'.dat',xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,weight1,xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2,weight2,const3)

   dxlib.write_out_dx_file(outfile+'.dx',xn1,yn1,zn1,dx1,dy1,dz1,origin1,new_values)

main()


