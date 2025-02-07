
import sys, os
import dx_gist_lib as dxlib

# this is written by Trent Balius in the Shoichet Lab
# written in 2014
#
# here we read in, manpulate and output dx file.
# the dx file are produed by GIST calculations.


def main():

   if len(sys.argv) != 5: # if no input
       print ("ERORR:")
       print ("syntex: dx-combine_grids.py gist-Eww-dens.dx gist-gO.dx 10.0 temp")
       print ("v[i] = {v1[i], if v2[i]>t or 0, if v2[i]<t")
       print ("v1 is energey, v2 is the density, t is the threshold that must be exceded.")
       return
 

   infile1   = sys.argv[1] # gist energy
   infile2   = sys.argv[2] # density
   threshold = float(sys.argv[3]) # density threshold
   outfile   = sys.argv[4]

   xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1 = dxlib.read_in_dx_file(infile1)
   xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2 = dxlib.read_in_dx_file(infile2)

   new_values = dxlib.combine_values_density_threshold(outfile+'.dat',xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2,threshold)

   dxlib.write_out_dx_file(outfile+'.dx',xn1,yn1,zn1,dx1,dy1,dz1,origin1,new_values)

main()


