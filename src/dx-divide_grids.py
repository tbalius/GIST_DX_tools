
import sys, os
import dx_gist_lib as dxlib
# this is written by Trent Balius in the Shoichet Lab
# written in 2015
#
# here we read in, manpulate and output dx file.
# the dx file are produed by GIST calculations.



def main():

   if len(sys.argv) != 5: # if no input
       print "ERORR:"
       print "syntex: dx-divide_grids.py gist-Eww-dens.dx gist-gO.dx 0.0329 temp"
       print "v[i] = v1[i] / (w2*v2[i]) "
       return
 

   infile1 = sys.argv[1]
   infile2 = sys.argv[2]
   weight2 = float(sys.argv[3])
   outfile = sys.argv[4]

   xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1 = dx.read_in_dx_file(infile1)
   xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2 = dx.read_in_dx_file(infile2)

   new_values = dx.divide_energy_density(outfile+'.dat',xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2,weight2)

   dx.write_out_dx_file(outfile+'.dx',xn1,yn1,zn1,dx1,dy1,dz1,origin1,new_values)

main()


