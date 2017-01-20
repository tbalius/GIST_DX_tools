
import sys, os
import dx_gist_lib as dxlib

# this is written by Trent Balius in the Shoichet Lab
# written in 2014
#
# here we read in, manpulate and output dx file.
# the dx file are produed by GIST calculations.

def main():

   if len(sys.argv) != 6: # if no input
       print "ERORR:"
       print "syntex: dx-combine_grids.py gist-Eww-dens.dx gist-gO.dx 1.0 10.0 temp"
       print "v[i] = {t1, if v1[i] >= t1 v2[i]>=t2 or 0, if v2[i]<t2"
       print "v1 is energey, v2 is the density, t1 is the energy cutoff t2 is the density threshold that must be exceded."
       return
 

   infile1   = sys.argv[1] # gist energy
   infile2   = sys.argv[2] # density
   t1        = float(sys.argv[3]) # energy threshold
   t2        = float(sys.argv[4]) # density threshold
   outfile   = sys.argv[5]

   xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1 = dxlib.read_in_dx_file(infile1)
   xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2 = dxlib.read_in_dx_file(infile2)

   #new_values = combine_values_density_threshold(outfile+'.dat',xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2,threshold)
   new_values = dxlib.combine_values_energy_density_thresholds(outfile+'.dat',xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,t1,xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2,t2)
   dxlib.write_out_dx_file(outfile+'.dx',xn1,yn1,zn1,dx1,dy1,dz1,origin1,new_values)

main()


