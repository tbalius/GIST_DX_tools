
import sys, os
import dx_gist_lib as dxlib

# this is written by Trent Balius in the Shoichet Lab
# written in 2014
#
# here we read in, manpulate and output dx file.
# the dx file are produed by GIST calculations.

def main():

   if len(sys.argv) != 4: # if no input
       print "ERORR:"
       print "syntex: dx-remove_extrema.py gistvalues.dx 5.0 gistvalues_cap"
       print "vnew[i] = {-t1, if v[i] <= -t1 "
       print "            t1, if v[i] >=  t1  "
       print "            v[i] otherwise      "
       return
 

   infile1   = sys.argv[1] # gist energy
   t1        = float(sys.argv[2]) # energy threshold
   outfile   = sys.argv[3]

   xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1 = dxlib.read_in_dx_file(infile1)

   new_values = []
   countn = 0
   countp = 0
   for val in values1:
       if (val < -1.0*t1):
           val = -1.0*t1
           countn = countn + 1
       elif (val > 1.0*t1):
           val =  1.0*t1
           countp = countp + 1
       #else:
       new_values.append(val)

   print "%d number of point are caped in positive direction" % countp
   print "%d number of point are caped in negative direction" % countn

   #new_values = combine_values_density_threshold(outfile+'.dat',xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2,threshold)
   dxlib.write_out_dx_file(outfile+'.dx',xn1,yn1,zn1,dx1,dy1,dz1,origin1,new_values)

main()


