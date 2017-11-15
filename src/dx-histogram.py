
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
       print "syntex: dx-histogram.py gistvalues.dx min max num output [vyes,vno]"
       return
 

   infile1   = sys.argv[1] # gist energy
   minv      = float(sys.argv[2]) # energy threshold
   maxv      = float(sys.argv[3])
   num       = int(sys.argv[4])
   output    = sys.argv[5]
   verbose   = sys.argv[6]
   #threshold = sys.argv[7]
   threshold  = 3.0

   if verbose == "vyes":
       boolverbose = True
   elif verbose == "vno":
       boolverbose = False
   else:
       print "error: must be vyes,vno"

   xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1 = dxlib.read_in_dx_file(infile1)


   histogram = []
   histogram_t = [] # trunc histogram
   histogram_nz = [] # histogram without zero values counted
   midbins    = []

   step = (maxv - minv)/float(num)
   print "infile = %s"%infile1
   print "minv = %f\nmaxv = %f\n num = %d\n step = %f\n"%(minv,maxv,num,step)
   midbin = minv + step/2.0
   for i in range(0,num):
       histogram.append(0)
       histogram_t.append(0)
       histogram_nz.append(0)
       if midbin > maxv: 
          print "error"
          exit()
       midbins.append(midbin)
       midbin = midbin + step
       
   countn = 0
   countp = 0
   min_ov = 1000 # minimum observed value
   max_ov = -1000 # maximum observed value

   ip3 = int( round( (+3.0 - minv) / (maxv-minv) * (num-1) )) # value of index for value of +3.0
   in3 = int( round( (-3.0 - minv) / (maxv-minv) * (num-1) ))# value of index for value of -3.0

   for val in values1:
       if (val < minv):
           countn = countn + 1
           continue
       if (val > maxv):
           countp = countp + 1
           continue
       i = int( round( (val - minv) / (maxv-minv) * (num-1) )) 
       if boolverbose:
          print "%f->%d->%f"%(val, i, midbins[i])
       histogram[i] = histogram[i] + 1
       if val != 0.0:
          histogram_nz[i]=histogram_nz[i]+1
       if val <= threshold and val >= -1.0*threshold:
          histogram_t[i]=histogram_t[i]+1
       elif (val > threshold):
          histogram_t[ip3] = histogram_t[ip3]+1
       elif (val < -1.0*threshold):
          histogram_t[in3] = histogram_t[in3]+1
   fh = open(output,'w')
   for i in range(0,num):
        print "%f, %d,%d,%d"%(midbins[i], histogram[i],histogram_nz[i],histogram_t[i])
        fh.write( "%f,%d,%d,%d\n"%(midbins[i], histogram[i],histogram_nz[i],histogram_t[i]))
   print "number less than min: %d\nnumber grater than max: %d"%(countn,countp)
main()


