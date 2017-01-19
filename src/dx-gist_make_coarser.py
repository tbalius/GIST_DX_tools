
import sys, os, math

# this is written by Trent Balius in the Shoichet Lab
# written in Oct, 2014
#
# here we read in, manpulate and output dx file.
# the dx file are produed by GIST calculations.

# here we want to make the grid courser. 



def read_in_dx_file(file):

  fileh = open(file,'r')

  flag_read_dx = False

  count = 0

  values = []
  for line in fileh:
      splitline = line.split()

      #print splitline
      #if len(splitline) < 2:
      if len(splitline) == 0:
          #print line
          continue

      ## this should be line 1
      if (splitline[0] == "object" and splitline[1] == "1"):
         print "count = ", count, " line = ", line
         xn = int(splitline[5])
         yn = int(splitline[6])
         zn = int(splitline[7])
         
      ## this should be line 2       
      if (splitline[0] == "origin"):
         #print line
         print "count = ", count, " line = ", line
         origin = [float(splitline[1]), float(splitline[2]), float(splitline[3])] 

      ## this should be lines 3-5
      if (splitline[0] == "delta"):
         #print line
         print "count = ", count, " line = ", line
         if (float(splitline[2]) == 0 and  float(splitline[3]) ==0):
            dx = float(splitline[1]) 
         elif (float(splitline[1]) == 0 and  float(splitline[3]) ==0):
            dy = float(splitline[2]) 
         elif (float(splitline[1])== 0 and  float(splitline[2])==0):
            dz = float(splitline[3]) 
            print dx, dy, dz 


      if (splitline[0] == "object" and splitline[1] == "2"):
         #print line
         print "count = ", count, " line = ", line
      if (splitline[0] == "object" and splitline[1] == "3"):
         #print line
         print "count = ", count, " line = ", line
         flag_read_dx = True
         continue # go to next line
      if (flag_read_dx):

         if (len(splitline) > 3): 
            print "Error: dx formate problem. more than 3 colums"
            exit()

         for value in splitline:
             values.append(float(value))  

      count = count + 1


  print len(values)
  fileh.close()
  return xn,yn,zn,dx,dy,dz,origin,values 

def write_out_dx_file(file,xn,yn,zn,dx,dy,dz,origin,values):

   
  fileh = open(file,'w')
#object 1 class gridpositions counts 40 40 40
#origin 35.31 27.576 18.265
#delta 0.5 0 0
#delta 0 0.5 0
#delta 0 0 0.5
#object 2 class gridconnections counts 40 40 40
#object 3 class array type float rank 0 items 64000 data follows

  fileh.write('object 1 class gridpositions counts %d %d %d\n' % (xn,yn,zn))
  fileh.write('origin %6.2f %6.2f %6.2f\n' % (origin[0],origin[1],origin[2]))
  fileh.write('delta %2.1f 0 0\n' % dx)
  fileh.write('delta 0 %2.1f 0\n' % dy)
  fileh.write('delta 0 0 %2.1f\n' % dz)
  fileh.write('object 2 class gridconnections counts %d %d %d\n' % (xn,yn,zn))
  fileh.write('object 3 class array type float rank 0 items %d data follows\n' % len(values))

  count = 1
  for value in values:
       if (value == 0.0): 
          fileh.write('%d' % 0)
       else:
          fileh.write('%f' % value)
       # print newline after 3rd number.
       if (count == 3): 
            fileh.write('\n')
            count = 0
       # print space after number but not at the end of the line.
       else:
            fileh.write(' ')
       count = count + 1

  # if the last line has less than 3 numbers then print the a newline.
  if (count < 3):
       fileh.write('\n')
  fileh.close()
       
######################################################################
# this function compute the atomic displacement grids.  
######################################################################

def make_courser(file,values,gridscale,xn,yn,zn,origin,number):
# file - log file
# values gist values 
# gridscale - the grid spacing
# number -- number of grid spacings to combine


    fileh = open(file,'w')

    new_values = []
    countx = 0
    county = 0
    countz = 0
    #yflag = False
    #zflag = False
    gridscale_new = gridscale*number

    ## here is how to calculate the new orgin: 
    ## x_new = x + gridscale / number * (1 + 2 + ... + number) 
    sum = 0
    for i in range (1,number): 
         sum = sum + i 
    print "num=",number, 'sum=',sum

    ox_n =  origin[0] + gridscale / number * sum
    oy_n =  origin[1] + gridscale / number * sum
    oz_n =  origin[2] + gridscale / number * sum
    origin_new = [ox_n,oy_n,oz_n]

    print origin_new

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
    
    # here we initialize the new grid. 
    grid = [] 
    xnew = int(math.floor(xn / number))
    ynew = int(math.floor(yn / number))
    znew = int(math.floor(zn / number))
    print xnew,ynew,znew
    for i in range(xnew):
        ydim = []
        for j in range(ynew):
            zdem = []
            for k in range(znew):
                zdem.append(0.0)
                count = count + 1
            ydim.append(zdem)
        grid.append(ydim) 

    # we fill the courser grid with the sumed values.
    for i in range(xn):
        i_n = int(math.ceil(i/number))
        if (i_n >= xnew): 
           print "continue"
           continue
        for j in range(yn):
            j_n = int(math.ceil(j/number))
            if (j_n >= ynew): 
               print "continue"
               continue
            for k in range(zn):
               k_n = int(math.ceil(k/number))
               if (k_n >= znew): 
                   print "continue"
                   continue
               grid[i_n][j_n][k_n] = grid[i_n][j_n][k_n] + grid_old[i][j][k]
               #fileh.write("%f = %f * %f +  %f * %f + %f\n" % (value,w1,values1[i],w2,values2[i],c3))

    values_new = []
    for i in range(xnew):
        for j in range(ynew):
            for k in range(znew):
               values_new.append(grid[i][j][k])

    fileh.close()
    return xnew,ynew,znew,origin_new,gridscale_new,values_new


######################################################################
# this function reads in energy values and densities as vectors.
# output vn[i] = v[i], if d[i] > cutoff
#     or vn[i] = 0 otherwise
######################################################################
#def apply_density_threshold_values(xn1,yn1,zn1,dx1,dy1,dz1,origin1,values,xn2,yn2,zn2,dx2,dy2,dz2,origin2,densities,threshold):

######################################################################
# this function apply a threshold to the dx file. 
# This will remove point that are close to zero 
# note the use of aboslute values.
# output vn[i] = v[i], if |v[i]| > cutoff
#     or vn[i] = 0 otherwise
#
######################################################################
#def apply_threshold_values(xn,yn,zn,dx,dy,dz,origin,values,threshold):


def main():

   if len(sys.argv) != 4: # if no input
       print "ERORR:"
       print "syntex: dx-gist_make_courser.py gist-Eww-dens.dx 1 temp"
       #print "syntex: dx-combine_grids.py gist-Eww-dens.dx 3 temp"
       print " "
       return
 

   infile1 = sys.argv[1]
   number  = int(sys.argv[2])
   outfile = sys.argv[3]

   xn,yn,zn,dx,dy,dz,origin,values = read_in_dx_file(infile1)

   gridscale = dx # assumes that they are all the same spaceing

   xnew,ynew,znew,origin_new,gridscale_new,values_new = make_courser(outfile+'.dat',values,gridscale,xn,yn,zn,origin,number)

   write_out_dx_file(outfile+'.dx',xnew,ynew,znew,gridscale_new,gridscale_new,gridscale_new,origin_new,values_new)

main()


