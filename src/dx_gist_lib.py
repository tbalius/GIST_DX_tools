
import sys, os

# this is written by Trent Balius in the Shoichet Lab
# written in 2014-2015
#
# here we read in, manpulate and output dx file.
# the dx file are produed by GIST calculations.


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
         print ("count = ", count, " line = ", line)
         xn = int(splitline[5])
         yn = int(splitline[6])
         zn = int(splitline[7])
         
      ## this should be line 2       
      if (splitline[0] == "origin"):
         #print line
         print ("count = ", count, " line = ", line)
         origin = [float(splitline[1]), float(splitline[2]), float(splitline[3])] 

      ## this should be lines 3-5
      if (splitline[0] == "delta"):
         #print line
         print ("count = ", count, " line = ", line)
         if (float(splitline[2]) == 0 and  float(splitline[3]) ==0):
            dx = float(splitline[1]) 
         elif (float(splitline[1]) == 0 and  float(splitline[3]) ==0):
            dy = float(splitline[2]) 
         elif (float(splitline[1])== 0 and  float(splitline[2])==0):
            dz = float(splitline[3]) 
            print (dx, dy, dz )


      if (splitline[0] == "object" and splitline[1] == "2"):
         #print line
         print ("count = ", count, " line = ", line)
      elif (splitline[0] == "object" and splitline[1] == "3"):
         #print line
         print ("count = ", count, " line = ", line)
         flag_read_dx = True
         continue # go to next line
      elif (splitline[0] == "object" ):
         print (" line = ", line)
         continue
      if (flag_read_dx):

         if (len(splitline) > 3): 
            print ("Error: dx formate problem. more than 3 colums")
            exit()

         for value in splitline:
             values.append(float(value))  

      count = count + 1


  print (len(values))
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
  #fileh.write('origin %6.2f %6.2f %6.2f\n' % (origin[0],origin[1],origin[2]))
  fileh.write('origin %7.4f %7.4f %7.4f\n' % (origin[0],origin[1],origin[2]))
  fileh.write('delta %6.3f 0 0\n' % dx)
  fileh.write('delta 0 %6.3f 0\n' % dy)
  fileh.write('delta 0 0 %6.3f\n' % dz)
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
# this function combines values file from the dx files.
# function v[i] = w1*v1[i] + w2*v2[i]+c3
######################################################################
def combine_values(file,xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,w1,xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2,w2,c3):
    if ( (xn1 != xn2) or (yn1 != yn2) or (zn1 != zn2) or (dx1 != dx2)
          or (dy1 != dy2) or (dz1 != dz2) or (origin1[0] != origin2[0])
          or (origin1[1] != origin2[1]) or (origin1[2] != origin2[2])):
       print ("header does not match")
    if (len(values1) != len(values2)):
       print ("the values don't have the same number of elements.")
       print (len(values1), len(values2))

    fileh = open(file,'w')

    new_values = []
    for i in range(len(values2)):
        #if (values1[i] == 0 or values2[i] == 0):
        if (values1[i] == 0 and values2[i] == 0): # change requested by Mossa Ghattas & Reed Stien. 
            value = 0.0
        else:
            value = w1*values1[i] + w2*values2[i]+c3

        #print "%f = %f * %f +  %f * %f + %f" % (value,w1,values1[i],w2,values2[i],c3)
        fileh.write("%f = %f * %f +  %f * %f + %f\n" % (value,w1,values1[i],w2,values2[i],c3))

        new_values.append(value)
    fileh.close()
    return new_values

       
######################################################################
# this function applies a density cutoff 
######################################################################
######################################################################
# this function reads in energy values and densities as vectors.
# output vn[i] = v[i], if d[i] > cutoff
#     or vn[i] = 0 otherwise
######################################################################
def combine_values_density_threshold(file,xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2,t):
    if ( (xn1 != xn2) or (yn1 != yn2) or (zn1 != zn2) or (dx1 != dx2)
          or (dy1 != dy2) or (dz1 != dz2) or (origin1[0] != origin2[0]) 
          or (origin1[1] != origin2[1]) or (origin1[2] != origin2[2])):
       print ("header does not match")
    if (len(values1) != len(values2)):
       print ("the values don't have the same number of elements.")
       print (len(values1), len(values2))

    fileh = open(file,'w')

    new_values = []
    for i in range(len(values2)):
        if (values2[i] < t):
            value = 0.0
        else: 
            value = values1[i]

        #print "%f = %f * %f +  %f * %f + %f" % (value,w1,values1[i],w2,values2[i],c3)
        fileh.write("%f ::  energy, %f; density, %f; threshold, %f.\n" % (value,values1[i],values2[i],t))

        new_values.append(value)
    fileh.close()
    return new_values

######################################################################
# this function applies a density and energy cutoff 
######################################################################
def combine_values_energy_density_thresholds(file,xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,t1,xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2,t2):
    if ( (xn1 != xn2) or (yn1 != yn2) or (zn1 != zn2) or (dx1 != dx2)
          or (dy1 != dy2) or (dz1 != dz2) or (origin1[0] != origin2[0])
          or (origin1[1] != origin2[1]) or (origin1[2] != origin2[2])):
       print ("header does not match")
    if (len(values1) != len(values2)):
       print ("the values don't have the same number of elements.")
       print (len(values1), len(values2))

    fileh = open(file,'w')

    new_values = []
    for i in range(len(values2)):
        if (values2[i] < t2 or values1[i] < t1): 
            value = 0.0
        else:   # high densitiy and high energy (unfavorable waters, red blobs)
            #value = values1[i]
            value = t1

        #print "%f = %f * %f +  %f * %f + %f" % (value,w1,values1[i],w2,values2[i],c3)
        fileh.write("%f ::  energy, %f; density, %f; thresholds, %f, %f.\n" % (value,values1[i],values2[i],t1,t2))

        new_values.append(value)
    fileh.close()
    return new_values


######################################################################
# this function will divide by the density. 
# values1 is the energy, values2 is the density, C_num_dens is the number denisty in bulk
# this fuction will return the normalized energy values. 
######################################################################
def divide_energy_density(file,xn1,yn1,zn1,dx1,dy1,dz1,origin1,values1,xn2,yn2,zn2,dx2,dy2,dz2,origin2,values2,C_num_dens):
    if ( (xn1 != xn2) or (yn1 != yn2) or (zn1 != zn2) or (dx1 != dx2)
          or (dy1 != dy2) or (dz1 != dz2) or (origin1[0] != origin2[0])
          or (origin1[1] != origin2[1]) or (origin1[2] != origin2[2])):
       print ("header does not match")
    if (len(values1) != len(values2)):
       print ("the values don't have the same number of elements.")
       print (len(values1), len(values2))

    #print C_num_dens

    fileh = open(file,'w')

    new_values = []
    for i in range(len(values2)):
        if (values2[i] != 0.0):
           value = values1[i] / (values2[i] * C_num_dens)
        else:
           value = 0.0

        #print "%f = %f * %f +  %f * %f + %f" % (value,w1,values1[i],w2,values2[i],c3)
        fileh.write("%f ::  energy, %f; density, %f; thresholds, %f.\n" % (value,values1[i],values2[i],C_num_dens))

        new_values.append(value)
    fileh.close()
    return new_values



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




