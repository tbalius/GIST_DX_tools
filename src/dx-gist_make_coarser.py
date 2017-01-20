
import sys, os, math
import dx_gist_lib as dxlib

# this is written by Trent Balius in the Shoichet Lab
# written in Oct, 2014
#
# here we read in, manpulate and output dx file.
# the dx file are produed by GIST calculations.

# here we want to make the grid courser. 


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

   xn,yn,zn,dx,dy,dz,origin,values = dxlib.read_in_dx_file(infile1)

   gridscale = dx # assumes that they are all the same spaceing

   xnew,ynew,znew,origin_new,gridscale_new,values_new = dxlib.make_courser(outfile+'.dat',values,gridscale,xn,yn,zn,origin,number)

   dxlib.write_out_dx_file(outfile+'.dx',xnew,ynew,znew,gridscale_new,gridscale_new,gridscale_new,origin_new,values_new)

main()


