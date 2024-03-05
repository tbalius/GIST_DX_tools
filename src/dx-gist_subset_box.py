
import sys, os, math
import mol2
import pdb_lib as pdb
import dx_gist_lib as dxlib

# This is written by Trent Balius in the Shoichet Lab
# written in march, 2016
# modifed in 2023/10 by Trent Balius at the FNLCR

def distance(v1,v2):
    if (len(v1)!=len(v2)):
       print "error" 
       exit()
    dist = 0.0
    for i in range(len(v1)):
        dist = dist + (v1[i]-v2[i])**2.0
    dist = math.sqrt(dist)
    return dist

def cal_grid(values,gridscale,xn,yn,zn,origin):
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
    return grid_old


def cal_box_min_max(box,pad):
    print len(box)
    #maxv = [box.atoms[0].X, box.atoms[0].Y, box.atoms[0].Z]
    #minv = [box.atoms[0].X, box.atoms[0].Y, box.atoms[0].Z]
    maxv = [box[0].X, box[0].Y, box[0].Z]
    minv = [box[0].X, box[0].Y, box[0].Z]

    #for a in box.atoms: 
    for a in box: 
         if maxv[0] < a.X:
             maxv[0] = a.X
         if minv[0] > a.X:
             minv[0] = a.X 

         if maxv[1] < a.Y:
             maxv[1] = a.Y 
         if minv[1] > a.Y:
             minv[1] = a.Y 

         if maxv[2] < a.Z:
             maxv[2] = a.Z 
         if minv[2] > a.Z:
             minv[2] = a.Z 

    maxv[0] = maxv[0] + pad
    maxv[1] = maxv[1] + pad
    maxv[2] = maxv[2] + pad

    minv[0] = minv[0] - pad
    minv[1] = minv[1] - pad
    minv[2] = minv[2] - pad

    return maxv, minv

######################################################################
# this function calculates gist score. 
# also calculates peratom controbutions distributing overlap 
# atoms equaly among all atoms contating the point
# also calculates volume displaced by molcule.
######################################################################

def calc_grid_in_box(values,gridscale,xn,yn,zn,origin,box_max,box_min):
# file - log file
# values gist values 
# gridscale - the grid spacing

     # we should check the new box is completely inside old box
     if (origin[0] > box_min[0] ) or (origin[1] > box_min[1] ) or (origin[2] > box_min[2] ):
         print("new box is not contained in the old box (min)")
         exit()

     oldmax_x = xn * gridscale + origin[0]
     oldmax_y = yn * gridscale + origin[1]
     oldmax_z = zn * gridscale + origin[2]
     
     if (oldmax_x < box_max[0] ) or (oldmax_y < box_max[1] ) or (oldmax_z < box_max[2] ): 
         print("new box is not contained in the old box (max)")
         exit() 
     
     min_grid_i = int(math.ceil((box_min[0] - origin[0] ) / gridscale))
     min_grid_j = int(math.ceil((box_min[1] - origin[1] ) / gridscale))
     min_grid_k = int(math.ceil((box_min[2] - origin[2] ) / gridscale))

     max_grid_i = int(math.floor((box_max[0] - origin[0] ) / gridscale))
     max_grid_j = int(math.floor((box_max[1] - origin[1] ) / gridscale))
     max_grid_k = int(math.floor((box_max[2] - origin[2] ) / gridscale))
  
     new_values = [] 
     for i in range(min_grid_i,max_grid_i):  
         for j in range(min_grid_j,max_grid_j):  
             for k in range(min_grid_k,max_grid_k):  
                    count = i*(yn*zn)+j*(zn) + k
                    new_values.append(values[count])

     xnn = max_grid_i - min_grid_i #+ 1
     ynn = max_grid_j - min_grid_j #+ 1
     znn = max_grid_k - min_grid_k #+ 1
     
     #max_grid_x = max_grid_i * gridscale + origin[0] 
     #max_grid_y = max_grid_j * gridscale + origin[1] 
     #max_grid_z = max_grid_k * gridscale + origin[2] 
     
     min_grid_x = min_grid_i * gridscale + origin[0] 
     min_grid_y = min_grid_j * gridscale + origin[1] 
     min_grid_z = min_grid_k * gridscale + origin[2] 
     
     originn = [min_grid_x, min_grid_y, min_grid_z]
     
     return new_values, xnn, ynn, znn, originn

def main():


   if len(sys.argv) != 5: # if no input
       print "ERORR:"
       print "syntex: dx-gist_bounding_box.py dx-file box.pdb pad output"
       print "dx-input-file input file in dx formate produed by gist, may be disities or energies"
       return
 

   infiledx     = sys.argv[1]
   infilebox    = sys.argv[2]
   pad          = float(sys.argv[3])
   outfile      = sys.argv[4]

   print (infiledx)
   print (infilebox)

   #vdwdict = intialize_vdw_parm('/nfs/home/tbalius/zzz.github/DOCK/proteins/defaults/vdw.parms.amb.mindock') 
   #vdwdict = intialize_vdw_parm('/home/baliuste/zzz.github/DOCK/ucsfdock/proteins/defaults/vdw.parms.amb.mindock') 
   #DOCKpath = os.getenv('DOCKBASE')
   #vdwdict = intialize_vdw_parm(DOCKpath+'/proteins/defaults/vdw.parms.amb.mindock') 

   xn,yn,zn,dx,dy,dz,origin,values = dxlib.read_in_dx_file(infiledx)

   gridscale = dx # assumes that they are all the same spaceing
   
   
   #box  = pdb.read_pdb(infilebox)[0]
   box  = pdb.read_pdb(infilebox)

   #xv, nv = cal_box_min_max(box) # get xv is maX Value and nv is miN Value. 
   box_max, box_min = cal_box_min_max(box,pad) # get xv is maX Value and nv is miN Value. 

   file1 = open(outfile+'gist_values.txt','w')

   #new_values = calc_score(outfile,values,gridscale,xn,yn,zn,origin, mol,vdwdict,file1,(N<3),Hchoice)
   new_values, new_xn,new_yn,new_zn,new_origin = calc_grid_in_box(values,gridscale,xn,yn,zn,origin,box_max,box_min)

   dxlib.write_out_dx_file(outfile+".dx",new_xn,new_yn,new_zn,dx,dy,dz,new_origin,new_values)

   file1.close()
main()


