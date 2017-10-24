import numpy as np
import math as m
#import matplotlib.pyplot as plt
from scipy import interpolate, ndimage
import sys
import os
from pkbinner import *



Base = './'  
#list_of_sims = ['LCDM','EXP001','EXP002','EXP003']
#snapshots = [92,56,41,35]
snapshots = [92,56]
list_of_sims = ['LCDM']
simu_set = len(list_of_sims)
counter = 0
list_snaps = [str(snap).zfill(3) for snap in snapshots ]

#list_particles = ['_type0_','_type1_','_']
list_particles = ['_']
#list_output =    ['_CoDECS_power_baryons_','_CoDECS_power_CDM_','_CoDECS_power_all_']
list_output =    ['_CoDECS_power_all_']

newBin = '/binned-new/'

bestparsfile='bestOverlapParams.txt'


fnames_list = [[[Base+simName+'/powerspec'+part+exts+'.txt' for part in list_particles] for exts in list_snaps] for simName in list_of_sims]

print("chosen names of files ")
print(fnames_list)

#__________________#
# Binning parameters
# These are the default values specified in the function rawPowers.
# If you want to change this, pass more arguments to the function rawPowers, such as:
# minModeCountA=25, minModeCountB=25, targetBinNumberA=125, targetBinNumberB=125
MinModeCount_A = 25 # counts the minimum number of modes per bin in order to reduce variance
TargetBinNumber_A = 125  # used to compute minimum logarithmic bin-size: k_range/TargetBinNumber_A
MinModeCount_B = 25 # counts the minimum number of modes per bin in order to reduce variance
TargetBinNumber_B = 125  # used to compute minimum logarithmic bin-size: k_range/TargetBinNumber_B

#__________________#

ndiscard_B = 4        # number of discarded bins at the end in k for toplevel mesh
ncut_A = 4           # number of saved  bins at the end of overlapping region in folded mesh 
swidth = 5            # smooth width: number of points around central to smooth with
thre_shot = 0.10      #  threshold in shot noise for discarding power spectra
intekind = 'linear'   #interpolation methods: 'linear','nearest', 'zero', 'slinear', 'quadratic', 'cubic'
#*******************#





for isi, simu in enumerate(list_of_sims):
    for j, snap in enumerate(list_snaps):
        for k, part in enumerate(list_particles):
            simuArr=[isi,j,k]

            Delta_A, K_A, Delta_B, K_B, Shotf = rawPowers(simuArr,fnames_list)
            PSArra = [Delta_A, K_A, Delta_B, K_B, Shotf] 
            
            
            
            fparstest = tryopen(bestparsfile,'best params file')
                       
            if fparstest==False:
	        print "Error reading best parameters file, continuing with default parameters."
            else:
                with open(bestparsfile) as fp:
                    countf=0
                    for line in fp:
                        paramslist=line.split()
                        if paramslist[0]==simu and paramslist[1]==snap:
                            ndiscard_B = int((paramslist[2].split('='))[1])
                            ncut_A     = int((paramslist[3].split('='))[1])
                            swidth     = int((paramslist[4].split('='))[1])
                            print "parameters adopted: \n"
                            print "ndiscardB (number of discarded bins from end of top mesh): "+str(ndiscard_B)
                            print "ncutA (number of discarded bins from beginning of folded mesh): "+str(ncut_A)
                            print "smoothwidth (width of points over to which apply smoothing): "+str(swidth)
                            countf=countf+1
                        elif (fparstest==True & countf==0):
                            print "Simulation numbers and/or snapshots in best parameters file do not match present settings."
            
            Delta2_all, K_all, Pk_all, K_list_A, Delta2_list_A, K_list_B, Delta2_list_B = overlapPS(PSArra,ndiscard_B,ncut_A,swidth,intekind,shotthreshold=thre_shot, fullout=True)
            
            

            #filling other arrays with zeros in order to write easily into file

            colsize=len(K_all)
            
            K_list_Afill = zerofill(K_list_A,colsize)
            Delta2_list_Afill = zerofill(Delta2_list_A,colsize)
            K_list_Bfill = zerofill(K_list_B,colsize)
            Delta2_list_Bfill = zerofill(Delta2_list_B,colsize)                              
            
            DataOutput = np.column_stack(( K_all, Pk_all, Delta2_all, K_list_Afill, Delta2_list_Afill, K_list_Bfill, Delta2_list_Bfill ))
            #DataOutput = np.column_stack(( K_list_all3, Delta2_list_all, Pk_all ))
            
            mkdirp(Base+simu+newBin)
            outfileSimu = Base+simu+newBin+simu+part+snap+'.txt'
            print '...printing: ', outfileSimu            
            np.savetxt(outfileSimu, DataOutput, fmt='%15.5e')
            with open(outfileSimu, 'a+') as fiE:
                fiE.write('#   '+'K in h/Mpc'.center(14)+'Power P(k)'.center(16)+'Delta2(k)'.center(16)+'K_A_Cut'.center(16)+'Delta2_A_Cut'.center(16)+'K_B_Cut'.center(16)+'Delta2_B_Cut'.center(16))

                fiE.write('\n'+'# params: '+'ndiscardB='+str(ndiscard_B)+" ncutA="+str(ncut_A)+" smoothwidth="+str(swidth))
            
                                      

            

print 'Finished'


        
            
            
      
      
      
      
      
      
      
    
    
    
  
  
  
