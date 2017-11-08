import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy import interpolate, ndimage
import sys
import os
from pkbinner import *



Base = './'  
#list_of_sims = ['LCDM','EXP001','EXP002','EXP003']
#snapshots = [92,56,41,35]
snapshots = [175,400]
list_of_sims = ['GR', 'F5']
simu_set = len(list_of_sims)
counter = 0
list_snaps = [str(snap).zfill(3) for snap in snapshots ]

#list_particles = ['_type0_','_type1_','_']
list_particles = ['_']
#list_output =    ['_CoDECS_power_baryons_','_CoDECS_power_CDM_','_CoDECS_power_all_']

#folder to store new binned Pk
newBin = '/binned-optim-ratio/'

#file where optimal binning and overlapping parameters are found
bestparamsfile="bestOverlapParams.txt"


fnames_list = [[[Base+simName+'/powerspec'+part+simName+part+exts+'.txt' for part in list_particles] for exts in list_snaps] for simName in list_of_sims]

print("chosen names of files ")
print(fnames_list)

#____#
#Units
unitsfact=1  # MGadget files have k in h/Mpc


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

# Interpolation, overlapping and smoothing parameters

ndiscard_B = 4        # number of discarded bins at the end in k for toplevel mesh
ncut_A = 5            # number of saved  bins at the end of overlapping region in folded mesh 
swidth = 3            # smooth width: number of points around central to smooth with
thre_shot = 0.10      #  threshold in shot noise for discarding power spectra
intekind = 'linear'   #interpolation methods: 'linear','nearest', 'zero', 'slinear', 'quadratic', 'cubic'
#*******************#

# list of parameters to explore in the loop. The parameters that minimize the oscillations of the derivatives of the ratio, are taken.
ndiscard_B_test_arra=[2,3,4,5]
ncut_A_test_arra=[3,4,5]
swidth_test_arra=[3,4,5]

bestTotalparams = []



fpars = open(bestparamsfile,'w+')

 

for j, snap in enumerate(list_snaps):
    simuL = [0,j,0]
    LDelta_A, LK_A, LDelta_B, LK_B, LShotf, LcountA, LcountB = rawPowers(simuL,fnames_list)
    LPSArra = [LDelta_A, LK_A, LDelta_B, LK_B, LShotf, LcountA, LcountB]
    #L is the baseline simulation. Spectra from simulation E will be divided by L and the ratio will be optimized. 
    for isi, simu in enumerate(list_of_sims[1:]):
        for k, part in enumerate(list_particles):
            
            simuE = [(isi+1),j,k]
            EDelta_A, EK_A, EDelta_B, EK_B, EShotf, EcountA, EcountB = rawPowers(simuE,fnames_list)
            EPSArra = [EDelta_A, EK_A, EDelta_B, EK_B, EShotf, EcountA, EcountB]       
            
            # parameters for loop minimization of estimator
            maxChi = 1000
            for ndiscard_B in ndiscard_B_test_arra:
                for ncut_A in ncut_A_test_arra:
                    for swidth in swidth_test_arra: 
                        
                        LD2_all, LK_all, LminMax = overlapPS(LPSArra,ndiscard_B,ncut_A,swidth,intekind,shotthreshold=thre_shot, fullout=False, unitsfactor=unitsfact)
                        ED2_all, EK_all, EminMax = overlapPS(EPSArra,ndiscard_B,ncut_A,swidth,intekind,shotthreshold=thre_shot, fullout=False, unitsfactor=unitsfact)
                        
                        LSimuData = np.transpose([LK_all,LD2_all])
                        ESimuData = np.transpose([EK_all,ED2_all])
                        
                        ratioData = calcRatioPS(LSimuData, ESimuData)
                        
                        ratio1stderivY, ratio1stderivX = calcDerivs(ratioData[:,0], ratioData[:,1])
                        ratio2ndderivY, ratio2ndderivX = calcDerivs(ratio1stderivX, ratio1stderivY)
                                                
                        ChiVec, Chi_kregion = ChiEstimator(ratio2ndderivY, ratio2ndderivX, LminMax)
                        
                        #plt.pause(0.5)
                        maxChiTemp = max(ChiVec)
                        print simu, snap, ndiscard_B,ncut_A, swidth, "maxChi= ", maxChiTemp, '\n'
                        if (maxChiTemp < maxChi) :
                            maxChi = maxChiTemp
                            bestSimuSnapParams = [simu,snap,ndiscard_B,ncut_A,swidth,maxChi]
                            ChiVecMax = np.copy(ChiVec)
                            
            
            #writing the best (optimal) overlapping parameters onto a file
            print "best= ", bestSimuSnapParams, '\n'
            print "saving best parameters to file: "+bestparamsfile
            bestTotalparams.append(bestSimuSnapParams)
            fpars.write(str(bestSimuSnapParams[0])+" "+str(bestSimuSnapParams[1])+" "+"ndiscardB="+str(bestSimuSnapParams[2])+" "+"ncutA="+str(bestSimuSnapParams[3])+" "+"smoothwidth="+str(bestSimuSnapParams[4])+" "+"maxChi="+str(bestSimuSnapParams[5])+"\n")
            #print parameters also for the baseline simulation L
            fpars.write(str(list_of_sims[0])+" "+str(bestSimuSnapParams[1])+" "+"ndiscardB="+str(bestSimuSnapParams[2])+" "+"ncutA="+str(bestSimuSnapParams[3])+" "+"smoothwidth="+str(bestSimuSnapParams[4])+" "+"maxChi="+str(bestSimuSnapParams[5])+"\n")
            
            #applying the best (optimal) parameters to the raw unbinned power spectra
            LD2_all_best, LK_all_best, LPk_all, Lcount_all, LK_list_A, LDelta2_list_A, Lcount_A, LK_list_B, LDelta2_list_B, Lcount_B, Lshot_list = overlapPS(LPSArra,bestSimuSnapParams[2],bestSimuSnapParams[3],bestSimuSnapParams[4],intekind,shotthreshold=thre_shot,fullout=True, unitsfactor=unitsfact)
            ED2_all_best, EK_all_best, EPk_all, Ecount_all, EK_list_A, EDelta2_list_A, Ecount_A, EK_list_B, EDelta2_list_B, Ecount_B, Eshot_list = overlapPS(EPSArra,bestSimuSnapParams[2],bestSimuSnapParams[3],bestSimuSnapParams[4],intekind,shotthreshold=thre_shot,fullout=True, unitsfactor=unitsfact)
            #filling arrays with zeros, in order to save all columns to a file
            
            colsize=len(LD2_all_best)
            LK_list_A = zerofill(LK_list_A,colsize)
            Lcount_A = zerofill(Lcount_A,colsize)
            LDelta2_list_A = zerofill(LDelta2_list_A,colsize)
            LK_list_B = zerofill(LK_list_B,colsize)
            Lcount_B = zerofill(Lcount_B,colsize)
            LDelta2_list_B = zerofill(LDelta2_list_B,colsize)
            
            colsize=len(ED2_all_best)
            EK_list_A = zerofill(EK_list_A,colsize)
            Ecount_A = zerofill(Ecount_A,colsize)
            EDelta2_list_A = zerofill(EDelta2_list_A,colsize)
            EK_list_B = zerofill(EK_list_B,colsize)
            Ecount_B = zerofill(Ecount_B,colsize)
            EDelta2_list_B = zerofill(EDelta2_list_B,colsize)
            
            DataOutputL = np.column_stack((LK_all_best, LPk_all, LD2_all_best, Lcount_all, LK_list_A, LDelta2_list_A, Lcount_A, LK_list_B, LDelta2_list_B, Lcount_B, Lshot_list))
            DataOutputE = np.column_stack((EK_all_best, EPk_all, ED2_all_best, Ecount_all, EK_list_A, EDelta2_list_A, Ecount_A, EK_list_B, EDelta2_list_B, Ecount_B, Eshot_list))
            
            mkdirp(Base+simu+newBin)
            
            #export the L binned spectra corresponding to the best parameters chosen, such that the ratio E/L is optimal.
            outfileL = Base+simu+newBin+list_of_sims[0]+part+snap+'.txt'
            #export the binned spectra for E corresponding to the best parameters chosen, such that the ratio E/L is optimal.
            outfileE = Base+simu+newBin+simu+part+snap+'.txt'
            
            print '...printing: ', outfileL, outfileE
            
            #saving Pk data to files
            np.savetxt(outfileL, DataOutputL, fmt='%15.5e')
            np.savetxt(outfileE, DataOutputE, fmt='%15.5e')
            
            #saving extra information to Pk files (for LCDM simulation)
            fiL = open(outfileL, 'a+')
            fiL.write('#   '+'K in h/Mpc'.center(14)+'Power P(k)'.center(16)+'Delta2(k)'.center(16)+'modeCount(k)'.center(16)+'K_A_Cut'.center(16)+'Delta2_A_Cut'.center(16)+'modeCountA(k)'.center(16)+'K_B_Cut'.center(16)+'Delta2_B_Cut'.center(16)+'modeCountB(k)'.center(16)+'ShotNoise(k)'.center(16))
            fiL.write('\n'+'# params: '+'ndiscardB='+str(bestSimuSnapParams[2])+" ncutA="+str(bestSimuSnapParams[3])+" smoothwidth="+str(bestSimuSnapParams[4]))
            fiL.close()
            
            #another pythonic  way of saving extra information to Pk files (for EXP simulation)
            with open(outfileE, 'a+') as fiE:
                fiE.write('#   '+'K in h/Mpc'.center(14)+'Power P(k)'.center(16)+'Delta2(k)'.center(16)+'modeCount(k)'.center(16)+'K_A_Cut'.center(16)+'Delta2_A_Cut'.center(16)+'modeCountA(k)'.center(16)+'K_B_Cut'.center(16)+'Delta2_B_Cut'.center(16)+'modeCountB(k)'.center(16)+'ShotNoise(k)'.center(16))
                fiE.write('\n'+'# params: '+'ndiscardB='+str(bestSimuSnapParams[2])+" ncutA="+str(bestSimuSnapParams[3])+" smoothwidth="+str(bestSimuSnapParams[4]))
            
    
            counter=counter+1
            
            
            
            
            
            
            


fpars.close()



print 'Finished'

 
