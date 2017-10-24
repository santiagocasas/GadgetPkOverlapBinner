import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy import interpolate, ndimage
import sys
import os


def mkdirp(what):
    try:
        os.makedirs(what)
    except OSError:
        pass


def tryopen(filenam, mssg="", mode='r+'):
    try:
        fop=open(filenam,mode)
        fop.close()
        return(True)
    except IOError as e:
        print(mssg+' file not found')
        return(False)

def zerofill(Arra,leng):
    if leng > len(Arra):
        Zeros = np.zeros(leng)
        for i in range(len(Arra)):
            Zeros[i]=Arra[i]
        return Zeros
    else:
        print "Error, array can not be filled with zeros"
        return Arra

def calcDerivs(Xvals, Yvals):
    deri = np.diff(Yvals)/np.diff(Xvals)
    xderi = (Xvals[1:]+Xvals[:-1])/2
    return deri, xderi


def calcRatioPS(dataSimu0, dataSimu1):
    ratiolist=np.copy(dataSimu0)
    dif = len(dataSimu0)-len(dataSimu1)
    if dif > 0:
        fInterp = interpolate.interp1d( dataSimu0[:,0], dataSimu0[:,1], kind='cubic')
        dataSimu0=dataSimu0[:-dif]
        ratiolist=ratiolist[:-dif]
        ratiolist[:,1]=dataSimu1[:,1]/fInterp(dataSimu1[:,0])
        #print str(dif)+' elements deleted for ratio'
    elif dif < 0:
        fInterp = interpolate.interp1d( dataSimu1[:,0], dataSimu1[:,1], kind='cubic')
        dataSimu1=dataSimu1[:dif]
        ratiolist[:,1]=fInterp(dataSimu0[:,0])/dataSimu0[:,1]
        #print str(dif)+' elements deleted for ratio'
    else:
        fInterp = interpolate.interp1d( dataSimu0[:,0], dataSimu0[:,1], kind='cubic')
        ratiolist[:,1]=dataSimu1[:,1]/fInterp(dataSimu1[:,0])
        #print 'equal size'
    return ratiolist



def ChiEstimator(deri2nd, kderi2nd, minMax):
    kregion = kderi2nd[(kderi2nd >= 1.0)]
    deri2region = deri2nd[(kderi2nd >= 1.0)]
    ChiVector = np.array([])
    for ik, ki in enumerate(deri2region[:-3]):
        quartet=[deri2region[ik],deri2region[ik+1],deri2region[ik+2],deri2region[ik+3]]
        aver=np.average(quartet)
        averd1=np.average((abs(quartet[1]-quartet[0]),abs(quartet[3]-quartet[2])))
        #(1+abs(aver))*#
        chi = (abs(quartet[2]-quartet[1]))/(1+abs(quartet[3]-quartet[0]))
        ChiVector = np.append(ChiVector,chi)
    
    return ChiVector, kregion[:-3]


def PowerBin(Bins, MinDlogK, MinModeCount, columns):
    istart = 0
    ind = [istart]
    K_list = []
    Delta2_list = []
    count_list = []
    ModeCount = columns[4]
    K         = columns[0]
    SumPower  = columns[8]
    ConvFac   = columns[9]
    Specshape = columns[7]
                
    while istart < Bins:
        count = sum(ModeCount[ind])
        deltak = m.log10(K[ind].max()) - m.log10(K[ind].min())
        
        if deltak >= MinDlogK and count >= MinModeCount:
            d2 = sum(SumPower[ind])/sum(ModeCount[ind])
            b = int(sum(ind*ModeCount[ind])/sum(ModeCount[ind]))
            kk = K[b]
            d2 = ConvFac[b]*d2*Specshape[b]
            K_list = np.append(K_list, kk)
            Delta2_list = np.append(Delta2_list, d2)
            count_list = np.append(count_list, sum(ModeCount[ind]))
            istart += 1
            ind = [istart]
        else:
            istart = istart + 1
            ind.append(istart)

    result = (K_list, Delta2_list, count_list)
    return result



def rawPowers(simuNums,fnames_list, minModeCountA=25, minModeCountB=25, targetBinNumberA=125, targetBinNumberB=125):
    isi=simuNums[0]
    j=simuNums[1]
    k=simuNums[2]
    
    fi = open(fnames_list[isi][j][k], "r")
    #print isi,j,k,"\n"
    #print fi.name
      
    time1 = float(fi.readline())
    binsA = int(fi.readline())
    dummy = float(fi.readline())
    dummyL = long(fi.readline())
    
    lines1 = [fi.readline() for i in range(binsA)]
    da1 = [[float(x) for x in lines1[i].split()] for i in range(len(lines1))]
    
    time2 = float(fi.readline())
    binsB = int(fi.readline())
    dummy = float(fi.readline())
    dummyL = long(fi.readline())
    
    lines2 = [fi.readline() for i in range(binsB)]
    da2 = [[float(x) for x in lines2[i].split()] for i in range(len(lines2))]
    
    fi.close()
    
    da1A = np.asarray(da1)
    da1B = np.asarray(da2)
    
    columnsA = [da1A[:,i] for i in range(10)]  #extract each column of the A array from text file
    K_A, Delta2_A, Shot_A, ModePow_A, ModeCount_A, Delta2Uncorrected_A, ModePowUncorrected_A, Specshape_A, SumPower_A, ConvFac_A  = columnsA
    
    columnsB = [da1B[:,i] for i in range(10)]  #extract each column of the B array from text file
    K_B, Delta2_B, Shot_B, ModePow_B, ModeCount_B, Delta2Uncorrected_B, ModePowUncorrected_B, Specshape_B, SumPower_B, ConvFac_B  = columnsB
    
  
    MinDlogK_A = (m.log10(max(K_A)) - m.log10(min(K_A)))/targetBinNumberA
    MinDlogK_B = (m.log10(max(K_B)) - m.log10(min(K_B)))/targetBinNumberB
  
    K_list_A, Delta2_list_A, count_list_A = PowerBin(binsA, MinDlogK_A, minModeCountA, columnsA)
  
    K_list_B, Delta2_list_B, count_list_B = PowerBin(binsB, MinDlogK_B, minModeCountB, columnsB)
    
    shotFunc = interpolate.interp1d(K_B, Shot_B, kind='linear')
    
    DA = np.copy(Delta2_list_A)
    KA = np.copy(K_list_A)
    KB = np.copy(K_list_B)
    DB = np.copy(Delta2_list_B)
    return DA, KA, DB, KB, shotFunc



def overlapPS(SimuPSArray,ndiscB,ncutA,swidth,intekind,shotthreshold=0.10,fullout='False'):
    threshot = shotthreshold     
    Delta2_list_A = np.copy(SimuPSArray[0])
    K_list_A = np.copy(SimuPSArray[1])    
    Delta2_list_B = np.copy(SimuPSArray[2])
    K_list_B = np.copy(SimuPSArray[3])
    shotFunc =  SimuPSArray[4]
    
    # discard last ndiscard_B entries of top level mesh
    K_list_B = K_list_B[:-ndiscB]
    Delta2_list_B = Delta2_list_B[:-ndiscB]    
    
    # leave only the desired entries from the folded mesh at the end of overlap region
    ind_A = np.where( K_list_A >= max(K_list_B) )[0]                 
    cutK = min(ind_A)-ncutA 
    K_list_A = K_list_A[cutK:]
    Delta2_list_A = Delta2_list_A[cutK:]           
    nA = len(K_list_A)
    nnA = len(ind_A)
    
    # range of overlapping power spectra    
    ind_AB = np.where( K_list_B >= min(K_list_A) )[0]
    nnB = len(ind_AB)            
    
    # k values of top level,overlap region and folded mesh
    K_list_AB = K_list_B[ind_AB]  
    K_list_B0 = K_list_B[(K_list_B <= min(K_list_A))]
    K_list_A3 = K_list_A[(K_list_A > max(K_list_B))]
    
    # we smooth the power with an uniform filter of size swidth
    Delta2_list_A = ndimage.filters.uniform_filter(Delta2_list_A, size=swidth)
    
    # we create an interpolated function of Delta2A in the range KA
    fInterpDA = interpolate.interp1d(K_list_A, Delta2_list_A, kind=intekind)
    
    # we overlap the top level and the folded powers according to a linear weight      
    linweight =     np.linspace(0,1,nnB)
    linweightcomp = 1-np.linspace(0,1,nnB)
    Delta2_list_AB1 = linweight*fInterpDA(K_list_AB)+linweightcomp*Delta2_list_B[ind_AB]
    
    # now we put together: the top level the overlapped and the folded mesh
    K_list_all = np.concatenate((K_list_B0,K_list_AB,K_list_A3))                        
    Delta2_list_all = np.concatenate((Delta2_list_B[(K_list_B <= min(K_list_A))],Delta2_list_AB1,Delta2_list_A[(K_list_A > max(K_list_B))]))
    
    # now we cut the power spectrum where the shot reaches a given threshold
    shotList = shotFunc(K_list_all)
    K_list_all = K_list_all[(shotList/Delta2_list_all) < threshot]
    Delta2_list_all = Delta2_list_all[(shotList/Delta2_list_all) < threshot]
    
    # convert to units of h/Mpc and calculate P(k) from D2(k)
    #K_list_all3 = 1.0e03*K_list_all
    K_list_all3 = K_list_all
    Pk_all      = 2.0*np.pi*np.pi*Delta2_list_all/np.power(K_list_all3, 3)
    krange = (max(K_list_AB) - min(K_list_AB))
    minmaxOverlap = (1000*min(K_list_AB),1000*max(K_list_AB),1000*krange)
    minOv = minmaxOverlap[0]
    maxOv = minmaxOverlap[1]
    rangOv = minmaxOverlap[2]
    avOv = np.average((minmaxOverlap[0],minmaxOverlap[1]))
    if fullout==False:
        return Delta2_list_all, K_list_all3, minmaxOverlap
    elif fullout==True:
        return Delta2_list_all, K_list_all3, Pk_all, K_list_A, Delta2_list_A, K_list_B, Delta2_list_B

