import matplotlib.pyplot as plt
import numpy as np
import sys
import copy
import matplotlib
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from scipy import interpolate
import tools

Base = '../' 
binfolderList=['/binned-optim-ratio/']
#simname=str(sys.argv[0])

#list_of_sims = ['LCDM', 'EXP001']
#snapshots = [92,56,41,35]
#list_of_reds = ['_z_000','_z_055','_z_100','_z_161']
list_of_sims = ['GR','F5']
snapshots = [175,400]
list_of_reds = ['z=1.0','z=0.']
list_particles = ['_']
#fill the snapshot numbers with zeros, to match filenames
list_snaps = [str(snap).zfill(3) for snap in snapshots ]

pklabels=['','$P(k)$', '$\Delta(k)$']
pkindex=1
#pkindex=1 -> dimensionful P(k),  pkindex=2 -> dimensionless Delta(k)

simu_set = len(list_of_sims)


#choose which file to plot
#simulation index
simindex=[0,1]
#snapshot (redshift) index
snaindex=[0,1]

#place where spectra for all sims are located
locpks=list_of_sims[1]

list_of_files_and_labels = [[(Base+locpks+binfolderList[0]+list_of_sims[sim]+list_particles[0]+str(list_snaps[sna])+'.txt', 
        'MGGadget '+pklabels[pkindex]+' '+list_of_sims[sim]+' '+list_of_reds[sna]) for sna in snaindex] for sim in simindex]

datalist = [[[np.loadtxt(list_of_files_and_labels[sim][sna][0])[:,(0,pkindex)], list_of_files_and_labels[sim][sna][1]] for sna in snaindex] for sim in simindex]



G = gridspec.GridSpec(4,4)

fig=plt.figure(1, figsize=(20,12), dpi=80,facecolor='w')

collist=[['b','g','r','c'], ['Indigo','Olive','OrangeRed','SkyBlue']]
markers=['-p', '-v', '-o','-s']

axes1 = fig.add_subplot(G[:,:])

minorLocator   = ticker.LogLocator(base=10,subs=np.linspace(1.0,10,10))

for sim in simindex:
    for sna in snaindex:
        axes1.loglog(datalist[sim][sna][0][:,0],datalist[sim][sna][0][:,1], markers[sna], label=datalist[sim][sna][1], color=collist[sna][sim], lw=3, ms=6, markevery=10, alpha=1.0 )
        axes1.text(0.02,0.1+0.1*sim,'Non-Linear power spectrum '+list_of_sims[sim], verticalalignment='top',horizontalalignment='left',transform=axes1.transAxes,fontsize=14,style='italic') 

xmini,xmaxi,ymini,ymaxi = axes1.axis('tight')
axes1.axis(ymax=ymaxi*1.2, xmin=xmini*0.9, xmax=xmaxi*1.5, ymin=ymini)
axes1.legend(loc='best',markerscale=1.5,prop={'size':12},numpoints=3,handlelength=3)
axes1.grid(True,which="major",ls=":")
axes1.tick_params(which='both',length=6, width=1, labeltop=True)
axes1.set_ylabel(pklabels[pkindex],size='x-large')
axes1.xaxis.set_minor_locator(minorLocator)
axes1.set_xlabel("k in $h/Mpc$",size='x-large')
minorLocator   = ticker.LogLocator(0.1,[0.01,0.02])

#axes2 = fig.add_subplot(G[2:,:])
sep="__"
joinreds=sep.join([list_of_reds[sna] for sna in snaindex])
joinreds=joinreds.replace('.','p')

joinsims=sep.join([list_of_sims[sim] for sim in simindex])

fig.savefig('nonlinear-Pk_'+joinsims+'_'+joinreds+'.png',dpi=200)

plt.show()
