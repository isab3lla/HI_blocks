import numpy as np
import HI_hybrid_lib as HIL
import readsnap
import Power_spectrum_library as PSL
import CIC_library as CIC
import sys


rho_crit=2.77536627e11 #h^2 Msun/Mpc^3
#############################################################################

## dims for calculating the Pk
dims=256

## info on the sim and snapshots

sim = '500_1024_171201'  #14 snapshots available
total_snap = 14
z_arr   = np.array([2.8,2.4,2.1,1.8,1.5,1.3,1.1,0.9,0.8,0.7,0.6,0.5,0.4,0.3])
z_str   = ['2p800','2p400','2p100','1p800','1p500','1p300','1p100',
					'0p900','0p800','0p700','0p600','0p500','0p400','0p300']

isnap=13

# sim = '120_256_180710'
#isnap = 0

#z_arr   = np.array([0.0])
#snap_ar = np.array(['00'])
#z_str   = ['0p000']



dumb1,dumb2,dumb3,snapshot_fname = HIL.getting_sim_paths(sim,isnap)

f_out = 'HI_Pk_'+sim[:-6]+'_z'+z_str[isnap]+'_fromSnap.dat'


#sim = 'super_new_snap_120'
#sim = 'final_snap_120'

dims3=dims**3
#############################################################################

#read snapshot head and obtain BoxSize, Omega_m and Omega_L
print '\nREADING SNAPSHOTS PROPERTIES'
head=readsnap.snapshot_header(snapshot_fname)
BoxSize=head.boxsize/1e3 #Mpc/h
Nall=head.nall
Masses=head.massarr*1e10 #Msun/h
Omega_m=head.omega_m
Omega_l=head.omega_l
redshift=head.redshift
Hubble=100.0*np.sqrt(Omega_m*(1.0+redshift)**3+Omega_l)  #h*km/s/Mpc
h=head.hubble

#find the total number of particles in the simulation
Ntotal=np.sum(Nall,dtype=np.uint64)
print 'Total number of particles in the simulation =',Ntotal

#sort the pos array
ID_unsort=readsnap.read_block(snapshot_fname,"ID  ",parttype=1)-1 #normalized
pos_unsort=readsnap.read_block(snapshot_fname,"POS ",parttype=1)/1e3 #Mpc/h
pos=np.empty((Ntotal,3),dtype=np.float32); pos[ID_unsort]=pos_unsort

print 'reading the HI mass of all the particles . . .'
M_HI_unsort=readsnap.read_block(snapshot_fname,"HIMS",parttype=1)*1e10 #Msun/h
print '\n done with the reading!'
M_HI=np.empty((Ntotal),dtype=np.float32); M_HI[ID_unsort]=M_HI_unsort

print 'assigned mass to particles:'
print '%e < M [Msun/h] < %e'%(np.min(M_HI),np.max(M_HI))
print len(M_HI), Ntotal

#compute the value of Omega_HI
Omega_HI=np.sum(M_HI,dtype=np.float64)/BoxSize**3/rho_crit
print '\nOmega_HI (recalculated) = %e'%Omega_HI

sys.exit()

#only keep particles with M_HI>0
indexes=np.where(M_HI>0.0)[0]; M_HI=M_HI[indexes]; pos=pos[indexes]
del indexes

#compute the mean neutral hydrogen mass per grid point
mean_M_HI=np.sum(M_HI,dtype=np.float64)/dims3
print 'mean HI mass per grid point = %2.2e Msun/h \n'%mean_M_HI


#compute the value of delta_HI = rho_HI / <rho_HI> - 1
delta_HI=np.zeros(dims3,dtype=np.float32)
print 'Computing the M_HI in the grid cells using CIC'
CIC.CIC_serial(pos,dims,BoxSize,delta_HI,M_HI) 
print '%e should be equal to:\n%e' %(np.sum(M_HI,dtype=np.float64),
                                       np.sum(delta_HI,dtype=np.float64))
delta_HI=delta_HI/mean_M_HI-1.0  #computes delta
print 'numbers may be equal:',np.sum(delta_HI,dtype=np.float64),0.0
print np.min(delta_HI),'< delta_HI <',np.max(delta_HI)

#compute the HI PS
Pk=PSL.power_spectrum_given_delta(delta_HI,dims,BoxSize)

#write total HI P(k) file
f=open(f_out,'w')
for i in range(len(Pk[0])):
    f.write(str(Pk[0][i])+' '+str(Pk[1][i])+'\n')
f.close()




