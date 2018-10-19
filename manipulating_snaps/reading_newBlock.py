import numpy as np
import readsnap
import sys,os


################################### INPUT #####################################

#snapshot_fname='/scratch/ipcarucci/lpicola_trial/output_files/60Mpc512Nall_z3p000'
#snapshot_fname='/scratch/ipcarucci/sims/snapshots/60_512_z3p000'

#snapshot_fname='marker_snap'
snapshot_fname='new_snap'



rho_crit=2.77536627e11 #h^2 Msun/Mpc^3
##############################################################################

#read snapshot head and obtain BoxSize, Omega_m and Omega_L
print '\n  READING SNAPSHOTS PROPERTIES\n'
head=readsnap.snapshot_header(snapshot_fname)
BoxSize=head.boxsize/1e3 #Mpc/h
Nall=head.nall
Masses=head.massarr*1e10 #Msun/h
Omega_m=head.omega_m
Omega_l=head.omega_l
redshift=head.redshift
Hubble=100.0*np.sqrt(Omega_m*(1.0+redshift)**3+Omega_l)  #h*km/s/Mpc
h=head.hubble
filenum=head.filenum

print Nall


print 'boxsize = %.1f Mpc/h, z = %.1f\n'%(BoxSize, redshift)
print 'DM particle mass = %.2e Msun/h\n'%Masses[1]
print 'number of particles = %d\n'%Nall[1]
print 'number of subfiles in the snapshot = %d\n'%filenum

#compute the values of Omega_CDM and Omega_B
Omega_cdm=Nall[1]*Masses[1]/BoxSize**3/rho_crit
Omega_b=Omega_m-Omega_cdm
print '\nOmega_CDM = %.3f\nOmega_B   = %0.3f\nOmega_M   = %.3f\n'\
    %(Omega_cdm,Omega_b,Omega_m)

#sys.exit()

#read the positions of all the particles
print 'reading the positions of all the particles . . .'
pos=readsnap.read_block(snapshot_fname,"POS ",parttype=1)/1e3 #Mpc/h
print '\n done with the reading!'
print '%.2f < X [Mpc/h] < %.2f'%(np.min(pos[:,0]),np.max(pos[:,0]))
print '%.2f < Y [Mpc/h] < %.2f'%(np.min(pos[:,1]),np.max(pos[:,1]))
print '%.2f < Z [Mpc/h] < %.2f\n'%(np.min(pos[:,2]),np.max(pos[:,2]))

#read the IDs of all the particles
print 'reading the IDs of all the particles . . .'
IDs=readsnap.read_block(snapshot_fname,"ID  ",parttype=1)
print '\n done with the reading!'

#read the marker of all the particles
print 'reading the markers of all the particles . . .'
marker=readsnap.read_block(snapshot_fname,"HIMS",parttype=1)
print '\n done with the reading!'
print '%.2f < marker < %.2f\n'%(np.min(marker[:]),np.max(marker[:]))

#for i in range(Nall[1]):
#	print 'ID: ',IDs[i], ' marker: ',marker[i]


## ordering the marker per ID
sorting_IDs_order = IDs.argsort()
sorted_marker = marker[sorting_IDs_order]

print sorted_marker[4999980:5000020]
print '\nciaoooo\n'
print sorted_marker[9999980:10000020]
