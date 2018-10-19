import numpy as np
import readsnap
import sys,os
import struct
import HI_hybrid_lib as HIL
import manipulate_snaps_lib as MSL


rho_crit=2.77536627e11 #h^2 Msun/Mpc^3

################################### INPUT #####################################

sim = '500_1024_171201'  #14 snapshots available
total_snap = 14
z_arr   = np.array([2.8,2.4,2.1,1.8,1.5,1.3,1.1,0.9,0.8,0.7,
                          0.6,0.5,0.4,0.3])

Np_halo = 64

print '\n\n    USING',Np_halo,'particles per halo\n\n'
# sim = '120_256_180710'
# total_snap = 1
# z_arr   = np.array([0.0])


##############################################################################


## looping over all snapshots
for i in range(total_snap):

	snapshot_fname,groups_fname,groups_number,out_fname = HIL.getting_sim_paths(sim,i)

	
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


	print 'boxsize = %.1f Mpc/h, z = %.1f\n'%(BoxSize, redshift)
	print 'DM particle mass = %.2e Msun/h\n'%Masses[1]
	print 'number of particles = %d\n'%Nall[1]
	print 'number of subfiles in the snapshot = %d\n'%filenum

	#compute the values of Omega_CDM and Omega_B
	Omega_cdm=Nall[1]*Masses[1]/BoxSize**3/rho_crit
	Omega_b=Omega_m-Omega_cdm
	print '\nOmega_CDM = %.3f\nOmega_B   = %0.3f\nOmega_M   = %.3f\n'\
	    %(Omega_cdm,Omega_b,Omega_m)

	##########################################################################

	# #read the IDs of all the particles
	# print 'reading the IDs of all the particles . . .'
	# IDs=readsnap.read_block(snapshot_fname,"ID  ",parttype=1)
	# print 'done with the reading!\n'

	print '\n\nNow passing to the HI distribution . . .' 
	## finding the HI-carrying particles
	OmegaHI,IDs_g,M_HI = HIL.CR15_HI_assignment(snapshot_fname,groups_fname,groups_number,Np_halo)

	del OmegaHI,IDs_g
	## the index of M_HI is the corresponding ID-1
	## i.e.  ID = index + 1

	print 'M_HI type: ',type(M_HI[0])

	print 'assigned HI mass to particles:'
	print '%e < M_HI [Msun/h] < %e'%(np.min(M_HI),np.max(M_HI))

	## sanity check
	print '\nsanity check: should be equal:'
	print Nall[1], len(M_HI)

	## setting units for M_HI as 1e10 Msun/h
	M_HI = M_HI/1.0e10
	MSL.new_block(snapshot_fname,out_fname,M_HI)


##############################################################################


