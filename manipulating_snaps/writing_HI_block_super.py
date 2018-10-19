import numpy as np
import readsnap
import sys,os
import struct
import HI_illustris as HIill


rho_crit=2.77536627e11 #h^2 Msun/Mpc^3

################################### INPUT #####################################

#sim = '500_1024_171201'  #14 snapshots available
#total_snap = 14
#z_arr   = np.array([2.8,2.4,2.1,1.8,1.5,1.3,1.1,0.9,0.8,0.7,
#                           0.6,0.5,0.4,0.3])

sim = '120_256_180710'
total_snap = 1
z_arr   = np.array([0.0])



out_fname='new_folder/super_new_snap'


##############################################################################


## looping over all snapshots
for i in range(total_snap):

	snapshot_fname,groups_fname,groups_number = HIill.getting_sim_paths(sim,i)


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
	OmegaHI,IDs_g,M_HI = HIill.illustris_HI_assignment(snapshot_fname,
												groups_fname,groups_number)

	del OmegaHI,IDs_g
	## the index of M_HI is the corresponding ID-1
	## i.e.  ID = index + 1

	print 'M_HI type: ',type(M_HI[0])

	print 'assigned HI mass to particles:'
	print '%e < M_HI [Msun/h] < %e'%(np.min(M_HI),np.max(M_HI))

	## sanity check
	print '\nsanity check: should be equal:'
	print Nall[1], len(M_HI)

	# ## ordering the HI mass per ID
	# sorting_IDs_order = IDs.argsort()
	# sorted_M_HI = M_HI[sorting_IDs_order]
	## we are ready:
	## 'sorted_M_HI' is the new block!
	print '\n\n###########################################'
	print 'Now we can add the new block to the snapshot.'
	print '###########################################\n'

	for i in range(filenum):

		curfilename = snapshot_fname+'.'+str(i)
		outputfname = out_fname+'.'+str(i)


		print 'modifyng ',curfilename

		## this is just to copy the snapshot subfile
		with open(outputfname, "wb") as newfile, open(curfilename, "rb") as snap:
		    newfile.write(snap.read())
		    newfile.close()

		head=readsnap.snapshot_header(curfilename)
		npart=head.npart; npart = npart[1]
		print 'which contains ',npart,' particles\n'

		swap=head.swap; format=head.format
		## read the IDs of SUB file
		offset,blocksize = readsnap.find_block(curfilename,format,swap,"ID  ",4)
		dt = np.uint64 # IDs are double integers
		add_offset = np.int32(0) # only one particle species (DM) in the sim

		f = open(curfilename,'rb')
		f.seek(offset + add_offset*np.dtype(dt).itemsize, os.SEEK_CUR)  
		subIDs = np.fromfile(f,dtype=dt,count=npart) # read data
		f.close()  
		if swap:
		    subIDs.byteswap(True)  

		new_data = M_HI[subIDs-1]

		## now we append the new data
		file = open(outputfname,'ab')
		file.write(struct.pack('<I',int(npart*4)))
		for j in range(npart):
		    file.write(struct.pack('<f',new_data[j]))    
		file.write(struct.pack('<I',int(npart*4)))
		file.close()
		del npart, head, subIDs, new_data



##############################################################################


