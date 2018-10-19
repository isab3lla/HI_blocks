## bundle of functions
## written for 21cm X CMBlens project
## isabella, May 2018


## routines available:

## 1) create a new snapshot with a new data block 
## new_block(snapshot_fname,out_fname,data)
## snapshot_fname --------> original snapshot
## out_fname -------------> name of new snapshot
## data ------------------> data to be appended
## DATA HAS TO BE SORTED IN IDs ORDER!





#########################################################

import numpy as np
import sys, os
import struct    # useful for writing binaries
import readsnap  # gadget format snapshot reader

#########################################################
####################   CONSTANTS  #######################
#########################################################

## gravitational constant ##
G = 4.302e-9 # Mpc (km/s)^2 Msun^-1





#########################################################
####   Gadget-like snapshot: appending a new block    ###
#########################################################
## 1) create a new snapshot with a new data block 
## new_block(snapshot_fname,out_fname,data)
## snapshot_fname --------> original snapshot
## out_fname -------------> name of new snapshot
## data ------------------> data to be appended

def new_block(snapshot_fname,out_fname,data):

	print '\n** Appending a new block to ',snapshot_fname,'\n'
	print '(data given should be sorted in IDs order)'

	print 'reading snapshot properties . . '
	head=readsnap.snapshot_header(snapshot_fname)
	Nall=head.nall
	filenum=head.filenum
	swap=head.swap
	format=head.format

	del head

	print 'number of particles = %d'%Nall[1]
	print 'number of subfiles in the snapshot = %d\n'%filenum


	################### ASSUMPTIONS ###################

	##  DATA HAS TO BE SORTED IN IDs ORDER! ##
	
	# snaps IDs are unsigned double integers
	dt = np.uint64           
	
	# sims w/ only 1 particle type (DM)
	add_offset = np.int32(0)
	
	# the data is an array with len of particles
	# i.e. a value for each particle
	if len(data)!=Nall[1]:
		print '\nDATA DOES NOT HAVE RIGHT LENGTH !!\n'
		sys.exit()

	# the data type is FLOAT
	if not isinstance(data[0],(np.float32,np.float64)):
		print '\nDATA MUST BE FLOAT !!\n'
		sys.exit()

	# the machine is INTEL, i.e. little-endian
	# (the appended data has to have same endianess
	#  as the rest of the snapshot) 
	if swap:
		endns = '>'		# to write w/ big-endian
	else:
		endns = '<'		# to write w/ little-endian

	###################################################	

	# checking size of data to append
	# i.e. 4 or 8 bytes?
	S = data.itemsize
	if S==8 :
		data_type = 'd'
	elif S==4 :
		data_type = 'f'
	else:
		print '\nCHECK DATA TYPE!!\n'
		sys.exit()


	## now looping on the subfiles of the snapshot
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
	    print 'which contains ',npart,' particles'
    
	    ## read the IDs of SUB file
	    offset,blocksize = readsnap.find_block(curfilename,format,swap,"ID  ",4)
	    
	    f = open(curfilename,'rb')
	    f.seek(offset + add_offset*np.dtype(dt).itemsize, os.SEEK_CUR)  
	    subIDs = np.fromfile(f,dtype=dt,count=npart) # read data
	    f.close()  
	    if swap:
	        subIDs.byteswap(True)  

	    ## indexes go as ID+1
	    sub_data = data[subIDs-1]
	    
	    ## now we append the new data
	    block_field = endns+'I'  # 4-bytes integer
	    block_data  = endns+data_type # i.e. either '<d' or '<f'
	    print '.. writing the new snapshot subfile ..\n'
	    file = open(outputfname,'ab')
	    file.write(struct.pack(block_field,int(npart*S)))
	    for j in range(npart):
	        file.write(struct.pack(block_data,sub_data[j]))    
	    file.write(struct.pack(block_field,int(npart*S)))
	    file.close()

	    del npart, head, subIDs, sub_data

	print '\nNew snapshot created!\n'


