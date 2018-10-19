import numpy as np
import readsnap
import sys,os


################################### INPUT #####################################

#snapshot_fname='/scratch/ipcarucci/lpicola_trial/output_files/60Mpc512Nall_z3p000'
#snapshot_fname='/scratch/ipcarucci/sims/snapshots/60_512_z3p000'
snapshot_fname='60_256_z3p000'

out_fname='new_snap'

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


## sanity check
print len(IDs), Nall[1]

marker = np.zeros((len(IDs)))

marker[np.where(IDs>5000000)[0]] = 1
marker[np.where(IDs>10000000)[0]] = 2

print 'marker type ',type(marker[3])

## ordering the marker per ID
sorting_IDs_order = IDs.argsort()
sorted_marker = marker[sorting_IDs_order]

print IDs[sorting_IDs_order]

print sorted_marker[4999980:5000020]
print '\nciaoooo\n'
print sorted_marker[9999980:10000020]




#sys.exit()
#for i in range(50):
#    print IDs[i], '  ', marker[i]


print 'and now let s see . . .'

import struct

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

    new_data = sorted_marker[subIDs-1]
    
    ## now we append the new data
    file = open(outputfname,'ab')
    file.write(struct.pack('<I',int(npart*8)))
    for j in range(npart):
        file.write(struct.pack('<d',new_data[j]))    
    file.write(struct.pack('<I',int(npart*8)))
    file.close()
    del npart, head, subIDs, new_data
