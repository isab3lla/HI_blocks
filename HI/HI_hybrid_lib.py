import numpy as np
import sys
import readsnap
import readfof
from scipy import interpolate

############################# routines available ##############################

# illustris_HI_assignment

# CR15_HI_assignment

# getting_sim_paths(sim,isnap)

##############################################################################

################################# UNITS #####################################
rho_crit=2.77536627e11 #h^2 Msun/Mpc^3

yr=3.15576e7  #seconds
km=1e5        #cm
Mpc=3.0856e24 #cm
kpc=3.0856e21 #cm
Msun=1.989e33 #g
Ymass=0.24   #helium mass fraction
mH=1.6726e-24  #proton mass in grams
gamma=5.0/3.0  #ideal gas
kB=1.3806e-26  #gr (km/s)^2 K^{-1}
nu0=1420.0     #21-cm frequency in MHz

pi=np.pi


#############################################################################

#This routine uses the Villaescusa-Navarro et al. 2018 method to assign HI to halos

# snapshot_fname --> name of the N-body snapshot
# groups_fname ----> name of the folder containing the FoF/Subfind halos
# groups_number ---> number of the FoF/Subfind file to read
# long_ids_flag ---> True if particle IDs are 64 bits. False otherwise
# Np_halo ---------> How many particles per halo (min is 32)

#the routine returns the IDs of the particles to whom HI has been assigned
#and its HI masses. Note that the HI masses array has a length equal to the 
#total number of particles in the simulation
#If the positions of the particles to which HI has been assigned is wanted,
#one should first sort the positions and then use the IDs to select them

# interpolating the parameters of Table 1 Villaescusa-Navarro et al. 2018
# i.e. HI-halo mass relation for FoF halos in Illustris

z_ill = np.array([0., 1., 2., 3., 4., 5.])
alpha = np.array([0.24,0.53,0.60,0.76,0.79,0.74])
M0    = np.array([43.,15.,13.,2.9,1.4,1.9])
M0    = M0*1e9  # Msun/h
Mmin  = np.array([200.,60.,36.,6.7,2.1,2.0])
Mmin  = Mmin*1e10  # Msun/h


def alpha_illustris(z):
	alpha_func = interpolate.interp1d(z_ill, alpha,1)
	return alpha_func(z)

def M0_illustris(z):
	M0_func = interpolate.interp1d(z_ill, M0,1)
	return M0_func(z)

def Mmin_illustris(z):
	Mmin_func = interpolate.interp1d(z_ill, Mmin,1)
	return Mmin_func(z)

# actual "soft cut-off" M_HI(M_halo) function
# i.e. eq.13 of Villaescusa-Navarro et al. 2018
def MHI_illustris(M_halo,z):
	M0z = M0_illustris(z)
	Mminz = Mmin_illustris(z)
	alphaz = alpha_illustris(z)
	return M0z*(M_halo/Mminz)**(alphaz)*np.exp(-(Mminz/M_halo)**0.35)

######## ACTUAL FUNCTION ########

def illustris_HI_assignment(snapshot_fname,groups_fname,groups_number,
							Np_halo=100,long_ids_flag=True):
					
	print '\n. . reading the snapshot header . .'
	#read snapshot header and obtain BoxSize, redshift and h
	head=readsnap.snapshot_header(snapshot_fname)
	BoxSize=head.boxsize/1e3 #Mpc/h
	Nall=head.nall
	redshift=head.redshift
	mass_DMparticle = head.massarr[1]*1e10 #Msun/h
	h=head.hubble; del head

	#find the total number of particles in the simulation
	Ntotal=np.sum(Nall,dtype=np.int64)
	print 'Total number of particles in the simulation: %d'%Ntotal
	
	#read FoF halos information
	halos=readfof.FoF_catalog(groups_fname,groups_number,
							  long_ids=long_ids_flag,swap=False)
	pos_FoF=halos.GroupPos/1e3   #Mpc/h
	M_FoF=halos.GroupMass*1e10   #Msun/h
	ID_FoF=halos.GroupIDs-1      #normalize IDs
	Len=halos.GroupLen           #number of particles in the halo
	Offset=halos.GroupOffset     #offset of the halo in the ID array
	del halos

	#some verbose
	print 'Number of FoF halos:',len(pos_FoF),len(M_FoF)
	print '%f < X [Mpc/h] < %f'%(np.min(pos_FoF[:,0]),np.max(pos_FoF[:,0]))
	print '%f < Y [Mpc/h] < %f'%(np.min(pos_FoF[:,1]),np.max(pos_FoF[:,1]))
	print '%f < Z [Mpc/h] < %f'%(np.min(pos_FoF[:,2]),np.max(pos_FoF[:,2]))
	print '%e < M [Msun/h] < %e\n'%(np.min(M_FoF),np.max(M_FoF))

	del pos_FoF

	#compute the total mass in halos from the halo catalogue
	print 'Total contributing mass from the catalogue = %e Msun/h'\
		%(np.sum(M_FoF,dtype=np.float64))

	print '. . considering halos with at least '+str(Np_halo)+' particles . . '
	# cutting catalogue to halos with at least 100 particles
	mass_100p = mass_DMparticle*Np_halo
	indeces = np.where(M_FoF>mass_100p)[0]
	M_FoF = M_FoF[indeces]; del indeces

	#compute the total mass in halos from the halo catalogue
	print 'Total contributing mass from the catalogue = %e Msun/h'\
		%(np.sum(M_FoF,dtype=np.float64))
	print 'Number of FoF halos:',len(M_FoF),'\n'
	print '%e < M [Msun/h] < %e'%(np.min(M_FoF),np.max(M_FoF))
	print str(Np_halo)+'* the DM particle mass = %e [Msun/h]\n'%(mass_DMparticle*Np_halo)


	# compute the values of M_HI(M_halo) function
	# see eq.13 of Villaescusa-Navarro et al. 2018
	print '\nParameters for M_HI(M_halo) function:'
	print '  M_0   (z=%2.1f) = %2.2e Msun/h'%(redshift,M0_illustris(redshift))
	print '  M_min (z=%2.1f) = %2.2e Msun/h'%(redshift,Mmin_illustris(redshift))
	print '  alpha (z=%2.1f) = %1.2f\n'%(redshift,alpha_illustris(redshift))


	#define the IDs array
	if long_ids_flag:
		IDs=np.empty(len(ID_FoF),dtype=np.uint64)
	else:
		IDs=np.empty(len(ID_FoF),dtype=np.uint32)

	#loop over the halos containing HI and assign the HI to the particles
	M_HI=np.zeros(Ntotal,dtype=np.float32); No_gas_halos=0; IDs_offset=0
	for index in range(len(M_FoF)):

		#select the IDs of all particles belonging to the halo
		indexes=ID_FoF[Offset[index]:Offset[index]+Len[index]]

		#fill the IDs array
		IDs[IDs_offset:IDs_offset+len(indexes)]=indexes
		IDs_offset+=len(indexes)

		#compute the total HI mass within the dark matter halo
		M_HI_halo = MHI_illustris(M_FoF[index],redshift)

		#if there are gas particles assign the HI to them
		Num_gas=len(indexes)
		if Num_gas>0:
			M_HI[indexes]+=(M_HI_halo*1.0/Num_gas)
		else:
			No_gas_halos+=1

	print '\nNumber of halos with no gas particles=',No_gas_halos

	#just keep the IDs of the particles to which HI has been assigned
	IDs=IDs[0:IDs_offset]

	#compute the value of OmegaHI
	OmegaHI = np.sum(M_HI,dtype=np.float64)/BoxSize**3/rho_crit
	print '\nOmega_HI (halos) = %e\n'%(OmegaHI)

	return OmegaHI,IDs,M_HI

#############################################################################


# This routine uses the Villaescusa-Navarro et al. 2018 method to assign HI to halos
# BUT imposing an overall normalization to recover Omega_HI
# as fitted by Crighton et al 2015

# snapshot_fname --> name of the N-body snapshot
# groups_fname ----> name of the folder containing the FoF/Subfind halos
# groups_number ---> number of the FoF/Subfind file to read
# long_ids_flag ---> True if particle IDs are 64 bits. False otherwise
# Np_halo ---------> How many particles per halo (min is 32)

#the routine returns the IDs of the particles to whom HI has been assigned
#and its HI masses. Note that the HI masses array has a length equal to the 
#total number of particles in the simulation
#If the positions of the particles to which HI has been assigned is wanted,
#one should first sort the positions and then use the IDs to select them

## Oemga_HI fit as in Crighton at al 2015
def Omega_CR15(z):
	A = 4.0;   gamma = 0.6
	return A*(1.+z)**gamma *1e-4

## an NON-NORMALIZED HI-Mhalo relation
def NN_MHI_illustris(M_halo,z):
	Mminz = Mmin_illustris(z)
	alphaz = alpha_illustris(z)
	return (M_halo/Mminz)**(alphaz)*np.exp(-(Mminz/M_halo)**0.35)

######## ACTUAL FUNCTION ########

def CR15_HI_assignment(snapshot_fname,groups_fname,groups_number,
							Np_halo=100,long_ids_flag=True):
					
	print '\n. . reading the snapshot header . .'
	#read snapshot header and obtain BoxSize, redshift and h
	head=readsnap.snapshot_header(snapshot_fname)
	BoxSize=head.boxsize/1e3 #Mpc/h
	Nall=head.nall
	redshift=head.redshift
	mass_DMparticle = head.massarr[1]*1e10 #Msun/h
	h=head.hubble; del head

	#find the total number of particles in the simulation
	Ntotal=np.sum(Nall,dtype=np.int64)
	print 'Total number of particles in the simulation: %d'%Ntotal
	
	#read FoF halos information
	halos=readfof.FoF_catalog(groups_fname,groups_number,
							  long_ids=long_ids_flag,swap=False)
	pos_FoF=halos.GroupPos/1e3   #Mpc/h
	M_FoF=halos.GroupMass*1e10   #Msun/h
	ID_FoF=halos.GroupIDs-1      #normalize IDs
	Len=halos.GroupLen           #number of particles in the halo
	Offset=halos.GroupOffset     #offset of the halo in the ID array
	del halos

	#some verbose
	print 'Number of FoF halos:',len(pos_FoF),len(M_FoF)
	print '%f < X [Mpc/h] < %f'%(np.min(pos_FoF[:,0]),np.max(pos_FoF[:,0]))
	print '%f < Y [Mpc/h] < %f'%(np.min(pos_FoF[:,1]),np.max(pos_FoF[:,1]))
	print '%f < Z [Mpc/h] < %f'%(np.min(pos_FoF[:,2]),np.max(pos_FoF[:,2]))
	print '%e < M [Msun/h] < %e\n'%(np.min(M_FoF),np.max(M_FoF))

	del pos_FoF

	#compute the total mass in halos from the halo catalogue
	print 'Total contributing mass from the catalogue = %e Msun/h'\
		%(np.sum(M_FoF,dtype=np.float64))

	print '. . considering halos with at least '+str(Np_halo)+' particles . . '
	# cutting catalogue to halos with at least 100 particles
	mass_100p = mass_DMparticle*Np_halo
	indeces = np.where(M_FoF>mass_100p)[0]
	M_FoF = M_FoF[indeces]; del indeces

	#compute the total mass in halos from the halo catalogue
	print 'Total contributing mass from the catalogue = %e Msun/h'\
		%(np.sum(M_FoF,dtype=np.float64))
	print 'Number of FoF halos:',len(M_FoF),'\n'
	print '%e < M [Msun/h] < %e'%(np.min(M_FoF),np.max(M_FoF))
	print str(Np_halo)+'* the DM particle mass = %e [Msun/h]\n'%(mass_DMparticle*Np_halo)

	# finding the normalization M0 
	# to reproduce Crighton15 Omega_HI values
	M0 = \
		(Omega_CR15(redshift)*rho_crit*BoxSize**3)/np.sum(NN_MHI_illustris(M_FoF,redshift),dtype=np.float64)

	# compute the values of M_HI(M_halo) function
	# see eq.13 of Villaescusa-Navarro et al. 2018
	print '\nParameters for M_HI(M_halo) function:'
	print '  M_min (z=%2.1f)   = %2.2e Msun/h'%(redshift,Mmin_illustris(redshift))
	print '  alpha (z=%2.1f)   = %1.2f'%(redshift,alpha_illustris(redshift))
	print '  normalization M_0 = %2.2e Msun/h'%M0
	print '        (instead of = %2.2e Msun/h )\n'%(M0_illustris(redshift))

	#define the IDs array
	if long_ids_flag:
		IDs=np.empty(len(ID_FoF),dtype=np.uint64)
	else:
		IDs=np.empty(len(ID_FoF),dtype=np.uint32)

	#loop over the halos containing HI and assign the HI to the particles
	M_HI=np.zeros(Ntotal,dtype=np.float32); No_gas_halos=0; IDs_offset=0
	for index in range(len(M_FoF)):

		#select the IDs of all particles belonging to the halo
		indexes=ID_FoF[Offset[index]:Offset[index]+Len[index]]

		#fill the IDs array
		IDs[IDs_offset:IDs_offset+len(indexes)]=indexes
		IDs_offset+=len(indexes)

		#compute the total HI mass within the dark matter halo
		M_HI_halo = M0*NN_MHI_illustris(M_FoF[index],redshift)

		#if there are gas particles assign the HI to them
		Num_gas=len(indexes)
		if Num_gas>0:
			M_HI[indexes]+=(M_HI_halo*1.0/Num_gas)
		else:
			No_gas_halos+=1

	print '\nNumber of halos with no gas particles=',No_gas_halos

	#just keep the IDs of the particles to which HI has been assigned
	IDs=IDs[0:IDs_offset]

	#compute the value of OmegaHI
	OmegaHI = np.sum(M_HI,dtype=np.float64)/BoxSize**3/rho_crit
	print '\nOmega_HI (halos)      = %e'%(OmegaHI)
	print 'Omega_HI (Crighton15) = %e\n'%(Omega_CR15(redshift))

	return OmegaHI,IDs,M_HI

#############################################################################

# This routine is for getting the snapshot path right
# Plus it works as a repository of the available sims

## info on the simulations available
parent_dic = '/scratch/ipcarucci/giulio/sims/'
parent_dir_HI = '/scratch/ipcarucci/giulio/sims/'
sim_avail = ['500_1024_171201','120_256_180710']

def getting_sim_paths(sim,isnap):
	## getting the right sim

	if sim=='500_1024_171201':
		
		z_arr   = np.array([2.8,2.4,2.1,1.8,1.5,1.3,1.1,0.9,0.8,0.7,
							0.6,0.5,0.4,0.3])

		snap_ar = np.array(['12','13','14','15','16','17','18','19','20','21',
							'22','23','24','25'])

		z_str   = ['2p800','2p400','2p100','1p800','1p500','1p300','1p100',
					'0p900','0p800','0p700','0p600','0p500','0p400','0p300']

		if isnap > len(z_arr):
			print 'specify an existing snapshot!'
			print '(isnap not valid)'
			sys.exit()

		print '\nSimulation: ',sim,' --> ',len(z_arr), 'redshifts available'
		print 'chosen isnap =',isnap,'  i.e. z =',z_arr[isnap]

	elif sim=='120_256_180710':

		z_arr   = np.array([0.0])

		snap_ar = np.array(['00'])

		z_str   = ['0p000']

		if isnap > len(z_arr):
			print 'specify an existing snapshot!'
			print '(isnap not valid)'
			sys.exit()

		print 'Simulation: ',sim,' --> ',len(z_arr), 'redshifts available'
		print 'chosen isnap =',isnap,'  i.e. z =',z_arr[isnap]


	else:
		print 'specify an existing simulation!!'
		print 'available ones: ',sim_avail
		sys.exit()


	# actual definition of the paths

	snap = snap_ar[isnap]

	groups_fname = parent_dic+sim+'/snapshots'
	groups_number = int(snap)

	final_piece = '/snapdir_0'+snap+'/'+sim[:-6]+'z'+z_str[isnap]+'_0'+snap
	
	snapshot_fname=groups_fname+final_piece
	snapshot_fname_HI=parent_dic+sim+'/snapshots_w_HI'+final_piece

	return [snapshot_fname,groups_fname,groups_number,snapshot_fname_HI]


