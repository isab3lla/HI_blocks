import numpy as np
import HI_hybrid_lib as HIL

#############################################################################

sim = '500_1024_171201'  #14 snapshots available
#total_snap = 14
#z_arr   = np.array([2.8,2.4,2.1,1.8,1.5,1.3,1.1,0.9,0.8,0.7,
#							0.6,0.5,0.4,0.3])
total_snap = 4
z_arr   = np.array([2.8,2.4,2.1,1.8,1.5,1.3,1.1,0.9,0.8,0.7,
							0.6,0.5,0.4,0.3])


#sim = '120_256_180710'
#total_snap = 1
#z_arr   = np.array([0.0])



Omega_HI_all = np.zeros((total_snap))
## looping over all snapshots
#for i in range(total_snap):

j=0
for i in [0,4,8,13]:

	snapshot_fname,groups_fname,groups_number,dumb0 = HIL.getting_sim_paths(sim,i)
	OmegaHI,dumb1,dumb2 = HIL.illustris_HI_assignment(snapshot_fname,groups_fname,groups_number,64)

	Omega_HI_all[j] = OmegaHI; j+=1
	del OmegaHI,dumb1,dumb2


#write total OmegHI file
f_out = 'Omega_HI_'+sim+'_64pHalo.dat'
f=open(f_out,'w')
#for i in range(len(Omega_HI_all)):
j=0
for i in [0,4,8,13]:
    f.write(str(z_arr[i])+' '+str(Omega_HI_all[j])+'\n')
    j+=1
f.close()
