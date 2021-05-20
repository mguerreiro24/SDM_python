#############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++Built by Miguel Fernandes Guerreiro+++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++02/03/2021++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#############################################################################
from data_prep import *
from SDM import *
import numpy as np
import matplotlib.pyplot as plt


def Plot_SDM(Z,A,xgrid,ygrid,species_name,land_ref,ppi=900,points=False,train=[],test=[]):
    if points:
        train = train[train['species'] == species_name]
        test = test[test['species'] == species_name]
    plt.subplot(1,1,1)
    X, Y = np.meshgrid(xgrid, ygrid)
##    ZZ = np.sum(Z,axis=0)
    ZZ = np.max(Z,axis=0)
    ZZ[:10,-10:] = 0
    ZZ[10:20,-20:-10] = 1
##    ZZ = Z[75]
##    l = np.linspace(0, ZZ.max(), 4)
    land_ref = land_ref[0]**2
    
    
    plt.contourf(X, Y, ZZ, levels=1, cmap=plt.cm.Reds)
    if points:
        plt.colorbar(format='%.2f')
        plt.scatter(train['dd long'], train['dd lat'],
                    s=2 ** 2, c='black',
                    marker='^', label='train')
        plt.scatter(test['dd long'], test['dd lat'],
                    s=2 ** 2, c='black',
                    marker='x', label='test')
        plt.legend()
    plt.contourf(X, Y, land_ref, colors=[(1,1,1,0),(0,0,0,1)])
##    plt.legend()
    plt.title(" ".join(species_name.split("_")))
    plt.axis('equal')
    plt.text(-16, 15.5, "AUC: %.3f" % A)#, ha="right")
    plt.savefig("{}_SDM.pdf".format(species_name),dpi = (ppi))
    plt.close()

if __name__=="__main__":
##    species = ['Abraliopsis_atlantica','Cranchia_scabra', 'Bathyteuthis_abyssicola']
    species = ['Abraliopsis_atlantica', 'Octopoteuthis_sicula', 'Helicocranchia_pfefferi', 'Liocranchia_reinhardti', 'Pterygioteuthis_gemmata', 'Bathyteuthis_abyssicola', 'Abralia_veranii', 'Stigmatoteuthis_arcturi', 'Bathothauma_lyromma', 'Abralia_redfieldi', 'Onychoteuthis_banksii', 'Leachia_atlantica', 'Heteroteuthis_dispar', 'Cranchia_scabra', 'Mastigopsis_hjorti']
##    species = ['Stigmatoteuthis_arcturi','Bathothauma_lyromma']#,'Abraliopsis_atlantica','Cranchia_scabra', 'Abralia_redfieldi', 'Abralia_veranii', 'Leachia_atlantica', 'Liocranchia_reinhardti', 'Heteroteuthis_dispar', 'Bathyteuthis_abyssicola']
##    Zs,Ls,As,xgrid,ygrid,train,test = calculated_MaxEnt_SDM(species)
##    Zs,As,xgrid,ygrid,train,test = calculate_RT_SDM(species)

##    Zs,As,xgrid,ygrid,train,test,land_reference = calculate_bioclim_SDM(species)
    import pickle
    with open('saveZBC.p','rb') as file:
        my_data = pickle.load(file)
    Zs = my_data['Zs']
    As = my_data['As']
    xgrid = my_data['x']
    ygrid = my_data['y']
    species = my_data['species']
    land_reference = my_data['mask']
    for i,specie in enumerate(species):
        Plot_SDM(Zs[i],As[i],xgrid,ygrid,specie,land_reference)#,points=True,train=train,test=test)
    my_data = {'species':species, 'Zs':Zs, 'As':As, 'x':xgrid, 'y':ygrid,'mask':land_reference}
##    pickle.dump(my_data, open('saveZRT.p','wb'))
    pickle.dump(my_data, open('saveZBC.p','wb'))



