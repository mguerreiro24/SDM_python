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


def Plot_SDM(Z,l,A,xgrid,ygrid,species_name,train,test,ppi=900):
    train = train[train['species'] == species_name]
    test = test[test['species'] == species_name]
    plt.subplot(1,1,1)
    X, Y = np.meshgrid(xgrid, ygrid)
    ZZ = np.max(Z,axis=0)
    plt.contourf(X, Y, ZZ, levels=l, cmap=plt.cm.Reds)
    plt.colorbar(format='%.2f')
    plt.scatter(train['dd long'], train['dd lat'],
                s=2 ** 2, c='black',
                marker='^', label='train')
    plt.scatter(test['dd long'], test['dd lat'],
                s=2 ** 2, c='black',
                marker='x', label='test')
    plt.legend()
    plt.title(" ".join(species_name.split("_")))
    plt.axis('equal')
    plt.text(-16, 15.5, "AUC: %.3f" % A)#, ha="right")
    plt.savefig("{}_SDM.pdf".format(species_name),dpi = (ppi))
    plt.close()

if __name__=="__main__":
    species = ['Pterygioteuthis_gemmata','Abraliopsis_atlantica','Cranchia_scabra', 'Abralia_redfieldi', 'Leachia_atlantica', 'Liocranchia_reinhardti', 'Heteroteuthis_dispar', 'Bathyteuthis_abyssicola', 'Helicocranchia_pfefferi', 'Helicocranchia_pfefferi', 'Mastigopsis_hjorti']
    Zs,Ls,As,xgrid,ygrid,train,test = calculated_MaxEnt_SDM(species)
    for i,specie in enumerate(species):
        Plot_SDM(Zs[i],Ls[i],As[i],xgrid,ygrid,specie,train,test)
