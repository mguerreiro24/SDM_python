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


def Plot_SDM(Z,l,A,xgrid,ygrid,species_name,data,ppi=900):
    train = data.train[data.train['species'] == species_name]
    test = data.test[data.test['species'] == species_name]
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
    species = ['Pterygioteuthis_gemmata','Abraliopsis_atlantica']
    d = Load_cephalopods_macaronesia()
    Zs,Ls,As = calculated_MaxEnt_SDM(d,species)
    for i,specie in enumerate(species):
        Plot_SDM(Zs[i],Ls[i],As[i],d.xgrid,d.ygrid,specie,d)
