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
import os

def Plot_SDM(Z,A,xgrid,ygrid,species_name,land,ppi=900,points=False,train=[],test=[]):
    if points:
        train = train[train['species'] == species_name]
        test = test[test['species'] == species_name]
    plt.subplot(1,1,1)
    X, Y = np.meshgrid(xgrid, ygrid)
    ZZ = np.max(Z,axis=0)

    l = np.linspace(0, 1, 25)

    plt.contourf(X, Y, ZZ, levels=l, cmap=plt.cm.Reds)
    plt.colorbar(format='%.2f')
    if points:
        plt.scatter(train['dd long'], train['dd lat'],
                    s=2 ** 2, c='black',
                    marker='^', label='train')
        plt.scatter(test['dd long'], test['dd lat'],
                    s=2 ** 2, c='black',
                    marker='x', label='test')
        plt.legend()
    plt.contourf(X, Y, land, colors=[(0,0,0,0),(0,0,0,1)])

    plt.title(" ".join(species_name.split("_")))
    plt.axis('equal')
    plt.text(-16, 15.5, "AUC: %.3f" % A)#, ha="right")
    plt.savefig("{}_SDM.pdf".format(species_name),dpi = (ppi))
    plt.close()


def land(xgrid,ygrid):
    X, Y = np.meshgrid(xgrid, ygrid)
    folder_g='D:\\PhD\\GIS-DBs\\GEBCO'
    filename_g='GEBCO_2019.nc'
    fileg = os.path.join(folder_g,filename_g)
    setD = read_gebco(fileg,"-2.0,50.0","-36,-0.1")
    grid = setD.grid
    xxgrid = setD.xgrid + grid/2
    yygrid = setD.ygrid + grid/2
    ix = np.searchsorted(xxgrid, X)
    iy = np.searchsorted(yygrid, Y)
    land_ref = setD.elevation[iy,ix]
    m = land_ref>-5
    land_ref = np.ones(land_ref.shape)
    land_ref[m] = 0.0
    return land_ref

if __name__=="__main__":
##    species = ['Abraliopsis_atlantica','Cranchia_scabra', 'Bathyteuthis_abyssicola']
    species = ['Abraliopsis_atlantica', 'Octopoteuthis_sicula', 'Helicocranchia_pfefferi', 'Liocranchia_reinhardti', 'Pterygioteuthis_gemmata', 'Bathyteuthis_abyssicola', 'Abralia_veranii', 'Stigmatoteuthis_arcturi', 'Bathothauma_lyromma', 'Abralia_redfieldi', 'Onychoteuthis_banksii', 'Leachia_atlantica', 'Heteroteuthis_dispar', 'Cranchia_scabra', 'Mastigopsis_hjorti', 'Vampyroteuthis_infernalis']
##    species = ['Stigmatoteuthis_arcturi','Bathothauma_lyromma']#,'Abraliopsis_atlantica','Cranchia_scabra', 'Abralia_redfieldi', 'Abralia_veranii', 'Leachia_atlantica', 'Liocranchia_reinhardti', 'Heteroteuthis_dispar', 'Bathyteuthis_abyssicola']
##    Zs,Ls,As,xgrid,ygrid,train,test = calculated_MaxEnt_SDM(species)
##    Zs,As,xgrid,ygrid,train,test = calculate_RT_SDM(species)

##    Zs,As,xgrid,ygrid,train,test,land_reference = calculate_bioclim_SDM(species)
    Zs,As,xgrid,ygrid,train,test = calculate_AquaMaps_SDM(species)
    import pickle
##    with open('saveAquamaps.p','rb') as file:
##        my_data = pickle.load(file)
##    Zs = my_data['Zs']
##    As = my_data['As']
##    xgrid = my_data['x']
##    ygrid = my_data['y']
##    species = my_data['species']
##    land = land(xgrid,ygrid)
    for i,specie in enumerate(species):
        print(specie)
        Plot_SDM(Zs[i],As[i],xgrid,ygrid,specie,land)#,points=True,train=train,test=test)
    my_data = {'species':species, 'Zs':Zs, 'As':As, 'x':xgrid, 'y':ygrid}
##    pickle.dump(my_data, open('saveZRT.p','wb'))
    pickle.dump(my_data, open('saveAquamaps.p','wb'))



