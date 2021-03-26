#############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++Built by Miguel Fernandes Guerreiro+++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++11/02/2021++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#############################################################################
#imports
import numpy as np
from sklearn.utils import Bunch#class Bunch(dict) keys==attributes
from sklearn import svm, metrics
from data_prep import *

#---------------------------------------------------------------------------
def build_polygon(lower_left_lon, lower_left_lat, size,data):#not under use
    """
requires:
    lower_left_lon    float()
    lower_left_lat    float()
    size              float()
    data              str()
ensures:
    Polygon entry for QGis (GeoJason style?)
"""
    return "POLYGON (( %f %f , %f %f , %f %f , %f %f ));%s\n" % (lower_left_lon, lower_left_lat,
                                                                 lower_left_lon+size, lower_left_lat,
                                                                 lower_left_lon+size, lower_left_lat+size,
                                                                 lower_left_lon, lower_left_lat+size,
                                                                 data)


def construct_grids(batch):
    """Construct the map grid from the batch object

    Parameters
    ----------
    batch : Bunch object
        
    Returns
    -------
    (xgrid, ygrid, zgrid) : 1-D arrays
        The grid corresponding to the values in batch.coverages
    """
    return (batch.xgrid,batch.ygrid,batch.zgrid)


def standardize_features(bunch):
    """
requires: bunch object
ensures: mean,std,train_cover_std
"""
    mean = bunch.cov_train.mean(axis=0)
    std = bunch.cov_train.std(axis=0)
    train_cover_std = (bunch.cov_train - mean) / std
    return mean, std, train_cover_std


def Fit_OneClassSVM(training_set):#<-model maximum entropy
    """ fit MaxEnt models into the train data
requires:
    training_set    array object with standardized observations
ensures:
    svm.OneClassSVM.fit()
"""
    clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.5)
    clf.fit(training_set)
    return clf


def calculated_MaxEnt_SDM(spec=["Abraliopsis_atlantica"]):
    test,train = train_test()


    # create a bunch for each species
    bunchs = []
    for sp in spec:
        bunchs.append(create_species_bunch(sp,
                                           train, test)
                      )

    # Load the compressed data
    data = Load_cephalopods_macaronesia()
    data.coverages[data.coverages<0] = -9999
    data.coverages[np.isnan(data.coverages)] = -9999
    data.coverages[np.isinf(data.coverages)] = -9999
    fmax = np.finfo(np.float64).max
    data.coverages[data.coverages>fmax] = fmax
    print("data prepped?")
    print(np.count_nonzero(np.isnan(data.coverages)),np.count_nonzero(np.isinf(data.coverages)))

    # Set up the data grid
    xgrid = data.xgrid
    ygrid = data.ygrid
    zgrid = data.zgrid

    # background points (grid coordinates) for evaluation
    np.random.seed(13)#
    background_points = np.c_[np.random.randint(low=0, high=data.Nz,
                                                size=10000),
                              np.random.randint(low=0, high=data.Ny,
                                                size=10000),
                              np.random.randint(low=0, high=data.Nx,
                                                size=10000)].T

    # water points
    land_reference = data.mask

    Zcoll = []
    levels = []
    AUCs = []
    # Fit, predict, and plot for each species.
    for i, species in enumerate(bunchs):
        print("_" * 80)
        print("Modeling distribution of species '%s'" % species.name)

        # Standardize features
        mean, std, train_cover_std = standardize_features(species)
##        print(mean)
        # Fit OneClassSVM | building maximum entropy model
        print(" - fit OneClassSVM ... ", end='')
        clf = Fit_OneClassSVM(train_cover_std)
        print("done.")
        print(" - predict species distribution ... ", end='')

        ## Predict species distribution using the training data | applying model
        Z = np.ones((data.Nz ,data.Ny, data.Nx), dtype=np.float64)

        # We'll predict only for the land points.
        idx = np.where(land_reference==True)
        coverages_land = data.coverages[:, idx[0], idx[1], idx[2]].T#T?
        
        pred = clf.decision_function((coverages_land - mean) / std)

        Z *= pred.min()#coastline 1
        Z[idx[0], idx[1], idx[2]] = pred

        levels.append(np.linspace(Z.min(), Z.max(), 25))
        Z[land_reference==False] = -9999#Land

        print("done.")
        # Compute AUC with regards to background points | testing AUC+ROC
        pred_background = Z[background_points[0], background_points[1], background_points[2]]
        pred_test = clf.decision_function((species.cov_test - mean) / std)
        scores = np.r_[pred_test, pred_background]
        y = np.r_[np.ones(pred_test.shape), np.zeros(pred_background.shape)]
        fpr, tpr, thresholds = metrics.roc_curve(y, scores)
        roc_auc = metrics.auc(fpr, tpr)
        
##        plt.text(-35, -70, "AUC: %.3f" % roc_auc, ha="right")
        print("\n Area under the ROC curve : %f" % roc_auc)
        Zcoll.append(Z)
        AUCs.append(roc_auc)
    return Zcoll,levels,AUCs,xgrid,ygrid,train,test


if __name__=="__main__":
    from data_prep import *
    # Load the compressed data
    Zs,Ls,As = calculated_MaxEnt_SDM(['Pterygioteuthis_gemmata','Abraliopsis_atlantica'])
