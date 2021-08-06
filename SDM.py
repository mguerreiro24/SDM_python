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
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import binom


#---------------------------------------------------------------------------
def lm_(p1,p2,v):
    """linear regression"""
    m =float(p2[0]-p1[0])/float(p2[1]-p1[1])
    b =float(p2[0]-p2[1]*m)
    return m*v+b
    

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


def standardize_features2(X_train):
    """standardizes data given and scaling object
requires:
    X_train   set used to build scaling object
ensures:
    scaler    object used to scale other sets using this fit
    X_train   original set used to build scaler object, but scaled

    E.g.: X_test = scaler.transform(X_test)
"""
    scaler = StandardScaler()
    scaler.fit(X_train)
    X_train = scaler.transform(X_train)
    return scaler,X_train


def set_Threshold(fpr,tpr,threshold,v):
    mm = tpr-fpr
    t = threshold[np.where(mm==mm.max())[0].tolist()[0]]
    out = np.zeros(v.shape)
    mask = v>=t
    out[mask] = 1.0
    return out


def background(data):
    """pseudo-absences models"""
    background_points = np.c_[np.random.randint(low=0, high=data.Nz,
                                                size=10000),
                              np.random.randint(low=0, high=data.Ny,
                                                size=10000),
                              np.random.randint(low=0, high=data.Nx,
                                                size=10000)].T
    return background_points


def Load_compressed_Data():
    # Load the compressed data
    data = Load_cephalopods_macaronesia()
    data.coverages[data.coverages<0] = -9999
    data.coverages[np.isnan(data.coverages)] = -9999
    data.coverages[np.isinf(data.coverages)] = -9999
    fmax = np.finfo(np.float64).max
    data.coverages[data.coverages>fmax] = fmax
    return data


def data_loading_Models():
    # Load the compressed data
    data = Load_compressed_Data()

    # Set up the data grid
    xgrid = data.xgrid
    ygrid = data.ygrid
    zgrid = data.zgrid

    # water points
    land_reference = data.mask
    
    return data, xgrid, ygrid, zgrid, land_reference


def get_bunches(spec):
    # create a bunch for each species
    data_avail = False
    try:
        bunchs = []
        for i,sp in enumerate(spec):
            with open('{}_.p'.format(sp),'rb') as file:
                b = pickle.load(file)
            bunchs.append(b)
            tst = b.pts_test
            tr = b.pts_train
            if i==0:
                test = tst
                train = tr
            else:
                test = np.r_[test,tst]
                train = np.r_[train,tr]
        data_avail = True
    except IOError:
        test,train = train_test()
        bunchs = create_community_bunch(spec,train, test)
    return data_avail, bunchs, train, test


def get_Xtra_bunches(spec):
    # create a bunch for each species[+presence absence data]
    data_avail = False
    try:
        bunchs = []
        for i,sp in enumerate(spec):
            with open('{}__.p'.format(sp),'rb') as file:
                b = pickle.load(file)
            bunchs.append(b)
            tst = b.pts_test
            tr = b.pts_train
##            tst_l = pts_test_label
##            tr_l = pts_train_label
            if i==0:
                test = tst
                train = tr
##                test_label = tst_l
##                train_label = tr_l
            else:
                test = np.r_[test,tst]
                train = np.r_[train,tr]
##                test_label = np.append(test_label,tst_l)
##                train_label = np.append(train_label,tr_l)
        data_avail = True
    except IOError:
        (test, test_label),(train, train_label) = train_test_PA()
        bunchs = create_community_bunch(spec, train, test)
        labels = dict(test_label=test_label, train_label=train_label)
        for i,sp in enumerate(spec):
            bunchs[i]["test_label"] = labels["test_label"][np.where(test['species']==sp)]
            bunchs[i]["train_label"] = labels["train_label"][np.where(train['species']==sp)]
    return data_avail, bunchs, test, train#, test_label, train_label

#------------------------MaxEnt--------------------------------------------#
def Fit_OneClassSVM(training_set):#<-model maximum entropy
    """ fit MaxEnt models into the train data
requires:
    training_set    array object with standardized observations
ensures:
    svm.OneClassSVM.fit()
"""
    clf = svm.OneClassSVM(nu=0.1,
                          kernel="rbf",
                          gamma=0.5)
    clf.fit(training_set)
    return clf


def calculated_MaxEnt_SDM(spec=["Abralia_redfieldi"]):
    # create a bunch for each species
    data_avail, bunchs, train, test = get_bunches(spec)

    data, xgrid, ygrid, zgrid, land_reference = data_loading_Models()
    # background points (grid coordinates) for evaluation
    background_points = background(data)

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
        Z = set_Threshold(fpr,tpr,thresholds,Z)
##        plt.text(-35, -70, "AUC: %.3f" % roc_auc, ha="right")
        print("\n Area under the ROC curve : %f" % roc_auc)
        Zcoll.append(Z)
        AUCs.append(roc_auc)
    if not data_avail:
        for i,sp in enumerate(spec):
            pickle.dump(bunchs[i], open('{}_.p'.format(sp),'wb'))
    return Zcoll,levels,AUCs,xgrid,ygrid,train,test


#------------------------Neural--------------------------------------------#
def Fit_MLPClassifier(training_set,labels):
    """ fit neural network into the train data
requires:
    training_set    array object with standardized observations
    labels          array with presence-absence (1 and 0 respectively)
ensures:
    MLPClassifier.fit()
"""
    clf = MLPClassifier(hidden_layer_sizes=(100,2),
                        solver='lbfgs',
                        max_iter=30000)#,
                        #random_state=1)
    clf.fit(training_set,labels)
    return clf


def calculate_NN_SDM(spec=["Abralia_redfieldi"]):
    test_t, train_t = train_test_PA()
    test, test_label = test_t
    train, train_label = train_t

    # create a bunch for each species
    bunchs = create_community_bunch(spec, train, test)

    # Load the compressed data
    data = Load_cephalopods_macaronesia()
    data.coverages[data.coverages<0] = -9999
    data.coverages[np.isnan(data.coverages)] = -9999
    data.coverages[np.isinf(data.coverages)] = -9999
    fmax = np.finfo(np.float64).max
    data.coverages[data.coverages>fmax] = fmax

    # Set up the data grid
    xgrid = data.xgrid
    ygrid = data.ygrid
    zgrid = data.zgrid

    # water points
    land_reference = data.mask

    Zcoll = []
    levels = []
    AUCs = []
    for i, species in enumerate(bunchs):
        print("_" * 80)
        print("Modeling distribution of species '%s'" % species.name)

        # Standardize features
        scaler, train_cover_std = standardize_features2(species.cov_train)
        # Fit MLPClassifier | building neural network model
        print(" - fit MLPClassifier ... ", end='')
##        clf = Fit_MLPClassifier(train_cover_std, train_label[np.where(train['species']==spec[i])])
        clf = Fit_MLPClassifier(train_cover_std, species.train_label)
        print("done.")
        print(" - predict species distribution ... ", end='')

        ## Predict species distribution using the training data | applying model
        Z = np.ones((data.Nz ,data.Ny, data.Nx), dtype=np.float64)

        # We'll predict only for the land points.
        idx = np.where(land_reference==True)
        coverages_land = data.coverages[:, idx[0], idx[1], idx[2]].T#T?
        
        pred = clf.predict_proba(scaler.transform(coverages_land))[:,1]

        Z *= pred.min()#coastline 1
        Z[idx[0], idx[1], idx[2]] = pred

        levels.append(np.linspace(Z.min(), Z.max(), 25))
        Z[land_reference==False] = -9999#Land

        print("done.")
        # Compute AUC with regards to background points | testing AUC+ROC
        pred_test = clf.predict_proba(scaler.transform(species.cov_test))[:,1]
##        fpr, tpr, thresholds = metrics.roc_curve(test_label[np.where(test['species']==spec[i])], pred_test)
        fpr, tpr, thresholds = metrics.roc_curve(species.test_label, pred_test)
        roc_auc = metrics.auc(fpr, tpr)
##        NN_score = clf.score(scaler.transform(species.cov_test), test_label[np.where(test['species']==spec[i])])
        NN_score = clf.score(scaler.transform(species.cov_test), species.test_label)
        Z = set_Threshold(fpr,tpr,thresholds,Z)
        print("\n Area under the ROC curve : %f" % roc_auc)
        print("\n Neural Network score : %f" % NN_score)
        Zcoll.append(Z)
        AUCs.append(roc_auc)
    if not data_avail:
        for i,sp in enumerate(spec):
            pickle.dump(bunchs[i], open('{}__.p'.format(sp),'wb'))
    return Zcoll,levels,AUCs,xgrid,ygrid,train,test


#------------------------Random Forest--------------------------------------------#
def Fit_randomForest(training_set,labels):
    """ fit random forest into the train data
requires:
    training_set    array object with observations
    labels          array with presence-absence (1 and 0 respectively)
ensures:
    RandomForestClassifier.fit()
"""
    clf = RandomForestClassifier(n_estimators=1000)
    clf.fit(training_set,labels)
    return clf


def calculate_RT_SDM(spec=["Abralia_redfieldi"]):
    test_t, train_t = train_test_PA()
    test, test_label = test_t
    train, train_label = train_t

    # create a bunch for each species
    bunchs = create_community_bunch(spec, train, test)

    # Load the compressed data
    data = Load_cephalopods_macaronesia()
    data.coverages[data.coverages<0] = -9999
    data.coverages[np.isnan(data.coverages)] = -9999
    data.coverages[np.isinf(data.coverages)] = -9999
    fmax = np.finfo(np.float64).max
    data.coverages[data.coverages>fmax] = fmax

    # Set up the data grid
    xgrid = data.xgrid
    ygrid = data.ygrid
    zgrid = data.zgrid



    # water points
    land_reference = data.mask

    Zcoll = []
    AUCs = []
    for i, species in enumerate(bunchs):
        print("_" * 80)
        print("Modeling distribution of species '%s'" % species.name)

        # Standardize features
        scaler, train_cover_std = standardize_features2(species.cov_train)
        # Fit Classifier | 
        print(" - fit Random Forest Classifier ... ", end='')
##        clf = Fit_randomForest(train_cover_std, train_label[np.where(train['species']==spec[i])])
        clf = Fit_randomForest(train_cover_std, species.train_label)

        print("done.")
        print(" - predict species distribution ... ", end='')

        ## Predict species distribution using the training data | applying model
        Z = np.ones((data.Nz ,data.Ny, data.Nx), dtype=np.float64)

        # We'll predict only for the land points.
        idx = np.where(land_reference==True)
        coverages_land = data.coverages[:, idx[0], idx[1], idx[2]].T#T?
        
        pred = clf.predict_proba(scaler.transform(coverages_land))[:,1]

        Z *= pred.min()#coastline 1
        Z[idx[0], idx[1], idx[2]] = pred


        Z[land_reference==False] = -9999#Land

        print("done.")
        # Compute AUC with regards to background points | testing AUC+ROC
        pred_test = clf.predict_proba(scaler.transform(species.cov_test))[:,1]
##        fpr, tpr, thresholds = metrics.roc_curve(test_label[np.where(test['species']==spec[i])], pred_test)
        fpr, tpr, thresholds = metrics.roc_curve(species.test_label, pred_test)
        roc_auc = metrics.auc(fpr, tpr)
##        NN_score = clf.score(scaler.transform(species.cov_test), test_label[np.where(test['species']==spec[i])])
        NN_score = clf.score(scaler.transform(species.cov_test), species.test_label)
        Z = set_Threshold(fpr,tpr,thresholds,Z)
        print("\n Area under the ROC curve : %f" % roc_auc)
        print("\n Random Forest score : %f" % NN_score)
        Zcoll.append(Z)
        AUCs.append(roc_auc)
    if not data_avail:
        for i,sp in enumerate(spec):
            pickle.dump(bunchs[i], open('{}__.p'.format(sp),'wb'))
    return Zcoll,AUCs,xgrid,ygrid,train,test


#----------------------------BIOCLIM----------------------------------------------#
"""bioclim envelop algorithm"""
class Fit_bioclim():
    def __init__(self,training_set):
        self.maximums = training_set.max(0)
        self.minimums = training_set.min(0)
    def decision_function(self,value):
        """value --> numpy array
requires: value, a numpy array
ensures: np.array with 1.0's for classified within envelope,
         0.0's for outside envelope
"""
        out = np.zeros(value.shape)
        for i,(mi,ma) in enumerate(zip(self.minimums,self.maximums)):
            mask = (value[:,i]<=ma) & (value[:,i]>=mi)
            out[:,i][mask] = 1
        out = out.sum(1)
        mask = out==len(self.minimums)
        return mask**2


def calculate_bioclim_SDM(spec=["Abralia_redfieldi"]):
    # create a bunch for each species
    data_avail, bunchs, train, test = get_bunches(spec)

    data, xgrid, ygrid, zgrid, land_reference = data_loading_Models()

    Zcoll = []
    levels = []
    AUCs = []
    # Fit, predict, and plot for each species.
    for i, species in enumerate(bunchs):
        print("_" * 80)
        print("Modeling distribution of species '%s'" % species.name)
        # Fit model
        print(" - fit BioClim ... ", end='')
        clf = Fit_bioclim(species.cov_train)
        print("done.")
        print(" - predict species distribution ... ", end='')

        ## Predict species distribution using the training data | applying model
        Z = np.ones((data.Nz ,data.Ny, data.Nx), dtype=np.float64)

        # We'll predict only for the land points.
        idx = np.where(land_reference==True)
        coverages_land = data.coverages[:, idx[0], idx[1], idx[2]].T#T?
        
        pred = clf.decision_function(coverages_land)

        Z *= pred.min()#coastline 1
        Z[idx[0], idx[1], idx[2]] = pred


        Z[land_reference==False] = -9999#Land

        print("done.")

        # Compute p-value
        pred_test = clf.decision_function(species.cov_test)
        p = np.count_nonzero(pred)/len(pred)
        k = np.count_nonzero(pred_test)
        n = len(pred_test)
        p_V = 1 - binom.cdf(k=k, n=n, p=p)
        
        print("\n p-value : %f" % p_V)
        Zcoll.append(Z)
        AUCs.append(p_V)
    if not data_avail:
        for i,sp in enumerate(spec):
            pickle.dump(bunchs[i], open('{}_.p'.format(sp),'wb'))
    return Zcoll,AUCs,xgrid,ygrid,train,test,land_reference

#----------------------------AquaMaps----------------------------------------------#
"""AquaMaps envelop algorithm"""

class Fit_AquaMaps():
    def __init__(self,training_set):
        #Quartiles and interquartile distance
        q75, q25 = np.percentile(training_set, [75 ,25])
        iqr = q75 - q25
        
        #core distribution
        self.p10 = np.percentile(training_set,10,axis=0)
        #between these two values == 1
        self.p90 = np.percentile(training_set,90,axis=0)

        #from here, gradient to == 0
        self.maximums = np.minimum(training_set.max(0), q25 - 1.5 * iqr)
        self.minimums = np.maximum(training_set.min(0), q75 + 1.5 * iqr)
    def predict_proba(self,value):
        """value --> numpy array
requires: value, a numpy array
ensures: np.array with probabilities for classified within envelope"""
        out = np.zeros(value.shape)
        #calculate prob for each dimension
        for i,(mi,p10,p90,ma) in enumerate(zip(self.minimums,self.p10,self.p90,self.maximums)):
            mask = (value[:,i]<=p90) & (value[:,i]>=p10)
            out[:,i][mask] = 1.0
            mask = (value[:,i]>p90) & (value[:,i]<=ma)
            out[:,i][mask] = lm_((1.0,p90),(.0,ma),out[:,i][mask])
            mask = (value[:,i]<p10) & (value[:,i]>=mi)
            out[:,i][mask] = lm_((.0,mi),(1.0,p10),out[:,i][mask])
        out = out.prod(1)
        return out


def calculate_AquaMaps_SDM(spec=["Abralia_redfieldi"]):
    # create a bunch for each species
    data_avail, bunchs, train, test = get_bunches(spec)

    data, xgrid, ygrid, zgrid, land_reference = data_loading_Models()

    Zcoll = []
    levels = []
    AUCs = []
    # Fit, predict, and plot for each species.
    for i, species in enumerate(bunchs):
        print("_" * 80)
        print("Modeling distribution of species '%s'" % species.name)
        # Fit model
        print(" - fit AquaMaps ... ", end='')
        clf = Fit_AquaMaps(species.cov_train)
        print("done.")
        print(" - predict species distribution ... ", end='')

        ## Predict species distribution using the training data | applying model
        Z = np.ones((data.Nz ,data.Ny, data.Nx), dtype=np.float64)

        # We'll predict only for the water points.
        idx = np.where(land_reference==True)
        coverages_land = data.coverages[:, idx[0], idx[1], idx[2]].T#T?

        pred = clf.predict_proba(coverages_land)

        Z *= pred.min()#coastline 1
        Z[idx[0], idx[1], idx[2]] = pred


        Z[land_reference==False] = -9999#Land and seabed

        print("done.")
        # Compute p-value
        pred_test = clf.predict_proba(species.cov_test)
        p = np.count_nonzero(pred)/len(pred)
        k = np.count_nonzero(pred_test)
        n = len(pred_test)
        p_V = 1 - binom.cdf(k=k, n=n, p=p)
        
        print("\n p-value : %f" % p_V)
        Zcoll.append(Z)
        AUCs.append(p_V)
    if not data_avail:
        for i,sp in enumerate(spec):
            pickle.dump(bunchs[i], open('{}__.p'.format(sp),'wb'))
    return Zcoll,AUCs,xgrid,ygrid,train,test




if __name__=="__main__":
    spec=["Abraliopsis_atlantica","Abralia_redfieldi"]











    
