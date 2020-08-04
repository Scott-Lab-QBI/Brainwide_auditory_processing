# -*- coding: utf-8 -*-
"""
Created on Feb 20 2020
Las modified Aug 04 2020
@author: leandroscholz
"""

def compute_MDS(mat_fpath, n_components=2):
    """
    compute_MDS receives the filepath to a .mat file 
    with multiple matrices containing x,y,z coordinates and a label (located in the last column of the matrix)
    and computes the 2-component Multi Dimensional Scaling (non-linear dimensionality reduction) 
    for the xyz coordinates
    """
    import numpy as np
    from scipy import io as sio 
    from sklearn import manifold
    
    mat = sio.loadmat(mat_fpath)
    metadata = ['__header__', '__version__', '__globals__']
    MDS_dict = {}
    
    #print(mat)
    for key, value in mat.items():
        if key not in metadata:
            print('Computing MDS for', key)
            X = value[:,:-1]
            mds = manifold.MDS(n_components, max_iter=100, n_init=1)
            MDS = mds.fit_transform(X)
            MDS = np.append(MDS,np.reshape(value[:,-1],(-1,1)), axis=1)
            MDS_dict.update({key+'_MDS': MDS})

    sio.savemat(mat_fpath[:-4]+'_MDS.mat', MDS_dict)
    print('Finished the computation of Multidimensional Scaling components.')
    return MDS_dict

def compute_KDE(mat_fpath, col=0, binned=False):
    """
    compute_KDE receives the file path to a .mat file
    with multiple matrices containing x,y coordinates and a label col 
    and computes the Kernel Density Estimates for each subset of 
    coordinates pertaining to each label value
    """
    from scipy import io as sio 
    from sklearn.neighbors import KernelDensity
    from sklearn.model_selection import GridSearchCV, LeaveOneOut
    import numpy as np 
    
    mat = sio.loadmat(mat_fpath)
    metadata = ['__header__', '__version__', '__globals__']
    
    density_dict = {}
    bandwidth_dict = {}
    
    for key, value in mat.items():
        if key not in metadata:
            labels = np.unique(value[:,-1])
            bandwidths = estimate_bandwidth(value[:,col], n_vals=500, order_magnitude=3)
            
            print('computing KDE for', key)
            for label in labels:
                
                #print('indices of observations with current label are', value[:,-1]==label)
                data = np.reshape(value[value[:,-1]==label,col],(-1,1))

                if len(data) > 1:
                    #print('data shape', data.shape)
                    grid = GridSearchCV(KernelDensity(kernel='gaussian'),
                                        {'bandwidth': bandwidths},
                                        cv=LeaveOneOut(), n_jobs=-1)
                    grid.fit(data)
                    print("Freq {0} best bandwidth: {1}".format(label, grid.best_params_['bandwidth']))
                    bandwidth_dict.update({key : grid.best_params_['bandwidth']})
                    kde = grid.best_estimator_ 
                else:
                    print('There is only one data point, using bandwidth = 1 for KDE computation')
                    kde = KernelDensity(bandwidth=1, kernel='gaussian')
                    kde.fit(data)
                    
                # construct a gaussian kernel density estimate of the distribution
                #kde = KernelDensity(bandwidth=10, kernel='gaussian')
                #kde.fit(data).
                if binned:
                    global_min = np.min(value[:,col])
                    global_max = np.max(value[:,col])
                    global_range = np.linspace((global_min - 0.1*(global_min)), (global_max + 0.1*(global_max)), 2000)
                    global_range = np.reshape(global_range,(-1,1))
                    #print(global_range)
                    #print('data shape: ', global_range.shape)
                    density_temp = np.reshape(np.exp(kde.score_samples(global_range)),(-1,1))
                    label_temp = np.reshape(np.repeat(label,len(global_range)),(-1,1))
                    temp = np.append(global_range, density_temp, axis=1)
                    temp = np.append(temp, label_temp, axis=1)
                else: 
                    #print(data)
                    #print('data shape: ', data.shape)
                    density_temp = np.reshape(np.exp(kde.score_samples(data)),(-1,1))
                    label_temp = np.reshape(np.repeat(label,len(data)),(-1,1))
                    temp = np.append(density_temp, label_temp, axis=1)
                
                if 'density' not in locals():
                    #print('density shape', temp.shape)
                    density = temp
                    print('density shape', density.shape)
                else:
                    density = np.append(density, temp, axis=0)
                    print('density shape', density.shape)
                
            density_dict.update({key+'_density' : density})
            del density
    if binned:        
        sio.savemat(mat_fpath[:-4]+'_global_range_density.mat', density_dict)
    else:
        sio.savemat(mat_fpath[:-4]+'_density.mat', density_dict)
    print('Finished the computation of Kernel Density Estimates!')
    return density_dict
    
def compute_KDE_2D(mat_fpath):
    """
    compute_KDE_2D receives the file path to a .mat file
    with multiple matrices containing x,y coordinates and a label 
    and computes the Kernel Density Estimates for each subset of 
    coordinates pertaining to each label value
    """
    from scipy import io as sio 
    from sklearn.neighbors import KernelDensity
    from sklearn.model_selection import GridSearchCV, LeaveOneOut
    import numpy as np 
    
    mat = sio.loadmat(mat_fpath)
    metadata = ['__header__', '__version__', '__globals__']
    
    density_dict = {}
    bandwidth_dict = {}
    
    for key, value in mat.items():
        if key not in metadata:
            labels = np.unique(value[:,-1])
            
            print('computing KDE for', key)
            for i, label in enumerate(labels):

                bandwidths = estimate_bandwidth(value[:,:-1], n_vals=300, order_magnitude=3)
                
                #print('indices of observations with current label are', value[:,-1]==label)
                row, col = value.shape
                # temporary hardcode to compute matrices with 3 and 4 columns
                if col == 3:
                    data = value[value[:,-1]==label,:-1]
                elif col == 4:
                    data = value[value[:,-1]==label,:-2]
                    
                if len(data) > 1:
                    if bandwidth_dict not in locals():
                        #print('data shape', data.shape)
                        grid = GridSearchCV(KernelDensity(kernel='gaussian'),
                                            {'bandwidth': bandwidths},
                                            cv=LeaveOneOut(), n_jobs=-1)
                        grid.fit(data)
                        print("Freq {0} best bandwidth: {1}".format(label, grid.best_params_['bandwidth']))
                        bandwidth_dict.update({key : grid.best_params_['bandwidth']})
                        kde = grid.best_estimator_ 
                    else:
                        # construct a gaussian kernel density estimate of the distribution
                        kde = KernelDensity(bandwidth=bandwidth_dict[key][i], kernel='gaussian')
                        kde.fit(data)
                else:
                    print('There is only one data point, using bandwidth = 1 for KDE computation')
                    kde = KernelDensity(bandwidth=1, kernel='gaussian')
                    kde.fit(data)
                    
                # construct a gaussian kernel density estimate of the distribution
                #kde = KernelDensity(bandwidth=10, kernel='gaussian')
                #kde.fit(data).
                density_temp = np.reshape(np.exp(kde.score_samples(data)),(-1,1))
                label_temp = np.reshape(np.repeat(label,len(data)),(-1,1))
                
                temp = np.append(density_temp, label_temp, axis=1)
                
                if 'density' not in locals():
                    #print('density shape', temp.shape)
                    density = temp
                    print('density shape', density.shape)
                else:
                    density = np.append(density, temp, axis=0)
                    print('density shape', density.shape)
                
            density_dict.update({key+'_density' : density})
            del density
            
    sio.savemat(mat_fpath[:-4]+'_density.mat', density_dict)
    print('Finished the computation of Kernel Density Estimates!')
    return density_dict
    
def compute_ks2d2s(X, labels=None, label_column=-1):
    """
    compute_ks2d2s calculates two-dimensional Kolmogorov-Smirnov test on
    two samples. Computation is done pair-wise using all combinations 
    of labels. 

    input
        X -> either m-by-3 or m-by-2 (if labels are given in X) matrix where 
        labels -> m-by-1 np.array. Default: None (verifies label_column) 
        label_column -> integer 0, 1 or 2. defines which column to consider as labels. Default: -1 (last column)
    """
    import itertools
    from ndtest.ndtest import ks2d2s

    p_values = {}
    
    if any(labels) is None:
        labels = np.unique(X[:,label_column])

        for i, j in itertools.combinations(labels,2):
            print('Testing 2D KS between', i,' and ', j)
            p = ks2d2s(X[labels == i, 0], \
                       X[labels == i, 1], \
                       X[labels == j, 0], \
                       X[labels == j, 1])
            print('p-value', p)
            p_values.update({str(i)+':'+str(j) : p})
    else:
        unique_labels = np.unique(labels)
        #print('unique labels', unique_labels)
        for i, j in itertools.combinations(unique_labels,2):
            print('Testing 2D KS between', i,' and ', j)
            p = ks2d2s(X[labels == i, 0], \
                       X[labels == i, 1], \
                       X[labels == j, 0], \
                       X[labels == j, 1])
            print('p-value', p)
            p_values.update({str(i)+':'+str(j) : p})
        
    return p_values
    
def p_values_to_matrix(p_values):
    """
    p_values_to_matrix transforms the dictionary output from compute_ks2d2s 
    into a dictionary with key 'p_val_matrix' and value as numpy array  
    
    input:
        p_values : dictionary output from compute_ks2d2s function 
    
    output 
        dictionary with key 'p_val_matrix' and value the numpy matrix with p_values (squareform)
    """
    import numpy as np 
    
    vars = []
    
    for key in p_values.keys():
        pair = key.split(':')
        
        for item in pair:
            if item not in vars:
                vars.append(item)
    
    vars.sort()
    print(vars)
    
    p_val_matrix = np.zeros((len(vars),len(vars)))

    for key, value in p_values.items():
    
        pair = key.split(':')
        row = vars.index(pair[0])
        col = vars.index(pair[1])
        p_val_matrix[row,col] = value
        p_val_matrix[col,row] = value
    
    print('p-values matrix:')
    print(p_val_matrix) 
    
    return {'matrix': p_val_matrix, 'vars' : vars} 
   
def estimate_bandwidth(coords, n_vals=100, order_magnitude = 3):
    """
    estimate_bandwidth receives a n-by-m numpy matrix and computes
    a range of values to be used in tests for the optimal 
    band(bin)-width for the input matrix. This is useful when used 
    in conjunction with sklearn.model_selectionGridSearchCV 
    when calculating densities and computing histograms
    """ 
    import numpy as np 
    
    if len(coords.shape)>1:
        row, col = coords.shape
        delta = np.zeros((col,))
        for c in range(col):
            delta[c] = np.max(coords[:,c]) - np.min(coords[:,c])
            
        max_delta = np.max(delta)/3
        exp = np.log10(max_delta)
    else:
        delta = np.max(coords) - np.min(coords)
        max_delta = delta/3
        exp = np.log10(max_delta)
    
    return 10 ** np.linspace(exp-order_magnitude, exp, n_vals)
  
def plot_KDE_2D(mat_fpath, columns=(0,1), bandwidth_dict=None):
    """
    from a mat file with variables of x,y,z and labels, saves png and svg files
    of contour plots of density values computed with KDE for each label

    inputs: 
    
        mat_fpath : filepath to a mat file (string)
        columns : pair of values to define columns of contour plot axes (tuple) 
        bandwidth_dict : in case you want a specific bandwith with which to 
                         compute KDE, give as a dictionary with keys being the same as variables in the mat file
    outputs: 
        indirectly, it outputs .png and .svg files in the same folder where the mat file is located
        
        bandwith_dict : dictionary with the optimal bandwidths estimates through GridSearchCV
                        if you wish to save xy, xz and yz plots, we recomment to use the same   
                        bandwith for all, to make plots look consistent across axes. 
    """
    from scipy import io as sio 
    from sklearn.neighbors import KernelDensity
    from sklearn.model_selection import GridSearchCV, LeaveOneOut
    import numpy as np 
    import matplotlib.pyplot as plt 
    
    mat = sio.loadmat(mat_fpath)
    metadata = ['__header__', '__version__', '__globals__']
    
    if bandwidth_dict is None:
        print('bandwidth dict does not exist, create one to start.')
        trigger = True
        bandwidth_dict = {}
    else: 
        print('using existing bandwidth dict.')
        trigger = False 
    density_dict = {}
    
    
    for key, value in mat.items():
        if key not in metadata:
            if trigger:
                bandwidth_dict.update({key : []})
            labels = np.unique(value[:,-1])
            # set a range of bandwidth values to test 
            bandwidths = estimate_bandwidth(value[:,:-1], n_vals=300, order_magnitude=3)
            
            col0, col1 = columns

            # create data needed to display contour plot 
            xx = np.linspace(min(value[:,col0]),max(value[:,col0]),1000)
            yy = np.linspace(min(value[:,col1]),max(value[:,col1]),1000)
            X, Y = np.meshgrid(xx, yy) 
            xy = np.vstack([X.ravel(), Y.ravel()]).T
            
            print('computing KDE for', key)
            for i, label in enumerate(labels):
                
                #print('indices of observations with current label are', value[:,-1]==label)
                row, col = value.shape
                # temporary hardcode to compute matrices with 3 and 4 columns
                if col == 3:
                    data = value[:,[col0,col1]]
                    data = data[value[:,-1]==label,:]
                    
                elif col == 4:
                    data = value[:,[col0,col1]]
                    data = data[value[:,-1]==label,:]
                    
                if len(data) > 1:
                    #print('data shape', data.shape)
                    if trigger:
                        #print('data shape', data.shape)
                        grid = GridSearchCV(KernelDensity(kernel='gaussian'),
                                            {'bandwidth': bandwidths},
                                            cv=LeaveOneOut(), n_jobs=-1)
                        grid.fit(data)
                        print("Freq {0} best bandwidth: {1}".format(label, grid.best_params_['bandwidth']))
                        bandwidth_dict[key].append(grid.best_params_['bandwidth'])
                        kde = grid.best_estimator_ 
                    elif not trigger:
                        # construct a gaussian kernel density estimate of the distribution
                        # with the existing optimal bandwidth
                        print('using existing optimal bandwidth of: ', str(bandwidth_dict[key][i]))
                        kde = KernelDensity(bandwidth=bandwidth_dict[key][i], kernel='gaussian')
                        kde.fit(data)
                else:
                    print('There is only one data point, using bandwidth = 1 for KDE computation')
                    kde = KernelDensity(bandwidth=1, kernel='gaussian')
                    kde.fit(data)
                
                # compute density estimates for all points in the grid 
                density_temp = np.exp(kde.score_samples(xy)).reshape(X.shape)
                density_temp2 = np.reshape(np.exp(kde.score_samples(data)),(-1,1))
                levels = np.linspace(0, density_temp.max(), 25)
                
                print('    checking sum of density points for the area: ')
                print('   ', np.sum(density_temp))
                print('    checking sum of density points for the ROIs only: ')
                print('   ', np.sum(density_temp2))
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.contourf(X, Y, density_temp, levels=levels, cmap='Blues')
                delta = (abs(xx[1]-xx[0]),abs(yy[1]-yy[0]))
                ax.set_xlim(xx[0]-delta[0]*0.1,xx[-1]+delta[0]*0.1)
                ax.set_ylim(yy[0]-delta[1]*0.1,yy[-1]+delta[1]*0.1)
                ax.set_aspect('equal')
                fig.savefig(fname = key+'_'+str(label)+'_'+str(col0)+str(col1)+'.svg')
                fig.savefig(fname = key+'_'+str(label)+'_'+str(col0)+str(col1)+'.png')
                plt.close()
            del density_temp
    return bandwidth_dict
    #print(bandwidth_dict)
    
if __name__ == '__main__':
    """
    main processes a specific .mat file through some operations of the workflow, this main function was created 
    based on our needs and you may adapt to what you wish. For example, there is plot_KDE_2D that function that
    is not used here but it is useful if you 
    
     Input: 
    
    -p : file path to a .mat file 
    
    Operations: 
        Compute Multidimensional scaling (MDS) with 2 components and the x,y and z coordinates as input (columns 0 to 2) 
        Kernel Density Estimates in 2D, also for the first coordinates of the matrix (columns 0 and 1) 
        Plot 3D plots where Z coordinates are the density of the datapoints in x,y
        Compute the Kolmogorov-Smirnov (KS) test for 2D dimensions and 2 samples (parwise frequency tests for this dataset) 
        Compute the p-value matrix from the KS test values
      
    """
    import argparse
    import numpy as np
    import plotly as py
    import plotly.express as px
    import plotly.graph_objects as go 
    from scipy import io as sio 
    from matplotlib.colors import Colormap
    
    AP = argparse.ArgumentParser()
    AP.add_argument("-p", "--file_path", required=True, help="Filepath where .mat file is")
    ARGS = vars(AP.parse_args())
    
    metadata = ['__header__', '__version__', '__globals__']
    
    mds_dict = compute_MDS(ARGS['file_path'])
    #mds_dict = sio.loadmat(ARGS['file_path'][:-4]+'_MDS.mat')
    mds_density_dict = compute_KDE_2D(ARGS['file_path'][:-4]+'_MDS.mat')
    #mds_density_dict = sio.loadmat(ARGS['file_path'][:-4]+'_MDS_density.mat')
    
    p_val_matrices = {}
    print(mds_dict)
    for key, value in mds_dict.items():
        if key not in metadata: 
            print('computing stats for', key)
            fig = go.Figure()
            #print('density:')
            #print(mds_density_dict[key+'_density'])
            fig.add_trace(go.Scatter3d(x=value[:,0], y=value[:,1], z=mds_density_dict[key+'_density'][:,0], 
                                       mode='markers', 
                                       marker=dict(size=10,
                                                   color= value[:,-1], # set color to an array/list of desired values
                                                   colorscale='Spectral',
                                                   opacity=0.8))) 
            fig.update_layout(title=key)
            fig.show()
            
            #print('calculate two-dimensional KS in \n ', value[:,:-1], 'labels', value[:,-1])
            temp_p_vals = compute_ks2d2s(value[:,:-1], labels=value[:,-1])
            temp_p_matrix = p_values_to_matrix(temp_p_vals)
            
            p_val_matrices.update({key : temp_p_matrix})
    
    raw_xy_dict = sio.loadmat(ARGS['file_path'])
    raw_xy_density = compute_KDE_2D(ARGS['file_path'])
    #raw_xy_density = sio.loadmat(ARGS['file_path'][:-4]+'_density.mat')
    #print(raw_xy_dict, raw_xy_density)
    
    for key, value in raw_xy_dict.items():
        if key not in metadata: 
            print('computing stats for', key)
            temp_p_vals = compute_ks2d2s(value[:,:-2], labels=value[:,-1])
            temp_p_matrix = p_values_to_matrix(temp_p_vals)   
            p_val_matrices.update({key : temp_p_matrix})
            
            fig = go.Figure()
            fig.add_trace(go.Scatter3d(x=value[:,0], y=value[:,1], z=raw_xy_density[key+'_density'][:,0], 
                           mode='markers', 
                           marker=dict(size=10,
                                       color= value[:,-1], # set color to an array/list of desired values
                                       colorscale='Spectral',
                                       opacity=0.8))) 
            fig.update_layout(title=key)
            fig.show()
    
    sio.savemat(ARGS['file_path'][:-4]+'_p_values_matrices.mat', p_val_matrices)
    
# How to use the main pipeline
# python tonotopy_computations.py -p ROI_per_brain_region_to_python.mat 
# how to use functions within another python script or online
# from density_computations import * 