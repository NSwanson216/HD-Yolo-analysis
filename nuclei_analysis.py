import os
import pandas as pd
import numpy as np
import scipy.spatial as scp
import matplotlib.pyplot as plt
import seaborn
import argparse
import sklearn.preprocessing as preprocessing
from sklearn.decomposition import FactorAnalysis, FastICA, PCA, KernelPCA, SparseCoder, DictionaryLearning
from sklearn import manifold
import scipy 
import torch

import nuclei_analysis_color as nac

os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:512"
checker = []

class Dictlist(dict):
    #subclass for multiple entries per dict key
    def __setitem__(self, key, value):
        try:
            self[key]
        except KeyError:
            super(Dictlist, self).__setitem__(key, [])
        self[key].append(value)


def nuclei_ratios(input_dir, with_parent=False):
    #check the average score of the nuclei segmentation 
    length_sum = 0
    total_nuclei = 0

    if with_parent == True:
        for item in os.listdir(input_dir):
            folder = os.path.join(input_dir, item)
            for sub_item in os.listdir(folder):
                if os.path.splitext(sub_item)[1] == '.csv':
                    df_path = os.path.join(folder, sub_item)
                    df = pd.read_csv(df_path, sep=',')
                    temp = df.loc[df['label'] != -100]
                    total_nuclei += len(temp)
                    scores = np.array(temp['score'])
                    length_sum += np.sum(scores)
    elif with_parent == False:
        for sub_item in os.listdir(input_dir):
                if os.path.splitext(sub_item)[1] == '.csv':
                    df_path = os.path.join(input_dir, sub_item)
                    df = pd.read_csv(df_path, sep=',')
                    temp = df.loc[df['label'] != -100]
                    total_nuclei += len(temp)
                    scores = np.array(temp['score'])
                    length_sum += np.sum(scores)
    else:
        ValueError
        
    return length_sum / total_nuclei


def calc_centroids(X, Y):
    if type(X) == float:
        print('problem here')

    X = np.array(X.split(','), dtype=float)
    Y = np.array(Y.split(','), dtype=float)
    
    #calc area
    A = 0.5 * np.abs(np.dot(X, np.roll(Y,1)) - np.dot(Y, np.roll(X,1)))
    
    #calc centroid coords
    length = len(X)
    sum_X = np.sum(X)
    sum_Y = np.sum(Y)

    return sum_X/length, sum_Y/length, A


def area(X, Y):
    X = np.array(X.split(','), dtype=float)
    Y = np.array(Y.split(','), dtype=float)
    return 0.5 * np.abs(np.dot(X, np.roll(Y,1)) - np.dot(Y, np.roll(X,1)))


def perimeter(X, Y):
    X = np.array(X.split(','), dtype=float)
    Y = np.array(Y.split(','), dtype=float)
    P_value = np.sum( np.sqrt((X - np.roll(X,-1))**2 + (Y - np.roll(Y,-1))**2) )

    return P_value
    

def axis_length(X, Y):
    X = np.array(X.split(','), dtype=float)
    Y = np.array(Y.split(','), dtype=float)
    
    distances = np.zeros((len(X),len(X)))
    for i in range(len(X)):        
        for j in range(len(X)):
            distances[i,j] = np.sqrt( (X[i] - X[j])**2 + (Y[i] - Y[j])**2 )
    
    max_dist = np.amax(distances)
    max_dist_coord = np.unravel_index(np.argmax(distances, axis=None), distances.shape)

    My = (Y[max_dist_coord[1]] + Y[max_dist_coord[0]]) / 2
    Mx = (X[max_dist_coord[1]] + X[max_dist_coord[0]]) / 2

    rho_M = max_dist / 2

    M = np.array([ X[max_dist_coord[0]] - Mx, Y[max_dist_coord[0]] - My])

    Θ_range = [(17 * np.pi) / 36, (19 * np.pi) / 36]
    Θ1, Θ2 = 1, 1

    while Θ1 > np.cos(Θ_range[0]) and Θ2 > np.cos(Θ_range[0]):
        for i in range(len(X)):
            if Θ1 > np.cos(Θ_range[0]):
                px1, py1 = X[max_dist_coord[0] + i], Y[max_dist_coord[0] + i]
                Θ1, rho_P1 = find_Θ(Mx, My, M, px1, py1, rho_M)
            if Θ2 > np.cos(Θ_range[0]):
                px2, py2 = X[max_dist_coord[0] - i], Y[max_dist_coord[0] - i]
                Θ2, rho_P2 = find_Θ(Mx, My, M, px2, py2, rho_M)

    min_dist = rho_P1 + rho_P2

    return max_dist, min_dist


def find_Θ(Mx, My, M, px, py, rho_M):
    rho_P = np.sqrt(( Mx - px)**2 + (My - py)**2 )
    P = np.array([px - Mx, py - My])
    Θ = (np.dot(M, P)) / (np.abs(rho_M) * np.abs(rho_P))
    return Θ, rho_P


def circularity(A, P):
    return (P**2) / (4 * np.pi * A)
    

def eccentricity(max_dist, min_dist):
    return (1 - min_dist/max_dist)


def solidity(X, Y, A): #the values of the maxes of solidity differ by ~10^-13
    X = np.array(X.split(','), dtype=float)
    Y = np.array(Y.split(','), dtype=float)
    XY = np.column_stack((X.T, Y.T))
    # A / area_of_convex_hull
    convex_hull = scp.ConvexHull(points=XY, qhull_options='QG4')
    return A / convex_hull.volume

    
def rel_distance(coord_0, coord_1):
    # calc distance between two points
    x0, y0 = coord_0[0], coord_0[1]
    x1, y1 = coord_1[0], coord_1[1]

    dist = np.sqrt( (x0 - x1)**2 + (y1 - y0)**2 )
    return dist


def delete_duplicates(dict1):
    check, TBD = [], []
    
    for key, value in dict1.items():
        for val in value:
            if val not in check:
                check.append(val)
        if key not in check:
            check.append(key)
        else:
            TBD.append(key)
    
    for item in TBD:
        del dict1[item]
    
    return dict1, check


def normalize(values, bounds):
    return [bounds['desired']['lower'] + (x - bounds['actual']['lower']) * 
            (bounds['desired']['upper'] - bounds['desired']['lower']) / 
            (bounds['actual']['upper'] - bounds['actual']['lower']) for x in values]


def area_distribution(nuclei_area):
    #create distribution for areas
    x = np.array(nuclei_area)
    """ mean = np.mean(x) 
    std = np.std(x) 
    y = 1/(std * np.sqrt(2 * np.pi)) * np.exp( - (x - mean)**2 / (2 * std**2))  """

    x_normed = normalize(x, {'actual': {'lower' : np.min(x), 'upper': np.max(x)}, 
                             'desired': {'lower' : 1, 'upper' : 4} })
    mean = np.mean(x_normed) 
    std = np.std(x_normed) 
    y = 1/(std * np.sqrt(2 * np.pi)) * np.exp( - (x_normed - mean)**2 / (2 * std**2)) 
    
    """ plt.figure(figsize = (6, 6)) 
    plt.plot(x_normed, y, color = 'black', linestyle = 'dashed') 
  
    plt.scatter( x_normed, y, marker = 'o', s = 25, color = 'red') 
    plt.show()   """

    return mean
    

def dataset_transformaiton(temp_matrix):
    return_list = []

    means = np.mean(temp_matrix, axis=0)
    return_list.append(means)

    maxes = np.max(temp_matrix, axis=0) 
    return_list.append(maxes)
    
    temp_matrix = torch.from_numpy(temp_matrix.T)
    temp_matrix = temp_matrix.cuda()
    U, sigma, V = torch.svd(temp_matrix)
    U = U.cpu().numpy()
    sigma = sigma.cpu().numpy()
    V = V.cpu().numpy()
    sigma = np.diag(sigma)
    V = V[0:15,0:15]
    A = U @ sigma @ V
    eigvals, eigenvector = np.linalg.eig(A)
    index = eigvals.argsort()[::-1]
    eigvals = eigvals[index]
    eigenvector = eigenvector[:,index]
    leading_eigenvector = eigenvector[:,0]
    leading_3_eigenvectors = np.concatenate( (eigenvector[:,0], eigenvector[:,1], eigenvector[:,2]), axis=0 )
    return_list.append(leading_eigenvector)
    return_list.append(leading_3_eigenvectors)

    del temp_matrix, U, sigma, V, A

    return return_list

    


def calculate_ploidy(input_dir, image_dir):
    c_ploidy, n_ploidy = 0, 0
    nuclei_area = pd.Series(dtype=object)
    raw_areas, perims, max_dist, min_dist = pd.Series(dtype=object), pd.Series(dtype=object), pd.Series(dtype=object), pd.Series(dtype=object)
    circuls, eccents, solids = pd.Series(dtype=object), pd.Series(dtype=object), pd.Series(dtype=object)
    color_feat_df_temp = pd.DataFrame(columns=['gray_mean', 'gray_std', 'saturation_mean', 'saturation_std', 
                                          'A_mean', 'A_std', 'B_mean', 'B_std'])
    

    for csv_file in os.listdir(input_dir):
        if os.path.splitext(csv_file)[1] == '.csv':
            for img_name in os.listdir(image_dir):
                if os.path.splitext(img_name)[0] == os.path.splitext(csv_file)[0][:-5]:
                    img_file = img_name
                    break
            df_path = os.path.join(input_dir, csv_file)
            df = pd.read_csv(df_path, sep=',')
            temp = df

            poly_x = list(temp['poly_x'])
            poly_y = list(temp['poly_y'])
            if len(poly_x) == 0:
                continue
            '''
            centroids = np.zeros((len(poly_x), 2), dtype=object)
            counter_name = []
            areas = np.zeros((len(poly_x)))

            #calc centroids of each nuclei
            for X, Y, counter in zip(poly_x, poly_y, range(len(poly_x))):
                centroids[counter, 0], centroids[counter, 1], areas[counter] = calc_centroids(X, Y) 
                counter_name.append('nuclei_' + str(counter))
            areas = pd.Series(areas)

            #calc relative distance for each nuclei
            centroids_arr = np.zeros((len(centroids), len(centroids)))
            for counter in range(len(centroids)):
                for counter2 in range(len(centroids)):
                    centroids_arr[counter, counter2] = rel_distance( centroids[counter, :], centroids[counter2,:] )
            
            centroids_df = pd.DataFrame(centroids_arr, columns=counter_name)
            centroids_df['nuclei_name'] = counter_name
            centroids_df.set_index('nuclei_name', inplace=True)

            grouping_coords = np.where( (centroids_arr > 0) & (centroids_arr < threshold) )
            multi_nuclei_df = centroids_df.iloc[grouping_coords]
                        
            multi_nuclei_dict = Dictlist()
            for index, column in zip(multi_nuclei_df.index, multi_nuclei_df.columns):
                multi_nuclei_dict[index] = column

            multi_nuclei_dict, multi_nuclei_list = delete_duplicates(multi_nuclei_dict)
            nuclei_list = list( centroids_df.columns.drop(multi_nuclei_list) )
            
            numer, denom = 0, 0
            for val in multi_nuclei_dict.values():
                numer += (len(val) + 1)
                denom += 1
            numer += len(nuclei_list)
            denom += len(nuclei_list)

            if denom == 0:
                continue
            
            c_ploidy += (numer / denom)

            nuclei_area = pd.concat((nuclei_area, areas))
            '''
            
            for X, Y, counter in zip(poly_x, poly_y, range(len(poly_x))):
                try:
                    raw_areas_temp = area(X, Y)
                    raw_areas = pd.concat((raw_areas, pd.Series(raw_areas_temp)))
                    perims_temp = perimeter(X, Y)
                    perims = pd.concat((perims, pd.Series(perims_temp)))
                    max_dist_temp, min_dist_temp = axis_length(X, Y)
                    max_dist = pd.concat((max_dist, pd.Series(max_dist_temp)))
                    min_dist = pd.concat((min_dist, pd.Series(min_dist_temp)))
                    circuls_temp = circularity(raw_areas_temp, perims_temp)
                    circuls = pd.concat((circuls, pd.Series(circuls_temp)))
                    eccents_temp = eccentricity(max_dist_temp, min_dist_temp)
                    eccents = pd.concat((eccents, pd.Series(eccents_temp)))
                    solids_temp = solidity(X, Y, raw_areas_temp)
                    solids = pd.concat((solids, pd.Series(solids_temp)))
                    
                except:
                    pass
                

            
            img_path = os.path.join(image_dir, img_file)
            color_feat_df_temp = nac.main(img_path, poly_x, poly_y, color_feat_df_temp)
           
            '''
            del denom, numer, multi_nuclei_dict, multi_nuclei_list, multi_nuclei_df, nuclei_list, 
            centroids, centroids_arr, centroids_df, 
            '''
            del raw_areas_temp, perims_temp, max_dist_temp, min_dist_temp,
            circuls_temp, eccents_temp, solids_temp
            
            
    temp_matrix = np.zeros((len(color_feat_df_temp), 15))
    temp_matrix[:,0] = raw_areas
    temp_matrix[:,1] = perims
    temp_matrix[:,2] = max_dist
    temp_matrix[:,3] = min_dist
    temp_matrix[:,4] = circuls
    temp_matrix[:,5] = eccents
    temp_matrix[:,6] = solids
    temp_matrix[:,7] = color_feat_df_temp['gray_mean']
    temp_matrix[:,8] = color_feat_df_temp['gray_std']
    temp_matrix[:,9] = color_feat_df_temp['saturation_mean']
    temp_matrix[:,10] = color_feat_df_temp['saturation_std']
    temp_matrix[:,11] = color_feat_df_temp['A_mean']
    temp_matrix[:,12] = color_feat_df_temp['A_std']
    temp_matrix[:,13] = color_feat_df_temp['B_mean']
    temp_matrix[:,14] = color_feat_df_temp['B_std']
    
    temp_matrix = temp_matrix[~np.isnan(temp_matrix).any(axis=1), :]

    datasets = []
    datasets = dataset_transformaiton(temp_matrix)
    
    del c_ploidy, n_ploidy, counter, raw_areas, perims, max_dist, min_dist, circuls, eccents, solids, color_feat_df_temp,
    nuclei_area, temp_matrix #areas, counter2, counter_name
    
    return datasets
            

def main(input_dir,image_dir, output_list):
     
    #ground_truth = pd.read_csv(ploidy_dir, sep=',')
    #ground_truth.set_index('Submitter ID', inplace=True)
    
    img_names = []
   
    for sub_folder in os.listdir(input_dir):
        for img_folder in os.listdir(image_dir):
            if img_folder == sub_folder:
                sub_img = img_folder
                break
        sub_input = os.path.join(input_dir, sub_folder)
        sub_img_input = os.path.join(image_dir, sub_img)
        datasets = calculate_ploidy(sub_input, sub_img_input)
        csv_names = []
        for output_dir in output_list:
            csv_names.append(os.path.join(output_dir, (sub_folder + '.csv')))
        
        for dataset, csv_name in zip(datasets, csv_names):
            if len(dataset) < 2:
                dataset = np.array([dataset]).T
            else:
                dataset = np.array([dataset])

            try:   
                dataset = pd.DataFrame(dataset)
            except:
                dataset = pd.DataFrame(dataset[:,:,0])
            dataset.to_csv(csv_name)


        sep='.'
        img_name = sub_folder.split(sep,1)[0]
        img_names.append(img_name)

    return img_names
    