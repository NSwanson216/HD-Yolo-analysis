import os
import pandas as pd
import numpy as np
import scipy.spatial as scp
import matplotlib.pyplot as plt
import seaborn
import argparse

#threshold = 15.54 #placeholder threshold from the liver cancer version
#threshold = 20 # this will be replaced as we work through it and find a suitable value

#nuclei_area_distribution =  #temp distribution - use mean ploidy of each slide(?

checker = []

class Dictlist(dict):
    #subclass for multiple entries per dict key
    def __setitem__(self, key, value):
        try:
            self[key]
        except KeyError:
            super(Dictlist, self).__setitem__(key, [])
        self[key].append(value)


class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y)
    def __add__(self, other):
        return Vector(self.x + other.x, self.y + other.y)
    def dot(self, other):
        return self.x * other.x + self.y * other.y
    def norm(self):
        return self.dot(self)**0.5
    def normalized(self):
        norm = self.norm()
        return Vector(self.x / norm, self.y / norm)
    def perp(self):
        return Vector(1, -self.x / self.y)
    def __mul__(self, scalar):
        return Vector(self.x * scalar, self.y * scalar)
    def __str__(self):
        return f'({self.x}, {self.y})'


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

    return np.sum( np.sqrt((X - np.roll(X,-1))**2 + (Y - np.roll(Y,-1))**2) )
    

def axis_length(X, Y):
    X = np.array(X.split(','), dtype=float)
    Y = np.array(Y.split(','), dtype=float)
    
    distances = np.zeros((len(X),len(X)))
    for i in range(len(X)):        
        for j in range(len(X)):
            distances[i,j] = np.sqrt( (X[i] - X[j])**2 + (Y[i] - Y[j])**2 )
    
    max_dist = np.amax(distances)
    
    ind = np.unravel_index(np.argmax(distances, axis=None), distances.shape)
    ind = np.array(ind)
    XY = Vector(X[ind[0]], Y[ind[0]]) - Vector(X[ind[1]], Y[ind[1]])
    XY_normed = XY.perp().normalized()
    P1 = Vector(X[ind[1]], Y[ind[1]]) + XY_normed * 3
    P2 = Vector(X[ind[1]], Y[ind[1]]) - XY_normed * 3
   
    min_dist = np.sqrt( (P1.x - P2.x)**2 + (P1.y - P2.y)**2 )
    
    return max_dist, min_dist
    
    
def circularity(A, P):
    return (P**2) / (4 * np.pi * A)
    

def eccentricity(max_dist, min_dist):
    return np.sqrt( 1 - ( max_dist**2 / min_dist**2 ))


def solidity(X, Y, A):
    X = np.array(X.split(','), dtype=float)
    Y = np.array(Y.split(','), dtype=float)
    XY = np.column_stack((X.T, Y.T))
    # A / area_of_convex_hull
    convex_hull = scp.ConvexHull(points=XY, qhull_options='QG4')
    return A / convex_hull.volume()

    
def rel_distance(coord_0, coord_1):
    # calc distance between two points
    x0, y0 = coord_0[0], coord_0[1]
    x1, y1 = coord_1[0], coord_1[1]

    dist = np.sqrt( (x0 - x1)**2 + (y1 - y0)**2 )
    return dist


def delete_duplicates(dict1):
    check, TBD = [], []

    #df = pd.DataFrame.from_dict(dict1, orient='index')
    #lengths = []
    #for val in dict1.values():
    #    lengths.append(len(val))
    #df['length'] = lengths
    #df.sort_values('length', inplace=True, ascending=False)
    ##continue to workshop this part
    #dict1 = df.T.to_dict(orient='list')
    
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
    

def calculate_ploidy(input_dir, ploidy_values, ground_truth, threshold):
    c_ploidy, n_ploidy = 0, 0
    nuclei_area = pd.Series(dtype=object)
    raw_areas, perims, max_dist, min_dist = pd.Series(dtype=object), pd.Series(dtype=object), pd.Series(dtype=object), pd.Series(dtype=object)
    circuls, eccents, solids = pd.Series(dtype=object), pd.Series(dtype=object), pd.Series(dtype=object)
        
    for csv_file in os.listdir(input_dir):
        if os.path.splitext(csv_file)[1] == '.csv':
            df_path = os.path.join(input_dir, csv_file)
            df = pd.read_csv(df_path, sep=',')
            temp = df.loc[df['label'] != -100]

            poly_x = list(temp['poly_x'])
            poly_y = list(temp['poly_y'])
             
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

            #calc area for each nuclei and compare it against a distribution to determine it's nuclear ploidy
            nuclei_area = pd.concat((nuclei_area, areas))

            # somehow take these and include them into the csv file
            # they are per tile, but we need them per sample - maybe save them as a list in a single df entry per sample (perimeter holds a comma
            # seperated list of all of the perimeter entries) and then expand that list when it comes to using it for ML
            # we could duplicate the rest of the stagnant values and then use these new values as a way to seperate each tile, all under the 
            # label of a single sample(? - but this brings back the issue of heterogeneity

            # take mean and stddev from each feature 
            
            for X, Y, counter in zip(poly_x, poly_y, range(len(poly_x))):
                raw_areas_temp = area(X, Y)
                raw_areas = pd.concat((raw_areas, raw_areas_temp))
                perims_temp = perimeter(X, Y)
                perims = pd.concat((perims, perims_temp))
                max_dist_temp, min_dist_temp = axis_length(X, Y)
                max_dist = pd.concat((max_dist, max_dist_temp))
                min_dist = pd.concat((min_dist, min_dist_temp))
                circuls_temp = circularity(raw_areas_temp, perims_temp)
                circuls = pd.concat((circuls, circuls_temp))
                eccents_temp = eccentricity(max_dist_temp, min_dist_temp)
                eccents = pd.concat((eccents, eccents_temp))
                solids_temp = solidity(X, Y, raw_areas_temp)
                solids = pd.concat((solids, solids_temp))


            del denom, numer, multi_nuclei_dict, multi_nuclei_list, multi_nuclei_df, nuclei_list, 
            centroids, centroids_arr, centroids_df
            

    n_ploidy = area_distribution(nuclei_area)
    areas = np.mean(raw_areas)
    areas_std = np.std(raw_areas)
    perimeters = np.mean(perims)
    perimeters_std = np.std(perims)
    major_axis = np.mean(max_dist)
    major_axis_std = np.std(max_dist)
    minor_axis = np.mean(min_dist)
    minor_axis_std = np.std(min_dist)
    circularities = np.mean(circuls)
    circularities_std = np.std(circuls)
    eccentricities = np.mean(eccents)
    eccentricities_std = np.std(eccents)
    solidities = np.mean(solids)
    solidities_std = np.std(solids)


    ploidy_values.loc[csv_file[0:9], 'ground_truth'] = ground_truth.loc[csv_file[0:9], 'Ploidy']
    ploidy_values.loc[csv_file[0:9], 'n_ploidy'] = n_ploidy
    ploidy_values.loc[csv_file[0:9], 'c_ploidy'] = c_ploidy / len(os.listdir(input_dir))

    ploidy_values.loc[csv_file[0:9], 'area'] = areas
    ploidy_values.loc[csv_file[0:9], 'area_std'] = areas_std
    ploidy_values.loc[csv_file[0:9], 'perimeter'] = perimeters
    ploidy_values.loc[csv_file[0:9], 'perimeter_std'] = perimeters_std
    ploidy_values.loc[csv_file[0:9], 'major_axis'] = major_axis
    ploidy_values.loc[csv_file[0:9], 'major_axis_std'] = major_axis_std
    ploidy_values.loc[csv_file[0:9], 'minor_axis'] = minor_axis
    ploidy_values.loc[csv_file[0:9], 'minor_axis_std'] = minor_axis_std
    ploidy_values.loc[csv_file[0:9], 'circularity'] = circularities
    ploidy_values.loc[csv_file[0:9], 'circularity_std'] = circularities_std
    ploidy_values.loc[csv_file[0:9], 'eccentricity'] = eccentricities
    ploidy_values.loc[csv_file[0:9], 'eccentricity_std'] = eccentricities_std
    ploidy_values.loc[csv_file[0:9], 'solidity'] = solidities
    ploidy_values.loc[csv_file[0:9], 'solidity_std'] = solidities_std

    
    del areas, c_ploidy, n_ploidy, counter, counter2, counter_name
    
    return ploidy_values
            

def main(input_dir, output_dir, ploidy_dir, threshold):#args, threshold):
    output_dir = os.path.join(output_dir, 'ploidy_values.csv')#args.output_dir, 'ploidy_values.csv')
    ploidy_values = pd.DataFrame( columns=['Sample_ID', 'c_ploidy', 'n_ploidy', 'ground_truth', 'area', 'area_std', 'perimeter', 'perimeter_std',
                                           'major_axis', 'major_axis_std', 'minor_axis', 'minor_axis_std', 'circularity', 'circularity_std',
                                           'eccentricity', 'eccentricity_std', 'solidity', 'solidity_std'] )
    ploidy_values.set_index('Sample_ID', inplace=True)

    ground_truth = pd.read_csv(ploidy_dir, sep=',')
    ground_truth.set_index('Submitter ID', inplace=True)

    for sub_folder in os.listdir(input_dir):#args.input_dir):
        sub_input = os.path.join(input_dir, sub_folder)#args.input_dir, sub_folder)
        ploidy_values = calculate_ploidy(sub_input, ploidy_values, ground_truth, threshold)

    return ploidy_values
    #os.makedirs(output_dir, mode=0o777, exist_ok=True)   
    #ploidy_values.to_csv(output_dir)
