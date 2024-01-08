from nuclei_analysis import *
import cv2

def crop_nuclei(img, X, Y):
    X = np.array(X.split(','), dtype=float)
    Y = np.array(Y.split(','), dtype=float)

    crop_img = img[X, Y]
    return crop_img

#from here, start to get the different nuclei features