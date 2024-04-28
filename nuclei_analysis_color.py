from nuclei_analysis import *
import cv2
from skimage.draw import polygon2mask


def crop_nuclei(img, X, Y):
    X = np.array(X.split(','), dtype=float)
    Y = np.array(Y.split(','), dtype=float)
    X = (np.rint(X)).astype(int)
    Y = (np.rint(Y)).astype(int)

    Xmax, Xmin = np.amax(X), np.amin(X)
    Ymax, Ymin = np.amax(Y), np.amin(Y)

    masked_img = np.array(img, dtype=float)[Xmin:Xmax, Ymin:Ymax, :]

    return masked_img


def grayscale(img):
    
    coeffs = np.array([0.114, 0.587, 0.229])
    images_gray = (img.astype(float) * coeffs).sum(axis=-1)
    
    gray_arr = np.array(images_gray, dtype=float)
    
    del images_gray
    return np.mean(gray_arr), np.std(gray_arr)


def saturation(img, X, Y):
    hsv_img = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    _,S,_ = cv2.split(hsv_img)
    #S_arr = np.array(S)
    X = np.array(X.split(','), dtype=float)
    Y = np.array(Y.split(','), dtype=float)
    X = (np.rint(X)).astype(int)
    Y = (np.rint(Y)).astype(int)

    Xmax, Xmin = np.amax(X), np.amin(X)
    Ymax, Ymin = np.amax(Y), np.amin(Y)

    S_arr = np.array(S)[Xmin:Xmax, Ymin:Ymax]

    del hsv_img, S
    return np.mean(S_arr), np.std(S_arr)


def LAB(img, X, Y):
    lab_img = cv2.cvtColor(img, cv2.COLOR_BGR2LAB)
    _,A,B = cv2.split(lab_img)
    
    X = np.array(X.split(','), dtype=float)
    Y = np.array(Y.split(','), dtype=float)
    X = (np.rint(X)).astype(int)
    Y = (np.rint(Y)).astype(int)

    Xmax, Xmin = np.amax(X), np.amin(X)
    Ymax, Ymin = np.amax(Y), np.amin(Y)

    A_arr = np.array(A)[Xmin:Xmax, Ymin:Ymax]
    B_arr = np.array(B)[Xmin:Xmax, Ymin:Ymax]

    del lab_img, A, B
    return np.mean(A_arr), np.std(A_arr), np.mean(B_arr), np.std(B_arr)


def main(img_dir, poly_x, poly_y, color_df):
    gray_mean, gray_std, sat_mean, sat_std = pd.Series(dtype=object), pd.Series(dtype=object), pd.Series(dtype=object), pd.Series(dtype=object)
    A_mean, A_std, B_mean, B_std = pd.Series(dtype=object), pd.Series(dtype=object), pd.Series(dtype=object), pd.Series(dtype=object)

    for X, Y in zip(poly_x, poly_y):
        try:
            X_check = np.array(X.split(','), dtype=float)
            Y_check = np.array(Y.split(','), dtype=float)
        except:
            continue
       
        img = cv2.imread(img_dir, cv2.IMREAD_COLOR)
        nuclei_img = crop_nuclei(img, X, Y)

        gmean_temp, gstd_temp = grayscale(nuclei_img)
        gray_mean = pd.concat((gray_mean, pd.Series(gmean_temp)))
        gray_std = pd.concat((gray_std, pd.Series(gstd_temp)))

        smean_temp, sstd_temp = saturation(img, X, Y)
        sat_mean = pd.concat((sat_mean, pd.Series(smean_temp)))
        sat_std = pd.concat((sat_std, pd.Series(sstd_temp)))

        AM_temp, AS_temp, BM_temp, BS_temp = LAB(img, X, Y)
        A_mean = pd.concat((A_mean, pd.Series(AM_temp)))
        A_std = pd.concat((A_std, pd.Series(AS_temp)))
        B_mean = pd.concat((B_mean, pd.Series(BM_temp)))
        B_std = pd.concat((B_std, pd.Series(BS_temp)))
        
            
    

    feat_vector = pd.DataFrame(np.column_stack((gray_mean, gray_std, sat_mean, sat_std, A_mean, A_std, B_mean, B_std)), columns=['gray_mean', 'gray_std', 'saturation_mean', 'saturation_std', 
                                          'A_mean', 'A_std', 'B_mean', 'B_std'])
    #print(np.shape(feat_vector))

    color_df = pd.concat( (color_df, feat_vector), axis=0)
    
    del feat_vector, img, nuclei_img, gmean_temp, gstd_temp, gray_mean, gray_std, smean_temp, sstd_temp, sat_mean, sat_std, AM_temp, AS_temp,
    BM_temp, BS_temp, A_mean, A_std, B_mean, B_std

    return color_df
