import hd_yolo_analysis
from hd_yolo_analysis import *
from time import time


if __name__ == '__main__':
    threshold = np.arange(15, 40, 5)
    input_dir = r"C:\Users\fiddl\IDC-GDC\CCG-MILD\test_patches_masks"
    output_dir = r"C:\Users\fiddl\IDC-GDC\CCG-MILD\Outputs"
    ploidy_dir = r"G:\My Drive\ISB - Work\Texture Analysis\Clinical Data\MILD\MILD_WGD_TRUNCATED.csv"
    i = 0
    time1 = time()
    
    for thresh in threshold:
        i+=1
        output_sub = os.path.join(output_dir, ('output_'+str(i)+'.csv'))
        ploidy_values = main(input_dir, output_sub, ploidy_dir, thresh)
        ploidy_values.to_csv(output_sub)
        elapsed = time() - time1
        print(elapsed)