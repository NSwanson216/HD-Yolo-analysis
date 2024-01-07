# Sorting through the files in glados in order remove the files we don't need, as well as doing the nuclei_ratio test

from nuclei_analysis import *
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Accuracy of nuclei segmentation scores.')
    parser.add_argument('--input_dir', required=True, type=str)
    args = parser.parse_args()
    input_dir = args.input_dir
    #input_dir = r'/titan/cancerregulome10/workspaces/users/btercan/ForBahar/Outputs/B587-TTP1'
    ratio = nuclei_ratios(input_dir, with_parent=False)
    print(ratio)
