import os
import shutil
import pandas as pd
from distutils.dir_util import copy_tree

# double check this before running it, just to be safe

if __name__ == '__main__':
    csv_path = r'/titan/cancerregulome10/workspaces/users/btercan/ForBahar/MILD/clinical_data/MILD_WGD_p-p.tsv'
    WGD_df = pd.read_csv(csv_path, sep=',')
    src = r'/titan/cancerregulome10/workspaces/users/btercan/ForBahar/MILD/tiles_FFPE'
    dest = r'/titan/cancerregulome10/workspaces/users/btercan/ForBahar/MILD/unused_tiles'

    WGD_df.set_index('Submitter ID', inplace=True)
    for name in os.listdir(src):
        full_path = os.path.join(src, name)
        if name not in WGD_df.index:
            shutil.move(full_path, dest)

    MILD_path = r'/titan/cancerregulome10/workspaces/users/btercan/ForBahar/MILD'
    subset_path = os.path.join(MILD_path, 'subsets')
    if os.path.exists(subset_path) != True:
        os.mkdir(subset_path)
    new_dir = []
    for i in range(1, 6):
        subset_name = 'subset_' + str(i)
        new_dir.append(os.path.join(subset_path, subset_name))
        if os.path.exists(new_dir[i-1]) != True:
            os.mkdir(new_dir[i-1])
    i = 0
    for name in os.listdir(src):
        full_path = os.path.join(src, name)
        new_path = os.path.join(new_dir[i], name)
        copy_tree(full_path, new_path)
        if i == 4:
            i = 0
        else:
            i += 1
        del full_path, new_path





        
