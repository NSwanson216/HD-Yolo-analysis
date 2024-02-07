import os
import shutil
import pandas as pd

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
