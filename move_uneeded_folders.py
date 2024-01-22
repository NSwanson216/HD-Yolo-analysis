import os
import shutil
import pandas as pd

# double check this before running it, just to be safe

if __name__ == '__main__':
    csv_path = r'path\to\WGD\csv'
    WGD_df = pd.read_csv(csv_path, sep=',')
    src = r'path\to\tile\folders'
    dest = r'new\path\to\uneeded\folders'

    WGD_df.set_index('Sample ID', inplace=True)
    for name in os.listdir(src):
        full_path = os.path.join(src, name)
        if name not in WGD_df.index:
            shutil.move(full_path, dest)