import os
import math
import numpy as np
import argparse

def convert_size(size_bytes):
   if size_bytes == 0:
       return "0B"
   size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
   i = int(math.floor(math.log(size_bytes, 1024)))
   p = math.pow(1024, i)
   s = round(size_bytes / p, 2)
   return "%s %s" % (s, size_name[i])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', required=True, type=str)
    args = parser.parse_args()
    current_dir = args.dir
    size = {}
    for folder in os.listdir(current_dir):
        full_folder = os.path.join(current_dir, folder)
        size[folder] = 0
        for ele in os.scandir(full_folder):
            size[folder] += os.path.getsize(ele)

    max_key = max(size, key=lambda k: size[k])
    max_size = convert_size(size[max_key])
    print(max_size)
    