import os
import glob
import argparse
from kpfunction import kp

# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('folder', type=str, help='Folder containing music files')
parser.add_argument('--tonics', type=str, default='tonics.txt', help='File containing tonics')
parser.add_argument('--pattern', type=str, default='*', help='File pattern to match')
args = parser.parse_args()

# initialize total degrees
total_degrees = [0 for _ in range(12)]

# load tonics from file
tonics = {}
with open(args.tonics, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            filename = parts[0]
            tonic = parts[1].lower()
            tonics[filename] = tonic

# process each file in the folder
file_pattern = os.path.join(args.folder, args.pattern)
for filepath in glob.glob(file_pattern):
    if os.path.isdir(filepath):
        continue
        
    filename = os.path.basename(filepath)
    if filename not in tonics:
        print(f"Skipping {filename}: No tonic defined")
        continue
    
    # analyze the file
    file_degrees = kp(filepath, tonics[filename])
    
    # add to total
    total_degrees = [total_degrees[i] + file_degrees[i] for i in range(12)]
    
    print(f"Processed {filename} (tonic: {tonics[filename]})")
    print(f"  Degrees: {file_degrees}")

# print the final total
print("\nTotal degree counts:")
print(total_degrees)
