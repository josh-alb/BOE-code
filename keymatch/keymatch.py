import argparse
import itertools
import re
import pathlib

# Get files
parser = argparse.ArgumentParser()
parser.add_argument('pred_file', type=pathlib.Path)
parser.add_argument('true_file', type=pathlib.Path)
args = parser.parse_args()

# Read keys
pred_keys = args.pred_file.read_text().strip().split()
true_keys = [
    key
    for key, repetitions in re.findall(
        r'([A-G][#b]?[m]?)\s*(\d+)',
        args.true_file.read_text()
    )
    for _ in range(int(repetitions))
]

n_pred_keys = len(pred_keys)
n_true_keys = len(true_keys)
if n_pred_keys != n_true_keys:
    raise ValueError(
        f'Key sequences are of inequal length: {n_pred_keys} vs {n_true_keys}'
    )

# Calculate accuracy
matching_keys = [
    pred_key == true_key
    for pred_key, true_key in itertools.zip_longest(pred_keys, true_keys)
]
print("correct", matching_keys.count(True))
print("total", len(matching_keys))
print(sum(matching_keys) / len(matching_keys))

