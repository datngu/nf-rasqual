#!/usr/bin/env python

import sys
import gzip
import argparse
import random

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract rasqual lead SNP - matrixQTL formated",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--normial", help="normial pval input - matrixQTL formated")
    parser.add_argument("--eigenMT", help="output of eigenMT")
    return parser.parse_args()

def process_eigenMT_file(filename):
    features = {}
    with open(filename, 'rt') as fi_BF:
        for line in fi_BF:
            if line.startswith("snp"):
                continue
            tokens = line.strip().split("\t")
            key = tokens[1]
            val = float(tokens[3])  # val is the p-value
            features[key] = (val, tokens[6], tokens[7].rstrip())
    return features


def process_normial_file(filename, features):
    results = {}
    with open(filename, 'rt') as fi_norm:
        # Read and process header line
        header = next(fi_norm).rstrip()
        print(f"{header}\tBF\tTESTS")
        
        # Process remaining lines
        for line in fi_norm:
            tokens = line.split("\t")
            key = tokens[1]
            
            try:
                val = float(tokens[3])
            except ValueError:
                continue
                
            smallest_pval, bf, tests = features.get(key, (-1, "NA", "NA"))
            
            if val <= smallest_pval:
                row = f'{line.rstrip()}\t{bf}\t{tests}'
                if key not in results:
                    results[key] = [row]
                else:
                    results[key].append(row)
    return results



def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Process input files
    features = process_eigenMT_file(args.eigenMT)
    results = process_normial_file(args.normial, features)
    
    # Set random seed and print results
    random.seed(2025)
    for key, val in results.items():
        if len(val) == 1:
            print(val[0])
        else:
            print(random.choice(val))

if __name__ == "__main__":
    main()


# python3 ties.py --eigenMT ALL_eigenMT_results.txt --normial nul_rasqual_normial_pval.txt > ties.txt
