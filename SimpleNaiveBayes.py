#%% This code contains functionality that read in probablity table defined by radiologist, patient feature.
# and calculates naive bayes probability for each disease
#%% Import libraries
# Read in the worksheet
import openpyxl
from openpyxl import Workbook, load_workbook
import csv
import numpy as np
import argparse
import os
#%% Utility functions
# Function that reads probability table
# Input:
#   table_file: the absolute path to the xlsx probability table
#   tableName: the worksheet name inside the file
# Output:
#
def readProbabilityTable(table_file, tableName):
    wb=openpyxl.load_workbook(table_file)
    worksheet = wb[tableName]

    #Read the keys
    rows = list(worksheet.rows)
    first_row = rows[0]
    second_row = rows[1]

    keys = [c.value for c in first_row]
    values = [c.value for c in second_row]

    values = [v for v in values if v is not None]
    keys = keys[:len(values)]

    values = values[1:]
    keys = keys[1:]

    key_value_dict = {}

    current_k = ""
    for k, v, i in zip(keys, values, range(2, len(values) + 2)):
        if k is not None:
            current_k = k.split(" ")[0].upper()
            key_value_dict[current_k] = {v.split(" ")[0].upper():i - 2}
        else:
            key_value_dict[current_k][v.split(" ")[0].upper()] = i - 2

    # Read the probabilty
    disease = []
    prob = []
    for row in rows[2:]:
        row_values = [c.value for c in row if c.value is not None]

        if len(row_values) ==0:
            break

        if row_values[0] is None:
            break
        disease.append(row_values[0])
        prob.append(row_values[1:])
    prob = np.array(prob)
    prob = prob/100.0

    return disease, prob, key_value_dict

# Function that reads patient info sheet
def ReadPatientFeatFile(feature_file):
    feat_arr = []
    with open(feature_file, 'rb') as csvfile:
        reader = csv.reader(csvfile)

        keys = next(reader)
        keys = [k.strip().upper().replace("\xef\xbb\xbf", "") for k in keys]

        for row in reader:
            row = [r.strip().upper() for r in row]
            if row[0] == "":
                continue
            feat_arr.append(dict(zip(keys, row)))

    return feat_arr

# The Naive bayes function
def NaiveBayes( feature, key_value_dict, prob):
    ndisease = prob.shape[0]

    pd = []
    value_chain = []
    # Calculate probability for each disease
    for d in range(ndisease):
        p = 1.0

        vc = []

        for k, v in feature.items():
            if k not in key_value_dict.keys():
                continue
            if v == "NA":
                continue
            c = key_value_dict[k][v]
            p = p * prob[d, c]
            vc.append((k, v, prob[d, c]))

        pd.append(p)
        value_chain.append(vc)

    pd = np.array(pd)
    pd = pd / np.sum(pd)
    return pd, value_chain

# Function that saves the output probability calculated by the naive bayes
def save_result(output_file, feat_arr, pd_arr, disease):
    with open(output_file, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(["ACC"] + disease)

        for feat, pd in zip(feat_arr, pd_arr):
            writer.writerow([feat["ACC"]] + list(pd))

def save_value_chain(vc_file, value_chain, disease):
    with open(vc_file, 'wb') as f:
        writer = csv.writer(f)
        for vc, d in zip(value_chain, disease):
            vc_to_print = []

            for v in vc:
                vc_to_print.extend(v)
            writer.writerow([d] + vc_to_print)

def save_top3(top3_file, feat_arr, pd_arr, disease):
    with open(top3_file, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(["ACC", "disease1", "disease2", "disease3", "prob1", "prob2", "prob3"])

        for feat, pd in zip(feat_arr, pd_arr):
            top3_idx = pd.argsort()[-3:][::-1]
            top3_pd = pd[top3_idx]
            top3_disease = [ disease[top3_idx[0]], disease[top3_idx[1]], disease[top3_idx[2]] ]
            writer.writerow([feat["ACC"], top3_disease[0], top3_disease[1], top3_disease[2], top3_pd[0], top3_pd[1], top3_pd[2] ] )

# This function adds trailing "/"
def addTrailingSlash(directory):
    if directory[-1] != "/":
        directory+="/"
    return directory


# The main function
def main():
    parser = argparse.ArgumentParser(description='Naive Bayes for the White matter/Basal ganglia project. ')
    parser.add_argument('table_file', type=str, help='The probability table in xlsx format.')
    parser.add_argument('tableName', type=str, help='The worksheet name the probability table file.')
    parser.add_argument('feat_file', type=str, help='The feature file for each patients.')
    parser.add_argument('output_file', type=str, help='The output file for each patients.')
    parser.add_argument('--vc_dir', type=str, help='The output directory of value chain file. Optional')
    parser.add_argument('--top3_file', type=str, help='The output file of top 3 disease. Optional')
    args = parser.parse_args()

    value_chain_arr = []
    pd_arr = []

    disease, prob, key_value_dict = readProbabilityTable(args.table_file, args.tableName)
    feat_arr = ReadPatientFeatFile(args.feat_file)

    for feat in feat_arr:
        pd, value_chain = NaiveBayes( feat, key_value_dict, prob)
        pd_arr.append(pd)
        value_chain_arr.append(value_chain)

    save_result(args.output_file, feat_arr, pd_arr, disease)

    if args.vc_dir:
        args.vc_dir = addTrailingSlash(args.vc_dir)
        for value_chain, feat in zip(value_chain_arr,feat_arr):
            acc = feat["ACC"]
            save_value_chain(args.vc_dir + acc + ".csv", value_chain, disease)

    if args.top3_file:
        save_top3(args.top3_file, feat_arr, pd_arr, disease)


if __name__ == "__main__":
    main()
