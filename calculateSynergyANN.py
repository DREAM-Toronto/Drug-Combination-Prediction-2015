import csv
import json
import sys

from pybrain.datasets import SupervisedDataSet
from pybrain.supervised import BackpropTrainer
from pybrain.tools.customxml import NetworkWriter
from pybrain.tools.shortcuts import buildNetwork


def read_drug_data(features):
    with open(features, 'r') as csvfile:
        feature_reader = csv.reader(csvfile, delimiter=",")
        next(feature_reader)
        drug_dict = {}
        for r in feature_reader:
            drug_dict[r[0]] = r[1:]
    return drug_dict


def read_synergy_data(synergy):
    with open(synergy, 'r') as csvfile:
        feature_reader = csv.reader(csvfile, delimiter=",")
        next(feature_reader)
        synergy_dict = {}
        for r in feature_reader:
            if not synergy_dict.get(r[1], None):
                synergy_dict[r[1]] = {}
            if not synergy_dict.get(r[2], None):
                synergy_dict[r[2]] = {}
            synergy_dict[r[1]][r[2]] = float(r[3])
            synergy_dict[r[2]][r[1]] = float(r[3])
    return synergy_dict


def dump_drug_dict_as_flat(drug_dict, out):
    pass


def build_training_input(drug_dict, synergy_dict):
    dim = 0
    training_input_dict = {}
    for d1 in synergy_dict:
        training_input_dict[d1] = {}
        for d2 in synergy_dict[d1]:
            training_input_dict[d1][d2] = {'DRUG_1': d1,
                                           'DRUG_2': d2,
                                           'INPUT': drug_dict[d1] + drug_dict[d2],
                                           'OUTPUT': synergy_dict[d1][d2]}
            dim = len(drug_dict[d1] + drug_dict[d2])
    return training_input_dict, dim


if __name__ == '__main__':

    features = sys.argv[1]
    synergy = sys.argv[2]

    out = "feature_sets.csv"

    drug_dict = read_drug_data(features)
    synergy_dict = read_synergy_data(synergy)
    dump_drug_dict_as_flat(drug_dict, out)
    training_input,input_len = build_training_input(drug_dict, synergy_dict)
    # input_len = training_input[list(training_input.keys())[0]]['INPUT']
    target_len = 1
    ds = SupervisedDataSet(input_len, target_len)
    for t1 in training_input:
        for t2 in training_input[t1]:
            print("Input Vector", training_input[t1][t2]['INPUT'], training_input[t1][t2]['OUTPUT'])
            ds.addSample(training_input[t1][t2]['INPUT'], training_input[t1][t2]['OUTPUT'])


    n = buildNetwork(ds.indim, 2, ds.outdim, bias=True)
    t = BackpropTrainer(n, learningrate=0.001, momentum=0.05, verbose=True)
    print("Training")
    t.trainUntilConvergence(ds,
                            verbose=True)
    NetworkWriter.writeToFile(n, 'trainedNetwork.xml')

    # n = NetworkReader.readFrom('trainedNetwork_2.xml')

    predictions = {}
    for d1 in drug_dict:
        if not predictions.get(d1, None):
            predictions[d1]={}
        for d2 in drug_dict:
            predictions[d1][d2] = n.activate(drug_dict[d1] + drug_dict[d2])[0]

    with open('predictions_4.json', 'w') as outfile:
        json.dump(predictions, outfile)


