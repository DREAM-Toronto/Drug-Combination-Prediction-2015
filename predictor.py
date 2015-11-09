import csv, os
import getopt
from numpy import random

#def main():
#    try:
#        args = getopt.getopt(sys.args[1:],'s:', ["file1", "file2"])
        

    

# Directory 
script_dir = os.getcwd() #<-- absolute dir the script is in
mono = "../Drug Synergy Data/Raw Data/Raw_Data_csv/ch1_ch2_monoTherapy"
comb = "../Drug Synergy Data/Raw Data/Raw_Data_csv/ch1_training_combinations"
abs_mono = os.path.join(script_dir, mono)
abs_comb = os.path.join(script_dir, comb)

# Training data file
path = os.path.join(script_dir, "../Drug Synergy Data")
file_lea = open(path+'/ch1_leaderBoard_monoTherapy.csv', "rt")
file_tra = open(path+'/ch1_train_combination_and_monoTherapy.csv', "rt")
#file_matrix = open(path+'/thisfile.csv', "rt")
data = csv.reader(file_tra, delimiter=' ', quotechar='|')
data_lea = csv.reader(file_lea, delimiter=' ', quotechar='|')

# Store data in dict
train_dict = {}
lea_dict = {}

#for word in  file_matrix:
#    print(word)

row_no = 1
for dat_a in data:
    values = dat_a[0].split(',')
    if (values[0][1:-1]!='CELL_LINE'): #and values[0]=='HCC1143'):
        row_no+=1
        train_dict[values[-1][1:-1]+"."+values[0][1:-1]] = [values[0][1:-1],
                                                            values[1][1:-1],values[2][1:-1],values[3][1:-1],
                                                            values[4][1:-1],values[5],values[6],values[7],
                                                            values[8],values[9],values[10],values[11],
                                                            values[12],values[13][1:-1],row_no]

row_no_lea = 1
b = 1
for dat_a_lea in data_lea:
    #print(dat_a_lea)
    if  dat_a_lea != []:
        values_lea = dat_a_lea[0].split(',')
    if (values_lea[0]!='CELL_LINE'): #and values_lea[0]=='HCC1143'):
        row_no_lea+=1
        lea_dict [values_lea[-1]+"."+values_lea[0]] = [values_lea[0],
                                                            values_lea[1],values_lea[2],values_lea[3],
                                                            values_lea[4],values_lea[5],values_lea[6],values_lea[7],
                                                            values_lea[8],values_lea[9],values_lea[10],values_lea[11],
                                                            values_lea[12],values_lea[13],row_no_lea]
        b = -1
    b = 1


def get_train(train_dict):
	#Populate train
	train = []
	for key in train_dict:
		train.append(([train_dict[key][0],train_dict[key][1],train_dict[key][2],train_dict[key][-4]]))
	return train

def get_target(lea_dict):
	#Populate target
	target = []
	for key in lea_dict:
		target.append(([lea_dict[key][0],lea_dict[key][1],lea_dict[key][2]]))
	return target

training = get_train(train_dict)
target = get_target(lea_dict)

def get_target_drugs(data):
    data_array = []
    for i in range(0, len(data)):    
        drug1 = data[i][0]
        drug2 = data[i][1]
        comb1 = [str(drug1), str(drug2)]
        comb2 = [str(drug2), str(drug1)]
        if comb1 and comb2 not in data_array:
            data_array.append(comb1)
    return data_array
    
def ddrug(a,b):
    return random.uniform(-1,1)

def dcell(a,b):
    return random.uniform(-1,1)

def get_top(psim_array):
    max_array = []
    for i in range(0, len(psim_array)):
        max_array.append(psim_array[i][3])
    return max(max_array)    

def get_psim(train, lea):
    top_array=[]
    for drg1 in target:
        psim_array = []
        for drg2 in training:
            psim = ddrug(drg1[1], drg2[1]) * ddrug(drg1[2], drg2[2]) * dcell(drg1[0], drg2[0])
            psim_array.append([drg1[0], drg1[1], drg1[2],   psim])
            
        top_array.append(get_top(psim_array))
    return top_array
    




