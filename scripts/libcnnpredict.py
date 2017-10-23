#!/usr/bin/python
# Badri Adhikari, 6-15-2017
# Subroutines for prediction

from keras.models import Sequential
from keras.layers import Activation, Flatten
from keras.layers import Convolution2D
from keras.layers.normalization import BatchNormalization
from keras.optimizers import Nadam
import numpy as np
import math
import os
import sys
import random
import keras.backend as K
epsilon = K.epsilon()


# Model architectures / Layers information
def read_model_arch(file_config):
	if not os.path.isfile(file_config):
		print ('Error! Could not find config file ' + file_config)
		sys.exit(1)
	layers = {}
	with open(file_config) as f:
		for line in f:
			if line.startswith('#'):
				continue
			if len(line) < 2:
				continue
			cols = line.strip().split()
			if len(cols) != 5:
				print ('Error! Config file ' + file_config + ' line ' + line + '??')
				sys.exit(1)
			layers[cols[0]] = cols[1] + ' ' + cols[2] + ' ' +  cols[3] + ' ' + cols[4] 
	print ('')
	print ('Read model architecture:')
	for k, v in sorted(layers.items()):
		print (k + ' : ' + v)
	print ('')
	return layers

# Feature file that has 0D, 1D, and 2D features (L is the first feature)
# Output size (a little >= L) to which all the features will be rolled up to as 2D features
def getX(feature_file, l_max):
	# calcualte the length of the protein (the first feature)
	reject_list = []
	reject_list.append('# PSSM')
	reject_list.append('# AA composition')
	L = 0
	with open(feature_file) as f:
		for line in f:
			if line.startswith('#'):
				continue
			L = line.strip().split()
			L = int(round(math.exp(float(L[0]))))
			break
	Data = []
	with open(feature_file) as f:
		accept_flag = 1
		for line in f:
			if line.startswith('#'):
				if line.strip() in reject_list:
					accept_flag = 0
				else:
					accept_flag = 1
				continue
			if accept_flag == 0:
				continue
			if line.startswith('#'):
				continue
			this_line = line.strip().split()
			if len(this_line) == 0:
				continue
			if len(this_line) == 1:
				# 0D feature
				feature2D = np.zeros((L, L))
				feature2D[:, :] = float(this_line[0])
				Data.append(feature2D)
			elif len(this_line) == L:
				# 1D feature
				feature2D1 = np.zeros((L, L))
				feature2D2 = np.zeros((L, L))
				for i in range (0, L):
					feature2D1[i, :] = float(this_line[i])
					feature2D2[:, i] = float(this_line[i])
				Data.append(feature2D1)
				Data.append(feature2D2)
			elif len(this_line) == L * L:
				# 2D feature
				feature2D = np.asarray(this_line).reshape(L, L)
				Data.append(feature2D)
			else:
				print (line)
				print ('Error!! Unknown length of feature in !!' + feature_file)
				print ('Expected length 0, ' + str(L) + ', or ' + str (L*L) + ' - Found ' + str(len(this_line)))
				sys.exit()
	F = len(Data)
	X = np.zeros((l_max, l_max, F))
	for i in range (0, F):
		X[0:L, 0:L, i] = Data[i]
	return X

def build_model_for_this_input_shape(model_arch, X):
	model = Sequential()
	for layer in range(1, 1000):
		if not model_arch.has_key('layer' + str(layer)):
			break
		parameters = model_arch['layer' + str(layer)]
		cols = parameters.split()
		num_kernels = int(cols[0])
		filter_size = int(cols[1])
		b_norm_flag = cols[2]
		activ_funct = cols[3]
		if layer == 1:
			model.add(Convolution2D(num_kernels, filter_size, filter_size, border_mode='same', input_shape=X[0, :, :, :].shape))
		else:
			model.add(Convolution2D(num_kernels, filter_size, filter_size, border_mode='same'))
		if b_norm_flag == '1':
			model.add(BatchNormalization())
		model.add(Activation(activ_funct))
	model.add(Flatten())
	return model

def make_prediction(model_arch, file_weights, X):
	model = build_model_for_this_input_shape(model_arch, X)
	model.load_weights(file_weights)
	P = model.predict(X)
	return P

def print_feature_summary(X):
	print 'FeatID         Avg        Med        Max        Sum        Avg[30]    Med[30]    Max[30]    Sum[30]'
	for ii in range(0, len(X[0, 0, 0, :])):
		(m,s,a,d) = (X[0, :, :, ii].flatten().max(), X[0, :, :, ii].flatten().sum(), X[0, :, :, ii].flatten().mean(), np.median(X[0, :, :, ii].flatten()))
		(m30,s30,a30, d30) = (X[0, 30, :, ii].flatten().max(), X[0, 30, :, ii].flatten().sum(), X[0, 30, :, ii].flatten().mean(), np.median(X[0, 30, :, ii].flatten()))
		print ' Feat%2s %10.4f %10.4f %10.4f %10.1f     %10.4f %10.4f %10.4f %10.4f' %(ii, a, d, m, s, a30, d30, m30, s30)

def get_x_from_this_file(feature_file):
	L = 0
	with open(feature_file) as f:
		for line in f:
			if line.startswith('#'):
				continue
			L = line.strip().split()
			L = int(round(math.exp(float(L[0]))))
			break
	x = getX(feature_file, L)
	F = len(x[0, 0, :])
	X = np.zeros((1, L, L, F))
	X[0, :, :, :] = x
	return X

def prediction2rr(P, fileRR):
	print 'Writing RR file ' + fileRR
	L  = int(math.sqrt(len(P)))
	PM = P.reshape(L, L)
	rr = open(fileRR, 'w')
	for i in range(0, L):
		for j in range(i, L):
			if abs(i - j) < 1:
				continue
			rr.write("%i %i 0 8 %.5f\n" %(i+1, j+1, PM[i][j]))
	rr.close()

def make_ensemble_prediction(weight_arch_dict, X):
	N = len(X[:, 0, 0, 0])
	L = len(X[0, :, 0, 0])
	P = np.zeros((N, int(L * L)))
	for weight in weight_arch_dict.keys():
		print ''
		print 'Running prediction using ' + weight + ' and ' + weight_arch_dict[weight]
		P0 = make_prediction(read_model_arch(weight_arch_dict[weight]), weight, X)
		for i in range (0, len(P0[:, 0])):
			P[i] = P[i] + P0[i]
	for i in range (0, len(P[:, 0])):
		P[i] = P[i] / len(weight_arch_dict)
	return P

