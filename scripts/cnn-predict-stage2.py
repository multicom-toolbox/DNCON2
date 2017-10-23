#!/usr/bin/python
# Badri Adhikari, 5-21-2017
# Make the prediction at stage2 with second-stage feature files

import sys
import os

from libcnnpredict import *
 
dir_config = sys.argv[1]
fileX      = sys.argv[2]
fileRR     = sys.argv[3]
 
print ''
print('SCRIPT     : ' + sys.argv[0])
print('dir_config : ' + dir_config )
print('fileX      : ' + fileX      )
print('fileRR     : ' + fileRR  )
print ''
 
weight_arch_dict = {}
for i in range(1, 21):
	print 'Reading weight file ' + dir_config + '/stage2-' + str(i) + '.hdf5 ..'
	weight_arch_dict[dir_config + '/stage2-' + str(i) + '.hdf5'] = dir_config + '/model-arch.config'

# Need to make X slightly bigger than L x L by padding zeros
# Building a model with L x L decreases performance
L = 0
with open(fileX) as f:
	for line in f:
		if line.startswith('#'):
			continue
		L = line.strip().split()
		L = int(round(math.exp(float(L[0]))))
		break

LMAX = L + 10
x = getX(fileX, LMAX)
F = len(x[0, 0, :])
X = np.zeros((1, LMAX, LMAX, F))
X[0, :, :, :] = x
print ''
print_feature_summary(X)

# Predict at (L+10) x (L+10) and trim it back to L x L
print ''
print 'Starting ensemble prediction.. '
P = make_ensemble_prediction(weight_arch_dict, X)
PF = ((P[0].reshape(LMAX, LMAX))[0:L, 0:L]).flatten()

# Convert prediction to CASP RR file without header
prediction2rr(PF, fileRR)
