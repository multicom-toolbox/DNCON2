#!/usr/bin/python
# Badri Adhikari, 5-21-2017
# Predict using the five models in first stage and append to feature file

import sys
import os

from libcnnpredict import *

dir_config     = sys.argv[1]
file_weights   = sys.argv[2]
string_header  = sys.argv[3]
fileX          = sys.argv[4]
fileX_stage2   = sys.argv[5]

file_weights  = dir_config + '/' + file_weights
model_arch = read_model_arch(dir_config + '/model-arch.config')

print ''
print('SCRIPT        : ' + sys.argv[0]   )
print('dir_config    : ' + dir_config    )
print('file_weights  : ' + file_weights  )
print('string_header : ' + string_header )
print('fileX         : ' + fileX         )
print('fileX_stage2  : ' + fileX_stage2  )
print ''

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

# Predict at (L+10) x (L+10) and trim it back to L x L
P = make_prediction(model_arch, file_weights, X)
PF = ((P[0].reshape(LMAX, LMAX))[0:L, 0:L]).flatten()

# Append this prediction as an additional feature to the existing feature file
file_tmp_this_feature = '/tmp/X-%s.txt' % os.getpid()
try:
	np.savetxt(file_tmp_this_feature, PF, fmt='%1.5f', newline = ' ', header = '# Prediction at ' + string_header + '\n', comments = '')
	# remove leading spaces
	os.system("sed -i 's/^ *//g' " + file_tmp_this_feature)
	# add a newline, so that it can be appended to original file
	os.system('echo  >> ' + file_tmp_this_feature)
	# this will be done by 'coevo-60A'
	if not os.path.isfile(fileX_stage2):
		os.system('cp ' + fileX + ' ' + fileX_stage2)
	# append this feature
	os.system('cat ' + file_tmp_this_feature + ' >> ' + fileX_stage2)
finally:
	os.remove(file_tmp_this_feature)
