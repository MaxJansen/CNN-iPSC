""" Extract test_set from hdf5 file """
import h5py
import numpy as np
Hdffile = '1CPM_islets.h5'
Fullread = h5py.File(Hdffile, 'r')

print("Keys: %s" % Fullread.keys())
a_group_key = list(Fullread.keys())[0]

act_file = list(Fullread.keys())[3]

# The test act file is 'test_out'
act_data = list(Fullread['test_out'])

np.savetxt('test_act.txt', act_data)
