""" Extract test_set from hdf5 file """
import h5py
import numpy as np
Hdffile = '1CPM_islets.h5'
Fullread = h5py.File(Hdffile, 'r')

print("Keys: %s" % Fullread.keys())

# The test act file is 'test_in'
hf = h5py.File('test_in.h5', 'w')
hf.create_dataset('test_seq', data = Fullread['test_in'])
hf.close()
