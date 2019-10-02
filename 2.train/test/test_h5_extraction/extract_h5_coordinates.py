import h5py
Hdffile = '1CPM_islets.h5'
Fullread = h5py.File(Hdffile, 'r')
a = Fullread['test_headers']
with open('test_coord.txt', 'w') as f:
    for item in a:
        f.write("%s\n" % item)
