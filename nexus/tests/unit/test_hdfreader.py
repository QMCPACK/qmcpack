
from versions import h5py_available
import testing
from testing import value_eq,object_eq


def test_import():
    import hdfreader
    from hdfreader import HDFreader,read_hdf
#end def test_import


if h5py_available:
    def test_read():
        import os
        import numpy as np
        import h5py
        from hdfreader import read_hdf

        ds = 'string value'
        di = 100
        df = np.pi
        das = np.array(tuple('abcdefghijklmnopqrstuvwxyz'),dtype=bytes)
        dai = np.arange(20,dtype=np.int64)
        daf = 0.1*np.arange(20,dtype=np.float64)

        path = testing.setup_unit_test_output_directory('hdfreader','test_read')

        def add_datasets(g):
            g.create_dataset('sdata',data=das)
            g.create_dataset('idata',data=dai)
            g.create_dataset('fdata',data=daf)
        #end def add_datasets

        def add_attrs(g):
            g.attrs['sval'] = ds
            g.attrs['ival'] = di
            g.attrs['fval'] = df
        #end def add_attrs

        def add_group(g,label=''):
            g = g.create_group('group'+str(label))
            add_attrs(g)
            add_datasets(g)
            return g
        #end def add_group

        testfile = os.path.join(path,'test.h5')
        f = h5py.File(testfile,'w')

        add_datasets(f)
        g1  = add_group(f,1)
        g2  = add_group(f,2)
        g11 = add_group(g1,1)
        g12 = add_group(g1,2)
        g21 = add_group(g2,1)
        g22 = add_group(g2,2)

        f.close()

        def check_datasets(g):
            assert(value_eq(g.sdata,das))
            assert(value_eq(g.idata,dai))
            assert(value_eq(g.fdata,daf))
        #end def check_datasets

        def check_groups(g):
            assert('group1' in g)
            assert('group2' in g)
            check_datasets(g.group1)
            check_datasets(g.group2)
        #end def check_groups

        h = read_hdf(testfile)

        check_datasets(h)
        check_groups(h)
        check_groups(h.group1)
        check_groups(h.group2)
    #end def test_read
#end if
