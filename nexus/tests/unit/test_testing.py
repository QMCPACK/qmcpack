

def test_import():
    from testing import value_diff,object_diff
    from testing import value_eq,object_eq
    from testing import value_neq,object_neq
#end def test_import



def test_value_checks():
    import numpy as np
    from testing import value_diff,value_eq,value_neq

    assert(id(value_diff)==id(value_neq))


    shift = 1e-8

    class Special(object):
        def __init__(self,*args):
            self.internal = list(args)
        #end def __init__

        def __len__(self):
            return len(self.internal)
        #end def __len__
    #end class Special


    def deep_dict(use_shift=False,int_as_float=False,overshift=False):
        vals = np.array([1.23,3.14])
        if use_shift:
            vals += shift
        #end if
        if overshift:
            vals += 0.1
        #end if
        v1,v2 = vals
        ivals = np.array([1,2],dtype=int)
        if int_as_float:
            ivals = np.array(ivals,dtype=float)
            if use_shift:
                ivals -= shift
            #end if
        #end if
        i1,i2 = ivals
        dd = dict(
            some_key = 'a',
            a_list = [
                1,
                (1,i2,'c'),
                [v1,v2,v1,v2],
                dict(a=i1,b=2,c=v2),
                dict(a=dict(b=dict(),c=dict(d=v2,e=4,f='abc'))),
                ['a',['a',['a',[v1,v2],i1],i1],i1],
                ],
            b = v1,
            d1 = {1:v2,
                  2:dict(a=1,b=['a1',dict(a=i1,b=['a',dict(a=1,b=[v1,v2],c=v1,d=Special()),1],c=v2),1],c=v1),
                  3:(1,2,3),
                  4:np.array([1,2,3],dtype=float),
                  },
            s1 = set([3,True,'abc',3.14])
            )
        return dd
    #end def deep_dict


    def deep_list(*args,**kwargs):
        return list(deep_dict(*args,**kwargs).values())
    #end def deep_list


    # identical, reshaped, and real value shifted objects should agree
    # also empty objects of the same type are assumed to be identical
    should_agree = [
        # simple types
        ( 3     , 3     ),
        ( True  , True  ),
        ( False , False ),
        ( ''    , ''    ),
        ( 'a'   , 'a'   ),
        ( 'abc' , 'abc' ),
        ( 3.14  , 3.14  ),
        ( 3.14  , 3.14+shift  ),
        ( 3.14  , 3.14-shift  ),
        # tuples of simple types
        ( tuple()                   , tuple()                   ),
        ( (3,True,'abc',3.14)       , (3,True,'abc',3.14)       ),
        ( (3,True,'abc',3.14)       , (3,True,'abc',3.14+shift) ),
        ( (3,True,'abc',3.14)       , (3,True,'abc',3.14-shift) ),
        # lists of simple types
        ( []                        , []                        ),
        ( [3,True,'abc',3.14]       , [3,True,'abc',3.14]       ),
        ( [3,True,'abc',3.14]       , [3,True,'abc',3.14+shift] ),
        ( [3,True,'abc',3.14]       , [3,True,'abc',3.14-shift] ),
        # nested lists of simple types
        ( [1,0,0,1,1,0,1,1,1]       , [[1,0,0],[1,1,0],[1,1,1]] ),
        ( [1,0,0,1,1,0,1,1,1]       , [[[1],[0],[0]],[[1],[1],[0]],[[1],[1],[1]]] ),
        ( [1,0,0,1,1,0,1,1,1.]      , [[1,0,0],[1,1,0],[1,1,1.+shift]] ),
        ( [1,0,0,1,1,0,1,1,1]       , [(1,0,0),(1,1,0),(1,1,1)] ),
        # arrays of simple types
        ( np.array([])              , np.array([])                    ),
        ( np.arange(10,dtype=int)   , np.arange(10,dtype=int)         ),
        ( np.arange(10,dtype=float) , np.arange(10,dtype=float)       ),
        ( np.arange(10,dtype=float) , np.arange(10,dtype=float)+shift ),
        ( np.arange(10,dtype=float) , np.arange(10,dtype=float)-shift ),
        ( np.eye(3,3,dtype=int)     , np.eye(3,3,dtype=int).ravel()   ),
        ( np.eye(3,3,dtype=float)   , np.eye(3,3,dtype=float).ravel() ),
        ( np.array([1,2,3.])        , np.array([1,2.,3])              ),
        ( np.array([1,2,3.])        , np.array([1,2,3],dtype=float)   ),
        ( np.array(tuple('abc'))    , np.array(tuple('abc'))          ),
        # sets of simple types
        ( set()                     , set()                     ),
        ( set([1,2,3])              , set([1,2,3])              ),
        ( set([3,True,'abc',3.14])  , set([3,True,'abc',3.14])  ),
        # nested sets of simple types
        ( set(['abc',(1,2,3)])      , set(['abc',(1,2,3)])      ),
        # dicts of simple types
        ( dict()                    , dict()                    ),
        ( dict(a=3,b=True,c='abc',d=3.14) , dict(a=3,b=True,c='abc',d=3.14) ),
        ( dict(a=3,b=True,c='abc',d=3.14) , dict(a=3,b=True,c='abc',d=3.14+shift) ),
        ( { 3:'a',True:'b','abc':'c',3.14:'d'} , { 3:'a',True:'b','abc':'c',3.14:'d'} ),
        # unknown types with zero length
        ( Special(), Special() ),
        # deeply nested types
        ( deep_list() , deep_list() ),
        ( deep_list() , deep_list(use_shift=True) ),
        ( deep_dict() , deep_dict() ),
        ( deep_dict() , deep_dict(use_shift=True) ),
        ]

    for v1,v2 in should_agree:
        assert(value_eq(v1,v2))
        assert(value_eq(v2,v1))
        assert(not value_neq(v1,v2))
        assert(not value_neq(v2,v1))
    #end for

    
    # checks using integer as float
    should_agree += [
        ( 3         , 3         ),
        ( 3         , 3.0       ),
        ( 3         , 3.0+shift ),
        ( 3         , 3.0-shift ),
        ( [1,2.,3]  , [1.,2,3.] ),
        ( np.array([1,2,3],dtype=int), np.array([1,2,3],dtype=float)       ),
        ( np.array([1,2,3],dtype=int), np.array([1,2,3],dtype=float)+shift ),
        ( np.array([1,2,3],dtype=int), np.array([1,2,3],dtype=float)-shift ),
        ( deep_list() , deep_list(int_as_float=True,) ),
        ( deep_list() , deep_list(int_as_float=True,use_shift=True) ),
        ( deep_dict() , deep_dict(int_as_float=True,) ),
        ( deep_dict() , deep_dict(int_as_float=True,use_shift=True) ),        
        ]

    for v1,v2 in should_agree:
        assert(value_eq(v1,v2,int_as_float=True))
        assert(value_eq(v2,v1,int_as_float=True))
        assert(not value_neq(v1,v2,int_as_float=True))
        assert(not value_neq(v2,v1,int_as_float=True))
    #end for


    # differing types should not agree
    values = [
        1,
        True,
        '1',
        1.0,
        (1,2,3),
        [1,2,3],
        np.array([1,2,3]),
        set([1,2,3]),
        {1:1,2:2,3:3},
        Special(),
        ]
    for i1,v1 in enumerate(values):
        for i2,v2 in enumerate(values):
            if i1!=i2:
                assert(not value_eq(v1,v2))
                assert(value_neq(v1,v2))
            #end if
        #end for
    #end for


    # differing values should not agree
    dr = 2e-6
    values_full = [
        # integers
        [-3,-2,-1,0,1,2,3],
        # floats
        [1.-3*dr,1.-2*dr,1.-dr,1.,1.+dr,1.+2*dr,1.+3*dr],
        # booleans
        [True,False],
        # strings
        ['','a','ab','abc'],
        # tuples
        [tuple(),(1,2),(1,3),(1,'a'),(1,2.0),(1,2,3),(1,2,3,4)],
        # lists
        [[],[1,2],[1,3],[1,'a'],[1,2.0],[1,2,3],[1,2,3,4]],
        # arrays
        [np.array([]),np.array([1,2]),np.array([1,3]),
         np.array([1,'a'],dtype=object),
         np.array([1,2,3]),np.array([1,2,3,4]),
         np.array([1.,2]),np.array([1,3.]),
         np.array([1.,'a'],dtype=object),
         np.array([1,2.,3]),np.array([1,2,3.,4]),
         ],
        # sets
        [set(),set([1]),set([1,2]),set([1,3]),set([1,'a']),
         set([1,2.0+shift]),set([1,2,3]),set([1,2,3,4])],
        # dicts
        [dict(),dict(a=1),dict(a=(1,2)),dict(a=1,b=2),dict(a=1,b=2.0)],
        # unknown types with differing lengths
        [Special(),Special(1),Special(1,2),Special(1,3),Special(1,2,3)],
        # deeply nested types with differing values inside
        [deep_dict(),deep_dict(int_as_float=True),deep_dict(overshift=True)],
        [deep_list(),deep_list(int_as_float=True),deep_list(overshift=True)],
        ]
    for values in values_full:
        for i1,v1 in enumerate(values):
            for i2,v2 in enumerate(values):
                if i1!=i2:
                    assert(not value_eq(v1,v2))
                    assert(value_neq(v1,v2))
                #end if
            #end for
        #end for
    #end for
    
#end def test_value_checks



def test_object_checks():
    import numpy as np
    from testing import object_diff,object_eq,object_neq

    assert(id(object_diff)==id(object_neq))


    shift = 1e-8

    distinct = dict()

    def check_agree(o1,o2=None):
        if id(o) not in distinct:
            distinct[id(o)] = o
        #end if
        if o2 is None:
            o2 = o.copy()
        #end if
        assert(object_eq(o1,o2,bypass=True))
        assert(not object_neq(o1,o2,bypass=True))
    #end def check_agree

    def check_disagree(o1,o2):
        assert(not object_eq(o1,o2,bypass=True))
        assert(object_neq(o1,o2,bypass=True))
    #end def check_agree


    o = dict()
    check_agree(o,o)
    o2 = dict()
    check_agree(o,o2)

    o = dict(a=1)
    check_agree(o)

    o = dict(a=1,b=2)
    check_agree(o)

    o = dict(a=1,b=2,c=3)
    check_agree(o)

    o = {
        'a/b/c' : 1.0,
        'a/e'   : (1,2,3),
        }
    o2 = {
        'a/b/c' : 1.0+shift,
        'a/e'   : (1,2,3),
        }
    check_agree(o,o2)

    o = {
        'a/b/c' : 1.0,
        'a/e'   : [1,2,3.],
        }
    o2 = {
        'a/b/c' : 1.0+shift,
        'a/e'   : [1,2,3.+shift],
        }
    check_agree(o,o2)

    distinct = list(distinct.values())
    for i1,o1 in enumerate(distinct):
        for i2,o2 in enumerate(distinct):
            if i1!=i2 and False:
                check_disagree(o1,o2)
            #end if
        #end for
    #end for

#end def test_object_checks
