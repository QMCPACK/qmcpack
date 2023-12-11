

def test_import():
    import versions
#end def test_import



def test_constants():
    from versions import nexus_version
    from versions import python_supported
    from versions import years_supported

    assert(len(nexus_version)==3)
    for n in nexus_version:
        assert(isinstance(n,int))
    #end for

    assert(python_supported=='python3')

    assert(isinstance(years_supported,int))
#end def test_constants



def test_time():
    from datetime import date
    from versions import time_ago

    long_ago = date(year=2017,month=9,day=1)

    years_ago = time_ago(long_ago)

    assert(years_ago>2-1e-6)
#end def test_time



def test_version_processing():
    from versions import process_version,version_to_string

    assert( (1,3,0) == process_version('1.3')     )
    assert( (1,3,0) == process_version(1,3)       )
    assert( (1,3,0) == process_version('1','3',0) )
    assert( (1,3,0) == process_version([1,3])     )
    assert( (1,0,0) == process_version(1)         )

    assert( '1.3.0' == version_to_string( (1,3,0) ) )
    assert( '1.3.0' == version_to_string( (1,3,0) ) )
    assert( '1.3.0' == version_to_string( (1,3,0) ) )
    assert( '1.3.0' == version_to_string( (1,3,0) ) )
    assert( '1.0.0' == version_to_string( (1,0,0) ) )

#end def test_version_processing



def test_versions_object():
    from versions import Versions,versions,raw_version_data
    from versions import currently_supported
    from versions import years_supported,time_ago
    from versions import check_versions

    # object construction successful
    assert(versions is not None)

    # object is singleton
    assert(len(Versions.instances)==1)
    assert(id(Versions.instances[0])==id(versions))

    # check integrity of raw version data
    deps = set('''
        python3
        numpy      
        scipy      
        h5py       
        matplotlib 
        pydot      
        spglib     
        seekpath   
        pycifrw    
        cif2cell
        '''.split())
    assert(set(raw_version_data.keys())==deps)

    # check integrity of object's internal data
    assert(versions.dependencies==deps)
    assert(set(versions.ordered_dependencies)==deps)
    req = versions.required_dependencies
    opt = versions.optional_dependencies
    assert( req | opt == deps)
    assert( req & opt == set())
    assert(set(versions.currently_supported.keys())==deps)
    assert(set(versions.dependency_available.keys())==deps)
    assert(set(versions.dependency_version.keys())==deps)
    assert(set(versions.dependency_supported.keys())==deps)
    for name,version in currently_supported:
        versions.is_dependency(name)
        versions.check_dependency(name)
        assert(versions.supported(name,version))
        assert(isinstance(versions.available(name),bool))
        supp = versions.supported(name)
        assert(isinstance(supp,bool) or supp is None)
        ver = versions.dependency_version[name]
        assert(isinstance(ver,tuple) or ver is None)
        assert(name in versions.currently_supported)
        assert(version == versions.currently_supported[name])
    #end for

    # current supported versions obey policy
    tol = 1e-6
    for name,version in currently_supported:
        vdate = versions.version_data[name][version]['date']
        assert(time_ago(vdate)>years_supported-tol)
    #end for

    # policy derived versions obey policy
    tol = 1e-6
    policy_versions = versions.policy_supported_version()
    assert(set(policy_versions.keys())==deps)
    for name,version in policy_versions.items():
        vdate = versions.version_data[name][version]['date']
        assert(time_ago(vdate)>years_supported-tol)
    #end for

    # writing abilities
    assert(isinstance(versions.write_current_versions(),str))
    assert(isinstance(versions.write_policy_versions(),str))
    assert(isinstance(versions.write_available_versions(),str))

    # check on version validity and report
    assert(isinstance(versions.check(write=False,exit=False,full=False),bool))
    assert(isinstance(check_versions(write=False,exit=False),bool))
#end def test_versions_object
