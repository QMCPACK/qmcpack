

def test_imports():
    from superstring import string2val
    from superstring import split_delims
    from superstring import valid_variable_name
    from superstring import find_matching_pair
    from superstring import remove_pair_sections
    from superstring import remove_empty_lines
#end def test_imports



def test_string2val():
    import numpy as np
    from testing import value_eq
    from superstring import string2val

    v = string2val('False')
    assert(isinstance(v,bool))
    assert(v==False)

    v = string2val('True')
    assert(isinstance(v,bool))
    assert(v==True)

    v = string2val('43')
    assert(isinstance(v,int))
    assert(v==43)

    v = string2val('3.14')
    assert(isinstance(v,float))
    assert(np.abs(v-3.14)<1e-6)

    v = string2val('1 2 3 4',delim=' ')
    assert(isinstance(v,np.ndarray))
    assert(value_eq(v,np.array([1,2,3,4],dtype=int)))

    v = string2val('1, 2, 3, 4',delim=',')
    assert(isinstance(v,np.ndarray))
    assert(value_eq(v,np.array([1,2,3,4],dtype=int)))

    v = string2val('1. 2 3 4',delim=' ')
    assert(isinstance(v,np.ndarray))
    assert(value_eq(v,np.array([1,2,3,4],dtype=float)))

    v = string2val('1., 2, 3, 4',delim=',')
    assert(isinstance(v,np.ndarray))
    assert(value_eq(v,np.array([1,2,3,4],dtype=float)))

#end def test_string2val



def test_split_delims():
    from superstring import split_delims

    s = '  ..a-._bc-d_.efg_.--h'
    assert(split_delims(s)==['a','bc','d','efg','h'])
#end def test_split_delims



def test_valid_variable_name():
    from superstring import valid_variable_name

    assert(valid_variable_name('valid_variable_name'))
    assert(valid_variable_name('_valid_variable_name'))
    assert(valid_variable_name('__valid_variable_name'))
    assert(valid_variable_name('valid_variable_name__'))
    assert(valid_variable_name('validVariableName'))
    assert(not valid_variable_name('invalid variable name'))
    assert(not valid_variable_name('invalid_variable-name'))
    assert(not valid_variable_name('inval!d_variable_name'))
    assert(not valid_variable_name('inv@lid_variable_name'))
    assert(not valid_variable_name('in>alid_variable_name'))
    assert(not valid_variable_name('invalid_variable_name '))

#end def test_valid_variable_name



def test_find_matching_pair():
    from superstring import find_matching_pair

    s = '''
        <simulation>
           <project id="vmc_hf_noj" series="0">
              <application name="qmcapp" role="molecu" class="serial" version="1.0"/>
           </project>
           <qmcsystem>
              <simulationcell>
                ...
              </simulationcell>
              <particleset>
                ...
              </particleset>
              <wavefunction>
                ...
              </wavefunction>
              <hamiltonian>
                ...
              </hamiltonian>
           </qmcsystem>
           <qmc method="vmc" move="pbyp">
              ...
           </qmc>
        </simulation>
        '''

    i1,i2 = find_matching_pair(s,'<>')
    assert(s[i1:i2]=='<simulation>')
    i1,i2 = find_matching_pair(s,'<>',i2)
    assert(s[i1:i2]=='<project id="vmc_hf_noj" series="0">')
    
    i1,i2 = find_matching_pair(s,('</','>'))
    assert(s[i1:i2]=='</project>')
    i1,i2 = find_matching_pair(s,('</','>'),i2)
    assert(s[i1:i2]=='</simulationcell>')
    i1,i2 = find_matching_pair(s,('</','>'),i2)
    assert(s[i1:i2]=='</particleset>')
    i1,i2 = find_matching_pair(s,('</','>'),i2)
    assert(s[i1:i2]=='</wavefunction>')
    i1,i2 = find_matching_pair(s,('</','>'),i2)
    assert(s[i1:i2]=='</hamiltonian>')
    i1,i2 = find_matching_pair(s,('</','>'),i2)
    assert(s[i1:i2]=='</qmcsystem>')
    i1,i2 = find_matching_pair(s,('</','>'),i2)
    assert(s[i1:i2]=='</qmc>')
    i1,i2 = find_matching_pair(s,('</','>'),i2)
    assert(s[i1:i2]=='</simulation>')

#end def test_find_matching_pair



def test_remove_pair_sections():
    from superstring import remove_pair_sections

    s = '''
        <simulation>
           <project id="vmc_hf_noj" series="0">
              <application name="qmcapp" role="molecu" class="serial" version="1.0"/>
           </project>
           <qmcsystem>
              <simulationcell>
                ...
              </simulationcell>
              <particleset>
                ...
              </particleset>
              <wavefunction>
                ...
              </wavefunction>
              <hamiltonian>
                ...
              </hamiltonian>
           </qmcsystem>
           <qmc method="vmc" move="pbyp">
              ...
           </qmc>
        </simulation>
        '''

    s = remove_pair_sections(s,('<hamiltonian>','</hamiltonian>'))

    s_no_h = '''
        <simulation>
           <project id="vmc_hf_noj" series="0">
              <application name="qmcapp" role="molecu" class="serial" version="1.0"/>
           </project>
           <qmcsystem>
              <simulationcell>
                ...
              </simulationcell>
              <particleset>
                ...
              </particleset>
              <wavefunction>
                ...
              </wavefunction>
              
           </qmcsystem>
           <qmc method="vmc" move="pbyp">
              ...
           </qmc>
        </simulation>
        '''
    assert(s==s_no_h)


    s = remove_pair_sections(s,('<qmcsystem>','</qmcsystem>'))

    s_no_sys = '''
        <simulation>
           <project id="vmc_hf_noj" series="0">
              <application name="qmcapp" role="molecu" class="serial" version="1.0"/>
           </project>
           
           <qmc method="vmc" move="pbyp">
              ...
           </qmc>
        </simulation>
        '''
    assert(s==s_no_sys)
    

    s = remove_pair_sections(s,('<simulation>','</simulation>'))

    assert(s.strip()=='')

#end def test_remove_pair_sections



def test_remove_empty_lines():
    from superstring import remove_empty_lines


    s = '''
        
        This string


        has a

        number of



        empty lines.
        '''

    sref = '''        This string
        has a
        number of
        empty lines.
'''


    assert(remove_empty_lines(s)==sref)

#end def test_remove_empty_lines
