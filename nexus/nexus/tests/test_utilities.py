import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.UTILITIES)

import os
from pathlib import Path


def test_relative_path():
    from nexus.utilities import relative_path

    for path_type in (str,Path):
        p1 = '/home/me'
        p2 = '/home/me/tmp'
        p1,p2 = [path_type(p) for p in (p1,p2)]
        assert relative_path(p1,p2)=='tmp'
        assert relative_path(p2,p1)=='..'

        p1 = './this/is/a/good/place'
        p2 = './this/is/../is/a/bad/place'
        p1,p2 = [path_type(p) for p in (p1,p2)]
        assert relative_path(p1,p2)=='../../bad/place'
        assert relative_path(p2,p1)=='../../good/place'

#end def test_relative_path





def test_is_relative_path():
    from nexus.utilities import is_relative_path,relative_path

    relative_paths = [
        '',
        '.',
        './',
        '..',
        '../',
        '../somewhere/..',
        './somewhere/above/',
        '../somewhere/below',
        'a/bald/or/bare/path',
        ]

    non_relative_paths = [
        '/',
        '/home',
        '/home/me',
        '/usr/bin/',
        '/lustre/proj/scratch',
        ] + [os.path.realpath(p) for p in relative_paths]

    for path_type in (str,Path):
        for p in relative_paths:
            p = path_type(p)
            print([p,str(p)])
            assert is_relative_path(p)

        for p in non_relative_paths:
            p = path_type(p)
            assert not is_relative_path(p)

    p1 = '/home/me'
    p2 = '/home/me/tmp'
    assert is_relative_path(relative_path(p1,p2))
    assert is_relative_path(relative_path(p2,p1))

    p1 = './this/is/a/good/place'
    p2 = './this/is/../is/a/bad/place'
    assert is_relative_path(relative_path(p1,p2))
    assert is_relative_path(relative_path(p2,p1))

#end def test_is_relative_path



def test_path_string():
    from nexus.utilities import path_string,path_object

    # path_string for str paths
    in_out_paths = [
        # in         out ps(str) out ps(Path)    
        ( ''        , ''        , '.'       ),
        ( '.'       , '.'       , '.'       ),
        ( './'      , '.'       , '.'       ),
        ( './..'    , './..'    , '..'      ),
        ( 'a/b'     , 'a/b'     , 'a/b'     ),
        ( 'a/b/'    , 'a/b'     , 'a/b'     ),
        ( './a/b'   , './a/b'   , 'a/b'     ),
        ( './a/b/'  , './a/b'   , 'a/b'     ),
        ( './a/./b' , './a/./b' , 'a/b'     ),
        ( '../a/..' , '../a/..' , '../a/..' ),
        ( 'a/../b'  , 'a/../b'  , 'a/../b'  ),
        ( '/'       , '/'       , '/'       ),
        ( '//'      , '//'      , '//'      ),
        ( '///'     , '///'     , '/'       ),
        ( '/a/'     , '/a'      , '/a'      ),
        ( '/a/b/'   , '/a/b'    , '/a/b'    ),
        ]

    for p_in, p_str_out, p_Path_out in in_out_paths:
        p1 = path_string(p_in)
        p2 = path_string(Path(p_in))
        p3 = str(Path(p_in))
        assert p1==p_str_out
        assert p2==p_Path_out
        assert p3==p_Path_out
        assert os.path.realpath(p1)==os.path.realpath(p2)

#end def test_path_string



def test_path_object():
    from nexus.utilities import path_object,path_string

    # fix Path information destruction
    assert str(Path(''   )) == '.'
    assert str(Path('./a')) == 'a'

    p,leading = path_object('')
    assert p==Path('')
    assert path_string(p,leading) == ''

    p,leading = path_object('./a')
    assert p==Path('./a')
    assert path_string(p,leading) == './a'

#end def test_path_object
