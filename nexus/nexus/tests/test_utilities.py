import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.UTILITIES)

import os
from pathlib import Path



def test_is_valid_path():
    from nexus.utilities import is_valid_path

    valid_paths = [
        '/home/me/file.txt',
        'back/../and/../forth/',
        '/tmp/session_data',
        './runs/qmc/vmc.s000.scalar.dat',
        'hyphens-are-healthy',
        'periods.produce.paradise',
        'underscores_unlock_understanding',
        'pi_is_about_3.14159265358979323846264338327950'
        ]

    invalid_paths = [
        'yes!',
        '@somewhere.com',
        '#comment',
        '$100',
        '20%',
        '2^3',
        'here&there',
        'height*width',
        'definitely;maybe',
        'this|that',
        'why?',
        's\ash',
        '`command`',
        "ain't",
        '"special"',
        'x,y',
        '(pa',
        'ren)',
        '[brac',
        'ket]',
        'lt<',
        'gt>',
        'spacious ',
        'pull\tab',
        '\newline',
        'ca\r\riage',
        'Do not go where the path may lead,' 
        ' go instead where there is no path'
        ' and leave a trail.',
        ]

    for path in valid_paths:
        assert is_valid_path(path)

    for path in invalid_paths:
        assert not is_valid_path(path)

#end def test_is_valid_path



def test_path_string():
    from nexus.utilities import path_string

    # path_string for str paths
    in_out_paths = [
        # in         out ps(str) out ps(Path)    
        ( ''        , ''        , '.'       ),
        ( '.'       , '.'       , '.'       ),
        ( './'      , './'      , '.'       ),
        ( './..'    , './..'    , '..'      ),
        ( 'a/b'     , 'a/b'     , 'a/b'     ),
        ( 'a/b/'    , 'a/b/'    , 'a/b'     ),
        ( './a/b'   , './a/b'   , 'a/b'     ),
        ( './a/b/'  , './a/b/'  , 'a/b'     ),
        ( './a/./b' , './a/./b' , 'a/b'     ),
        ( '../a/..' , '../a/..' , '../a/..' ),
        ( 'a/../b'  , 'a/../b'  , 'a/../b'  ),
        ( '/'       , '/'       , '/'       ),
        ( '//'      , '//'      , '//'      ),
        ( '///'     , '///'     , '/'       ),
        ( '/a/'     , '/a/'     , '/a'      ),
        ( '/a/b/'   , '/a/b/'   , '/a/b'    ),
        ]

    for p_in, p_str_out, p_Path_out in in_out_paths:
        p1 = path_string(p_in)
        p2 = path_string(Path(p_in))
        p3 = str(Path(p_in))
        assert p1==p_in
        assert p1==p_str_out
        assert p2==p_Path_out
        assert p3==p_Path_out
        assert os.path.realpath(p1)==os.path.realpath(p2)


    # test bytes
    in_out_paths = [
        # in         out ps(str))    
        ( b''        , ''        ),
        ( b'.'       , '.'       ),
        ( b'./'      , './'      ),
        ( b'./..'    , './..'    ),
        ( b'a/b'     , 'a/b'     ),
        ( b'a/b/'    , 'a/b/'    ),
        ( b'./a/b'   , './a/b'   ),
        ( b'./a/b/'  , './a/b/'  ),
        ( b'./a/./b' , './a/./b' ),
        ( b'../a/..' , '../a/..' ),
        ( b'a/../b'  , 'a/../b'  ),
        ( b'/'       , '/'       ),
        ( b'//'      , '//'      ),
        ( b'///'     , '///'     ),
        ( b'/a/'     , '/a/'     ),
        ( b'/a/b/'   , '/a/b/'   ),
        ]

    for p_in, p_str_out in in_out_paths:
        assert path_string(p_in)==p_str_out

#end def test_path_string



def test_is_relative_path():
    from nexus.utilities import is_relative_path

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

#end def test_is_relative_path



def test_path_dual_typing():
    '''Illustrate issues w/ mixing pathlib.Path with os.path'''

    # joining paths 
    #   os.system('ls') and os.system('./ls') are very different!
    assert os.path.join('','ls')                       == 'ls'   # right
    assert str(Path('')/Path('ls'))                    == 'ls'   # right
    assert os.path.join(Path(''),Path('ls'))           == './ls' # wrong!!
    assert os.path.join(str(Path('')),str(Path('ls'))) == './ls' # wrong!!


    # splitting paths
    assert os.path.split('dir/name')       == ('dir','name') # right
    assert os.path.split(Path('dir/name')) == ('dir','name') # right
    assert str(Path('dir/name').parent)    == 'dir'          # right
    assert Path('dir/name').name           == 'name'         # right

    assert os.path.split('name')       == ('','name') # right
    assert os.path.split(Path('name')) == ('','name') # right
    assert str(Path('name').parent)    == '.'         # wrong!!
    assert Path('name').name           == 'name'      # right

    assert os.path.split('dir/')       == ('dir',''   )  # right
    assert os.path.split(Path('dir/')) == (''   ,'dir')  # wrong!!
    assert str(Path('dir/').parent)    == '.'            # wrong!!
    assert Path('dir/').name           == 'dir'          # wrong!!

    #   or in short
    assert os.path.dirname( 'dir/name') == str(Path('dir/name').parent)
    assert os.path.basename('dir/name') ==     Path('dir/name').name
    
    assert os.path.dirname( 'name'    ) != str(Path('name').parent)
    assert os.path.basename('name'    ) ==     Path('name').name
    
    assert os.path.dirname( 'dir/'    ) != str(Path('dir/').parent)
    assert os.path.basename('dir/'    ) !=     Path('dir/').name


    # relative paths
    assert os.path.relpath('/home/me/dir','/home/me')        == 'dir'
    assert str(Path('/home/me/dir').relative_to('/home/me')) == 'dir'

    assert os.path.relpath('/home/me','/home/me/dir') == '..'
    try:
        Path('/home/me').relative_to('/home/me/dir')
        works = True
    except ValueError:
        works = False
    assert not works

    assert os.path.relpath('./up/..', '.'      ) == '.'
    assert os.path.relpath('.'      , './up/..') == '.'

    assert str(Path('./up/..').relative_to('.')) == 'up/..'
    try:
        Path('.').relative_to('./up/..')
        works = True
    except ValueError:
        works = False
    assert not works

#end def test_path_dual_typing
