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



bare_relative_paths = [
    '',
    'a/bald/or/bare/path',
    ]

relative_paths = [
    '.',
    './',
    '..',
    '../',
    '../somewhere/..',
    './somewhere/above/',
    '../somewhere/below',
    ] + bare_relative_paths

non_relative_paths = [
    '/',
    '/home',
    '/home/me',
    '/usr/bin/',
    '/lustre/proj/scratch',
    ] + [os.path.realpath(p) for p in relative_paths]



def test_is_relative_path():
    from nexus.utilities import is_relative_path,relative_path

    for path_type in (str,Path):
        for p in relative_paths:
            p = path_type(p)
            assert is_relative_path(p)

        for p in bare_relative_paths:
            if path_type==Path and p=='':
                # str(Path(''))=='.'
                continue
            p = path_type(p)
            is_rel,is_bare = is_relative_path(p,ret_bare=True)
            assert is_rel
            assert is_bare

        for p in non_relative_paths:
            p = path_type(p)
            assert not is_relative_path(p)

        for p in non_relative_paths:
            p = path_type(p)
            is_rel,is_bare = is_relative_path(p,ret_bare=True)
            assert not is_rel
            assert not is_bare

    p1 = '/home/me'
    p2 = '/home/me/tmp'
    assert is_relative_path(relative_path(p1,p2))
    assert is_relative_path(relative_path(p2,p1))

    p1 = './this/is/a/good/place'
    p2 = './this/is/../is/a/bad/place'
    assert is_relative_path(relative_path(p1,p2))
    assert is_relative_path(relative_path(p2,p1))

#end def test_is_relative_path



def test_path_object():
    from nexus.utilities import path_object

    # dual-type tests
    for path_type in (str,Path):
        # all paths
        for p in relative_paths+non_relative_paths:
            p = path_type(p)
            po = path_object(p)
            assert isinstance(po,Path)
            assert po==Path(p)

        # relative paths
        for p in relative_paths:
            p = path_type(p)
            p,is_rel = path_object(p,ret_rel=True)
            assert is_rel

        for p in relative_paths:
            # no exception is raised
            p = path_type(p)
            p = path_object(p,relative=True)

        for p in relative_paths:
            p = path_type(p)
            try:
                p = path_object(p,relative=False)
                expected = False
            except ValueError:
                expected = True
            assert expected


        # non-relative paths
        for p in non_relative_paths:
            p = path_type(p)
            p,is_rel = path_object(p,ret_rel=True)
            assert not is_rel

        for p in non_relative_paths:
            p = path_type(p)
            try:
                p = path_object(Path(p),strict=True)
                expected = False
            except ValueError:
                expected = True
            assert expected

        for p in non_relative_paths:
            # no exception is raised
            p = path_type(p)
            p = path_object(p,relative=False)

        for p in non_relative_paths:
            p = path_type(p)
            try:
                p = path_object(p,relative=True)
                expected = False
            except ValueError:
                expected = True
            assert expected

#end def test_path_object



def test_path_string():
    from nexus.utilities import path_string,path_object

    # path_string with str preserves trailing /
    assert path_string('')     == '.'
    assert path_string('.')    == '.'
    assert path_string('./')   == './'
    assert path_string('..')   == '..'
    assert path_string('../')  == '../'
    assert path_string('/')    == '/'
    assert path_string('//')   == '//'
    assert path_string('a')    == './a'
    assert path_string('./a/') == './a/'

    # path_string with Path does not preserve trailing /
    # this is because Path itself does not preserve trailing /
    assert path_string(Path(''))     == '.'
    assert path_string(Path('.'))    == '.'
    assert path_string(Path('./'))   == '.'
    assert path_string(Path('..'))   == '..'
    assert path_string(Path('../'))  == '..'
    assert path_string(Path('/'))    == '/'
    assert path_string(Path('//'))   == '//'
    assert path_string(Path('a'))    == './a'
    assert path_string(Path('./a/')) == './a'

    # leave non-bare path strings alone
    for p in relative_paths+non_relative_paths:
        if p in bare_relative_paths:
            continue
        assert path_string(p)==p

    # dual-type tests
    for path_type in (str,Path):
        for p in relative_paths+non_relative_paths:
            p = path_type(p)
            p = path_string(p)
            assert isinstance(p,str)

        # prepend ./ to bare paths
        for p in bare_relative_paths:
            pin = p
            p = path_type(p)
            ps = path_string(p)
            if pin=='':
                assert ps=='.'
            else:
                assert ps=='./'+pin

        # invert path_object operation
        for p in relative_paths+non_relative_paths:
            if p in bare_relative_paths:
                continue
            pin = p
            p = path_type(p)
            po = path_object(p)
            ps = path_string(po)
            if pin!='/':
                assert ps==pin.rstrip('/')
            else:
                assert ps==pin


        # relative paths
        for p in relative_paths:
            p = path_type(p)
            p,is_rel = path_string(p,ret_rel=True)
            assert is_rel

        for p in relative_paths:
            # no exception is raised
            p = path_type(p)
            p = path_string(p,relative=True)

        for p in relative_paths:
            p = path_type(p)
            try:
                p = path_string(p,relative=False)
                expected = False
            except ValueError:
                expected = True
            assert expected


        # non-relative paths
        for p in non_relative_paths:
            p = path_type(p)
            p,is_rel = path_string(p,ret_rel=True)
            assert not is_rel

        for p in non_relative_paths:
            # no exception is raised
            p = path_type(p)
            p = path_string(p,relative=False)

        for p in non_relative_paths:
            p = path_type(p)
            try:
                p = path_string(p,relative=True)
                expected = False
            except ValueError:
                expected = True
            assert expected

#end def test_path_string
