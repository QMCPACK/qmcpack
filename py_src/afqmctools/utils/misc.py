import sys
import subprocess


class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


def get_git_hash():
    """ Return git revision.

    Adapted from:
        http://stackoverflow.com/questions/14989858/get-the-current-git-hash-in-a-python-script

    Returns
    -------
    sha1 : string
        git hash with -dirty appended if uncommitted changes.
    """

    try:
        src = [s for s in sys.path if 'afqmctools' in s][0]

        sha1 = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                       cwd=src).strip()
        suffix = subprocess.check_output(['git', 'status',
                                          '--porcelain',
                                          '../'],
                                         cwd=src).strip()
    except:
        suffix = False
        sha1 = 'none'.encode()
    if suffix:
        return sha1.decode('utf-8') + '-dirty'
    else:
        return sha1.decode('utf-8')
