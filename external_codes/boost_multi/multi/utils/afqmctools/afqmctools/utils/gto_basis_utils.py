import numpy
import os.path

def extend_gto(fname,sym='X',expo=None):
    if os.path.exists(fname):
        base_gto = open(fname,"r").read()
        if expo is None:
            return base_gto
    else:
        base_gto = ''
    assert( len(expo) > 0 )
    assert( len(expo) <= 7 )
    lname = numpy.array(['S','P','D','F','G','H','I'])
    for i in range(len(expo)):
        base_gto += sym + "   " + lname[i] + " \n  {}  1.00\n".format(expo[i])
    return base_gto

def extend_gto_id(fname,sym,expo,ids):
    if os.path.exists(fname):
        base_gto = open(fname,"r").read()
        if expo is None:
            return base_gto
    else:
        base_gto = ''
    assert( len(expo) <= len(ids) )
    assert( len(expo) > 0 )
    for i in range(len(expo)):
        base_gto += sym + "   " + ids[i] + " \n  {}  1.00\n".format(expo[i])
    return base_gto

def default_basis_map(Lmax,atoms):
    assert( Lmax >= 2 )
    assert( Lmax <= 6 )
    basis_map = {}
    x_ = []
    def add_shell(L,bmap,x):
        labels = ['S','P','D','F','G','H','I']
        for i in range(0,L+1):
            bmap.append(labels[i])
        if L==2:
            x.append(0.5)
            x.append(0.5)
            x.append(0.5)
        elif L==3:
            x.append(2.0)
            x.append(2.0)
            x.append(2.0)
            x.append(0.5)
        elif L==4:
            x.append(4.0)
            x.append(4.0)
            x.append(4.0)
            x.append(2.0)
            x.append(0.5)
        elif L==5:
            x.append(7.0)
            x.append(7.0)
            x.append(7.0)
            x.append(4.0)
            x.append(2.0)
            x.append(0.5)
        elif L==6:
            x.append(9.0)
            x.append(9.0)
            x.append(9.0)
            x.append(7.0)
            x.append(4.0)
            x.append(2.0)
            x.append(0.5)
        return bmap,x
    for a in atoms:
        bmap = []
        for l_ in range(2,Lmax+1):
            bmap, x_ = add_shell(l_,bmap,x_)
        basis_map.update({a:bmap})
    return numpy.array(x_),basis_map

def double_basis_map(Lmax,atoms):
    assert( Lmax >= 2 )
    assert( Lmax <= 6 )
    basis_map = {}
    x_ = []
    def add_shell(L,bmap,x):
        labels = ['S','P','D','F','G','H','I']
        for i in range(0,L+1):
            bmap.append(labels[i])
            bmap.append(labels[i])
        if L==2:
            x.append(0.2) #S
            x.append(3.0) #S
            x.append(0.2) #P
            x.append(3.0) #P
            x.append(0.2) #D
            x.append(3.0) #D
        elif L==3:
            x.append(1.0) #S
            x.append(2.0) #S
            x.append(1.0) #P
            x.append(2.0) #P
            x.append(1.0) #D
            x.append(2.0) #D
            x.append(0.2) #F
            x.append(1.0) #F
        elif L==4:
            x.append(4.0) #S
            x.append(5.0) #S
            x.append(4.0) #P
            x.append(5.0) #P
            x.append(4.0) #D
            x.append(5.0) #D
            x.append(2.0) #F
            x.append(3.0) #F
            x.append(0.2) #G
            x.append(1.0) #G
        elif L==5:
            x.append(7.0) #S
            x.append(8.0) #S
            x.append(7.0) #P
            x.append(8.0) #P
            x.append(7.0) #D
            x.append(8.0) #D
            x.append(4.0) #F
            x.append(5.0) #F
            x.append(2.0) #G
            x.append(3.0) #G
            x.append(0.2) #H
            x.append(1.0) #H
        elif L==6:
            x.append(9.0) #S
            x.append(10.0) #S
            x.append(9.0) #P
            x.append(10.0) #P
            x.append(9.0) #D
            x.append(10.0) #D
            x.append(7.0) #F
            x.append(8.0) #F
            x.append(4.0) #G
            x.append(5.0) #G
            x.append(2.0) #H
            x.append(3.0) #H
            x.append(0.2) #I
            x.append(1.0) #I
        return bmap,x
    for a in atoms:
        bmap = []
        for l_ in range(2,Lmax+1):
            bmap, x_ = add_shell(l_,bmap,x_)
        basis_map.update({a:bmap})
    return numpy.array(x_),basis_map
