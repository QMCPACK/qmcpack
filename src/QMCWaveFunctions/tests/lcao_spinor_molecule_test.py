import h5py
import numpy as np
from scipy.special import sph_harm, factorial2

def write_h5_file():
    hf = h5py.File('lcao_spinor_molecule.h5','w')
    
    #atoms
    atoms = hf.create_group('atoms')
    nat = np.array([2])
    nsp = np.array([1])
    pos = np.array([[0.1,0.2,0.3],[-0.3,-0.2,-0.1]])
    ids = np.array([0,0])
    atoms.create_dataset('number_of_atoms', data=nat)
    atoms.create_dataset('number_of_species', data=nsp)
    atoms.create_dataset('positions', data=pos)
    atoms.create_dataset('species_ids', data=ids)
    sp = atoms.create_group('species_0')
    
    atnum  = np.array([1])
    charge = np.array([1])
    core   = np.array([1])
    name = "H"
    mylen = "S"+str(len(name))
    strList = [name]
    asciiList = [n.encode("ascii", "ignore") for n in strList]
    sp.create_dataset("atomic_number", data=atnum)
    sp.create_dataset("charge", data=charge)
    sp.create_dataset("core", data=core)
    sp.create_dataset("name", (1,), mylen, asciiList)
    
    #PBC
    pbc = hf.create_group("PBC")
    pbc.create_dataset("PBC",(1,), dtype="b1", data=False)
    
    #application
    app = hf.create_group("application")
    code = "generic"
    mylen = "S"+str(len(code))
    strList = [code]
    asciiList = [n.encode("ascii", "ignore") for n in strList]
    app.create_dataset("code",(1,), mylen, asciiList)
    
    #basisset
    bs = hf.create_group("basisset")
    bs.create_dataset("NbElements", data=np.array([1]))
    name="LCAOBSet"
    mylen="S"+str(len(name))
    strList=[name]
    asciiList=[n.encode("ascii","ignore") for n in strList]
    bs.create_dataset("name", (1,), mylen, asciiList)
    atbs = bs.create_group("atomicBasisSet0")
    
    atbs.create_dataset("NbBasisGroups", data=np.array([1]))
    mystr = "cartesian"
    mylen = "S"+str(len(mystr))
    strList = [mystr]
    asciiList = [n.encode("ascii","ignore") for n in strList]
    atbs.create_dataset("angular",(1,), mylen, asciiList)
    mystr = "H"
    mylen = "S"+str(len(mystr))
    strList = [mystr]
    asciiList = [n.encode("ascii","ignore") for n in strList]
    atbs.create_dataset("elementType",(1,), mylen, asciiList)
    mystr = "Gamess"
    mylen = "S"+str(len(mystr))
    strList = [mystr]
    asciiList = [n.encode("ascii","ignore") for n in strList]
    atbs.create_dataset("expandYlm",(1,), mylen, asciiList)
    atbs.create_dataset("grid_npts", data=np.array([1001]))
    atbs.create_dataset("grid_rf", data=np.array([100]))
    atbs.create_dataset("grid_ri", data=np.array([1e-06]))
    mystr = "log"
    mylen = "S"+str(len(mystr))
    strList = [mystr]
    asciiList = [n.encode("ascii","ignore") for n in strList]
    atbs.create_dataset("grid_type",(1,), mylen, asciiList)
    mystr = "Gaussian"
    mylen = "S"+str(len(mystr))
    strList = [mystr]
    asciiList = [n.encode("ascii","ignore") for n in strList]
    atbs.create_dataset("name",(1,), mylen, asciiList)
    mystr = "no"
    mylen = "S"+str(len(mystr))
    strList = [mystr]
    asciiList = [n.encode("ascii","ignore") for n in strList]
    atbs.create_dataset("normalized",(1,), mylen, asciiList)
    
    bg = atbs.create_group("basisGroup0")
    bg.create_dataset("NbRadFunc", data=np.array([1]))
    bg.create_dataset("l", data=np.array([0]))
    bg.create_dataset("n", data=np.array([0]))
    mystr = "H00"
    mylen = "S"+str(len(mystr))
    strList = [mystr]
    asciiList = [n.encode("ascii","ignore") for n in strList]
    bg.create_dataset("rid",(1,), mylen, asciiList)
    mystr = "Gaussian"
    mylen = "S"+str(len(mystr))
    strList = [mystr]
    asciiList = [n.encode("ascii","ignore") for n in strList]
    bg.create_dataset("type",(1,), mylen, asciiList)
    rf = bg.create_group("radfunctions")
    dr = rf.create_group("DataRad0")
    dr.create_dataset("contraction", data=np.array([1.0]))
    dr.create_dataset("exponent", data=np.array([2.5]))
    
    
    kpts = hf.create_group("Super_Twist")
    kpts.create_dataset("eigenset_0", data=np.array([[0.075, 0.15]]))
    kpts.create_dataset("eigenset_0_imag", data=np.array([[0.225, 0.45]]))
    kpts.create_dataset("eigenset_1", data=np.array([[-0.12, -0.06]]))
    kpts.create_dataset("eigenset_1_imag", data=np.array([[0.48, 0.24]]))
    hf.close()

class cartGauss:
    def __init__(self,expt,l=0,i=0,j=0,k=0):
        self.expt = expt
        self.l = l
        self.i = i
        self.j = j
        self.k = k
        assert(i+j+k == l)
    def norm(self):
        n = (2*self.expt / np.pi)**(3./4.)
        n *= np.sqrt(2.**(self.l) / factorial2(2*self.i - 1) / factorial2(2*self.j - 1) / factorial2(2*self.k - 1)) * np.sqrt(2*self.expt)**self.l
        return n
    def val(self,pos):
        r = np.linalg.norm(pos)
        norm = self.norm()
        return norm *pos[0]**self.i * pos[1]**self.j * pos[2]**self.k * np.exp(-self.expt * r * r)

def get_reference_values(pos, s):
    cs = np.cos(s)
    ss = np.sin(s)
    eis = cs + 1.j*ss
    emis = cs - 1.j*ss

    print("Position: {}".format(pos))
    print("Spin: {}".format(s))

    g0 = cartGauss(2.5, 0, 0, 0, 0)
    g1 = cartGauss(2.5, 0, 0, 0, 0)

    R0 = np.array([0.1,0.2,0.3])
    R1 = np.array([-0.3,-0.2,-0.1])

    c0 = 0.3
    c1 = 0.6

    upcoef = (0.25 + 0.75j)
    dncoef = (-0.2 + 0.8j)

    dr = 1e-7

    g0val = g0.val(pos-R0)
    g0px  = g0.val(pos-(R0 + np.array([dr,0,0])))
    g0mx  = g0.val(pos-(R0 - np.array([dr,0,0])))
    g0py  = g0.val(pos-(R0 + np.array([0,dr,0])))
    g0my  = g0.val(pos-(R0 - np.array([0,dr,0])))
    g0pz  = g0.val(pos-(R0 + np.array([0,0,dr])))
    g0mz  = g0.val(pos-(R0 - np.array([0,0,dr])))

    g1val = g1.val(pos-R1)
    g1px  = g1.val(pos-(R1 + np.array([dr,0,0])))
    g1mx  = g1.val(pos-(R1 - np.array([dr,0,0])))
    g1py  = g1.val(pos-(R1 + np.array([0,dr,0])))
    g1my  = g1.val(pos-(R1 - np.array([0,dr,0])))
    g1pz  = g1.val(pos-(R1 + np.array([0,0,dr])))
    g1mz  = g1.val(pos-(R1 - np.array([0,0,dr])))

    #atom 0
    uppx = c0*g0px + c1*g1val
    upmx = c0*g0mx + c1*g1val
    updx = (uppx - upmx) / (2*dr)
    dnpx = c1*g0px + c0*g1val
    dnmx = c1*g0mx + c0*g1val
    dndx = (dnpx - dnmx) / (2*dr)
    uppy = c0*g0py + c1*g1val
    upmy = c0*g0my + c1*g1val
    updy = (uppy - upmy) / (2*dr)
    dnpy = c1*g0py + c0*g1val
    dnmy = c1*g0my + c0*g1val
    dndy = (dnpy - dnmy) / (2*dr)
    uppz = c0*g0pz + c1*g1val
    upmz = c0*g0mz + c1*g1val
    updz = (uppz - upmz) / (2*dr)
    dnpz = c1*g0pz + c0*g1val
    dnmz = c1*g0mz + c0*g1val
    dndz = (dnpz - dnmz) / (2*dr)

    spdx = upcoef * updx * eis + dncoef * dndx * emis
    spdy = upcoef * updy * eis + dncoef * dndy * emis
    spdz = upcoef * updz * eis + dncoef * dndz * emis

    print("grad atom 0: {}, {}, {}".format(spdx, spdy, spdz))

    #atom 1
    uppx = c0*g0val + c1*g1px
    upmx = c0*g0val + c1*g1mx
    updx = (uppx - upmx) / (2*dr)
    dnpx = c1*g0val + c0*g1px
    dnmx = c1*g0val + c0*g1mx
    dndx = (dnpx - dnmx) / (2*dr)
    uppy = c0*g0val + c1*g1py
    upmy = c0*g0val + c1*g1my
    updy = (uppy - upmy) / (2*dr)
    dnpy = c1*g0val + c0*g1py
    dnmy = c1*g0val + c0*g1my
    dndy = (dnpy - dnmy) / (2*dr)
    uppz = c0*g0val + c1*g1pz
    upmz = c0*g0val + c1*g1mz
    updz = (uppz - upmz) / (2*dr)
    dnpz = c1*g0val + c0*g1pz
    dnmz = c1*g0val + c0*g1mz
    dndz = (dnpz - dnmz) / (2*dr)

    spdx = upcoef * updx * eis + dncoef * dndx * emis
    spdy = upcoef * updy * eis + dncoef * dndy * emis
    spdz = upcoef * updz * eis + dncoef * dndz * emis

    print("grad atom 1: {}, {}, {}".format(spdx, spdy, spdz))



if __name__ == "__main__":
   write_h5_file()
   pos = np.array([0.01, -0.02, 0.03])
   s = 0.6
   get_reference_values(pos, s)
    
