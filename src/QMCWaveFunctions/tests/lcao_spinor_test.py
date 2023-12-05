import h5py
import numpy as np
from scipy.special import sph_harm, factorial2

def write_h5_file():
    hf = h5py.File('lcao_spinor.h5','w')
    
    #atoms
    atoms = hf.create_group('atoms')
    nat = np.array([1])
    nsp = np.array([1])
    pos = np.array([[0.0,0.0,0.0]])
    ids = np.array([0])
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
    kpts.create_dataset("eigenset_0", data=np.array([[0.25],[0.75]]))
    kpts.create_dataset("eigenset_0_imag", data=np.array([[0.75],[0.25]]))
    kpts.create_dataset("eigenset_1", data=np.array([[-0.2],[0.8]]))
    kpts.create_dataset("eigenset_1_imag", data=np.array([[0.8],[-0.2]]))
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

    print("Position: {}".format(pos))
    print("Spin: {}".format(s))

    #gaussian basis function values
    g = cartGauss(2.5,0,0,0,0) #s function
    dr = 1e-6
    val   = g.val(pos)
    valpx = g.val(pos + np.array([dr,0,0]))
    valmx = g.val(pos - np.array([dr,0,0]))
    valpy = g.val(pos + np.array([0,dr,0]))
    valmy = g.val(pos - np.array([0,dr,0]))
    valpz = g.val(pos + np.array([0,0,dr]))
    valmz = g.val(pos - np.array([0,0,dr]))
    dx = (valpx - valmx) / (2*dr)
    dy = (valpy - valmy) / (2*dr)
    dz = (valpz - valmz) / (2*dr)
    ddx = (valpx - 2*val + valmx) / (dr*dr)
    ddy = (valpy - 2*val + valmy) / (dr*dr)
    ddz = (valpz - 2*val + valmz) / (dr*dr)
    lap = ddx+ddy+ddz

    print("Basis ")
    print("  Val : {}".format(val))
    print("  Grad: {}  {}  {}".format(dx,dy,dz))
    print("  Lap : {}".format(lap))

    #build spinor info
    upcoef = (0.25+0.75j)
    dncoef = (-0.2+0.8j)

    upval = upcoef*val
    updx = upcoef*dx
    updy = upcoef*dy
    updz = upcoef*dz
    uplap = upcoef*lap
    dnval = dncoef*val
    dndx = dncoef*dx
    dndy = dncoef*dy
    dndz = dncoef*dz
    dnlap = dncoef*lap

    spval = (cs + 1j*ss)*upval + (cs - 1j*ss)*dnval
    spdx  = (cs + 1j*ss)*updx  + (cs - 1j*ss)*dndx
    spdy  = (cs + 1j*ss)*updy  + (cs - 1j*ss)*dndy
    spdz  = (cs + 1j*ss)*updz  + (cs - 1j*ss)*dndz
    splap = (cs + 1j*ss)*uplap + (cs - 1j*ss)*dnlap
    spds  = (-ss + 1j*cs)*upval + (-ss - 1j*cs)*dnval

    print(" 1st Spinor:")
    print("  Val     : {}".format(spval))
    print("  Grad    : {} {} {}".format(spdx,spdy,spdz))
    print("  Lap     : {}".format(splap))
    print("  SpinGrad: {}".format(spds))

    upcoef = (0.75+0.25j)
    dncoef = (0.8-0.2j)

    upval = upcoef*val
    updx = upcoef*dx
    updy = upcoef*dy
    updz = upcoef*dz
    uplap = upcoef*lap
    dnval = dncoef*val
    dndx = dncoef*dx
    dndy = dncoef*dy
    dndz = dncoef*dz
    dnlap = dncoef*lap

    spval = (cs + 1j*ss)*upval + (cs - 1j*ss)*dnval
    spdx  = (cs + 1j*ss)*updx  + (cs - 1j*ss)*dndx
    spdy  = (cs + 1j*ss)*updy  + (cs - 1j*ss)*dndy
    spdz  = (cs + 1j*ss)*updz  + (cs - 1j*ss)*dndz
    splap = (cs + 1j*ss)*uplap + (cs - 1j*ss)*dnlap
    spds  = (-ss + 1j*cs)*upval + (-ss - 1j*cs)*dnval

    print(" 2nd Spinor:")
    print("  Val     : {}".format(spval))
    print("  Grad    : {} {} {}".format(spdx,spdy,spdz))
    print("  Lap     : {}".format(splap))
    print("  SpinGrad: {}".format(spds))

    print()
    print()

if __name__ == "__main__":
    write_h5_file()
    
    pos = np.array([0.1,-0.3, 1.7])
    s   = 0.6
    get_reference_values(pos,s)
    
    pos = np.array([-0.4,1.5,-0.2])
    s   = -1.3
    get_reference_values(pos,s)
