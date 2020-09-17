import numpy
import os.path


# Simple class 
class BasisBlock():
    """Objects of this class hold a list of exponents and contraction coefficients
       defining a block of a set angular momentum in the basis set. 
    """

    def __init__(self,sym,L,n_gauss=1,n_fun=1,exponents=None):
        """Constructor
    Parameters
    ----------
    sym: string.
        Symbol/label of species for this basis set block.
    L: integer.
        Angular momentum of the basis block.
    n_gauss: integer. Default: 1
        Number of gaussian basis functions in the basis set block.
    n_fun: integer. Default: 1
        Number of functions in the basis set block.
    exponents: Array of floating point numbers. Default: None
        If not None, must contain a list of n_gauss exponents.
        When exponents are provided, they are kept fixed and not included in the
        list of optimizable parameters.
        """
        
        if exponents is not None:
            assert(len(exponents) == n_gauss)
        assert( L == 'S' or L == 'P' or L == 'D' or L == 'F' or 
                L == 'G' or L == 'H' or L == 'I' )
        self.L = L
        self.sym = sym
        self.num_params = 0
        self.n_gauss = n_gauss
        self.n_fun = n_fun
        self.exponents = exponents
        self.fix_exponents = (exponents is not None)
        if n_gauss == 1: 
            assert(n_fun==1)
            if not self.fix_exponents: 
                self.num_params = 1
        else:        
            self.num_params = n_gauss * n_fun 
            if not self.fix_exponents: 
                self.num_params += n_gauss 

    def basis_str(self,x):
        """ Returns the string associated with this basis block.

    Parameters
    ----------
    x : array 
        Parameters of the basis set in linearized array. 
        The expected order is:
        [ exp1, c11, c12, ..., exp2, c21, c22, ..., exp3, c31, ...] 
        where expn is the exponent of the nth basis function and
        cni is the nth contraction coefficient of the ith basis. 

    Returns
    -------
    string ::
        Basis block. 
        """
        assert( len(x) == self.num_params )
        if self.n_gauss == 1:
            if self.fix_exponents:
                return self.sym + "   " + self.L + " \n  {}  1.00\n".format(self.exponents[0]) 
            else:
                return self.sym + "   " + self.L + " \n  {}  1.00\n".format(x[0]) 
        else:
            gto_basis = self.sym + "   " + self.L + " \n" 
            n=0
            for i in range(self.n_gauss):
                if self.fix_exponents:
                    gto_basis += "  {}  ".format(self.exponents[i]) 
                else:
                    gto_basis += "  {}  ".format(x[n]) 
                    n+=1
                for j in range(self.n_fun):
                    gto_basis += "  {}  ".format(x[n]) 
                    n+=1
                gto_basis += "\n"
# Simple class 
class EvenTemperedBasisBlock():
    """Objects of this class hold a list of exponents and contraction coefficients
       defining a block of a set angular momentum in the basis set. 
       Exponents are generated from even tempered set.
         expo[i] = alpha*beta**(i-1.0) 
       For now, only beta and contraction coefficients are mutable. 
    """

    def __init__(self,sym,L,alpha=0.1,n_gauss=5,n_fun=1,beta=None):
        """Constructor
    Parameters
    ----------
    sym: string.
        Symbol/label of species for this basis set block.
    L: integer.
        Angular momentum of the basis block.
    alpha: floating point. Default: 0.1
        Base in even temperted formula.
    n_gauss: integer. Default: 1
        Number of gaussian basis functions in the basis set block.
    n_fun: integer. Default: 1
        Number of functions in the basis set block.
    beta: floating point. Default: None
        Power factor in even temperted formula.
        If not None, it is fixed and not included in the list of parameters.
        """
        
        assert n_gauss > 1, "Use BasisBlock if n_gauss==1"
        assert alpha > 0.0, "Error: alpha < 0.0"
        assert( L == 'S' or L == 'P' or L == 'D' or L == 'F' or 
                L == 'G' or L == 'H' or L == 'I' )
        self.L = L
        self.sym = sym
        self.alpha = alpha
        self.fix_exponents = beta is not None 
        if beta is not None:
            assert beta > 1.0, "Error: beta <= 1.0" 
            self.beta = beta
        else:
            self.beta = 0.0 
        self.n_gauss = n_gauss
        self.n_fun = n_fun
        self.num_params = n_gauss * n_fun 
        if not self.fix_exponents: 
            self.num_params += 1 

    def basis_str(self,x):
        """ Returns the string associated with this basis block.

    Parameters
    ----------
    x : array 
        Parameters of the basis set in linearized array. 
        The expected order is:
        [ exp1, c11, c12, ..., exp2, c21, c22, ..., exp3, c31, ...] 
        where expn is the exponent of the nth basis function and
        cni is the nth contraction coefficient of the ith basis. 

    Returns
    -------
    string ::
        Basis block. 
        """
        assert( len(x) == self.num_params )
        gto_basis = self.sym + "   " + self.L + " \n" 
        if self.fix_exponents:
            beta = self.beta
            n=0
        else:
            beta = x[0]
            n=1
        for i in range(self.n_gauss):
            gto_basis += "  {}  ".format(self.alpha*beta**(i)) 
            for j in range(self.n_fun):
                gto_basis += "  {}  ".format(x[n]) 
                n+=1
            gto_basis += "\n"
        return gto_basis

class OptimizableBasisSet(): 

    def __init__(self, atom_id, def_basis={}):
        """
        Parameters
        ----------
        atom_id: array of strings
            Array with the species labels. Must be consistent and in the same order as the
            QE calculation.
        def_basis: Python Dictionary. Default: {}
            Default (unoptimizable) basis set. If a species symbol is defined in the dictionary,
            the mapped value in the dictionary must correspond to the loocation of a file with 
            a basis set for this element. The basis set in the file will be added to the basis 
            set of that element, but this part of the basis set is unmutable.
        """
        # names of the species
        assert(len(atom_id) > 0)
        self.atom_id = atom_id 

        # default_basis is a dictionary that contains the default basis
        # set (if it exists in def_basis) for each species.
        # If no default basis is provided it will contain an empty string
        # The actual basis for each species will be the default basis, followed
        # by any extensions included with "add".
        self.default_basis = {}
        for I,atm in enumerate(self.atom_id):
            assert( not (atm in self.default_basis) )
            if atm in def_basis:
                assert(os.path.exists(def_basis[atm]))
                self.default_basis.update({atm:def_basis[atm]})
            else:
                self.default_basis.update({atm:""})
            
        # Number of optimizable parameters 
        self.number_of_params = 0

        # Contains the bounds, within the parameter vector, of the parameters
        # associated with a given basis set block
        self.__p0 = {}
        self.__pN = {}

        # Contains the blocks that define the extended basis
        # each block is an element of the BasisBlock class.
        self.basis_blocks = {}  

        # setup empty structures
        for I,atm in enumerate(self.atom_id):
            self.basis_blocks.update({atm:[]})
            self.__p0.update({atm:[]})
            self.__pN.update({atm:[]})
        
            
    # Extends the basis set of species "atm" with an new basis block
    def add(self, bblk):
        """Extends the basis set of a species by adding a new (possibly optimizable)block.

        Parameters:
        ----------
        bblk: Object of type BasisBlock or EvenTemperedBasisBlock.
            A basis set block. 
        """
        assert( isinstance(bblk,BasisBlock) or isinstance(bblk,EvenTemperedBasisBlock) )
        atm = bblk.sym
        assert(atm in self.basis_blocks)
        assert(atm in self.__p0)
        assert(atm in self.__pN)
        atm_bblocks = self.basis_blocks[atm]
        atm_p0 = self.__p0[atm]
        atm_pN = self.__pN[atm]
        assert( len(atm_bblocks) == len(atm_p0) )
        assert( len(atm_p0) == len(atm_pN) )
        atm_bblocks.append(bblk)
        atm_p0.append(self.number_of_params)
        self.number_of_params += bblk.num_params
        atm_pN.append(self.number_of_params)

    # Returns the string for thhe extended part of the basis set of a species.
    #  atm: id of species
    #  x: parameter vector
    def extended_basis(self, atm, x):
        """Returns a string that corresponds to the extended part of the basis set for
           species atm. The extended part of the basis corresponds to the blocks
            included with self.add(...)

        Parameters
        ----------
            atm: string
                String with the label of a species.
            x: 1-D floating point array.
                Array with variational parameters in the basis.
                Dimension must be equal to self.number_of_params.
        Returns
        -------
            extra_basis: String            
                A string with the extended basis of species atm.
        """
        assert( atm in self.basis_blocks )
        assert( len(x) == self.number_of_params )
        extra_basis = ""
        for I,bblk in enumerate(self.basis_blocks[atm]):
            extra_basis += bblk.basis_str(x[self.__p0[atm][I]:self.__pN[atm][I]])
        return extra_basis 

    # Generates the string that contains the basis set for a species.
    #  atm: id of species
    #  x: parameter vector
    def basis_str(self, atm, x):
        """Returns a string with the basis set of a given species. 

        Parameters
        ----------
            atm: string
                String with the label of a species.
            x: 1-D floating point array.
                Array with variational parameters in the basis.
                Dimension must be equal to self.number_of_params.
        Returns
        -------
            extra_basis: String            
                A string with the basis set of species atm.
        """
        assert( atm in self.atom_id )
        assert( atm in self.default_basis )
        assert( len(x) == self.number_of_params )
        if os.path.exists(self.default_basis[atm]):
            return open(self.default_basis[atm],"r").read() + self.extended_basis(atm,x)
        else:
            return self.extended_basis(atm,x)

    # Prints the basis set of a species to a file.
    #   fname: file name
    #   atm: id of a species
    #   x: parameter vector
    def print_to_file(self,fname,atm,x):
        """Prints the basis set of a species to a file.

        Parameters
        ----------
            fname: string
                Name of file.
            atm: string
                String with the label of a species.
            x: 1-D floating point array.
                Array with variational parameters in the basis.
                Dimension must be equal to self.number_of_params.
        """
        with open(fname, "w") as f:
            f.write(self.basis_str(atm,x))

    # Prints all basis sets for species in fnames to file
    #   fnames: dictionary that maps species ids to file names.
    #   x: parameter vector
    def print_all_to_file(self,fnames,x):
        """Prints all the requested basis sets to separate file.

        Parameters
        ----------
            fnames: Python Dictionary 
                Dictionary with file names. Only those species contained in the object
                and defined in the dictionary will be printed.
            x: 1-D floating point array.
                Array with variational parameters in the basis.
                Dimension must be equal to self.number_of_params.
        """
        for I,atm in enumerate(self.atom_id):
            if atm in fnames:
                self.print_to_file(fnames[atm],atm,x)

    # Prints the basis set of a species to a file.
    #   fname: file name
    #   atm: id of a species
    #   x: parameter vector
    def print_to_stdout(self,atm,x):
        """Prints the basis set of a species to a stdout.

        Parameters
        ----------
            atm: string
                String with the label of a species.
            x: 1-D floating point array.
                Array with variational parameters in the basis.
                Dimension must be equal to self.number_of_params.
        """
        print("\n")
        print(self.basis_str(atm,x))
        print("\n")

    # Prints all basis sets for species in fnames to file
    #   fnames: dictionary that maps species ids to file names.
    #   x: parameter vector
    def print_all_to_stdout(self,x):        
        """Prints all the basis sets contained by the object to stdout. 

        Parameters
        ----------
            x: 1-D floating point array.
                Array with variational parameters in the basis.
                Dimension must be equal to self.number_of_params.
        """
        print("\n")
        for I,atm in enumerate(self.atom_id):
            print(self.basis_str(atm,x))
            print("\n")
        print("\n")


def default_basis_set(Lmax,atoms):
    """ Generates a OptimizableBasisSet object with default parameters.
        A single uncontacted gaussian per angular momentum will be used.
        Parameters
        ----------
            Lmax: integer. Requirement: 2 <= L <= 6
                Maximum angular momentum. For each l=2,3,...,Lmax,
                shells will be added with angular momentum: 0,1,...,l.
            atoms: array of strrings
                For each symbol in the array, a new basis set is created.
        Returns
        -------
            Object of type OptimizableBasisSet.
    """
    assert( Lmax >= 2 )
    assert( Lmax <= 6 )
    bset = OptimizableBasisSet(atoms)
    x = [] 
    def get_shell(sym,L,x):
        labels = ['S','P','D','F','G','H','I']
        bblk_list = []
        for i in range(0,L+1):
            # single uncontracted gaussians 
            bblk_list.append(BasisBlock(sym,labels[i]))
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
        return bblk_list, x 
    for a in atoms:
        for l_ in range(2,Lmax+1):
            bblk_list, x = get_shell(a,l_,x)
            for bblk in bblk_list: 
                bset.add(bblk)
    return numpy.array(x),bset

