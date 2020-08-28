import numpy as np
from scipy.special import factorial2

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

def cartesian_ordering_test(order):
    #using default cartesian ordering...namely GAMESS
    
    if order == "Gamess":
        ordering = {0: [(0,0,0)],
                    1: [(1,0,0),(0,1,0),(0,0,1)],
                    2: [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)],
                    3: [(3,0,0),(0,3,0),(0,0,3),(2,1,0),(2,0,1),(1,2,0),(0,2,1),(1,0,2),(0,1,2),(1,1,1)],
                    4: [(4,0,0),(0,4,0),(0,0,4),(3,1,0),(3,0,1),(1,3,0),(0,3,1),(1,0,3),(0,1,3),(2,2,0),(2,0,2),(0,2,2),(2,1,1),(1,2,1),(1,1,2)]
                   }
        print("using Gamess")
        print(ordering)
    elif order == "Dirac":
        ordering = {0: [(0,0,0)],
                    1: [(1,0,0),(0,1,0),(0,0,1)],
                    2: [(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)],
                    3: [(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),(0,3,0),(0,2,1),(0,1,2),(0,0,3)],
                    4: [(4,0,0),(3,1,0),(3,0,1),(2,2,0),(2,1,1),(2,0,2),(1,3,0),(1,2,1),(1,1,2),(1,0,3),(0,4,0),(0,3,1),(0,2,2),(0,1,3),(0,0,4)]
                   }
        print("using Dirac")
        print(ordering)

    coeffs = {0: [0.5],
              1: [0.3,-0.2,0.1],
              2: [0.6,-0.5,0.4,-0.3,0.2,-0.1],
              3: [1.0,-0.9,0.8,-0.7,0.6,-0.5,0.4,-0.3,0.2,-0.1],
              4: [1.5,-1.4,1.3,-1.2,1.1,-1.0,0.9,-0.8,0.7,-0.6,0.5,-0.4,0.3,-0.2,0.1]
             }

    pos = np.array([0.1,-0.3,0.2])

    orbitalValue = 0
    for l in range(len(ordering)):
        for n in range(len(ordering[l])):
            i,j,k = ordering[l][n]
            gbf = cartGauss(1.0,l,i,j,k)
            print(i,j,k,gbf.val(pos))
            orbitalValue += coeffs[l][n] * gbf.val(pos)
    print("Orbital Value: {}".format(orbitalValue))

if __name__ == "__main__":
    cartesian_ordering_test("Gamess")
    cartesian_ordering_test("Dirac")

