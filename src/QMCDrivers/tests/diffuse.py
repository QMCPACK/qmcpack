
# Compute diffusion using fake RNG
# See qmc_algorithms/Diffusion/DMC_Propagator.ipynb

#  Generate values for test_vmc.cpp and test_vmc_omp.cpp
#
import numpy as np
import math
import sys

# Box-Muller, copies Utilties/RandomGenerator.h
def gaussian_rng_list(n):
    input_rng = [0.5]*(n+1) # Fake RNG currently outputs 0.5 every time
    slightly_less_than_one = 1.0 - sys.float_info.epsilon
    vals = []
    for i in range(0,n,2):
        temp1 = math.sqrt(-2.0 * math.log(1.0- slightly_less_than_one*input_rng[i]))
        temp2 = 2*math.pi*input_rng[i+1]
        vals.append(temp1*math.cos(temp2))
        vals.append(temp2*math.sin(temp2))
    if n%2 == 1:
        temp1 = math.sqrt(-2.0 * math.log(1.0- slightly_less_than_one*input_rng[n-1]))
        temp2 = 2*math.pi*input_rng[n]
        vals.append(temp1*math.cos(temp2))
    return vals


#  These functions taken from DMC_Propagator Jupyter notebook
def scaled_drift_func(tau, Fmag2, F_i):
  return (F_i*((tau) if (Fmag2 < 2.22044604925031e-16) else ((((math.sqrt(2*Fmag2*tau + 1) - 1)/Fmag2) if (True) else None))))

def drift_diffuse_func(r,F_i,chi,dt):
  return chi + r + F_i*dt

tau = 0.1

chi_vals = np.array(gaussian_rng_list(6)).reshape((2,3))

R = np.array([[1.0, 0.0, 0.0],
              [0.0, 0.0, 1.0]])

scaled_chi_vals = chi_vals * math.sqrt(tau)

print 'One step'
for r_val, chi_val in zip(R, scaled_chi_vals):
  rp_val = drift_diffuse_func(r_val, np.zeros(3), chi_val, tau)
  print ['%.15g'%v for v in rp_val]

# In test_vmc_omp, there are two steps taken - one for the initial 'randomize' step,
#  and one for for the actual VMC step.
print 'Two steps'
for r_val, chi_val in zip(R, scaled_chi_vals):
  rp_val = drift_diffuse_func(r_val, np.zeros(3), chi_val, tau)
  rp2_val = drift_diffuse_func(rp_val, np.zeros(3), chi_val, tau)
  print ['%.15g'%v for v in rp2_val]

