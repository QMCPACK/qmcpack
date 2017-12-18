#!/usr/bin/env python
import re

def check_next_particle_grad_lap(mm,rel_tot,val_thr=1e-2,header='For particle'):
  """ check if particle gradient and laplacian errors are within tolerance
   assume wftest output contains the blocks such as:
For particle #0 at         2.77662296       2.025233638       3.660986223
Gradient      =         6.68693298       -3.95727398       7.476365008
  Finite diff =        6.686934158      -3.957273925       7.476365603
  Error       =   -1.178022464e-06  -5.538599357e-08  -5.949851918e-07
  Relative Error = 1.761678288e-07 1.39959967e-08 7.95821487e-08

Laplacian     = -93.12287222
  Finite diff =  -93.1228783
  Error       = 6.076049388e-06  Relative Error = 6.524765449e-08

  Args:
   mm (mmap): memory map of wftest output file e.g. wftest.000
   rel_tot (float,optional): tolerance for relative error.
   val_thr (float,optional): threshold for significance of derivative, default is 1e-2. Derivatives below this threshold are treated as zero. This threhold exists because small error in a tiny value can result in a large relative error.
   header (str,optional): default is 'For particle'
  Returns:
   bool: success
  """

  idx = mm.find(header)
  if idx == -1:
    raise RuntimeError('failed to find %s'%header)
  # end if

  mm.seek(idx)
  mm.readline() # skip header

  success = True # determine success with the following checks

  # 1. gradient error
  idx = mm.find('Gradient')
  mm.seek(idx)
  grad_line = mm.readline()
  grad_xyz  = grad_line.split('=')[-1]
  try: # real values
    grad_val = map(float,grad_xyz.split())
  except: # complex values
    xyzl = re.split(r'[(,)]',grad_xyz.strip('\n'))
    grad_real = map(float,xyzl[1::3])
    grad_imag = map(float,xyzl[2::3])
    grad_val  = [grad_real[i] + 1j*grad_imag[i] for i in range(3)]
  # end try
  idx = mm.find('Relative Error')
  mm.seek(idx)
  grad_line = mm.readline()
  grad_re   = map(float,grad_line.split()[-3:]) # relative error
  if (sum(grad_re)>rel_tot): # check if error is significant
    if (abs(sum(grad_val))>val_thr): # ignore small absolute errors
      success = False
    # end if
  # end if

  # 2. laplacian error
  idx = mm.find('Laplacian')
  mm.seek(idx)
  lap_valt= mm.readline().split('=')[-1]
  try: # real
    lap_val = float(lap_valt)
  except: # complex
    lapl = re.split(r'[(,)]',lap_valt.strip('\n'))
    lap_val = float(lapl[1]) + 1j*float(lapl[2])
  # end try
  idx = mm.find('Relative Error')
  mm.seek(idx)
  lap_line = mm.readline()
  tokens   = lap_line.split()
  lap_re   = float(tokens[-1]) # relative error
  if (lap_re>rel_tot):
    if (abs(lap_val)>val_thr):
      success = False
    # end if
  # end if

  return success
# end def

def all_lines_with_tag(mm,tag,nline_max):
    """ return a list of memory indices pointing to the start of tag """
    mm.seek(0) # rewind file
    all_idx = []
    for iline in range(nline_max):
        idx = mm.find(tag)
        if idx == -1:
            break
        # end if
        mm.seek(idx)
        all_idx.append(idx)
        mm.readline()
    # end for iline

    # guard
    if iline >= nline_max-1:
        raise NotImplementedError('may need to increase nline_max')
    # end if
    return all_idx
# end def all_lines_with_tag

if __name__ == '__main__':

  fname = 'wftest.000'
  nline_max = int(1e6) # assume fname has at most 1 million lines
  ratio_tol = 1e-6 # tolerance for ratio test
  rel_tol   = 1e-3 # tolerance for relative error (0.1%)

  from mmap import mmap
  with open(fname,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with

  # 1. grade finite-difference test
  plocs = all_lines_with_tag(mm,'For particle',nline_max=nline_max)
  fsuccess = True # finite-difference test success
  for ploc in plocs:
    mm.seek(ploc)
    success = check_next_particle_grad_lap(mm,rel_tol)
    fsuccess= fsuccess&success
  # end for plocs

  # 2. grade ratio test
  """ use function designed for: 
  Deriv  Numeric Analytic Diff
      to parse:
  Particle       Ratio of Ratios     Computed Ratio   Internal Ratio
  """
  from check_deriv import parse_deriv_block
  data = parse_deriv_block(mm,'Particle       Ratio of Ratios')
  # parse_deriv_block was not written for grad_lap, need to rename
  gl_name_map = {'Particle':'iparam','ratio':'numeric','computed':'analytic','internal':'diff'}
  rr = data[gl_name_map['ratio']]
  rsuccess = abs(sum(rr)/len(rr) - 1.) < ratio_tol

  if fsuccess:
    print('Finite difference test: PASS')
  else:
    print('Finite difference test: FAIL')
  if rsuccess:
    print('Ratio test: PASS')
  else:
    print('Ratio test: FAIL')

  if fsuccess and rsuccess:
    exit(0)
  else:
    exit(1)

# end __main__
