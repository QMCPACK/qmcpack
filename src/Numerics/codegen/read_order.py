
# Generate order of i,j,k indices from Gamess output

# Also see the SETLAB routine in inputb.src for the ordering

# Example notation:
# XYZ
# 2X
def count_vals(s):
    x,y,z = 0,0,0
    mult = 1
    for c,c1 in zip(s,s[1:]+' '):
        #print 'c,c1',c,c1
        if c.isdigit():
            continue
        if c1.isdigit():
            mult = int(c1)
        if c == 'X':
            x += mult
        elif c == 'Y':
            y += mult
        elif c == 'Z':
             z += mult
        mult = 1
        if not c.isdigit() and c not in ['S','X','Y','Z']:
             print 'Error, unknown character',c
    return x,y,z

# order by max number of repeats, then x,y,z
def to_string(x,y,z):
    order = [('x',x),('y',y),('z',z)]
    order.sort(key=lambda x:x[1],reverse=True)
    s = ''
    for c,n in order:
        s += c*n
    return s
        


# order by x,y,z
def to_string_simple(x,y,z):
    s = 'x'*x
    s += 'y'*y
    s += 'z'*z
    return s

# used in Python script to autogen
def create_get_ijk(order_list):
    out_str = """
def get_ijk():
  ijk = []
%s
  return ijk"""

    shell_name = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
    current_sum = -1
    body_str = ""
    for x,y,z,s in order_list:
        new_sum = x+y+z
        if new_sum != current_sum:
            body_str +=  '  # '+shell_name[new_sum] + '\n'
            current_sum = new_sum
        body_str += '  ijk.append( (%d,%d,%d,"%s") )\n'%(x,y,z,s)

    return out_str%body_str

# generate C++ used in CartesianTensor.h
def create_getABC(order_list):
    out_str = """
template<class T, class Point_t, class Tensor_t, class GGG_t>
void CartesianTensor<T,Point_t, Tensor_t, GGG_t>::getABC(int n, int& a, int& b, int& c)
{
// following Gamess notation
  switch(n)
  {
%s
  default:
    std::cerr <<"CartesianTensor::getABC() - Incorrect index." << std::endl;
    APP_ABORT("");
    break;
  }
}
"""
    shell_name = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
    idx = 0
    body_str = ''
    current_sum = -1
    for x,y,z,s in order_list:
        new_sum = x+y+z
        if new_sum != current_sum:
            body_str +=  '  // '+shell_name[new_sum] + '\n'
            current_sum = new_sum
        body_str += '  case %d: // %s\n'%(idx,s)
        body_str += '    a = %d; b = %d; c = %d;\n'%(x,y,z)
        body_str += '    break;\n'
        idx += 1
    return out_str%body_str

# generate C++ used in SoaCartesianTensor.h
def create_getABC_SoA(order_list):
    out_str = """
template<class T>
void SoaCartesianTensor<T>::getABC(int n, int& a, int& b, int& c)
{
// following Gamess notation
  switch(n)
  {
%s
  default:
    std::cerr <<"CartesianTensor::getABC() - Incorrect index." << std::endl;
    APP_ABORT("");
    break;
  }
}
"""
    shell_name = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
    idx = 0
    body_str = ''
    current_sum = -1
    for x,y,z,s in order_list:
        new_sum = x+y+z
        if new_sum != current_sum:
            body_str +=  '  // '+shell_name[new_sum] + '\n'
            current_sum = new_sum
        body_str += '  case %d: // %s\n'%(idx,s)
        body_str += '    a = %d; b = %d; c = %d;\n'%(x,y,z)
        body_str += '    break;\n'
        idx += 1
    return out_str%body_str

def create_gms_order_for_pyscf(order_list):
    pad = ' '*4
    out_str = pad + "d_gms_order = {\n%s%s}"
    body_str = ''
    current_sum = -1
    shell_list = list()
    for x,y,z,s in order_list:
      new_sum = x + y + z
      if new_sum != current_sum:
        if current_sum >= 0:
          body_str += 2*pad + '%d:[%s],\n'%(current_sum, ','.join(shell_list))
        shell_list = list()
        current_sum = new_sum
      expanded_shell_name = to_string(*count_vals(s))
      shell_list.append('"%s"'%expanded_shell_name.lower())

    body_str += 2*pad + '%d:[%s],\n'%(current_sum, ','.join(shell_list))

    return out_str%(body_str,pad)


def read_order(fname):
    order_list = []
    already_seen = dict()
    with open(fname, 'r') as f:
        for line in f:
            #line = line.rstrip()
            if not line:
                continue
            order = line[9:16].strip()
            if order.startswith('1'):
                order = order[1:]
            if order not in already_seen:
                x,y,z = count_vals(order.strip())
                order_list.append((x,y,z,order.strip()))
                already_seen[order] = 1
            #new_sum = x + y + z
            #if new_sum != current_sum:
            #    print '  # ',shell_name[new_sum]
            #    current_sum = new_sum
    
            #print order,x,y,z
            #print "  ijk.append('%s')"%to_string(x,y,z)
    return order_list

if __name__ == "__main__":
    # The file 'order.txt' is taken from Gamess output
    order_list = read_order('order.txt')

    # Create get_ijk for gen_cartesian_tensor.py
    #print create_get_ijk(order_list)

    # Create getABC for CartesianTensor.h.in
    #print create_getABC(order_list)

    # Create getABC for SoaCartesianTensor.h.in
    #print create_getABC_SoA(order_list)

    # Create order dictionary for PyscfToQmcpack.py in QMCTools
    print create_gms_order_for_pyscf(order_list)
