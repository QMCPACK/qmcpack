
import numpy as np

# Script to generate data for inverse and inverse update tests

# Output for copy-and-paste to C++
def output_for_cpp(A,var_name='a'):
  N = A.shape[0]
  for i in range(N):
    for j in range(N):
      print '%s(%d,%d) = %12.10g;'%(var_name,i,j,A[i,j])


A = np.array([
[2.3, 4.5, 2.6],
[0.5, 8.5, 3.3],
[1.8, 4.4, 4.9]])

#A = np.array([
#[2.3, 4.5, 2.6, 1.2],
#[0.5, 8.5, 3.3, 0.3],
#[1.8, 4.4, 4.9, 2.8],
#[0.8, 4.1, 3.2, 1.1]
#])

print A

print 'det A',np.abs(np.linalg.det(A))
print 'log det A',np.log(np.abs(np.linalg.det(A)))
Ainv = np.linalg.inv(A)
print Ainv

output_for_cpp(Ainv)

print ''
row = np.array([1.9, 2.0, 3.1])
#row  = np.array([[1.9, 2.0, 3.1],
#                 [0.1, 4.2, 1.4]])

#row = np.array([3.2, 0.5, 5.9, 3.7])

B = A.copy()
# update row
#B[0,:] = row
# update column
B[:,0] = row
#B[:,0] = row[0]
#B[:,1] = row[1]
print 'Updated A with column to get matrix B:'
print B
print 'det B',np.abs(np.linalg.det(B))
print 'log det B',np.log(np.abs(np.linalg.det(B)))
print 'det ratio',np.linalg.det(B)/np.linalg.det(A);
Binv = np.linalg.inv(B)

print ''
output_for_cpp(Binv,var_name='b')
