
# Generate sparse matrix multiply test cases

import numpy as np
import scipy.sparse as sps


def flatten(a):
  '''Flatten a nested list'''
  return [item for d in a for item in d]

def string_elements(a):
  '''Convert elements of a list to strings'''
  return [str(item) for item in a]

def matrix_vector_mult():
  A = np.eye(3)
  X = np.ones((3,1))
  Y = np.ones((3,1))
  alpha = 1.0
  beta = 2.0

  output = generate_case(A, X, Y, alpha, beta, type='double')
  output += generate_case(A, X, Y, alpha, beta, type='complex<double>')

  A = np.eye(3)
  A[1,2] = 2.0
  X = np.arange(0,6).reshape(3,2)
  Y = np.arange(0,6).reshape(3,2)
  alpha = 1.0
  beta = 2.0

  output += '\n'
  output += generate_case(A, X, Y, alpha, beta, type='double')
  output += generate_case(A, X, Y, alpha, beta, type='complex<double>')

  A = np.eye(3)
  A[1,1] = 0.0
  A[1,2] = 1.0
  X = np.arange(0,9).reshape(3,3)
  Y = np.arange(0,9).reshape(3,3)
  alpha = 1.0
  beta = 2.0

  output += '\n'
  output += generate_case(A, X, Y, alpha, beta, type='double')
  output += generate_case(A, X, Y, alpha, beta, type='complex<double>')

  A = np.eye(2)
  A[0,1] = 1.0
  X = np.arange(0,6).reshape(2,3)
  Y = np.arange(0,6).reshape(2,3)
  alpha = 1.0
  beta = 2.0

  output += '\n'
  output += generate_case(A, X, Y, alpha, beta, type='double')
  output += generate_case(A, X, Y, alpha, beta, type='complex<double>')

  with open('sparse_mult_cases.cpp', 'w') as f:
    f.write(output)


case_id = 1
def generate_case(A, X, Y, alpha, beta, type='double'):
  global case_id

  b = alpha*(A.dot(X)) + beta*Y

  init_list = []
  for k,v in sps.dok_matrix(A).iteritems():
    init_list.append("s2Dd({i},{j},{v})".format(i=k[0],j=k[1],v=v))

  x_init_list = flatten(X.tolist())
  y_init_list = flatten(Y.tolist())
  ans_init_list = flatten(b.tolist())

  out = '''
TEST_CASE("sparse_matrix_{test_case}", "[sparse_matrix]")
{{
  typedef s2D<{type}> s2Dd;

  const int M = {M};
  const int N = {N};
  const int K = {K};

  std::vector<s2Dd> I = {{{init_list}}};
  SparseMatrix<{type}> A(M,K);
  A.initFroms2D(I,false);

  {type} B[K*N] = {{{x_init_list}}};
  {type} C[K*N] = {{{y_init_list}}};
  {type} alpha = {alpha};
  {type} beta = {beta};

  {type} finalC[K*N] = {{{ans_init_list}}};

  int ldb = N;
  int ldc = N;
  SparseMatrixOperators::product_SpMatM(M, N, K, alpha, A, B, ldb, beta, C, ldc);

  for (int i = 0; i < M; i++)
  {{
    for (int j = 0; j < N; j++)
    {{
      //printf("%d %d %g %g\\n",i,j,C[i*N + j],finalC[i*N + j]);
      REQUIRE(realPart(C[i*N + j]) == Approx(realPart(finalC[i*N + j])));
    }}
  }}
}}
  '''.format(M=A.shape[0],
             N=Y.shape[1],
             K=A.shape[1],
             init_list = ',\n    '.join(init_list),
             test_case = str(case_id),
             x_init_list = ','.join(string_elements(x_init_list)),
             y_init_list = ','.join(string_elements(y_init_list)),
             alpha = 1.0,
             beta = 2.0,
             ans_init_list = ','.join(string_elements(ans_init_list)),
             type = type
           )

  case_id += 1

  return out


if False:
  # Some experiments with CSR matrix format in scipy.sparse
  # Could potentially generate tests for generation of CSR matrices
  z = np.zeros((1,3))
  print z
  z[0,0] = 1.0
  z[0,1] = 2.0
  z[0,2] = 3.0

  M = sps.csr_matrix(z)

  print 'nnz = ', M.nnz
  print 'colums = ',M.indices
  print 'indptr = ',M.indptr
  print 'values = ',M.data

matrix_vector_mult()

