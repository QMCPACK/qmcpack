ParticleAttribSoA
=================
Experimental pacakage implements classes to facilitate SIMD algorithms.
* SoaContainer.h
* VectorSoaContainer.h : replacing container<TinyVector<T,3>>
* TensorSoaContainer.h : replacing container<Tensor<T,3>>

VectorSoaContainer
---------------
\code
VectorSoaContainer<double,3> R;
auto r=R[i];                       //get the value of the i-th position
R(i)=0;                            //assign to the i-th position
R(i)=TinyVector<double,3>(-1,2,3); //assign  to the i-th position
\endcode

Access operators to each compoenent are provided.

\code
R.data(0); //return the starting address of X component
R.data(1); //return the starting address of Y component
R.data(2); //return the starting address of Z component
\endode

TensorSoaContainer
------------------
\code
TensorSoaContainer<double,3> H;
auto h=H[i];                             //get the value of the i-th hessian as Tensor<T,3>
H(i)=0;                                  //assign to the i-th hessian
H(i)=Tensor<double,3>(0,1,2,1,3,4,2,4,5);//assign  to the i-th hessian
\endcode

Access operators to each compoenent are provided.

\code
H.data(int i, int j); //return the starting address of (i,j) component
\endode

