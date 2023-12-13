#include <multi/adaptors/blas/gemv.hpp>
#include <multi/array.hpp>

#include <iostream>

int main() {
    float matA[3][3] = {
        {1.1, 2.2, 3.3},
        {4.4, 5.5, 6.6},
        {7.7, 8.8, 9.9},
    };
    float vecB[3] = {1.0, 2.0, 3.0};
    float vecC[3] = {0.0, 0.0, 0.0};

    namespace multi = boost::multi;

    {  // make references to c-arrays
        multi::array_ref<float, 2> A{matA};
        multi::array_ref<float, 1> B{vecB};
        multi::array_ref<float, 1> C{vecC};

        multi::blas::gemv(1.0, A, B, 0.0, C);  // C is output
    }
    {  // make references to c-arrays
        auto const& A = multi::ref(matA);  // deduce element type and dimensionality
        auto const& B = multi::ref(vecB);
        auto&&      vecCref = multi::ref(vecC);

        multi::blas::gemv(1.0, A, B, 0.0, vecCref);  // vecC holds the result
    }
    {  // one liner
        multi::blas::gemv(1.0, multi::ref(matA), multi::ref(vecB), 0.0, multi::ref(vecC));  // vecC holds the result
    }
    {  // one liner with output
        multi::ref(vecC) = multi::blas::gemv(1.0, multi::ref(matA), multi::ref(vecB));
    }
    {  // using the library, not references to c-arrays
        multi::array<float, 2> A = {
            {1.1, 2.2, 3.3},
            {4.4, 5.5, 6.6},
            {7.7, 8.8, 9.9},
        };
        multi::array<float, 1> B = {1.0, 2.0, 3.0};
 
        multi::array<float, 1> C = multi::blas::gemv(1.0, A, B);
    }
    {
        multi::array<float, 2> A = {
            {1.1, 2.2, 3.3},
            {4.4, 5.5, 6.6},
            {7.7, 8.8, 9.9},
        };
        multi::array<float, 1> B = {1.0, 2.0, 3.0};

        auto C =+ multi::blas::gemv(1.0, A, B);  // create (allocate) the result in C
    }
}
