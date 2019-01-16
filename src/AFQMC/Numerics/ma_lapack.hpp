//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Alfredo Correa, correaa@llnl.gov 
//    Lawrence Livermore National Laboratory 
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Alfredo Correa, correaa@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef MA_LAPACK_HPP
#define MA_LAPACK_HPP

#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Numerics/detail/lapack.hpp"
#include "multi/array_ref.hpp"
#include<cassert>

namespace ma{

template<class MultiArray2D>
int getrf_optimal_workspace_size(MultiArray2D && A){
        assert(A.strides()[0] > 0);
        assert(A.strides()[1] == 1);

        int res;
        getrf_bufferSize(A.shape()[1], A.shape()[0],A.origin(),A.strides()[0],res);
        return res;
}

template<class MultiArray2D, class Array1D, class Buffer>
MultiArray2D getrf(MultiArray2D&& m, Array1D& pivot, Buffer&& WORK){
        assert(m.strides()[0] >= std::max(std::size_t(1), std::size_t(m.shape()[1])));
        assert(m.strides()[1] == 1);
        assert(pivot.size() >= std::min(m.shape()[1], m.shape()[0]+1));

        int status = -1;
        getrf(
                m.shape()[1], m.shape()[0], m.origin(), m.strides()[0],
                pivot.data(),
                status,
                WORK.data()
        );
        assert(status==0);
        return std::forward<MultiArray2D>(m);
}

template<class MultiArray2D>
int getri_optimal_workspace_size(MultiArray2D && A){
        assert(A.strides()[1] == 1);
        assert(A.shape()[0] == A.shape()[1]);
        int lwork = -1;
        getri_bufferSize(A.shape()[0], A.origin(), A.strides()[0],lwork);
        return lwork;
}

template<class MultiArray2D, class MultiArray1D, class Buffer>
MultiArray2D getri(MultiArray2D&& A, MultiArray1D const& IPIV, Buffer&& WORK){
//	assert(A.strides()[0] > std::max(std::size_t(1), A.shape()[1]));
	assert(A.strides()[1] == 1);
	assert(IPIV.size() >= size_t(A.shape()[0]));
	assert(WORK.size() >= std::max(std::size_t(1), size_t(A.shape()[0])));
	
	int status = -1;
	getri(
		A.shape()[0], A.origin(), A.strides()[0], 
		IPIV.data(), 
		WORK.data(), WORK.size(), 
		status
	);
	assert(status == 0);
	return std::forward<MultiArray2D>(A);
}

template<class MultiArray2D>
int geqrf_optimal_workspace_size(MultiArray2D const& A){
	assert(A.strides()[0] > 0);
	assert(A.strides()[1] == 1);

        int res;
        geqrf_bufferSize(A.shape()[1], A.shape()[0],A.origin(),A.strides()[0],res);
        return res;
}

template<class MultiArray2D, class Array1D, class Buffer>
MultiArray2D geqrf(MultiArray2D&& A, Array1D&& TAU, Buffer&& WORK){
        // why was this here???
	//assert(A.strides()[0] > std::max(std::size_t(1), A.shape()[0]));
	assert(A.strides()[1] == 1);
	assert(TAU.strides()[0] == 1);
	assert(TAU.size() >= std::max(std::size_t(1), size_t(std::min(A.shape()[0], A.shape()[1]))));
	assert(WORK.size() >= std::max(std::size_t(1), size_t(A.shape()[0])));
	
	int status = -1;
	geqrf(
		A.shape()[1], A.shape()[0], A.origin(), A.strides()[0], 
		TAU.data(), 
		WORK.data(), WORK.size(),
		status
	);
	assert(status==0);
	return std::forward<MultiArray2D>(A);
}

template<class MultiArray2D>
int gelqf_optimal_workspace_size(MultiArray2D const& A){
	assert(A.strides()[0] > 0);
	assert(A.strides()[1] == 1);

        int res;
        gelqf_bufferSize(A.shape()[1], A.shape()[0],A.origin(),A.strides()[0],res);
	return res;
}

template<class MultiArray2D, class Array1D, class Buffer>
MultiArray2D gelqf(MultiArray2D&& A, Array1D&& TAU, Buffer&& WORK){
	assert(A.strides()[1] > 0);
	assert(A.strides()[1] == 1);
	assert(TAU.strides()[0] == 1);
	assert(TAU.size() >= std::max(std::size_t(1), size_t(std::min(A.shape()[0], A.shape()[1]))));
	assert(WORK.size() >= std::max(std::size_t(1), size_t(A.shape()[1])));

	int status = -1;
	gelqf(
		A.shape()[1], A.shape()[0], A.origin(), A.strides()[0], TAU.data(),
		WORK.data(), WORK.size(), 
		status
	);
	assert(status==0);
	return std::forward<MultiArray2D>(A);
}


template<class MultiArray2D>
int gqr_optimal_workspace_size(MultiArray2D const& A){
	assert(A.strides()[0] > 0);
	assert(A.strides()[1] == 1);

        int res;
        gqr_bufferSize(A.shape()[1], A.shape()[0],
                       std::max(std::size_t(1), size_t(std::min(A.shape()[0], A.shape()[1]))),
                       A.origin(),A.strides()[0],res);
        return res;
}

template<class MultiArray2D, class Array1D, class Buffer>
MultiArray2D gqr(MultiArray2D&& A, Array1D&& TAU, Buffer&& WORK){
	assert(A.strides()[1] == 1);
	assert(TAU.strides()[0] == 1);
	assert(TAU.size() >= std::max(std::size_t(1), size_t(std::min(A.shape()[0], A.shape()[1]))));
	assert(WORK.size() >= std::max(std::size_t(1), size_t(A.shape()[0])));

	int status = -1;
	gqr(
		A.shape()[1], A.shape()[0], std::max(std::size_t(1), size_t(std::min(A.shape()[0], A.shape()[1]))), 
		A.origin(), A.strides()[0], TAU.data(), 
		WORK.data(), WORK.size(), 
		status
	);
	assert(status==0);
	return std::forward<MultiArray2D>(A);
}

template<class MultiArray2D>
int glq_optimal_workspace_size(MultiArray2D const& A){
	assert(A.strides()[0] > 0);
	assert(A.strides()[1] == 1);

        int res;
        glq_bufferSize(A.shape()[1], A.shape()[0],
                       std::max(std::size_t(1), size_t(std::min(A.shape()[0], A.shape()[1]))),
                       A.origin(),A.strides()[0],res);
        return res;
}

template<class MultiArray2D, class Array1D, class Buffer>
MultiArray2D glq(MultiArray2D&& A, Array1D&& TAU, Buffer&& WORK){
	assert(A.strides()[1] == 1);
	assert(TAU.strides()[0] == 1);
	assert(TAU.size() >= std::max(std::size_t(1), size_t(std::min(A.shape()[0], A.shape()[1]))));
	assert(WORK.size() >= std::max(std::size_t(1), size_t(A.shape()[1])));

	int status = -1;
	glq(
		A.shape()[1], A.shape()[0], std::max(std::size_t(1), size_t(std::min(A.shape()[0], A.shape()[1]))), 
		A.origin(), A.strides()[0], TAU.data(), 
		WORK.data(), WORK.size(), 
		status
	);
	assert(status==0);
	return std::forward<MultiArray2D>(A);
}

template<class MultiArray2D,
         typename = typename std::enable_if_t<MultiArray2D::dimensionality == 2>
        >
MultiArray2D potrf(MultiArray2D&& A) {
        assert(A.shape()[0]==A.shape()[1]);
        int INFO;
        potrf('U',A.shape()[0],A.origin(),A.strides()[0],INFO);
        if(INFO != 0) throw std::runtime_error(" error in ma::potrf: Error code != 0");
}

template<class MultiArray1D,
         class MultiArray2D,
         typename = typename std::enable_if_t<MultiArray1D::dimensionality == 1>,
         typename = typename std::enable_if_t<MultiArray2D::dimensionality == 2>
        >
std::pair<MultiArray1D,MultiArray2D> symEig(MultiArray2D const& A) {
        using eigSys = std::pair<MultiArray1D,MultiArray2D>;
        using Type = typename MultiArray2D::element; 
        using RealType = typename qmcplusplus::afqmc::remove_complex<Type>::value_type; 
        using extensions = typename boost::multi::layout_t<1u>::extensions_type;
        assert(A.shape()[0]==A.shape()[1]);
        assert(A.strides()[1]==1);
        assert(A.shape()[0]>0);
        int N = A.shape()[0];
        int LDA = A.strides()[0];
        
            MultiArray1D eigVal(extensions{N});
            MultiArray2D eigVec({N,N});
            MultiArray2D A_({N,N});
            for(int i=0; i<N; i++)
              for(int j=0; j<N; j++) 
                A_[i][j] = conj(A[i][j]);                
            char JOBZ('V');
            char RANGE('A');
            char UPLO('U');
            RealType VL=0;
            RealType VU=0;
            int IL=0;
            int IU=0;
            RealType ABSTOL=0;//DLAMCH( 'Safe minimum' );
            int M; // output: total number of eigenvalues found
            std::vector<int> ISUPPZ(2*N);
            std::vector<Type> WORK(1); // set with workspace query
            int LWORK=-1;
            std::vector<RealType> RWORK(1); // set with workspace query
            int LRWORK=-1;
            std::vector<int> IWORK(1);
            int LIWORK=-1;
            int INFO;

            hevr (JOBZ, RANGE, UPLO, N, A_.origin(), LDA, VL, VU, IL, IU, ABSTOL, 
                          M,eigVal.origin(), eigVec.origin(), N, ISUPPZ.data(), WORK.data(), LWORK, 
                          RWORK.data(), LRWORK, IWORK.data(), LIWORK, INFO);

            LWORK = int(real(WORK[0]));
            WORK.resize(LWORK);
            LRWORK = int(RWORK[0]);
            RWORK.resize(LRWORK);
            LIWORK = int(IWORK[0]);
            IWORK.resize(LIWORK);

            hevr (JOBZ, RANGE, UPLO, N, A_.origin(), LDA, VL, VU, IL, IU, ABSTOL, 
                          M,eigVal.origin(), eigVec.origin(), N, ISUPPZ.data(), WORK.data(), LWORK, 
                          RWORK.data(), LRWORK, IWORK.data(), LIWORK, INFO);
            if(INFO != 0) throw std::runtime_error(" error in ma::eig: Error code != 0");
            if(M != N) throw std::runtime_error(" error in ma::eig: Not enough eigenvalues"); 
            for(int i=0; i<N; i++)
              for(int j=i+1; j<N; j++) 
                std::swap(eigVec[i][j],eigVec[j][i]);                

            return std::pair<MultiArray1D,MultiArray2D>{eigVal,eigVec};

}

}

#endif

