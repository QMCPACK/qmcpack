#ifndef AFQMC_READWFN_HPP
#define AFQMC_READWFN_HPP

#include<cstdlib>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<ctype.h>

#include "Utilities/SimpleParser.h"

#include "AFQMC/config.h"
#include "boost/multi_array.hpp"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Matrix/csr_matrix_construct.hpp"

namespace qmcplusplus
{

namespace 
{

void read_header(std::ifstream& in, std::string& type, int& wfn_type,  bool& fullMOMat, bool& Cstyle, int& ndet)
{

  std::vector<std::string> words;
  getwords(words,in);
  do {
    if(words.size() == 0)
      app_error()<<"Format error in ASCII integral file. End of file in header. \n";
    for(std::vector<std::string>::iterator it=words.begin(); it!=words.end(); it++) {
      if(*it == "&FCI") {
        // do nothing 
      } else if(*it == "Type" || *it == "TYPE" || *it == "type") {
        if( it+1 == words.end() )  {
          app_error()<<"Format error in ASCII integral file. NAEB \n";
          APP_ABORT("Format error in ASCII integral file. NAEB \n");
        }
        type = *(it+1);
        it++;
      } else if(*it == "NCI" || *it == "nci") {
        if( it+1 == words.end() )  {
          app_error()<<"Format error in ASCII integral file. NETOT \n";
          APP_ABORT("Format error in ASCII integral file. NETOT \n");
        }
        ndet = atoi((++it)->c_str());
      } else if(*it == "UHF" || *it == "GHF") {
        if( it+1 == words.end() ) {
          app_error()<<"Format error in ASCII integral file. UHF/GHF \n";
          APP_ABORT("Format error in ASCII integral file. UHF/GHF \n");
        }
        wfn_type = atoi((++it)->c_str());
        switch(wfn_type) {
          case 0:
          {
            app_log()<<"Reading a RHF-type trial wave-function. \n";
            break;
          }
          case 1:
          {
            app_log()<<"Reading a UHF-type trial wave-function. \n";
            break;
          }
          case 2:
          {
            app_log()<<"Reading a GHF-type trial wave-function. \n";
            app_error()<<" GHF type not implemented. \n";
            APP_ABORT(" GHF type not implemented in readWfn.. \n");
            break;
          }
          default:
          {
            app_error()<<"Unknown wave-function type in AFQMC/Utilities/readWfn.hpp: " <<wfn_type <<std::endl;
            APP_ABORT("Unknown wave-function type in AFQMC/Utilities/readWfn.hpp \n"); 
          }
        }
      } else if(*it == "FullMO" || *it == "FULLMO") {
        fullMOMat = true;
      } else if(*it == "CMajor") {
        Cstyle = false;
      } 
    }
    getwords(words,in);
    if(words.size() == 0)
      app_error()<<"Format error in ASCII integral file. End of file in header. \n";
  } while((words[0].find(std::string("/"))==std::string::npos && words[0].find(std::string("&END"))==std::string::npos));

} 

void skip_determinant(std::ifstream& in, bool Cstyle, bool fullMOMat, int NMO, int NAEA)
{
      int nread;
      ComplexType dummy;
      if(Cstyle) {
        nread = fullMOMat?NMO:NAEA;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in readWfn.hpp.  \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              APP_ABORT("Problems reading ASCII file in readWfn.hpp.  \n");  
            }
          }
      } else {
        nread = fullMOMat?NMO:NAEA;
        for(int j=0; j<nread; j++)
          for(int i=0; i<NMO; i++) {
            in>>dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in readWfn.hpp. \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              APP_ABORT("Problems reading wfn in readWfn.hpp.\n"); 
            }
          }
      }
}

template<class Mat>
void read_mat(std::ifstream& in, Mat&& OrbMat, bool Cstyle, bool fullMOMat, int NMO, int NAEA)
{
      int nread;
      ComplexType dummy;
      if(Cstyle) {
        nread = fullMOMat?NMO:NAEA;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEA) OrbMat[i][j] = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in readWfn.hpp.  \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              APP_ABORT("Problems reading ASCII file in readWfn.hpp.  \n");
            }
          }
      } else {
        nread = fullMOMat?NMO:NAEA;
        for(int j=0; j<nread; j++)
          for(int i=0; i<NMO; i++) {
            in>>dummy;
            if(j<NAEA) OrbMat[i][j] = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in readWfn.hpp. \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              APP_ABORT("Problems reading wfn in readWfn.hpp.\n");
            }
          }
      }
}

}

namespace afqmc
{

WALKER_TYPES getWalkerType(std::string filename)
{
  std::ifstream in;
  in.open(filename.c_str());
  if(in.fail()) {
     app_error()<<"Problems opening file:  " <<filename <<std::endl;
     APP_ABORT("Problems opening ASCII integral file. \n");
  }

  bool fullMOMat = false;
  bool Cstyle = true;
  int wfn_type=0;
  int ndet_in_file=-1;
  std::string type;

  read_header(in,type,wfn_type,fullMOMat,Cstyle,ndet_in_file);
  in.close();

  if(wfn_type == 0) return CLOSED;
  else if(wfn_type == 1) return COLLINEAR;
  else if(wfn_type == 2) return NONCOLLINEAR;
  else return UNDEFINED_WALKER_TYPE;
}

/*
 * Reads ndets from the ascii file. 
 * If pureSD == false, PsiT contains the Slater Matrices of all the terms in the expansion.
 *      For walker_type==1, PsiT contains 2*ndets terms including both Alpha/Beta components.
 * If pureSD == true, PsiT contains only the reference determinant and excitations contains
 *      the occupation strings of all the determinants in the expansion, including the reference.  
 */ 
void read_wavefunction(std::string filename, int& ndets, std::string& type, WALKER_TYPES walker_type,
        boost::mpi3::shared_communicator& comm, int NMO, int NAEA, int NAEB,
        std::vector<PsiT_Matrix>& PsiT, std::vector<ComplexType>& ci,
        std::vector<int>& excitations) 
{
  assert(walker_type!=UNDEFINED_WALKER_TYPE);
  std::ifstream in;
  in.open(filename.c_str());
  if(in.fail()) {
     app_error()<<"Problems opening file:  " <<filename <<std::endl;
     APP_ABORT("Problems opening file. \n");
  }

  bool fullMOMat = false;
  bool Cstyle = true;
  int wfn_type=0;
  int ndet_in_file=-1;
  int NEL = NAEA;
  if(walker_type!=CLOSED) NEL+=NAEB;

  /*
   * type:
   *   - occ: All determinants are specified with occupation numbers 
   *   - matrix: Slater matrix for all terms in the expansion
   *   - mixed: mixed representation. Reference determinant in matrix form, 
   *            determinant list (including reference) in occupation numbers. 
   */ 
  read_header(in,type,wfn_type,fullMOMat,Cstyle,ndet_in_file);

  std::transform(type.begin(),type.end(),type.begin(),(int (*)(int)) tolower);

  /*
   * Expected order of inputs and tags:
   * Coefficients: 
   * Occupations: 
   * Determinant: 
   */

  if(ndets <= 0) ndets = ndet_in_file; 
  
  if(ndet_in_file < ndets) 
    APP_ABORT("Error: Requesting too many determinants from wfn file.\n");
  if(type != "occ" && type != "matrix" && type != "mixed") 
    APP_ABORT("Error: Unknown wavefunction type in file.\n"); 
  if(wfn_type == 1 && walker_type == CLOSED) 
    APP_ABORT("Error in read_wavefunction: walker_type < wfn_type. \n");
  if(wfn_type == 2 && (walker_type == CLOSED || walker_type == COLLINEAR)) 
    APP_ABORT("Error in read_wavefunction: walker_type < wfn_type. \n");

  if( walker_type == COLLINEAR && (NAEA==NAEB && wfn_type == 0) ) 
    app_log()<<"  MESSAGE: Using walker_type=colinear with a closed-shell wfn with NAEA==NAEB. \n"
             <<"           Consider doing a closed shell calculation ( walker_type=closed in WalkerSet)\n"; 

  ci.reserve(ndets);
  excitations.clear();

  std::string tag;
  in>>tag;
  if(tag != "Coefficients:")
    APP_ABORT(" Error: Expecting Coefficients: tag in wavefunction file. \n");
  ComplexType dummy;
  for(int i=0; i<ndet_in_file; i++) {
    in>>dummy;
    if(i<ndets) ci.emplace_back(dummy);
  }
  
  if(type == "occ" || type == "mixed") {
    in>>tag;
    if(tag != "Occupations:")  
      APP_ABORT(" Error: Expecting Occupations: tag in wavefunction file. \n");
    excitations.reserve(NEL*ndets);
    if(walker_type == NONCOLLINEAR)
      APP_ABORT(" Error: walker_type==NONCOLLINEAR not yet implemented in read_wavefunction.\n");
    for(int i=0; i<ndet_in_file; i++) {
      if(i<ndets) {
        for(int k=0, q=0; k<NAEA; k++) {
          in>>q;
          if(q < 1 || q > NMO)
            APP_ABORT("Error: Bad occupation number in wavefunction file. \n");
          excitations.emplace_back(q-1);
        }  
        if(walker_type == COLLINEAR) {
          if(wfn_type==0) {
            for(int k=0, q=0; k<NAEB; k++) {
              q = *(excitations.rbegin()+NAEA-1)+NMO;
              excitations.emplace_back(q);
            }
          } else {
            for(int k=0, q=0; k<NAEB; k++) {
              in>>q;
              if(q <= NMO || q > 2*NMO)
                APP_ABORT("Error: Bad occupation number in wavefunction file. \n");
              excitations.emplace_back(q-1);
            } 
          }    
        }
      } else { 
        int nel_ = NAEA;
        if(wfn_type > 0) nel_+=NAEB; 
        for(int k=0, q=0; k<nel_; k++) in>>q; 
      }
    }
    // sort excitations
    auto it = excitations.begin();
    for(int i=0; i<ndets; i++) {
      std::sort(it,(it+=NAEA));
      if(walker_type == COLLINEAR) std::sort(it,(it+=NAEB));
    }
  }
  
  if(type == "occ")  {
    // empty PsiT in this case
  } else if(type == "mixed") {
    // what you read into PsiT is the orbital set from which  
    // slater matrices are built
    //PsiT.reserve( (walker_type!=1)?1:2 );
    APP_ABORT(" FINISH implementation of read_wavefunction \n\n\n");
  } else if(type == "matrix") {
    PsiT.reserve( (walker_type!=COLLINEAR)?ndets:2*ndets );

    if(wfn_type == 0) {

      boost::multi_array<ComplexType,2> OrbMat(extents[NMO][NAEA]);
      for(int i=0,q=0; i<ndets; i++) {
        if(comm.rank()==0) {
          in>>tag >>q;
          if(tag != "Determinant:" || q!=i+1)  
            APP_ABORT(" Error: Expecting Determinant: # tag in wavefunction file. \n");
          read_mat(in,OrbMat,Cstyle,fullMOMat,NMO,NAEA);
        }
        PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat,1e-8,'H',comm));
        if(walker_type==COLLINEAR)
          PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[indices[range_t()][range_t(0,NAEB)]],
                                        1e-8,'H',comm));
      }  

    } else if(wfn_type == 1) {

      boost::multi_array<ComplexType,2> OrbMat(extents[NMO][NAEA]);
      for(int i=0,q=0; i<ndets; i++) {
        if(comm.rank()==0) {
          in>>tag >>q;
          if(tag != "Determinant:" || q!=i+1)
            APP_ABORT(" Error: Expecting Determinant: # tag in wavefunction file. \n");
          read_mat(in,OrbMat,Cstyle,fullMOMat,NMO,NAEA);
        }
        PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat,1e-8,'H',comm));
        if(comm.rank()==0)  read_mat(in,OrbMat[indices[range_t()][range_t(0,NAEB)]],
                Cstyle,fullMOMat,NMO,NAEB);
        PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[indices[range_t()][range_t(0,NAEB)]],
                                        1e-8,'H',comm));
      }

    } else if(wfn_type == 2) {

      boost::multi_array<ComplexType,2> OrbMat(extents[2*NMO][NAEA]);
      for(int i=0,q=0; i<ndets; i++) {
        if(comm.rank()==0) {
          in>>tag >>q;
          if(tag != "Determinant:" || q!=i+1)
            APP_ABORT(" Error: Expecting Determinant: # tag in wavefunction file. \n");
          read_mat(in,OrbMat,Cstyle,fullMOMat,2*NMO,NAEA);
        }
        PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat,1e-8,'H',comm));
        if(comm.rank()==0)  read_mat(in,OrbMat[indices[range_t()][range_t(0,NAEB)]],
                                        Cstyle,fullMOMat,NMO,NAEB);
        PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat[indices[range_t()][range_t(0,NAEB)]],
                                        1e-8,'H',comm));
      }

    } //type 
 
  } else {
    APP_ABORT("Error: Unknown wavefunction type in file.\n");
  }

  in.close();
}


// old routine
bool readWfn( std::string fileName, ComplexMatrix& OrbMat, int NMO, int NAEA, int NAEB)
{

  std::ifstream in;
  in.open(fileName.c_str());
  if(in.fail()) {
     app_error()<<"Problems opening ASCII integral file:  " <<fileName <<std::endl;
     APP_ABORT("Problems opening ASCII file in readWfn.hpp.  \n");  
  }

  bool fullMOMat = false;
  bool Cstyle = true;
  int wfn_type;
  int ndet;
  std::string type;

  read_header(in,type,wfn_type,fullMOMat,Cstyle,ndet);

    int ncols = NAEA;
    //int nrows = NMO;
    if(wfn_type == 2)
      ncols = NAEA+NAEB;

    ComplexType dummy;
    int nread;

    if (wfn_type == 0 ) {

      OrbMat.resize(2*NMO,ncols);
      if(Cstyle) {
        nread = fullMOMat?NMO:NAEA;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEA) OrbMat(i,j) = dummy;
            if(j<NAEB) OrbMat(i+NMO,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in readWfn.hpp.  \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              APP_ABORT("Problems reading wfn in readWfn.hpp.\n"); 
            }
          }
      } else {
        nread = fullMOMat?NMO:NAEA;
        for(int j=0; j<nread; j++)
          for(int i=0; i<NMO; i++) {
            in>>dummy;
            if(j<NAEA) OrbMat(i,j) = dummy;
            if(j<NAEB) OrbMat(i+NMO,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in readWfn.hpp. \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              APP_ABORT("Problems reading wfn in readWfn.hpp.\n");
            }
          }
      }

    } else if(wfn_type == 1) {

      OrbMat.resize(2*NMO,ncols);
      if(Cstyle) {
        nread = fullMOMat?NMO:NAEA;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEA) OrbMat(i,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. (alpha) \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              APP_ABORT("Problems reading ASCII file in PureSingleDeterminant. (alpha) \n");
            }
          }
        nread = fullMOMat?NMO:NAEB;
        for(int i=0; i<NMO; i++)
          for(int j=0; j<nread; j++) {
            in>>dummy;
            if(j<NAEB) OrbMat(i+NMO,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. (beta) \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              return false;
            }
          }
      } else {
        nread = fullMOMat?NMO:NAEA;
        for(int j=0; j<nread; j++)
          for(int i=0; i<NMO; i++) {
            in>>dummy;
            if(j<NAEA) OrbMat(i,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. (alpha) \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              return false;
            }
          }
        nread = fullMOMat?NMO:NAEB;
        for(int j=0; j<nread; j++)
          for(int i=0; i<NMO; i++) {
            in>>dummy;
            if(j<NAEB) OrbMat(i+NMO,j) = dummy;
            if(in.fail()) {
              app_error()<<"Problems reading ASCII file in PureSingleDeterminant. (beta) \n";
              app_error()<<i <<" " <<j <<std::endl;
              in.close();
              return false;
            }
          }
      }
    } else if(wfn_type == 2) {

      if(Cstyle) {

       nread = fullMOMat?NMO:NAEA;
       int nread2 = fullMOMat?NMO:NAEB;
       for(int i=0; i<NMO; i++) {
        for(int j=0; j<nread; j++) {
          in>>dummy;
          if(j<NAEA) OrbMat(i,j) = dummy;
          if(in.fail()) {
            app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
            in.close();
            return false;
          }
        }
        for(int j=0; j<nread2; j++) {
          in>>dummy;
          if(j<NAEB) OrbMat(i,j+NAEA) = dummy;
          if(in.fail()) {
            app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
            in.close();
            return false;
          }
        }
       }
       for(int i=0; i<NMO; i++) {
        for(int j=0; j<nread; j++) {
          in>>dummy;
          if(j<NAEA) OrbMat(i+NMO,j) = dummy;
          if(in.fail()) {
            app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
            in.close();
            return false;
          }
        }
        for(int j=0; j<nread2; j++) {
          in>>dummy;
          if(j<NAEB) OrbMat(i+NMO,j+NAEA) = dummy;
          if(in.fail()) {
            app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
            in.close();
            return false;
          }
        }
       }

     } else {

      nread = fullMOMat?2*NMO:NAEA+NAEB;
      for(int j=0; j<nread; j++) {
        for(int i=0; i<2*NMO; i++) {
         in>>dummy;
         if(j<NAEA+NAEB) OrbMat(i,j) = dummy;
         if(in.fail()) {
           app_error()<<"Problems reading ASCII file in PureSingleDeterminant. \n";
           in.close();
           return false;
         }
       }
      }

     } // Cstyle
    }

  in.close();

}

// modify for multideterminant case based on type
int readWfn( std::string fileName, boost::multi_array<ComplexType,3>& OrbMat, int NMO, int NAEA, int NAEB, int det = 0)
{
  std::ifstream in;
  in.open(fileName.c_str());
  if(in.fail()) {
     app_error()<<"Problems opening ASCII integral file:  " <<fileName <<std::endl;
     APP_ABORT("Problems opening ASCII integral file. \n");
  }

  bool fullMOMat = false;
  bool Cstyle = true;
  int wfn_type;
  int ndet=1;
  std::string type;

  read_header(in,type,wfn_type,fullMOMat,Cstyle,ndet);

  if(ndet != 1) 
    APP_ABORT("Error: readWfn is for single determinant wave functions. \n");
  if(type != "matrix")
    APP_ABORT("Error: Only type=matrix accepted in readWfn. \n"); 

  ComplexType dummy;
  std::string tag;
  in>>tag >>dummy;
  if(tag != "Coefficients:")
    APP_ABORT(" Error: Expecting Coefficients: tag in wavefunction file. \n");

  int q;
  in>>tag >>q;
  if(tag != "Determinant:" || q!=1)
    APP_ABORT(" Error: Expecting Determinant: 1 tag in wavefunction file. \n");

  if (wfn_type == 0 ) {

     OrbMat.resize(extents[1][NMO][NAEA]);
     read_mat(in,OrbMat[0],Cstyle,fullMOMat,NMO,NAEA);

  } else if(wfn_type == 1) {

    OrbMat.resize(extents[2][NMO][NAEA]);
    read_mat(in,OrbMat[0],Cstyle,fullMOMat,NMO,NAEA);
    read_mat(in,OrbMat[1],Cstyle,fullMOMat,NMO,NAEB);

  } else if(wfn_type == 2) {

    OrbMat.resize(extents[1][2*NMO][NAEA+NAEB]);
    read_mat(in,OrbMat[0][indices[range_t()][range_t(0,NAEA)]],Cstyle,fullMOMat,2*NMO,NAEA);
    read_mat(in,OrbMat[0][indices[range_t()][range_t(NAEA,NAEA+NAEB)]],Cstyle,fullMOMat,2*NMO,NAEB);

  } //type 

  in.close();

  return wfn_type;  

}

}

}
#endif

