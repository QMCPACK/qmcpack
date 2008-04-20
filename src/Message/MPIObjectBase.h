//////////////////////////////////////////////////////////////////
// (c) Copyright 2008- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#ifndef QMCPLUSPLUS_MPIOBJECTBASE_H
#define QMCPLUSPLUS_MPIOBJECTBASE_H
#include "Message/Communicate.h"

namespace APPNAMESPACE
{

  class MPIObjectBase
  {
    public:

      typedef Communicate::mpi_comm_type mpi_comm_type;
      typedef Communicate::intra_comm_type intra_comm_type;

      ///default constructor
      MPIObjectBase(Communicate* c=0);
      
      ///virtual destructor
      virtual ~MPIObjectBase();

      /** initialize myComm
       * @param c communicator
       */
      void initCommunicator(Communicate* c);

      ///return the rank of this node
      inline int rank() const {return myComm->rank();}

      inline Communicate* getCommunicator() const 
      {
        return myComm;
      }

      inline mpi_comm_type getMPI() const 
      { 
        return myComm->getMPI();
      }

      /** return true if the rank == 0
      */
      inline bool is_manager() const 
      {
        return !myComm->rank();
      }

    protected:
      Communicate* myComm;

    private:
      //disable copy constructor for now
      MPIObjectBase(const MPIObjectBase& a){}
  };

}
#endif // QMCPLUSPLUS_MPIOBJECTBASE_H
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2468 $   $Date: 2008-02-22 09:27:30 -0500 (Fri, 22 Feb 2008) $
 * $Id: Communicate.h 2468 2008-02-22 14:27:30Z jnkim $ 
 ***************************************************************************/
