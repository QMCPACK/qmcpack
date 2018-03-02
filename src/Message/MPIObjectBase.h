//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file MPIObjectBase.h
 * @brief declaration of MPIObjectBase
 */

#ifndef QMCPLUSPLUS_MPIOBJECTBASE_H
#define QMCPLUSPLUS_MPIOBJECTBASE_H
#include "Message/Communicate.h"

namespace qmcplusplus
{

/** Base class for any object which needs to know about a MPI communicator.
 */
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

  ///return the rank of the communicator
  inline int rank() const
  {
    return myComm->rank();
  }

  ///return the group id of the communicator
  inline int getGroupID() const
  {
    return myComm->getGroupID();
  }

  ///return myComm
  inline Communicate* getCommunicator() const
  {
    return myComm;
  }

  ///return MPI communicator if one wants to use MPI directly
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

  /** default is 1 minal
   */
  void setReportLevel(int level=1);

  ///return the name
  inline const std::string& getName() const
  {
    return myName;
  }

  inline void setName(const std::string& aname)
  {
    myName = aname;
  }

protected:
  ///level of report
  int ReportLevel;
  /** pointer to Communicate
   * @todo use smart pointer
   */
  Communicate* myComm;
  /** class Name
   */
  std::string ClassName;
  /** name of the object */
  std::string myName;

  //private:
  //  //disable copy constructor for now
  //  MPIObjectBase(const MPIObjectBase& a);
};

}
#endif // QMCPLUSPLUS_MPIOBJECTBASE_H
