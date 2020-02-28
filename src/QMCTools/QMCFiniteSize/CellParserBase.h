//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCFINITESIZE_CELL_PARSER_H
#define QMCFINITESIZE_CELL_PARSER_H
using namespace qmcplusplus;

class CellParserBase
{
  public:
    CellParserBase();
    ~CellParserBase();
    
    virtual bool put(char* fname);
    
  private:
};
#endif 
