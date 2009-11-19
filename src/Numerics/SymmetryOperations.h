//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_SYMMETRYOPERATIONS_H
#define QMCPLUSPLUS_SYMMETRYOPERATIONS_H

#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus {
  
    struct SymmetryGroup
    {
      public:
      SymmetryGroup( string nm="invalid"): name(nm), nClasses(0), nSymmetries(0) 
      {
      }
      
      ~SymmetryGroup()
      {
      }
      
      void addOperator(Matrix<double> op, vector<double> characterlist, int cls)
      {
        SymOps.push_back(op); nSymmetries++;
        Characters.push_back(characterlist);
        Classes.push_back(cls); nClasses=std::max(nClasses,cls+1);
        if (nClasses>characterlist.size())
        {
          app_log()<<" Character table size or class number is wrong."<<endl;
          APP_ABORT("SymmetryGroup::addOperator");
        }
      }
      
      void putClassCharacterTable(vector<vector<double> > cct){ CharacterTableByClass=cct; }
      
      double getsymmetryCharacter(int symmetryOperation, int irrep)
      {
        return Characters[symmetryOperation][irrep];
      }
      
      double getclassCharacter(int clss, int irrep)
      {
        return CharacterTableByClass[clss][irrep];
      }
      
      int getClass(int symmetryOperation) {return Classes[symmetryOperation];}
      
      int getSymmetriesSize() {return nSymmetries;}
      int getClassesSize() {return nClasses;}
      
      void TransformSinglePosition(ParticleSet::ParticlePos_t& oldPos, int symNumber, int el=0)
      {
        ParticleSet::ParticlePos_t rv(1);
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) rv[0][i] += SymOps[symNumber][i][j]*oldPos[el][j];
        oldPos[el]=rv[0];
      }
      
      void TransformAllPositions(ParticleSet::ParticlePos_t& oldPos, int symNumber)
      {
        ParticleSet::ParticlePos_t rv(oldPos.size());
        for (int k=0;k<oldPos.size();k++) for (int i=0;i<3;i++) for (int j=0;j<3;j++) rv[k][i] += SymOps[symNumber][i][j]*oldPos[k][j];
        oldPos=rv;
      }
      private:
      vector<Matrix<double> > SymOps;
      vector<vector<double> > Characters;
      vector<vector<double> > CharacterTableByClass;
      vector<int> Classes;
      int nClasses;
      int nSymmetries;
      string name;
      
    };
    

  /**Builds the symmetry class
  */
  class SymmetryBuilder
  {
  public:
    /// Constructor.
    
    SymmetryBuilder(string symclass, xmlNodePtr q): symname(symclass)
    {
      put(q);
    }
    
    ~SymmetryBuilder(){}; 
    
    SymmetryGroup* getSymmetryGroup()
    {
      return &symgrp;
    }
    
    void put(xmlNodePtr q)
    {
//       string symgrp("invalid");
//       
//       ParameterSet aAttrib;
//       aAttrib.add(symgrp,"group","string");
//       aAttrib.put(q);
      
      //doesn't actually work yet
//       vector<double> principleAxis(3,0);
//       while(kids != NULL) 
//       {
//         string cname((const char*)(kids->name));
//         if(cname == "axis") 
//         {
//           putContent(principleAxis,kids);
//         }
//         kids=kids->next;
//       }
      
      if (symname=="D2H") buildD2H();
      else  if (symname=="C2V") buildC2V();
      else 
      {
        app_log()<<"Symmetry Class "<< symname <<" not yet implemented"<<endl;
        APP_ABORT("SymmetryClass::put");
      }
    }

    
  private:
    SymmetryGroup symgrp;
    string symname;
    
    void buildI(SymmetryGroup& I, vector<double> ctable, int cls)
    {
      Matrix<double> matrix_i(3,3);
      for(int i=0;i<3;i++) matrix_i[i][i]=1.0;
      I.addOperator(matrix_i, ctable, cls);
    };
//    C2V
    void buildC2Vx(SymmetryGroup& C2, vector<double> ctable, int cls)
    {
      Matrix<double> matrix_c2x(3,3);
      matrix_c2x[0][0]=1; matrix_c2x[1][1]=-1; matrix_c2x[2][2]=-1;
      C2.addOperator(matrix_c2x, ctable, cls);
    };
    void buildC2Vy(SymmetryGroup& C2, vector<double> ctable, int cls)
    {
      Matrix<double> matrix_c2y(3,3);
      matrix_c2y[0][0]=1; matrix_c2y[1][2]=1; matrix_c2y[2][1]=-1;
      C2.addOperator(matrix_c2y, ctable, cls);
    };
    void buildC2Vz(SymmetryGroup& C2, vector<double> ctable, int cls)
    {
      Matrix<double> matrix_c2z(3,3);
      matrix_c2z[0][0]=1; matrix_c2z[1][2]=-1; matrix_c2z[2][1]=1;
      C2.addOperator(matrix_c2z, ctable, cls);
    };
//  D2H
    void buildD2Hx(SymmetryGroup& C2, vector<double> ctable, int cls)
    {
      Matrix<double> matrix_c2x(3,3);
      matrix_c2x[0][0]=1; matrix_c2x[1][1]=-1; matrix_c2x[2][2]=-1;
      C2.addOperator(matrix_c2x, ctable, cls);
    };
    void buildD2Hy(SymmetryGroup& C2, vector<double> ctable, int cls)
    {
      Matrix<double> matrix_c2y(3,3);
      matrix_c2y[0][0]=-1; matrix_c2y[1][1]=1; matrix_c2y[2][2]=-1;
      C2.addOperator(matrix_c2y, ctable, cls);
    };
    void buildD2Hz(SymmetryGroup& C2, vector<double> ctable, int cls)
    {
      Matrix<double> matrix_c2z(3,3);
      matrix_c2z[0][0]=-1; matrix_c2z[1][1]=-1; matrix_c2z[2][2]=1;
      C2.addOperator(matrix_c2z, ctable, cls);
    };
    
//     void buildD2(SymmetryGroup& C2, vector<double> ctable, int cls)
//     {
//       Matrix<double> matrix_c2x(3,3), matrix_c2y(3,3), matrix_c2z(3,3);
//       matrix_c2z[0][0]=-1; matrix_c2z[1][1]=-1; matrix_c2z[2][2]=1;
//       matrix_c2y[0][0]=-1; matrix_c2y[1][1]=1; matrix_c2y[2][2]=-1;
//       matrix_c2x[0][0]=1; matrix_c2x[1][1]=-1; matrix_c2x[2][2]=-1;
//       C2.addOperator(matrix_c2x, ctable, cls);
//       C2.addOperator(matrix_c2y, ctable, cls);
//       C2.addOperator(matrix_c2z, ctable, cls);
//       return C2;
//     };
    
//     void buildSigmaD(SymmetryGroup& SD, vector<double> ctable, int cls)
//     {
//       Matrix<double> matrix_sd1(3,3), matrix_sd2(3,3), matrix_sd3(3,3), matrix_sd4(3,3), matrix_sd5(3,3), matrix_sd6(3,3);
//       
//       matrix_sd1[0][0]=1; matrix_sd1[2][1]=1; matrix_sd1[1][2]=1;
//       matrix_sd2[0][0]=1; matrix_sd2[2][1]=-1; matrix_sd2[1][2]=-1;
//       matrix_sd3[0][1]=1; matrix_sd3[1][0]=1; matrix_sd3[2][2]=1;
//       matrix_sd4[0][1]=-1; matrix_sd4[1][0]=-1; matrix_sd4[2][2]=1;
//       matrix_sd5[0][2]=1; matrix_sd5[1][1]=1; matrix_sd5[2][0]=1;
//       matrix_sd6[0][2]=-1; matrix_sd6[1][1]=1; matrix_sd6[2][0]=-1;
//       
//       SD.addOperator(matrix_sd1, ctable, cls);
//       SD.addOperator(matrix_sd2, ctable, cls);
//       SD.addOperator(matrix_sd3, ctable, cls);
//       SD.addOperator(matrix_sd4, ctable, cls);
//       SD.addOperator(matrix_sd5, ctable, cls);
//       SD.addOperator(matrix_sd6, ctable, cls);
//       
//       return SD;
//     };
    
//     void buildS_4(SymmetryGroup& SD, vector<double> ctable, int cls)
//     {
//       Matrix<double> matrix_sd1(3,3), matrix_sd2(3,3), matrix_sd3(3,3), matrix_sd4(3,3), matrix_sd5(3,3), matrix_sd6(3,3);
//       
//       matrix_sd1[0][0]=1; matrix_sd1[2][1]=-1; matrix_sd1[1][2]=1;
//       matrix_sd2[0][0]=1; matrix_sd2[2][1]=1; matrix_sd2[1][2]=-1;
//       matrix_sd3[0][1]=1; matrix_sd3[1][0]=-1; matrix_sd3[2][2]=1;
//       matrix_sd4[0][1]=-1; matrix_sd4[1][0]=1; matrix_sd4[2][2]=1;
//       matrix_sd5[0][2]=1; matrix_sd5[1][1]=1; matrix_sd5[2][0]=-1;
//       matrix_sd6[0][2]=-1; matrix_sd6[1][1]=1; matrix_sd6[2][0]=1;
//       
//       SD.addOperator(matrix_sd1, ctable, cls);
//       SD.addOperator(matrix_sd2, ctable, cls);
//       SD.addOperator(matrix_sd3, ctable, cls);
//       SD.addOperator(matrix_sd4, ctable, cls);
//       SD.addOperator(matrix_sd5, ctable, cls);
//       SD.addOperator(matrix_sd6, ctable, cls);
//       
//       return SD;
//     };
    
    void buildD2H()
    {
      //Character table
      vector<vector<double> > CT(4,vector<double>(4,0));
      CT[0][0]=1; CT[0][1]=1; CT[0][2]=1; CT[0][3]=1;
      CT[1][0]=1; CT[1][1]=1; CT[1][2]=-1; CT[1][3]=-1;
      CT[2][0]=1; CT[2][1]=-1; CT[2][2]=1; CT[2][3]=-1;
      CT[3][0]=1; CT[3][1]=-1; CT[3][2]=-1; CT[3][3]=1;
     
      buildI(symgrp, CT[0],0);
      buildD2Hx(symgrp, CT[1],1);
      buildD2Hy(symgrp, CT[2],2);
      buildD2Hz(symgrp, CT[3],3);
      symgrp.putClassCharacterTable(CT);
    };
    
    void buildC2V()
    {
      //Character table
      vector<vector<double> > CT(4,vector<double>(4,0));
      CT[0][0]=1; CT[0][1]=1; CT[0][2]=1; CT[0][3]=1;
      CT[1][0]=1; CT[1][1]=1; CT[1][2]=-1; CT[1][3]=-1;
      CT[2][0]=1; CT[2][1]=-1; CT[2][2]=1; CT[2][3]=-1;
      CT[3][0]=1; CT[3][1]=-1; CT[3][2]=-1; CT[3][3]=1;
      
      buildI(symgrp, CT[0],0);
      buildC2Vx(symgrp, CT[1],1);
      buildC2Vy(symgrp, CT[2],2);
      buildC2Vz(symgrp, CT[3],3);
      symgrp.putClassCharacterTable(CT);
    };
  };
  
  
  
  
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: kesler $
 * $Revision: 3814 $   $Date: 2009-05-04 13:58:01 -0500 (Mon, 04 May 2009) $
 * $Id: WaveFunctionTester.h 3814 2009-05-04 18:58:01Z kesler $ 
 ***************************************************************************/
