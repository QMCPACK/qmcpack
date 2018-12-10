//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_SYMMETRYOPERATIONS_H
#define QMCPLUSPLUS_SYMMETRYOPERATIONS_H

#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

struct SymmetryGroup
{
public:
  SymmetryGroup( std::string nm="invalid"): name(nm), nClasses(0), nSymmetries(0)
  {
  }

  ~SymmetryGroup()
  {
  }

  void addOperator(Matrix<double> op, std::vector<double> characterlist, int cls)
  {
    SymOps.push_back(op);
    nSymmetries++;
    Characters.push_back(characterlist);
    Classes.push_back(cls);
    nClasses=std::max(nClasses,cls);
    if (nClasses>characterlist.size())
    {
      app_log()<<" Character table size or class number is wrong."<< std::endl;
      APP_ABORT("SymmetryGroup::addOperator");
    }
  }

  void putClassCharacterTable(std::vector<std::vector<double> > cct)
  {
    CharacterTableByClass=cct;
  }

  double getsymmetryCharacter(int symmetryOperation, int irrep)
  {
    return Characters[symmetryOperation][irrep];
  }

  double getclassCharacter(int clss, int irrep)
  {
    return CharacterTableByClass[clss][irrep];
  }

  int getClass(int symmetryOperation)
  {
    return Classes[symmetryOperation];
  }

  int getSymmetriesSize()
  {
    return nSymmetries;
  }
  int getClassesSize()
  {
    return nClasses;
  }

  void TransformSinglePosition(ParticleSet::ParticlePos_t& oldPos, int symNumber, int el=0)
  {
    ParticleSet::ParticlePos_t rv(1);
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        rv[0][i] += SymOps[symNumber][i][j]*oldPos[el][j];
    oldPos[el]=rv[0];
  }

  void TransformAllPositions(ParticleSet::ParticlePos_t& oldPos, int symNumber)
  {
    ParticleSet::ParticlePos_t rv(oldPos.size());
    for (int k=0; k<oldPos.size(); k++)
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
          rv[k][i] += SymOps[symNumber][i][j]*oldPos[k][j];
    oldPos=rv;
  }
private:
  std::vector<Matrix<double> > SymOps;
  std::vector<std::vector<double> > Characters;
  std::vector<std::vector<double> > CharacterTableByClass;
  std::vector<int> Classes;
  int nClasses;
  int nSymmetries;
  std::string name;

};


/**Builds the symmetry class
*/
class SymmetryBuilder
{
public:
  /// Constructor.

  SymmetryBuilder()
  {
  }

  ~SymmetryBuilder() {};

  SymmetryGroup* getSymmetryGroup()
  {
    return &symgrp;
  }

  void put(xmlNodePtr q)
  {
    symname="";
    xmlNodePtr kids = q->children;
    while(kids != NULL)
    {
      std::string cname((const char*)(kids->name));
      if(cname == "symmetryclass")
      {
        ParameterSet aAttrib;
        aAttrib.add(symname,"name","string");
        aAttrib.put(kids);
      }
      kids=kids->next;
    }
    if (symname=="D2H")
      buildD2H();
    else
      if (symname=="C2V")
        buildC2V();
      else
      {
        buildByHand(q);
//         app_log()<<"Symmetry Class "<< symname <<" not yet implemented"<< std::endl;
//         APP_ABORT("SymmetryClass::put");
      }
  }


private:
  SymmetryGroup symgrp;
  std::string symname;

  void buildI(SymmetryGroup& I, std::vector<double> ctable, int cls)
  {
    Matrix<double> matrix_i(3,3);
    for(int i=0; i<3; i++)
      matrix_i[i][i]=1.0;
    I.addOperator(matrix_i, ctable, cls);
  };
//    C2V
  void buildC2Vx(SymmetryGroup& C2, std::vector<double> ctable, int cls)
  {
    Matrix<double> matrix_c2x(3,3);
    matrix_c2x[0][0]=-1;
    matrix_c2x[1][1]=1;
    matrix_c2x[2][2]=-1;
    C2.addOperator(matrix_c2x, ctable, cls);
  };
  void buildC2Vy(SymmetryGroup& C2, std::vector<double> ctable, int cls)
  {
    Matrix<double> matrix_c2y(3,3);
    matrix_c2y[0][2]=1;
    matrix_c2y[1][1]=1;
    matrix_c2y[2][0]=1;
    C2.addOperator(matrix_c2y, ctable, cls);
  };
  void buildC2Vz(SymmetryGroup& C2, std::vector<double> ctable, int cls)
  {
    Matrix<double> matrix_c2z(3,3);
    matrix_c2z[0][2]=-1;
    matrix_c2z[1][1]=1;
    matrix_c2z[2][0]=-1;
    C2.addOperator(matrix_c2z, ctable, cls);
  };
//  D2H
  void buildD2Hx(SymmetryGroup& C2, std::vector<double> ctable, int cls)
  {
    Matrix<double> matrix_c2x(3,3);
    matrix_c2x[0][0]=1;
    matrix_c2x[1][1]=-1;
    matrix_c2x[2][2]=-1;
    C2.addOperator(matrix_c2x, ctable, cls);
  };
  void buildD2Hy(SymmetryGroup& C2, std::vector<double> ctable, int cls)
  {
    Matrix<double> matrix_c2y(3,3);
    matrix_c2y[0][0]=-1;
    matrix_c2y[1][1]=1;
    matrix_c2y[2][2]=-1;
    C2.addOperator(matrix_c2y, ctable, cls);
  };
  void buildD2Hz(SymmetryGroup& C2, std::vector<double> ctable, int cls)
  {
    Matrix<double> matrix_c2z(3,3);
    matrix_c2z[0][0]=-1;
    matrix_c2z[1][1]=-1;
    matrix_c2z[2][2]=1;
    C2.addOperator(matrix_c2z, ctable, cls);
  };

//     void buildD2(SymmetryGroup& C2, std::vector<double> ctable, int cls)
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

//     void buildSigmaD(SymmetryGroup& SD, std::vector<double> ctable, int cls)
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

//     void buildS_4(SymmetryGroup& SD, std::vector<double> ctable, int cls)
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
    std::vector<std::vector<double> > CT(4,std::vector<double>(4,0));
    CT[0][0]=1;
    CT[0][1]=1;
    CT[0][2]=1;
    CT[0][3]=1;
    CT[1][0]=1;
    CT[1][1]=1;
    CT[1][2]=-1;
    CT[1][3]=-1;
    CT[2][0]=1;
    CT[2][1]=-1;
    CT[2][2]=1;
    CT[2][3]=-1;
    CT[3][0]=1;
    CT[3][1]=-1;
    CT[3][2]=-1;
    CT[3][3]=1;
    buildI(symgrp, CT[0],0);
    buildD2Hx(symgrp, CT[1],1);
    buildD2Hy(symgrp, CT[2],2);
    buildD2Hz(symgrp, CT[3],3);
    symgrp.putClassCharacterTable(CT);
  };

  void buildC2V()
  {
    //Character table
    std::vector<std::vector<double> > CT(4,std::vector<double>(4,0));
    CT[0][0]=1;
    CT[0][1]=1;
    CT[0][2]=1;
    CT[0][3]=1;
    CT[1][0]=1;
    CT[1][1]=1;
    CT[1][2]=-1;
    CT[1][3]=-1;
    CT[2][0]=1;
    CT[2][1]=-1;
    CT[2][2]=1;
    CT[2][3]=-1;
    CT[3][0]=1;
    CT[3][1]=-1;
    CT[3][2]=-1;
    CT[3][3]=1;
    buildI(symgrp, CT[0],0);
    buildC2Vx(symgrp, CT[1],1);
    buildC2Vy(symgrp, CT[2],2);
    buildC2Vz(symgrp, CT[3],3);
    symgrp.putClassCharacterTable(CT);
  };

  void buildByHand(xmlNodePtr q)
  {
    int nClasses(0);
    int nSymmetries(0);
    xmlNodePtr kids = q->children;
    xmlNodePtr symclssptr(NULL);
    while(kids != NULL)
    {
      std::string cname((const char*)(kids->name));
      if(cname == "symmetryclass")
      {
        symclssptr=kids;
        OhmmsAttributeSet aAttrib;
        aAttrib.add(nClasses,"classes" );
        aAttrib.add(nSymmetries,"operators" );
        aAttrib.add(nSymmetries,"symmetries" );
        aAttrib.put(kids);
      }
      kids=kids->next;
    }
    std::vector<std::vector<double> > CT;
    std::vector<std::vector<double> > OPS;
    std::vector<int> clCT;
    kids = symclssptr->children;
    while(kids != NULL)
    {
      std::string cname((const char*)(kids->name));
      if(cname == "charactertable")
      {
        xmlNodePtr kids2 = kids->children;
        while(kids2 != NULL)
        {
          std::string cname2((const char*)(kids2->name));
          if(cname2 == "class")
          {
            std::string clss;
            OhmmsAttributeSet oAttrib;
            oAttrib.add(clss,"name");
            oAttrib.put(kids2);
            std::vector<double> c;
            putContent(c,kids2);
            CT.push_back(c);
          }
          kids2=kids2->next;
        }
      }
      kids=kids->next;
    }
    kids = symclssptr->children;
    while(kids != NULL)
    {
      std::string cname((const char*)(kids->name));
      if (cname=="symmetries")
      {
        xmlNodePtr kids2 = kids->children;
        while(kids2 != NULL)
        {
          std::string cname2((const char*)(kids2->name));
          if(cname2 == "operator")
          {
            int clss(-1);
            OhmmsAttributeSet oAttrib;
            oAttrib.add(clss,"class");
            oAttrib.put(kids2);
            if(clss==-1)
            {
              app_log()<<"  Must label class of operators"<< std::endl;
              APP_ABORT("SymmetryClass::put");
            }
            clCT.push_back(clss);
            std::vector<double> c;
            putContent(c,kids2);
            OPS.push_back(c);
          }
          kids2=kids2->next;
        }
      }
      kids=kids->next;
    }
    //make operator matrices, associate them with the class and add them to the symmetry group.
    for(int i=0; i<OPS.size(); i++)
    {
      Matrix<double> op(3,3);
      for(int j=0; j<3; j++)
        for(int k=0; k<3; k++)
          op(j,k)=OPS[i][j*3 + k];
      symgrp.addOperator(op, CT[clCT[i]-1], clCT[i]);
    }
    symgrp.putClassCharacterTable(CT);
//       for(int j=0; j<CT.size();j++) for(int k=0; k<CT[0].size();k++) app_log()<<CT[j][k]<<" ";
  }


};




}
#endif
