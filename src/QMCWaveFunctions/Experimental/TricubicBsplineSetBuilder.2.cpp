//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/PlaneWave/PWParameterSet.h"
#include "QMCWaveFunctions/TricubicBsplineSetBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/OhmmsBlas.h"
#include "Message/OpenMP.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
namespace qmcplusplus
{

void TricubicBsplineSetBuilder::readComplex2RealDataWithTruncation(const char* hroot,
    const std::vector<int>& occSet, int spinIndex, int degeneracy)
{
  ReportEngine PRE(ClassName,"readComplex2RealDataWithTruncation");
  if(BigDataSet.size())
  {
    int norb=occSet.size()/degeneracy;
    for(int iorb=0; iorb<norb; iorb++)
    {
      std::ostringstream wnshort;
      wnshort<<curH5Fname << "#"<<occSet[iorb] << "#" << spinIndex;
      std::map<std::string,RSOType*>::iterator it(BigDataSet.find(wnshort.str()));
      RSOType& phi(*((*it).second));
      activeBasis->Centers[iorb]=phi.Center;
      activeBasis->Origins[iorb]=phi.Origin;
      activeBasis->add(iorb,phi.Coeffs);
      PRE << "Reusing spline function " << wnshort.str()  <<  '\n';
      PRE << "  center=" << activeBasis->Centers[iorb] << '\n';
      PRE << "  origin=" << activeBasis->Origins[iorb] << '\n';
    }
    return;
  }
  PosType center(0.0),origin(0.0);
  int norb=occSet.size()/degeneracy;
  //input centers
  std::vector<PosType> centerIn(norb);
  //read the centers first
  for(int iorb=0; iorb<norb; iorb++)
  {
    std::string centerName=myParam->getCenterName(hroot,occSet[iorb]);
    HDFAttribIO<PosType > cdummy(centerIn[iorb]);
    cdummy.read(hfileID,centerName.c_str());
  }
  TinyVector<int,DIM> pbc;
  if(OpenEndGrid)
    pbc=BoxGrid;
  else
    pbc=TinyVector<int,DIM>(BoxGrid[0]-1,BoxGrid[1]-1,BoxGrid[2]-1);
  //default grid for the truncated orbitals is the original grid
  TinyVector<int,DIM> bc(pbc);
  //determine the small box: input cutoff + buffer
  if(myParam->BufferRadius>0.0)
  {
    RealType rc=myParam->Rcut+myParam->BufferRadius;
    bc=TinyVector<int,DIM>(2*static_cast<int>(rc/dataKnot.dx)+1,
                           2*static_cast<int>(rc/dataKnot.dy)+1,
                           2*static_cast<int>(rc/dataKnot.dz)+1);
    //need to overwrite the grid: cannot be bigger than the original grid
    for(int idim=0; idim<DIM; idim++)
      bc[idim]=(bc[idim]<pbc[idim])?bc[idim]:pbc[idim];
  }
  //reset the grid with [0,bc*delta] and zero-boudary conditions
  activeBasis->setGrid(0.0,bc[0]*dataKnot.dx, 0.0,bc[1]*dataKnot.dy, 0.0,bc[2]*dataKnot.dz,
                       bc[0],bc[1],bc[2],false,false,false,true);
  app_log() << "    Truncated grid for localized orbitals " << bc << "\n";
  //app_log() << "    Grid-spacing of the input grid " << dataKnot.dx << " " << dataKnot.dy << " " << dataKnot.dz << std::endl;
  //app_log() << "    Grid-spacing of the truncated grid " << activeBasis->bKnots.dx << " "
  //  << activeBasis->bKnots.dy << " " << activeBasis->bKnots.dz << std::endl;
  std::vector<std::vector<int>* > gIndex;
  for(int idim=0; idim<DIM; idim++)
    gIndex.push_back(new std::vector<int>(bc[idim]));
  StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
  Array<ComplexType,3> inTemp(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
  Array<ValueType,3> inCopy(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
  Array<ValueType,3> inTrunc(bc[0],bc[1],bc[2]);
  std::string fname(curH5Fname);
  fname.insert(fname.size()-2,"tr.");
  hid_t h5out = H5Fcreate(fname.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  myParam->writeParameters(h5out);
  hid_t heig=H5Gcreate(h5out,"eigenstates",0);
  //write grid information
  hid_t hg=H5Gcreate(heig,"grid",0);
  int yes=1;
  HDFAttribIO<int> tg(yes);
  tg.write(hg,"translated");
  yes=0;
  tg.write(hg,"closed");
  PosType spacing(dataKnot.dx,dataKnot.dy,dataKnot.dz);
  HDFAttribIO<PosType> spacing_out(spacing);
  spacing_out.write(hg,"spacing");
  HDFAttribIO<TinyVector<int,DIM> > bc_out(bc);
  bc_out.write(hg,"dimensions");
  H5Gclose(hg);
  hid_t htwist=H5Gcreate(heig,"twist0",0);
  HDFAttribIO<PosType> hta(TwistAngle);
  hta.write(htwist,"twist_angle");
  inTrunc = ValueType();
  for(int iorb=0; iorb<norb; iorb++)
  {
    std::ostringstream wnshort;
    std::string eigvName=myParam->getEigVectorName(hroot,occSet[iorb],spinIndex);
    wnshort<<curH5Fname << "#"<<occSet[iorb] << "#" << spinIndex;
    std::map<std::string,RSOType*>::iterator it(BigDataSet.find(wnshort.str()));
    if(it == BigDataSet.end())
    {
      HDFAttribIO<Array<ComplexType,3> > dummy(inTemp);
      dummy.read(hfileID,eigvName.c_str());
      dataKnot.Find(centerIn[iorb][0],centerIn[iorb][1],centerIn[iorb][2]);
      TinyVector<int,DIM> n0(dataKnot.ix0,dataKnot.iy0,dataKnot.iz0);
      TinyVector<int,DIM> nL(n0-bc/2);
      PosType corner;
      for(int idim=0; idim<DIM; idim++)
      {
        std::vector<int>& cIndex(*gIndex[idim]);
        for(int ii=nL[idim], ic=0; ic<bc[idim]; ii++)
        {
          int iloc=ii;
          if(ii<0)
            iloc += pbc[idim];
          else
            if(ii>=pbc[idim])
              iloc -= pbc[idim];
          if(iloc<pbc[idim])
            cIndex[ic++]=iloc;
        }
        //for(int kkk=0; kkk<cIndex.size(); kkk++) std::cout << std::setw(3) << cIndex[kkk];
        //cout << std::endl;
      }
      //Centers[iorb]=activeBasis->Centers[iorb]=PosType(n0[0]*dataKnot.dx,n0[1]*dataKnot.dy,n0[2]*dataKnot.dz);
      activeBasis->Centers[iorb]=center=centerIn[iorb];
      activeBasis->Origins[iorb]=origin=PosType(nL[0]*dataKnot.dx,nL[1]*dataKnot.dy,nL[2]*dataKnot.dz);
      //cout << "  center index "<< n0 << std::endl;
      //  << n0[0]*dataKnot.dx << " " << n0[1]*dataKnot.dy << " " << n0[2]*dataKnot.dz << std::endl
      //  << nL[0]*dataKnot.dx << " " << nL[1]*dataKnot.dy << " " << nL[2]*dataKnot.dz << std::endl
      //  << std::endl;
      for(int ic=1; ic<bc[0]-1; ic++)
        //for(int ic=0; ic<bc[0]; ic++)
      {
        int ic_in=(*gIndex[0])[ic];
        for(int jc=1; jc<bc[1]-1; jc++)
          //for(int jc=0; jc<bc[1]; jc++)
        {
          int jc_in=(*gIndex[1])[jc];
          const ComplexType* restrict inPtr = inTemp.data()+BoxGrid[2]*(jc_in+ic_in*BoxGrid[1]);
          ValueType* restrict outPtr = inTrunc.data()+bc[2]*(jc+ic*bc[1]);
          const std::vector<int>& zIndex(*gIndex[2]);
          for(int kc=1; kc<bc[2]-1; kc++)
            //for(int kc=0; kc<bc[2]; kc++)
          {
            convert(inPtr[zIndex[kc]],outPtr[kc]);
          }
        }
      }
      ///skip writing
      //char newoname[8];
      //sprintf(newoname,"band%d",iorb);
      //hid_t h1= H5Gcreate(htwist,newoname,0);
      //sprintf(newoname,"spin%d",spinIndex);
      //hid_t s1= H5Gcreate(h1,newoname,0);
      //HDFAttribIO<PosType> w1(centerIn[iorb]);
      //w1.write(s1,"center");
      //HDFAttribIO<PosType> w2(activeBasis->Origins[iorb]);
      //w2.write(s1,"origin");
      //HDFAttribIO<Array<ValueType,3> > w3(inTrunc);
      //w3.write(s1,"eigenvector");
      //HDFAttribIO<Array<ValueType,3> > w4(inCopy);
      ////BLAS::copy(inTemp.size(),inTemp.data(),inCopy.data());
      ////w4.write(h1,"eigenvector_full");
      //H5Gclose(s1);
      //H5Gclose(h1);
      StorageType* newP =new StorageType;
      BigDataSet[wnshort.str()]=new RSOType(center,origin,newP);
      activeBasis->add(iorb,inTrunc,newP);
      app_log() << "Reading spline function " << eigvName << " (" << wnshort.str()  << ")"  << std::endl;
      app_log() << "  center=" << activeBasis->Centers[iorb] << std::endl;
      app_log() << "  origin=" << activeBasis->Origins[iorb] << std::endl;
    }
    else
    {
      RSOType& phi(*((*it).second));
      activeBasis->Centers[iorb]=phi.Center;
      activeBasis->Origins[iorb]=phi.Origin;
      activeBasis->add(iorb,phi.Coeffs);
      app_log() << "Reusing spline function " << eigvName << " (" << wnshort.str()  << ")\n";
      app_log() << "  center=" << activeBasis->Centers[iorb] << "\n";
      app_log() << "  origin=" << activeBasis->Origins[iorb] << "\n";
    }
  }
  H5Gclose(htwist);
  H5Gclose(heig);
  H5Fclose(h5out);
  delete_iter(gIndex.begin(),gIndex.end());
  //abort();
}

}
