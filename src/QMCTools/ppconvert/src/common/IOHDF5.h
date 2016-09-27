//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef INPUT_OUTPUT_HDF5_H
#define INPUT_OUTPUT_HDF5_H
#include "IOBase.h"
#include <iostream>
#include <stack>
#include <hdf5.h>


namespace IO {

  /// This class stores a section of an HDF5 file.  The boolean value,
  /// IsRoot, store whether this particular section is a the root node
  /// of an HDF5 file.
  class IOTreeHDF5Class : public IOTreeClass
  {
  private:
    bool IsOpen;
    hid_t BoolType, ComplexType;
    /// ReadGroup reads a HDF5 group, given by name, from the file.
    /// It reads in all variables and groups within the file, calling
    /// itself recursively for groups within itself.
    void ReadGroup (hid_t parentGroupID, std::string name, IOTreeClass *parent);
    /// StripName strips the trailing ".#" from a string.  These were
    /// added by the HDF5 writer in order to have multiples sections
    /// with the same name.

    void PrintTree(int numIndent );
    void StripName ( std::string str,std::string &newString,
		    int &myInt);
  public:
    /// This is the HDF5 handle for the group.
    hid_t GroupID;
    IOTreeClass* NewSection( std::string name);
    /// This prints the variables and sections below me, mostly for
    /// debugging purposes.
  
    IOFileType GetFileType();
    void PrintTree();
    void GroupIterator ( std::string member_name);
    bool OpenFile ( std::string fileName, std::string mySectionName,
		   IOTreeClass *parent);
    bool NewFile( std::string fileName,std::string myName,IOTreeClass* parent);
    void IncludeSection (IOTreeClass *);
    void CloseFile();
    void FlushFile();
    hid_t GetBoolType()    { return BoolType; }
    hid_t GetComplexType() { return ComplexType; }
    IOTreeHDF5Class() : IOTreeClass()
    {
      IsOpen=false;
      CurrSecNum=0;
    }
  };

}


#endif
