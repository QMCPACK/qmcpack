/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
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
    void ReadGroup (hid_t parentGroupID, string name, IOTreeClass *parent);
    /// StripName strips the trailing ".#" from a string.  These were
    /// added by the HDF5 writer in order to have multiples sections
    /// with the same name.

    void PrintTree(int numIndent );
    void StripName (string str,string &newString,
		    int &myInt);
  public:
    /// This is the HDF5 handle for the group.
    hid_t GroupID;
    IOTreeClass* NewSection(string name);
    /// This prints the variables and sections below me, mostly for
    /// debugging purposes.
  
    IOFileType GetFileType();
    void PrintTree();
    void GroupIterator (string member_name);
    bool OpenFile (string fileName, string mySectionName,
		   IOTreeClass *parent);
    bool NewFile(string fileName,string myName,IOTreeClass* parent);
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
