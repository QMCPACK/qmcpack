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

#ifndef IO_ASCII_H
#define IO_ASCII_H
#include "IOBase.h"
#include <iostream>
#include <stack>
#include <string>
#include <list>


namespace IO {

  /// This class holds an ASCII token, which is just a std::string and the
  /// line number in which it appeared in the file.
  class TokenClass
  {
  public:
    std::string Str;
    int LineNumber;
  };


  /// This is the ASCII specialization of IOTreeClass for ASCII text
  /// files.  It's syntax is as follows:
  /// Section (SectionName)
  /// {
  ///   double x = 3;
  ///   blitz::Array<int,1> y(3) = [1, 2, 3];
  ///   blitz::Array<int,3> z(2,2,1) = [ 1, 2, 
  ///                             3, 4 ];
  ///   Section (Species, "species1.h5");
  /// }
  class IOTreeASCIIClass : public IOTreeClass
  {
    /// Reads a text file into a buffer eliminating c++ and c-style
    /// comments.  
    bool ReadWithoutComments( std::string fileName, blitz::Array<char,1> &buffer);
    /// Reads a section from a list of TokenClass objects.  iter should
    /// refer to the current place in the list that we should start
    /// reading at.  iter should point to a place just after the '{'.
    /// If wantEndBrace is true, it will look for an ending '}'.
    /// Otherwise it will read until the list of Tokens runs out.  
    bool ReadSection (IOTreeClass *parent, std::string name,
		      std::list<TokenClass>::iterator &iter,
		      std::list<TokenClass> &tokenList,
		      bool wantEndBrace);
  public:
    void WriteSection(ofstream &outFile,int indent);
    IOFileType GetFileType();
    /// Print an indented tree of section variable names.
    void PrintTree(int level);
    /// Same thing, just calls above with level 0;
    void PrintTree();

    IOTreeClass* NewSection( std::string name);
    void IncludeSection (IOTreeClass *);
    /// Takes the name of a file to read, the name of my section and a
    /// pointer to my parent.  Reads the file into a tree of
    /// IOTreeClass's.
    bool OpenFile ( std::string filename, std::string myName, 
		   IOTreeClass *parent);
    bool NewFile ( std::string fileName, std::string mySectionName,
		  IOTreeClass *parent);
    /// Do any file handling necessary and delete the whole tree of data.
    void CloseFile();
    void FlushFile();

    IOTreeASCIIClass()
    { IsModified = false; }
  };

}
#endif
