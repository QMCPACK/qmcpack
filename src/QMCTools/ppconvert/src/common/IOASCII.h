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

#ifndef IO_ASCII_H
#define IO_ASCII_H
#include "IOBase.h"
#include <iostream>
#include <stack>
#include <string>
#include <list>


namespace IO {

  /// This class holds an ASCII token, which is just a string and the
  /// line number in which it appeared in the file.
  class TokenClass
  {
  public:
    string Str;
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
    bool ReadWithoutComments(string fileName, blitz::Array<char,1> &buffer);
    /// Reads a section from a list of TokenClass objects.  iter should
    /// refer to the current place in the list that we should start
    /// reading at.  iter should point to a place just after the '{'.
    /// If wantEndBrace is true, it will look for an ending '}'.
    /// Otherwise it will read until the list of Tokens runs out.  
    bool ReadSection (IOTreeClass *parent, string name,
		      list<TokenClass>::iterator &iter,
		      list<TokenClass> &tokenList,
		      bool wantEndBrace);
  public:
    void WriteSection(ofstream &outFile,int indent);
    IOFileType GetFileType();
    /// Print an indented tree of section variable names.
    void PrintTree(int level);
    /// Same thing, just calls above with level 0;
    void PrintTree();

    IOTreeClass* NewSection(string name);
    void IncludeSection (IOTreeClass *);
    /// Takes the name of a file to read, the name of my section and a
    /// pointer to my parent.  Reads the file into a tree of
    /// IOTreeClass's.
    bool OpenFile (string filename, string myName, 
		   IOTreeClass *parent);
    bool NewFile (string fileName, string mySectionName,
		  IOTreeClass *parent);
    /// Do any file handling necessary and delete the whole tree of data.
    void CloseFile();
    void FlushFile();

    IOTreeASCIIClass()
    { IsModified = false; }
  };

}
#endif
