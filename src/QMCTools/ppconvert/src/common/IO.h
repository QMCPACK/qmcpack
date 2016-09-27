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

#ifndef IO_H
#define IO_H

#include "IOBase.h"
#include "IOVar.h"
#include "IOASCII.h"

namespace IO {
  // This output stream is for verbose output from programs.  It is
  // connected to stderr if SetVerbose(true) is called.
  extern std::ostream verr;
  void SetVerbose(bool verb);

  template<> inline bool 
  IOTreeClass::WriteVar ( std::string name, const char* val)
  {
    return WriteVar(name, std::string(val));
  }
  template<typename T> bool 
  IOTreeClass::WriteVar ( std::string name, T val)
  {
    if (GetFileType() == ASCII_TYPE) {
      VarList.push_back(new IOVarASCII<T,0>(name,val));
    }
    else {
      std::cerr << "Unknown file type in WriteVar.\n";
      abort();
    }
    return true;
  }

  template<typename T, int LEN> bool
  IOTreeClass::WriteVar ( std::string name, const TinyVector<T,LEN> &val)
  {
    Array<T,1> aVal(LEN);
    for (int i=0; i<LEN; i++)
      aVal(i) = val[i];
    WriteVar (name, aVal);
  }

    

  template<typename T, int RANK> bool 
  IOTreeClass::WriteVar ( std::string name, const blitz::Array<T,RANK> &val)
  {
    if (GetFileType() == ASCII_TYPE) {
      VarList.push_back(new IOVarASCII<T,RANK>(name, val));
    }
    else {
      std::cerr << "Unknown file type in WriteVar.\n";
      abort();
    }
    return true;
  }

  template<typename T, int RANK, int LEN> bool
  IOTreeClass::WriteVar ( std::string name, const blitz::Array<TinyVector<T,LEN>,RANK> &val)
  {
    TinyVector<int,RANK+1> shape;
    for (int dim=0; dim<RANK; dim++)
      shape[dim] = val.extent(dim);
    shape[RANK] = LEN;

    Array<T,RANK+1> aval((T*)&(val(0)[0]), shape, blitz::neverDeleteData);
    return WriteVar (name, aval);
  }



  /// In the file name format name.extn, returns the extension.
  /// Actually returns everything after the trailing.
  inline std::string Extension ( std::string fileName);


  /// This function takes a filename, determines it extension, creates a
  /// new IOTreeASCIIClass or IOTreeHDF5Class based on the
  /// extension, and calls OpenFile on the new object.
  /// Extensions:  
  /// .h5:            HDF5
  /// .xml:           XML
  /// .anything_else  ASCII
  IOTreeClass *ReadTree ( std::string fileName, std::string myName, IOTreeClass *parent);

  IOTreeClass *NewTree ( std::string fileName, std::string myName, IOTreeClass *parent);



  ///  Wrapper class for IOTreeClass that gives a nearly identical
  ///  interface as the OutputSectionClass.
  class IOSectionClass
  {
  private:
    IOTreeClass *CurrentSection;
  public:

    /// Opens the file reference by fileName and reads the contents into
    /// the tree in CurrentSection.  Creates a new object based on the
    /// extnesion of the filename.  For ".h5", it creates an
    /// IOTreeHDF5Class.  For ".xml" it creaes an IOTreeXMLClass.
    /// After creating the object, it calls the objects virtual OpenFile
    /// function, reading the contents of the file into the tree.
    bool OpenFile ( std::string fileName);
    std::string GetName(){ return CurrentSection->Name;}
    std::string GetFileName();
    std::string GetVarName(int num){ return GetVarPtr(num)->GetName();}
    /// Creates a file at the top level, choosing the appropriate type
    /// based on the file extension.
    bool NewFile ( std::string fileName);

    /// Calls CurrentSections close file and then deletes the
    /// CurrentSection.  
    void CloseFile ();

    /// Flush all buffers to disk for safety
    void FlushFile();

    /// Opens the num'th section with the given name.  The default
    /// value for num is 0.
    bool OpenSection ( std::string name, int num=0);

    /// Opens the num'th section below CurrentSection.
    bool OpenSection (int num);

    /// This mounts a file in the current tree under CurrentSection at
    /// the end of CurrentsSection's SectionList.  It does not change
    /// what CurrentSection points to, ie. it does not descend to the
    /// newly-opened section.
    bool IncludeSection ( std::string name, std::string fileName);

    /// Creates a new section of the same type as currentSection under
    /// currentSection.  Pushes the new section to the end of the
    /// section list.
    inline void NewSection ( std::string name)
    {  CurrentSection = CurrentSection->NewSection(name); }

    /// This function creates a new file of the appropriate type as
    /// determined by the extension of fileName and mounts it at the end
    /// of the list under CurrentSection.  Returns false if the file
    /// couldn't be created.
    bool NewSection ( std::string name, std::string fileName);

    /// Closes the current section.  That is, CurrentSection becomes
    /// CurrentSection's parent.
    void CloseSection ();

    /// Template function which reads a variable in the present section
    /// into the passed-by-reference T variable.
    template<class T>
    bool ReadVar( std::string name, T &var)
    {  return (CurrentSection->ReadVar(name, var)); }

    template<class T>
    bool ReadVar( std::string name, T &var, T Default)
    { 
      bool success = ReadVar(name, var);
      if (!success)
	var = Default;
      return (success);
    }

    /// Writes a variable under the current section.
    template<typename T> bool
    WriteVar ( std::string name, T val)
    { return CurrentSection->WriteVar(name, val); }

    template<typename T, int RANK> bool
    WriteVar ( std::string name, const blitz::Array<T,RANK>& val)
    { return CurrentSection->WriteVar(name, val); }
  
    template<class T> bool
    AppendVar( std::string name, T val)
    { return CurrentSection->AppendVar(name, val); }

    template<typename T, int RANK> bool
    AppendVar ( std::string name, const blitz::Array<T,RANK>& val)
    { return CurrentSection->AppendVar(name, val); }
  
    inline IOVarBase *GetVarPtr( std::string name)
    {    return (CurrentSection->GetVarPtr(name)); }

    inline IOVarBase *GetVarPtr(int num)
    { return (CurrentSection->GetVarPtr(num)); }

    inline void SetUnderscores(bool use)
    { CurrentSection->SetUnderscores(use); }

    /// Returns the number of subsections within the present section
    /// which have the name name.  If called without a name, it returns
    /// the total number of sections.
    inline int CountSections( std::string name="")
    { return (CurrentSection->CountSections(name)); }
    inline int CountVars()
    {return (CurrentSection->CountVars());}
    /// Calls CurrentSections virtual PrintTree() function.  This is for
    /// debugging purposes.  It spits out a hierarchy of the sections
    /// and variable names.
    void PrintTree()
    { CurrentSection->PrintTree(); }

    IOSectionClass(IOSectionClass &io)
    {
      CurrentSection = io.CurrentSection;
    }

    IOSectionClass() 
    {
      CurrentSection=NULL;
    }
	
  };


} // Ends namespace IO

#endif
