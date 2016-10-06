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

#ifndef INPUT_OUTPUT_BASE_H
#define INPUT_OUTPUT_BASE_H

#include <string>
#include <list>
#include <stack>
#include "blitz/array.h"
#include <fstream>

#include "IOVar.h"

namespace IO 
{
  /// This class stores a tree of input file sections.  Each member of
  /// the tree contains a list of tree nodes below it and a list of
  /// variables contained in the present node.
  class IOTreeClass
  {
  protected:
    // USE ME!  I'm not being used yet.
    bool IsModified;
    bool UseUnderscores;
  public:
    std::list<IOVarBase*> VarList;
    std::list<IOTreeClass*> SectionList;

    inline void MarkModified();
    /// This is used to ensure proper ordering of sections in the HDF
    /// version in which there is no guarantee that the sections will
    /// come out of the file in the same order you put them in.
    int MyNumber, CurrSecNum;
    virtual void PrintTree()=0;
    virtual void PrintTree(int numIndent)=0;
    virtual IOFileType GetFileType() = 0;

    IOTreeClass* Parent;
    /// This is the empty std::string unless I'm the root node of some file. 
    std::string FileName;
    std::string Name;
    inline void InsertSection (IOTreeClass *newSec);
    inline bool FindSection ( std::string name, IOTreeClass * &sectionPtr, 
			     int num=0);
    inline int CountSections( std::string name);
    inline int CountVars();

    template<class T>
    bool ReadVar( std::string name, T &var)
    {
      bool readVarSuccess;
      std::list<IOVarBase*>::iterator varIter=VarList.begin();
      while ((varIter!=VarList.end()) && ((*varIter)->GetName()!=name)) 
	varIter++;

      bool found = varIter != VarList.end();
      if (found)
	readVarSuccess=(*varIter)->Read(var);
      else if (Parent!=NULL)
	readVarSuccess=Parent->ReadVar(name,var);
      else 
	return false;
    
      return readVarSuccess;	 
    }

    inline IOVarBase* GetVarPtr( std::string name)
    {
      std::list<IOVarBase *>::iterator iter = VarList.begin();
      while ((iter != VarList.end()) && ((*iter)->GetName() != name)){
	iter++;
      }
      if (iter == VarList.end())
	return NULL;
      else {
	MarkModified();  // If we're getting the pointer, we're probably 
	// gonna modify the variable.
	return *iter;
      }
    }
    inline IOVarBase* GetVarPtr(int num)
    {
      std::list<IOVarBase *>::iterator iter = VarList.begin();
      int i=0;
      while ((iter != VarList.end()) && (i != num)){
	iter++;
	i++;
      }
      if (iter == VarList.end())
	return NULL;
      else {
	MarkModified();  // If we're getting the pointer, we're probably 
	// gonna modify the variable.
	return *iter;
      }
    }




    /// Write me!
    virtual IOTreeClass* NewSection( std::string name)=0;
  
    virtual bool OpenFile ( std::string fileName, 
			   std::string mySectionName, 
			   IOTreeClass *parent) = 0;
    virtual bool NewFile ( std::string fileName,
			  std::string mySectionName,
			  IOTreeClass *parent) = 0;
    /// Inserts a new Include directive in the present section.
    virtual void IncludeSection (IOTreeClass *) = 0;
    virtual void CloseFile() = 0;
    virtual void FlushFile() = 0;

    /// These create a new variable with the given name and value:
    template<typename T> bool WriteVar ( std::string name, T val);
    template<typename T, int  LEN> bool  WriteVar ( std::string name, const TinyVector<T,LEN> &val);
    template<typename T, int RANK> bool WriteVar ( std::string name, const blitz::Array<T,RANK> &val);
    template<typename T, int RANK, int LEN> bool WriteVar 
    ( std::string name, const blitz::Array<TinyVector<T,LEN>,RANK> &val);


    /// Append a value to a variable of dimension of 1 higher than val.
    /// i.e. Add a double to an blitz::Array<double,1> or add blitz::Array<double,1>
    /// to an blitz::Array<double,2>, etc.
    template<class T>
    inline bool AppendVar( std::string name, T val);
    inline void SetUnderscores (bool use)
    { UseUnderscores = use; }

    inline IOTreeClass() : UseUnderscores(false), FileName("")
    { 
      // Nothing for now
    }
  };

  void IOTreeClass::MarkModified()
  {
    if (FileName != "")
      IsModified = true;
    else if (Parent!=NULL)
      Parent->MarkModified();
  }


  template<class T>
  inline bool IOTreeClass::AppendVar( std::string name, T val)
  { 
    IOVarBase *var = GetVarPtr(name);
    if (var == NULL)
      return false;
    MarkModified();
    int dim0 = var->GetExtent(0);
    var->Resize(dim0+1);
    return var->Write(val, 0);
    //    return var->Append(val);
  }



  /// Returns the number of subsections with the given name within the
  /// present section.
  inline int IOTreeClass::CountSections( std::string name)
  {
    std::list<IOTreeClass*>::iterator sectionIter;
    sectionIter=SectionList.begin();
    int numSections=0;
    while (sectionIter!=SectionList.end()){
      if ((name==(*sectionIter)->Name) || (name == "")){
	numSections++;
      }
      sectionIter++;
    }
    return numSections;
  }

  inline int IOTreeClass::CountVars()
  {
    std::list<IOVarBase*>::iterator varIter;
    varIter=VarList.begin();
    int numVars=0;
    while (varIter!=VarList.end()){
      numVars++;
      varIter++;
    }
    return numVars;
  }





  /// FindSection locates a subsection with the given name within the
  /// section in contains and returns it in the pointer, sectionPtr,
  /// which is passed a reference.  Returns true if the section is
  /// found.  The final parameter, which default value "true",
  /// optionally resets the section iterator to the beginning of the
  /// section.  Thus, one may control whether or not order is
  /// significant.  
  inline bool IOTreeClass::FindSection ( std::string name, 
					IOTreeClass* &sectionPtr,
					int num)
  {
  
    std::list<IOTreeClass*>::iterator Iter=SectionList.begin();
    int counter=0;
    while(counter<=num && Iter!=SectionList.end()){
      if ((*Iter)->Name==name){
	counter++;
      }
      if (counter<=num)
	Iter++;
    }
    bool found = Iter != SectionList.end(); 
    if (found){
      sectionPtr = *Iter;
    }
    return (found); 
  }



  inline void IOTreeClass::InsertSection(IOTreeClass *newSec)
  {
    std::list<IOTreeClass *>::iterator iter;
  
    if (SectionList.empty())
      SectionList.push_back(newSec);
    else
      {
	iter = SectionList.begin();
	while ((iter != SectionList.end()) && 
	       ((*iter)->MyNumber < newSec->MyNumber))
	  iter++;
	if (iter!=SectionList.end())
	  SectionList.insert(iter, newSec);
	else
	  SectionList.push_back(newSec);
      }
  }


}

#endif
