//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "IO.h"
#include <fstream>

namespace IO
{
  std::filebuf vbuf;
  std::ostream verr(&vbuf);

  void SetVerbose (bool verb) {
    if (verb)
      verr.rdbuf(std::cerr.rdbuf());
  }

  /// In the file name format name.extn, returns the extension.
  /// Actually returns everything after the trailing.
  std::string Extension (std::string fileName)
  {
    std::string extn;
    std::stack<char> bwExtn;
    int pos = fileName.length()-1;
    while ((pos >= 0) && fileName[pos]!='.') {
      bwExtn.push(fileName[pos]);
      pos--;
    }
    
    if (fileName[pos] == '.') 
      while (!bwExtn.empty()) {
	extn += bwExtn.top();
	bwExtn.pop();
      }
    else
      extn = "";
    return (extn);
  }

  /// This function takes a filename, determines it extension, creates a
  /// new IOTreeASCIIClass or IOTreeHDF5Class based on the
  /// extension, and calls OpenFile on the new object.
  /// Extensions:  
  /// .h5:            HDF5
  /// .xml:           XML
  /// .anything_else  ASCII
  IOTreeClass *ReadTree (std::string fileName, std::string myName, IOTreeClass *parent)
  {
    IOTreeClass *newTree;
    std::string extn = Extension (fileName);
    newTree = new IOTreeASCIIClass;
    
    newTree->FileName = fileName;
    bool success = newTree->OpenFile (fileName, myName, parent);
    if (success)
      return (newTree);
    else{
      delete newTree;
      return (NULL);
    }
  }
  
  IOTreeClass *NewTree (std::string fileName, std::string myName, IOTreeClass *parent)
  {
    IOTreeClass *newTree;
    std::string extn = Extension (fileName);
    newTree = new IOTreeASCIIClass;
    
    bool success = newTree->NewFile (fileName, myName, parent);
    if (success)
      return (newTree);
    else{
      delete newTree;
      return (NULL);
    }
  }

  bool 
  IOSectionClass::OpenFile (std::string fileName)
  {
    CurrentSection = ReadTree (fileName, "Root", NULL);
    if (CurrentSection == NULL)
      return (false);
    else
      return (true);
  }

  bool 
  IOSectionClass::NewFile (std::string fileName)
  {
    CurrentSection = NewTree (fileName, "Root", NULL);
    if (CurrentSection == NULL)
      return (false);
    else
      return (true);
  }


  void 
  IOSectionClass::CloseFile ()
  {
    while (CurrentSection->Parent != NULL)
      CloseSection();
    CurrentSection->FlushFile();
    CurrentSection->CloseFile();
    delete (CurrentSection);
  }


  void 
  IOSectionClass::FlushFile()
  {
    IOTreeClass *tree = CurrentSection;
    while (tree->Parent != NULL)
      tree = tree->Parent;
    tree->FlushFile();
  }


  bool IOSectionClass::OpenSection (std::string name, int num)
  {
    IOTreeClass *newSection = NULL;
    bool success;
    success = CurrentSection->FindSection(name, newSection, num);
    if (success)
      CurrentSection=newSection;
    return success;
  }


  bool 
  IOSectionClass::OpenSection (int num)
  {
    IOTreeClass *newSection;
    std::list<IOTreeClass*>::iterator Iter=CurrentSection->SectionList.begin();
    int i = 0;
    while ((i<num) && 
	   (Iter != CurrentSection->SectionList.end())){
      i++;
      Iter++;
    }
    if (i<num)
      return false;
    else {
      CurrentSection = *Iter;
      return true;
    }
  }


  bool 
  IOSectionClass::IncludeSection (std::string name, std::string fileName)
  {
    IOTreeClass *newSection;
    newSection = ReadTree (fileName, name, CurrentSection);
    if (newSection == NULL)
      return false;
    else {
      CurrentSection->IncludeSection(newSection);
      return true;
    }
  }

  ///Don't think this pushes to back of list like it should nor does newfile
  bool 
  IOSectionClass::NewSection (std::string name, std::string fileName)
  {
    IOTreeClass *newSection;
    newSection = NewTree (fileName, name, CurrentSection);
    if (newSection == NULL)
      return false;
    else {
      CurrentSection->IncludeSection(newSection);
      CurrentSection = newSection;
      return true;
    }
  }

  void 
  IOSectionClass::CloseSection ()
  {
    //cerr << "Closing Section " << CurrentSection->Name << endl;
    assert (CurrentSection->Parent != NULL);
    CurrentSection = CurrentSection->Parent;
  }


  std::string
  IOSectionClass::GetFileName()
  {
    IOTreeClass *tree = CurrentSection;
    while (tree->Parent != NULL)
      tree = tree->Parent;
    return tree->FileName;
  }





}
