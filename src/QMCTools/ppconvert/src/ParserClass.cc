//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "ParserClass.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cassert>

std::streamsize ParserClass::FileSize(std::string fname)
{
  std::ifstream infile;
  infile.open(fname.c_str(), std::ifstream::in);
  if (!infile.is_open())
    return false;
  std::streampos fileSize = 0;
  infile.seekg(0, std::ios_base::end);
  fileSize = infile.tellg();
  infile.close();
  return fileSize;
}

bool 
MemParserClass::OpenFile(std::string fname)
{
  std::ifstream infile;
  infile.open(fname.c_str(), std::ifstream::in);
  if (!infile.is_open())
    return false;
  std::streampos fileSize = 0;
  infile.seekg(fileSize, std::ios_base::end);
  fileSize = infile.tellg();
  infile.seekg(0, std::ios_base::beg);
  if (fileSize != (std::streampos)-1) {
    Buffer.resize(fileSize);
    infile.read(&(Buffer[0]), fileSize);
    infile.close();
  }
  else {
    int fd = open (fname.c_str(), O_RDONLY);
    char ch[100000];
    int len = read (fd, ch, 99999);
    ch[len] = '\0';
    close (fd);
    Buffer = ch;
  }
  Pos = 0;
  return true;
}

void
MemParserClass::CloseFile()
{
  Buffer.clear();
}

bool
MemParserClass::FindToken(std::string token)
{
  bool found = false;
  int toklen = token.size();
  int tokenPos = Buffer.find(token, Pos);
  if (tokenPos == -1)
    return false;
  Pos = tokenPos+token.size();
  return true;
}


bool 
MemParserClass::ReadInt (int &val)
{
//   int numChars;
//   int success = sscanf (&(Buffer[Pos]), " %d %n", &val, &numChars);
//   if (success) 
//     Pos += numChars;
//   return (success == 1);
  char * endptr;
  val = strtol (&(Buffer[Pos]), &endptr, 10);
  if (endptr == &(Buffer[Pos]))
    return false;
  Pos += endptr - &(Buffer[Pos]);
  return (true);  
}

bool 
MemParserClass::ReadLong (long &val)
{
  char * endptr;
  val = strtol (&(Buffer[Pos]), &endptr, 10);
  if (endptr == &(Buffer[Pos]))
    return false;
  Pos += endptr - &(Buffer[Pos]);
  return (true);  
}


bool 
MemParserClass::ReadDouble (double &val)
{
  char *endptr;
  val =strtod (&(Buffer[Pos]), &endptr);
  if (endptr == &(Buffer[Pos]))
    return false;
  Pos += endptr - &(Buffer[Pos]);
  return (true);
}

void
MemParserClass::SavePos()
{  saved = Pos;  }


void
MemParserClass::RestorePos()
{  Pos = saved;  }

void
FileParserClass::SavePos()
{  saved = Pos;  }

void
FileParserClass::RestorePos()
{  
  Pos = saved;
  Infile.seekg(Pos, std::ios_base::end);
}


bool
ParserClass::ReadComplex (std::complex<double> &val)
{
  double re, im;
  if (FindToken ("("))
    if (ReadDouble(re))
      if (FindToken(","))
	if (ReadDouble(im))
	  if (FindToken(")")) {
	    val = std::complex<double>(re,im);
	    return true;
	  }
  return false;
}

bool isWhiteSpace (char c)
{
  return ((c==' ') || (c=='\t') || (c=='\n') || (c=='\r'));
}

bool
MemParserClass::ReadWord (std::string &word)
{
  word = "";
  char str[2];
  str[1] = '\0';
  while (isWhiteSpace (Buffer[Pos]) && (Pos<(Buffer.size()-1)))
    Pos++;
  while (!isWhiteSpace(Buffer[Pos]) && (Pos<Buffer.size()-1)) {
    str[0] = Buffer[Pos];
    word.append(str);
    Pos++;
  }
  return true;
}

bool
MemParserClass::ReadLine (std::string &word)
{
  word = "";
  char str[2];
  str[1] = '\0';
  while ((Buffer[Pos]!='\n') && (Pos<Buffer.size()-1)) {
    str[0] = Buffer[Pos];
    word.append(str);
    Pos++;
  }
  return true;
}

bool
MemParserClass::NextLine ()
{
  while ((Buffer[Pos]!='\n') && (Pos<Buffer.size()-1))
    Pos++;
  if (Pos < Buffer.size()) {
    Pos++;
    return true;
  }
  else
    return false;
}

void
MemParserClass::Reset()
{
  Pos = 0;
}



bool
FileParserClass::OpenFile (std::string fname)
{
  Infile.open(fname.c_str(), std::ifstream::in);
  if (!Infile.is_open())
    return false;
  FileSize = 0;
  Infile.seekg(FileSize, std::ios_base::end);
  FileSize = Infile.tellg();
  Infile.seekg((std::streampos) 0, std::ios_base::beg);
  Pos = 0;
  return true;
}

void
FileParserClass::CloseFile()
{
  if (!Infile.is_open()) 
    std::cerr << "Tried to close a FileParserClass that's not open.\n";
  else
    Infile.close();
}

bool
FileParserClass::FindToken(std::string token)
{
  assert (Infile.is_open());
  char compare[token.size()+1];
  Pos = Infile.tellg();
  bool found = false;
  while (!found) {
    Infile.seekg(Pos, std::ios_base::beg);
    Infile.get(compare, (std::streamsize)token.size()+1, '\0');
    if (token == compare) {
      Pos += token.size();
      Infile.seekg(Pos, std::ios_base::beg);
      return true;
    }
    if (Pos >= FileSize)
      return false;
    Pos+=1;
  }
  return false;
}

bool
FileParserClass::ReadInt (int &val)
{
  Infile >> val;
  return !Infile.fail();
}

bool
FileParserClass::ReadLong(long &val)
{
  Infile >> val;
  return !Infile.fail();
}

bool
FileParserClass::ReadDouble(double &val)
{
  Infile >> val;
  return !Infile.fail();
}

bool
FileParserClass::ReadComplex (std::complex<double> &val)
{
  double re, im;
  if (FindToken ("("))
    if (ReadDouble(re))
      if (FindToken(","))
	if (ReadDouble(im))
	  if (FindToken(")")) {
	    val = std::complex<double>(re,im);
	    return true;
	  }
  return false;
}

bool
FileParserClass::ReadWord (std::string &word)
{
  word = "";
  char ch;
  char str[2];
  str[1] = '\0';
  while (isWhiteSpace(ch = Infile.get()) && !Infile.eof());
  if (Infile.eof())
    return false;
  while (!isWhiteSpace(ch = Infile.get()) && !Infile.eof()) {
    str[0] = ch;
    word.append (str);
  }
  if (isWhiteSpace(ch))
    Infile.unget();
  return true;
}

bool
FileParserClass::ReadLine (std::string &line)
{
  line = "";
  char ch;
  char str[2];
  str[1] = '\0';
  if (Infile.eof())
    return false;
  while (((ch = Infile.get())!='\n') && !Infile.eof()) {
    str[0] = ch;
    line.append (str);
  }
  return true;
}

bool
FileParserClass::NextLine ()
{
  if (Infile.eof())
    return false;
  while (Infile.get()!='\n' && !Infile.eof());

  return true;
}

void
FileParserClass::Reset()
{
  assert (Infile.is_open());
  Pos = 0;
  Infile.seekg (Pos, std::ios_base::beg);
}




////////////////////////////////////////////////////////////
//                   FileParserClass2                     //
////////////////////////////////////////////////////////////
bool
FileParserClass2::OpenFile(std::string fname)
{
  Infile.open(fname.c_str(), std::ifstream::in);
  if (!Infile.is_open())
    return false;
  FileSize = 0;
  Infile.seekg(FileSize, std::ios_base::end);
  FileSize = Infile.tellg();
  Infile.seekg((std::streampos) 0, std::ios_base::beg);
  Pos = 0;
  ReadChunk (0);
  return true;
}

void
FileParserClass2::Reset()
{
  assert (Infile.is_open());
  Pos = 0;
  Infile.seekg (Pos, std::ios_base::beg);
  ReadChunk(0);
}


void 
FileParserClass2::CloseFile()
{
  if (!Infile.is_open()) 
    std::cerr << "Tried to close a FileParserClass that's not open.\n";
  else
    Infile.close();
  Buffer.resize(0);
}



void
FileParserClass2::ReadChunk (long start)
{
  long n = std::min(MaxBufferSize, FileSize-start);
  if (Buffer.size() != n)
    Buffer.resize(n);
  Infile.seekg(start, std::ios_base::beg);
  Infile.read(&Buffer[0], n);
  BufferStart = start;
  Pos = start;
}

bool
FileParserClass2::FindToken (std::string token)
{
  bool fileEnd = false;
  if (Pos < BufferStart)
    ReadChunk (Pos);
  do {
    long tokenPos = Buffer.find (token, Pos-BufferStart);
    if (tokenPos != -1) {
      Pos = tokenPos + BufferStart + token.size();
      return true;
    }
    else if (BufferStart + Buffer.size() >= FileSize) {
      return false;
    }
    else {
      ReadChunk (BufferEnd()- (token.size()+1));
      Pos = BufferStart;
    }
  } while (true);
}

void
FileParserClass2::SavePos()
{  saved = Pos;  }

void
FileParserClass2::RestorePos()
{  
  Pos = saved;
  if (Pos < BufferStart)
    ReadChunk (Pos);
  else if (BufferStart + Buffer.size() < Pos)
    ReadChunk (Pos);
}



bool 
FileParserClass2::ReadInt(int &val)
{
  if ((BufferEnd() - Pos) < 100)
    ReadChunk (Pos);
  char *endptr;
  long offset = Pos - BufferStart;
  val = strtol (&(Buffer[offset]), &endptr, 10);
  if (endptr == &(Buffer[offset]))
    return false;
  Pos += endptr - &(Buffer[offset]);
  return (true);  
}

bool 
FileParserClass2::ReadLong(long &val)
{
  if (BufferEnd() - Pos < 100)
    ReadChunk (Pos);
  char *endptr;
  long offset = Pos - BufferStart;
  val = strtol (&(Buffer[offset]), &endptr, 10);
  if (endptr == &(Buffer[offset]))
    return false;
  Pos += endptr - &(Buffer[offset]);
  return (true);  
}


bool 
FileParserClass2::ReadDouble(double &val)
{
  if (BufferEnd() - Pos < 100)
    ReadChunk (Pos);
  char *endptr;
  long offset = Pos - BufferStart;
  val = strtod (&(Buffer[offset]), &endptr);
  if (endptr == &(Buffer[offset])) {
    std::cerr << "Couldn't file double.\n";
    return false;
  }
  Pos += endptr - &(Buffer[offset]);
  return (true);  
}

bool
FileParserClass2::ReadWord (std::string &word)
{
  if (BufferEnd() - Pos < 100)
    ReadChunk (Pos);
   char str[2];
  str[1] = '\0';
  long offset = Pos - BufferStart;
  while (isWhiteSpace (Buffer[offset]) && (offset<(Buffer.size()-1))) {
    offset++; 
    Pos++;
  }

  while (!isWhiteSpace(Buffer[offset]) && (offset<Buffer.size()-1)) {
    str[0] = Buffer[offset];
    word.append(str);
    offset++;
    Pos++;
  }
  return true;
}

bool
FileParserClass2::ReadLine (std::string &line)
{
  if (BufferEnd() - Pos < 100)
    ReadChunk (Pos);
  char str[2];
  str[1] = '\0';
  long offset = Pos - BufferStart;

  while ((Buffer[offset] != '\n') && (offset<Buffer.size()-1)) {
    str[0] = Buffer[offset];
    line.append(str);
    offset++;
    Pos++;
  }
  return true;
}

bool
FileParserClass2::NextLine ()
{
  if (BufferEnd() - Pos < 100)
    ReadChunk (Pos);
  long offset = Pos - BufferStart;

  while ((Buffer[offset] != '\n') && (offset<Buffer.size()-1)) {
    offset++;
    Pos++;
  }
  return true;
}
