#ifndef PARSER_CLASS_H
#define PARSER_CLASS_H

#include<vector>
#include<string>
#include<complex>
#include<cstdio>
#include<fstream>

using namespace std;

class ParserClass
{
public:
  streamsize FileSize (string fname);
  virtual bool OpenFile (string fname)  = 0;
  virtual void CloseFile()              = 0;
  virtual bool FindToken (string token) = 0;
  virtual bool ReadInt (int &val)       = 0;
  virtual bool ReadLong (long &val)     = 0;
  virtual bool ReadDouble(double &val)  = 0;
  virtual bool ReadWord (string &word)  = 0;
  virtual bool ReadLine (string &line)  = 0;
  virtual void Reset()                  = 0;
  virtual void SavePos()                = 0;
  virtual void RestorePos()             = 0;
  bool ReadComplex(complex<double> &val);  
};

class MemParserClass : public ParserClass 
{
private:
  string Buffer;
  int Pos, saved;
public:
  bool OpenFile (string fname);
  void CloseFile ();
  bool FindToken (string token);
  bool ReadInt (int &val);
  bool ReadLong (long &val);
  bool ReadDouble(double &val);
  bool ReadComplex(complex<double> &val);
  bool ReadWord (string &word);
  bool ReadLine (string &line);
  void SavePos();
  void RestorePos();
  inline void Reset();

  MemParserClass() {
    // do nothing for now
  }
};




class FileParserClass : public ParserClass
{
private:
  ifstream Infile;
  streampos FileSize;
  streampos Pos, saved;
public:
  bool OpenFile (string fname);
  void CloseFile ();
  bool FindToken (string token);
  bool ReadInt (int &val);
  bool ReadLong (long &val);
  bool ReadDouble(double &val);
  bool ReadComplex(complex<double> &val);
  bool ReadWord (string &word);
  bool ReadLine (string &line);
  void SavePos();
  void RestorePos();
  void Reset();

  FileParserClass() {
    // do nothing for now
  }
};

class FileParserClass2 : public ParserClass
{
private:
  ifstream Infile;
  long FileSize;
  long Pos, saved;
  string Buffer;
  long BufferStart, MaxBufferSize;
  void ReadChunk (long startpos);
  inline long BufferEnd() 
  { return BufferStart + (long)Buffer.size(); }
public:
  bool OpenFile (string fname);
  void CloseFile ();
  bool FindToken (string token);
  bool ReadInt (int &val);
  bool ReadLong (long &val);
  bool ReadDouble(double &val);
  bool ReadComplex(complex<double> &val);
  bool ReadWord (string &word);
  bool ReadLine (string &line);
  void SavePos();
  void RestorePos();
  void Reset();

  FileParserClass2(int buffSize=16777216) {
    MaxBufferSize = buffSize;
  }
};


#endif
