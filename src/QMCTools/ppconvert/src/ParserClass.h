#ifndef PARSER_CLASS_H
#define PARSER_CLASS_H

#include<vector>
#include<string>
#include<complex>
#include<cstdio>
#include<fstream>


class ParserClass
{
public:
  std::streamsize FileSize ( std::string fname);
  virtual bool OpenFile ( std::string fname)  = 0;
  virtual void CloseFile()              = 0;
  virtual bool FindToken ( std::string token) = 0;
  virtual bool ReadInt (int &val)       = 0;
  virtual bool ReadLong (long &val)     = 0;
  virtual bool ReadDouble(double &val)  = 0;
  virtual bool ReadWord ( std::string &word)  = 0;
  virtual bool ReadLine ( std::string &line)  = 0;
  virtual bool NextLine ()              = 0;
  virtual void Reset()                  = 0;
  virtual void SavePos()                = 0;
  virtual void RestorePos()             = 0;
  bool ReadComplex(std::complex<double> &val);  
};

class MemParserClass : public ParserClass 
{
private:
  std::string Buffer;
  int Pos, saved;
public:
  bool OpenFile ( std::string fname);
  void CloseFile ();
  bool FindToken ( std::string token);
  bool ReadInt (int &val);
  bool ReadLong (long &val);
  bool ReadDouble(double &val);
  bool ReadComplex(std::complex<double> &val);
  bool ReadWord ( std::string &word);
  bool ReadLine ( std::string &line);
  bool NextLine ();
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
  std::ifstream Infile;
  std::streampos FileSize;
  std::streampos Pos, saved;
public:
  bool OpenFile ( std::string fname);
  void CloseFile ();
  bool FindToken ( std::string token);
  bool ReadInt (int &val);
  bool ReadLong (long &val);
  bool ReadDouble(double &val);
  bool ReadComplex(std::complex<double> &val);
  bool ReadWord ( std::string &word);
  bool ReadLine ( std::string &line);
  bool NextLine ();
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
  std::ifstream Infile;
  long FileSize;
  long Pos, saved;
  std::string Buffer;
  long BufferStart, MaxBufferSize;
  void ReadChunk (long startpos);
  inline long BufferEnd() 
  { return BufferStart + (long)Buffer.size(); }
public:
  bool OpenFile ( std::string fname);
  void CloseFile ();
  bool FindToken ( std::string token);
  bool ReadInt (int &val);
  bool ReadLong (long &val);
  bool ReadDouble(double &val);
  bool ReadComplex(std::complex<double> &val);
  bool ReadWord ( std::string &word);
  bool ReadLine ( std::string &line);
  bool NextLine ();
  void SavePos();
  void RestorePos();
  void Reset();

  FileParserClass2(int buffSize=16777216) {
    MaxBufferSize = buffSize;
  }
};


#endif
