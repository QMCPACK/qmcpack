#ifndef PARSE_COMMAND_H
#define PARSE_COMMAND_H

#include <map>
#include <string>
#include <iostream>
#include <assert.h>
#include <list>
#include <vector>

using namespace std;

class CommandLineParserClass;

class ParamClass
{
private:
  string Name;
  vector<string> Args;
  bool Found;
  int NumArgsNeeded;
  friend class CommandLineParserClass;

public:
  ///////////////
  // Accessors //
  ///////////////
  string GetName () { return Name; }

  string GetArg  () 
  { 
    assert (NumArgsNeeded==1); 
    return Args[0];
  }

  string GetArg (int arg)
  {
    return Args[arg];
  }

  void SetArg (string arg) {
    assert (Args.size() < NumArgsNeeded);
    Args.push_back(arg);
  }
  
  //////////////////
  // Constructors //
  //////////////////
  ParamClass (string name, int numArgs) {
    Name = name;
    NumArgsNeeded = numArgs;
    Found = false;
  }

  ParamClass (string name, bool needsArg) {
    Name     = name;
    NumArgsNeeded = needsArg ? 1 : 0;
    Found    = false;
  }

  ParamClass() {
    NumArgsNeeded = 0;
    Found    = false;
  }
};

class
CommandLineParserClass
{
private:
  map<string, ParamClass> ArgMap;
  vector<string> Files;
public:
  bool Parse (int argc, char **argv);
  inline bool Found (string name) {
    return ArgMap[name].Found;
  }
  inline string GetArg (string name) {
    return ArgMap[name].GetArg();
  }
  inline string GetArg (string name, int num) {
    return ArgMap[name].GetArg(num);
  }
  inline int NumFiles() { return Files.size(); }
  string GetFile(int i) { return Files[i]; }

  CommandLineParserClass (list<ParamClass> &argList);
};


#endif
