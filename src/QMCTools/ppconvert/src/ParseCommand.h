#ifndef PARSE_COMMAND_H
#define PARSE_COMMAND_H

#include <map>
#include <string>
#include <iostream>
#include <assert.h>
#include <list>
#include <vector>


class CommandLineParserClass;

class ParamClass
{
private:
  std::string Name;
  std::vector<std::string> Args;
  bool Found;
  int NumArgsNeeded;
  friend class CommandLineParserClass;

public:
  ///////////////
  // Accessors //
  ///////////////
  std::string GetName () { return Name; }

  std::string GetArg  () 
  { 
    assert (NumArgsNeeded==1); 
    return Args[0];
  }

  std::string GetArg (int arg)
  {
    return Args[arg];
  }

  void SetArg ( std::string arg) {
    assert (Args.size() < NumArgsNeeded);
    Args.push_back(arg);
  }
  
  //////////////////
  // Constructors //
  //////////////////
  ParamClass ( std::string name, int numArgs) {
    Name = name;
    NumArgsNeeded = numArgs;
    Found = false;
  }

  ParamClass ( std::string name, bool needsArg) {
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
  std::map<std::string, ParamClass> ArgMap;
  std::vector<std::string> Files;
public:
  bool Parse (int argc, char **argv);
  inline bool Found ( std::string name) {
    return ArgMap[name].Found;
  }
  inline std::string GetArg ( std::string name) {
    return ArgMap[name].GetArg();
  }
  inline std::string GetArg ( std::string name, int num) {
    return ArgMap[name].GetArg(num);
  }
  inline int NumFiles() { return Files.size(); }
  std::string GetFile(int i) { return Files[i]; }

  CommandLineParserClass (list<ParamClass> &argList);
};


#endif
