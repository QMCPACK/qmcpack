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

namespace IO {

  IOFileType IOTreeASCIIClass::GetFileType()
  {
    return ASCII_TYPE;
  }

  /// Simply prints 3*num spaces
  inline void ASCIIPrintIndent(int num)
  {
    for (int counter=0;counter<num*2;counter++)
      cout<<' ';
  }

  /// Simply prints 3*num spaces
  inline void ASCIIPrintIndent(int num,ofstream &outFile)
  {
    for (int counter=0;counter<num*2;counter++)
      outFile<<' ';
  }


  /// Prints an indented hierarchy of sections and variable names to
  /// cout. 
  void IOTreeASCIIClass::PrintTree(int indentNum)
  {
    ASCIIPrintIndent(indentNum);
    cout<<"Section: "<<Name<<endl;
    list<IOVarBase*>::iterator varIter=VarList.begin();
    while (varIter!=VarList.end()){
      ASCIIPrintIndent(indentNum+1);
      cout<<"Variable: "<<(*varIter)->GetName()<<" "<<endl;
      varIter++;
    }
    list<IOTreeClass*>::iterator secIter=SectionList.begin();
    while (secIter!=SectionList.end()){
      //    cout<<"Section: "<<(*secIter)->Name<<endl;
      (*secIter)->PrintTree(indentNum+1);
      secIter++;
    }
  }

  /// Calls PrintTree(0)
  void IOTreeASCIIClass::PrintTree()
  {
    PrintTree(0);
  }


  /// Returns true if theChar is a special character that should be
  /// parsed into its own token.
  bool isSpecial(char theChar)
  {
    return ( (theChar=='(') ||
	     (theChar==')') ||
	     (theChar=='{') ||
	     (theChar=='}') ||
	     (theChar=='[') ||
	     (theChar==']') ||
	     (theChar=='<') ||
	     (theChar=='>') ||	
	     (theChar=='=') ||
	     (theChar==';') ||
	     (theChar==','));
  }
	   
  /// Returns true if theChar is a space, newline, tab, or carriage return.      
  bool isWhiteSpace(char theChar)
  {
    return ( (theChar=='\n') ||
	     (theChar==' ' ) ||
	     (theChar=='\t') ||
	     (theChar=='\r'));
  }
      

  /// Returns true if theChar is a letter or underscore		      
  bool isAlpha(char theChar)
  {
    return ((theChar>='a' && theChar<='z') || (theChar>='A' && theChar<='Z')
	    ||theChar=='_');
  }

  /// Returns true if theChar is a digit
  bool isDigit(char theChar)
  {
    return (theChar>='0' && theChar<='9');
  }

  /// Returns true if theChar is the a valid character for starting a
  /// number.  Includes a digit, a '.' or a '-'
  bool isNumStart(char theChar)
  {
    return ((isDigit(theChar)) || (theChar=='.') || (theChar=='-'));
  }

  /// Returns true if ch is a valid character comprising a number.
  bool isNumChar (char ch)
  {
    return (isDigit(ch) || (ch =='.') || (ch=='e') 
	    || (ch=='-') || (ch=='+'));
  }


  /// Tokenize takes an array of characters and constructs a list of
  /// TokenClass objects.  Each token has a string and a line number.
  /// Valid tokens are special characters: "(){}[]<>,", quoted strings,
  /// words, or numbers.  White space is not significant, except in
  /// separating tokens.
  void Tokenize(blitz::Array<char,1> buffer, list<TokenClass>& tokenList)
  {
    int pos=0;
    int lineNum=1;
    while (pos<buffer.size()){
      if (isSpecial(buffer(pos))){
	TokenClass tempToken;
	tempToken.Str+=buffer(pos);
	tempToken.LineNumber = lineNum;
	tokenList.push_back(tempToken);
	pos++;
      }
      else if (buffer(pos)=='\"'){
	TokenClass tempToken;
	tempToken.Str="\"";
	pos++;
	while (buffer(pos)!='\"'){
	  tempToken.Str+=buffer(pos);
	  pos++;
	}
	pos++;
	tempToken.Str+='\"';
	tempToken.LineNumber = lineNum;
	tokenList.push_back(tempToken);
      }
      else if (isAlpha(buffer(pos))){
	TokenClass tempToken;
	while ((isAlpha(buffer(pos)) || isDigit(buffer(pos)))){
	  tempToken.Str+=buffer(pos);
	  pos++;
	}
	tempToken.LineNumber = lineNum;
	tokenList.push_back(tempToken);
      }
      else if (isWhiteSpace(buffer(pos))){
	if (buffer(pos)=='\n')
	  lineNum++;
	pos++;
      }
      else if (isNumStart(buffer(pos))){
	TokenClass tempToken;
	while (isNumChar(buffer(pos))){
	  tempToken.Str+=buffer(pos);
	  pos++;
	}
	tempToken.LineNumber = lineNum;
	tokenList.push_back(tempToken);
      }
      else if (buffer(pos)=='\0')
	break;
      else {
	cerr<<"There was a token we do not recognize in line "<<lineNum<<endl;
	cerr <<"The rest of the file is as follows:\n";
	while (pos<buffer.size()) {
	  cerr << (int)buffer(pos);
	  pos++;
	}
	exit(1);
      }
    }
  }	     
	  

  /// Just a shortcut to look at two characters at a time.
  bool checkPair(blitz::Array<char,1> &buffer,int counter,char* toSee)
  {
    if (counter+1>=buffer.size()){
      return false;
    }
    if (buffer(counter)==toSee[0] && buffer(counter+1)==toSee[1]){
      return true;
    }
    else return false;

  }

  /// Reads a file into a character array, removing C and C++ style
  /// comments. 
  bool IOTreeASCIIClass::ReadWithoutComments(string fileName,
					     blitz::Array<char,1> 
					     &buffer)
  {
    ifstream infile;
    infile.open(fileName.c_str());
    if (!infile.is_open()) 
      return false;
    blitz::Array<char,1> tmpBuffer;
    int counter=0;
    bool inQuote=false;
    char dummyChar;
    while (!infile.eof()){
      infile.get(dummyChar);    
      counter++;
    }

    tmpBuffer.resize(counter);
    buffer.resize(counter);
    counter=-1;
    infile.close();
    ifstream infile2;
    infile2.open(fileName.c_str());
    while (!infile2.eof()){
      counter++;
      infile2.get(dummyChar);
      tmpBuffer(counter)=dummyChar;    
    }
    //  cout<<tmpBuffer;
    int bufferLoc=0;
    for (int counter=0;counter<tmpBuffer.size();counter++){
      if (inQuote){
	while ( counter<tmpBuffer.size() && tmpBuffer(counter) != '\"'){
	  buffer(bufferLoc)=tmpBuffer(counter);
	  counter++;
	  bufferLoc++;
	}
	buffer(bufferLoc)=tmpBuffer(counter);
	bufferLoc++;
	inQuote=false;
      }
      else {
	if (checkPair(tmpBuffer,counter,"//")){
	  while (tmpBuffer(counter)!='\n' && counter<tmpBuffer.size()){
	    counter++;
	  }
	  buffer(bufferLoc)=tmpBuffer(counter); //copy the \n over
	  bufferLoc++;
	}
	else if (checkPair(tmpBuffer,counter,"/*")){
	  while (!checkPair(tmpBuffer,counter,"*/") && counter<tmpBuffer.size()){
	    counter++;
	    if (tmpBuffer(counter)=='\n'){
	      buffer(bufferLoc)=tmpBuffer(counter);
	      bufferLoc++;
	    }		   
	  }
	  counter++; //end up in the / of comment
	}
	else if (tmpBuffer(counter)=='\"'){
	  inQuote=true;
	  buffer(bufferLoc)=tmpBuffer(counter);
	  bufferLoc++;
	}
	else {
	  buffer(bufferLoc)=tmpBuffer(counter);
	  bufferLoc++;
	}
      }
    }
    buffer.resizeAndPreserve(bufferLoc);
    infile2.close();
    return (true);
  }



  /// If isError is true, then we print out an error message giving the
  /// line number and the string passed to us in ErrorStr.
  inline void ReadAbort (bool isError, int lineNumber, string ErrorStr)
  {
    if (isError) {
      cerr << "Error in input file at line number " << lineNumber 
	   << ":\n";
      cerr << ErrorStr;
      exit(1);
    }
  }

  /// Removes all double quotes from the input string and return it.
  string StripQuote(string str)
  {
    string newString;
    int i=0;
    assert (str[0] == '\"');
    assert (str[str.length()-1] == '\"');
    while (i<str.length())
      {
	if (str[i] != '\"')
	  newString += str[i];
	i++;
      }
    return newString;
  }


  /// Looks at the string passed to it and returns the corresponding
  /// enumerated type.  If the type is not recognized, it returns
  /// INVALID.  
  IODataType GetType (string typeString)
  {
    if (typeString=="double")
      return DOUBLE_TYPE;
    else if (typeString=="int")
      return INT_TYPE;
    else if (typeString=="string")
      return STRING_TYPE;
    else if (typeString=="bool")
      return BOOL_TYPE;
    else if (typeString=="complex")
      return COMPLEX_TYPE;
    else return INVALID;
  }


  /// Takes a token and reads its value into a double, aborting if there
  /// is a problem.
  void ReadAtomicVar(TokenClass token,double &d)
  {

    char* endPtr;
    d=strtod(token.Str.c_str(),&endPtr);
    ReadAbort(*endPtr!='\0',token.LineNumber,"Expected Double\n");
  }

  /// Takes a token and reads its value into an int, aborting if there
  /// is a problem.
  void ReadAtomicVar(TokenClass token,int &d)
  {

    char* endPtr;
    d=strtol(token.Str.c_str(),&endPtr,10);
    ReadAbort(*endPtr!='\0',token.LineNumber,"Expected Int\n");
  }

  /// Takes a token and reads its value into a string, aborting if there
  /// is a problem.
  void ReadAtomicVar(TokenClass token,string &d)
  {

    ReadAbort (token.Str[0] != '\"', token.LineNumber, 
	       "Expected '\"'.");
    ReadAbort (token.Str[token.Str.length()-1] != '\"', token.LineNumber, 
	       "Expected '\"'.");
    d=StripQuote(token.Str);
  }

  /// Takes a token and reads its value into a bool, aborting if there
  /// is a problem.
  void ReadAtomicVar(TokenClass token, bool &b)
  {
    if (token.Str=="true"){
      b=true;
    }
    else if (token.Str=="false"){
      b=false;
    }
    else ReadAbort(true,token.LineNumber,"Expected true or false\n");
  }

  /// Takes a token and reads its value into a bool, aborting if there
  /// is a problem.
  void ReadAtomicVar(TokenClass token, complex<double> &a)
  {
    cerr << "Reading complex not yet implemented.";
  }

  /// This template function reads a 1-D array from a token list into
  /// the array.  The syntax requires an opening '[' the a
  /// comma-separated list of values, then a closing ']'.
  template <class T>
  void ReadArrayData(list<TokenClass>::iterator &iter,
		     list<TokenClass> &tokenList,
		     blitz::Array<T,1> valArray)
  {
    ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
    iter++;
    for (int counter=0;counter<valArray.extent(0)-1;counter++){
      ReadAtomicVar(*iter,valArray(counter));
      iter++;
      ReadAbort(iter->Str != ",", iter->LineNumber, "Expected , not found\n");
      iter++;
    }
    //Read last value
    ReadAtomicVar(*iter,valArray(valArray.extent(0)-1));
    iter++;
    ReadAbort(iter->Str != "]", iter->LineNumber, "Expected ] not found\n");
    iter++;
    ReadAbort(iter->Str != ";", iter->LineNumber, "Expected ; not found\n");
    iter++;
  }


  /// This template function reads a 2-D array from a token list into
  /// the array.  The syntax requires an opening '[' the a
  /// comma-separated list of values, then a closing ']'.  The data is
  /// read row-ordered, i.e. the first index changes fastests as we read
  /// in the values.
  template <class T>
  void ReadArrayData(list<TokenClass>::iterator &iter,
		     list<TokenClass> &tokenList,
		     blitz::Array<T,2> valArray)
  {
    ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
    iter++;
    for (int j=0; j<valArray.extent(0); j++)
      for (int i=0;i<valArray.extent(1);i++) {
	ReadAtomicVar(*iter,valArray(j,i));
	iter++;
	// Read comma if this isn't the last value.
	if ((i!=valArray.extent(1)-1) || 
	    (j!=valArray.extent(0)-1)) {
	  ReadAbort(iter->Str != ",", iter->LineNumber, 
		    "Expected , not found\n");
	  iter++;
	}
      }
    ReadAbort(iter->Str != "]", iter->LineNumber, "Expected ] not found\n");
    iter++;
    ReadAbort(iter->Str != ";", iter->LineNumber, "Expected ; not found\n");
    iter++;
  }


  /// This template function reads a 3-D array from a token list into
  /// the array.  The syntax requires an opening '[' the a
  /// comma-separated list of values, then a closing ']'.  The data is
  /// read row-ordered, i.e. the first index changes fastests as we read
  /// in the values.
  template <class T>
  void ReadArrayData(list<TokenClass>::iterator &iter,
		     list<TokenClass> &tokenList,
		     blitz::Array<T,3> valArray)
  {
    ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
    iter++;
    for (int k=0; k<valArray.extent(0); k++)
      for (int j=0; j<valArray.extent(1); j++)
	for (int i=0;i<valArray.extent(2);i++) {
	  ReadAtomicVar(*iter,valArray(k,j,i));
	  iter++;
	  // Read comma if this isn't the last value.
	  if ((i!=valArray.extent(2)-1) || 
	      (j!=valArray.extent(1)-1) ||
	      (k!=valArray.extent(0)-1)) {
	    ReadAbort(iter->Str != ",", iter->LineNumber, 
		      "Expected , not found\n");
	    iter++;
	  }
	}
    ReadAbort(iter->Str != "]", iter->LineNumber, "Expected ] not found\n");
    iter++;
    ReadAbort(iter->Str != ";", iter->LineNumber, "Expected ; not found\n");
    iter++;
  }


  /// This template function reads a 4-D array from a token list into
  /// the array.  The syntax requires an opening '[' the a
  /// comma-separated list of values, then a closing ']'.  The data is
  /// read row-ordered, i.e. the first index changes fastests as we read
  /// in the values.
  template <class T>
  void ReadArrayData(list<TokenClass>::iterator &iter,
		     list<TokenClass> &tokenList,
		     blitz::Array<T,4> valArray)
  {
    ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
    iter++;
    for (int k=0; k<valArray.extent(0); k++) 
      for (int j=0; j<valArray.extent(1); j++)
	for (int i=0;i<valArray.extent(2);i++) 
	  for (int h=0;h<valArray.extent(3);h++) {
	    ReadAtomicVar(*iter,valArray(k,j,i,h));
	    iter++;
	    // Read comma if this isn't the last value.
	    if ((h!=valArray.extent(3)-1) || 
		(i!=valArray.extent(2)-1) || 
		(j!=valArray.extent(1)-1) ||
		(k!=valArray.extent(0)-1)) {
	      ReadAbort(iter->Str != ",", iter->LineNumber, 
			"Expected , not found\n");
	      iter++;
	    }
	  }
    ReadAbort(iter->Str != "]", iter->LineNumber, "Expected ] not found\n");
    iter++;
    ReadAbort(iter->Str != ";", iter->LineNumber, "Expected ; not found\n");
    iter++;
  }


  IOVarBase *NewASCIIVar (string name, IODataType newType, int ndim,
			  blitz::Array<int,1> dims)
  {
    if (ndim == 0) {
      if (newType == DOUBLE_TYPE)
	return new IOVarASCII<double,0>(name);
      else if (newType == INT_TYPE)
	return new IOVarASCII<int,0>(name);
      else if (newType == STRING_TYPE)
	return new IOVarASCII<string,0>(name);
      else if (newType == BOOL_TYPE)
	return new IOVarASCII<bool,0>(name);
      else if (newType == COMPLEX_TYPE)
	return new IOVarASCII<complex<double>,0>(name);
    }
    else if (ndim == 1) {
      if (newType == DOUBLE_TYPE) {
	IOVarASCII<double,1> *newVar = new IOVarASCII<double,1>(name);
	newVar->ArrayValue.resize(dims(0));
	return newVar;
      }
      else if (newType == INT_TYPE) {
	IOVarASCII<int,1> *newVar = new IOVarASCII<int,1>(name);
	newVar->ArrayValue.resize(dims(0));
	return newVar;
      }
      else if (newType == STRING_TYPE) {
	IOVarASCII<string,1> *newVar = new IOVarASCII<string,1>(name);
	newVar->ArrayValue.resize(dims(0));
	return newVar;
      }
      else if (newType == BOOL_TYPE) {
	IOVarASCII<bool,1> *newVar = new IOVarASCII<bool,1>(name);
	newVar->ArrayValue.resize(dims(0));
	return newVar;
      }
      else if (newType == COMPLEX_TYPE) {
	IOVarASCII<complex<double>,1> *newVar = new IOVarASCII<complex<double>,1>(name);
	newVar->ArrayValue.resize(dims(0));
	return newVar;
      }
    }
    else if (ndim == 2) {
      if (newType == DOUBLE_TYPE) {
	IOVarASCII<double,2> *newVar = new IOVarASCII<double,2>(name);
	newVar->ArrayValue.resize(dims(0), dims(1));
	return newVar;
      }
      else if (newType == INT_TYPE) {
	IOVarASCII<int,2> *newVar = new IOVarASCII<int,2>(name);
	newVar->ArrayValue.resize(dims(0), dims(1));
	return newVar;
      }
      else if (newType == STRING_TYPE) {
	IOVarASCII<string,2> *newVar = new IOVarASCII<string,2>(name);
	newVar->ArrayValue.resize(dims(0), dims(1));
	return newVar;
      }
      else if (newType == BOOL_TYPE) {
	IOVarASCII<bool,2> *newVar = new IOVarASCII<bool,2>(name);
	newVar->ArrayValue.resize(dims(0), dims(1));
	return newVar;
      }
      else if (newType == COMPLEX_TYPE) {
	IOVarASCII<complex<double>,2> *newVar = new IOVarASCII<complex<double>,2>(name);
	newVar->ArrayValue.resize(dims(0), dims(1));
	return newVar;
      }
    }  
    else if (ndim == 3) {
      if (newType == DOUBLE_TYPE) {
	IOVarASCII<double,3> *newVar = new IOVarASCII<double,3>(name);
	newVar->ArrayValue.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
      else if (newType == INT_TYPE) {
	IOVarASCII<int,3> *newVar = new IOVarASCII<int,3>(name);
	newVar->ArrayValue.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
      else if (newType == STRING_TYPE) {
	IOVarASCII<string,3> *newVar = new IOVarASCII<string,3>(name);
	newVar->ArrayValue.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
      else if (newType == BOOL_TYPE) {
	IOVarASCII<bool,3> *newVar = new IOVarASCII<bool,3>(name);
	newVar->ArrayValue.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
      else if (newType == COMPLEX_TYPE) {
	IOVarASCII<complex<double>,3> *newVar = new IOVarASCII<complex<double>,3>(name);
	newVar->ArrayValue.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
    }
    else if (ndim == 4) {
      if (newType == DOUBLE_TYPE) {
	IOVarASCII<double,4> *newVar = new IOVarASCII<double,4>(name);
	newVar->ArrayValue.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
      else if (newType == INT_TYPE) {
	IOVarASCII<int,4> *newVar = new IOVarASCII<int,4>(name);
	newVar->ArrayValue.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
      else if (newType == STRING_TYPE) {
	IOVarASCII<string,4> *newVar = new IOVarASCII<string,4>(name);
	newVar->ArrayValue.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
      else if (newType == BOOL_TYPE) {
	IOVarASCII<bool,4> *newVar = new IOVarASCII<bool,4>(name);
	newVar->ArrayValue.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
      else if (newType == COMPLEX_TYPE) {
	IOVarASCII<complex<double>,4> *newVar = new IOVarASCII<complex<double>,4>(name);
	newVar->ArrayValue.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
    }
    else
      return NULL;
    return NULL;
  }



  /// Reads an array from a list of tokens, starting at the token
  /// pointed to by iter.  It places the array into the newVar object.
  /// It expects to begin reading after the word "Array".
  IOVarBase * ReadArray(list<TokenClass>::iterator &iter,
			list<TokenClass> &tokenList)
  {
    ReadAbort(iter->Str != "<", iter->LineNumber, "Expected < not found\n");
    iter++;
    IODataType myType=GetType(iter->Str);
    ReadAbort(myType==INVALID,iter->LineNumber,
	      "Array does not have atomic type\n");
    iter++;
    ReadAbort(iter->Str != ",", iter->LineNumber, "Expected , not found\n");
    iter++;
    int numDim;
    ReadAtomicVar(*iter,numDim);
    iter++;
    ReadAbort(iter->Str != ">", iter->LineNumber, "Expected , not found\n");
    iter++;
  
    blitz::Array<int,1> dimSize(numDim);
  
    string myName=iter->Str;

    iter++;
    ReadAbort(iter->Str != "(", iter->LineNumber, "Expected ( not found\n");
    iter++;
    for (int counter=0;counter<numDim-1;counter++){
      ReadAtomicVar(*iter,dimSize(counter));
      iter++;
      ReadAbort(iter->Str != ",", iter->LineNumber, "Expected , not found\n");
      iter++;
    }

    //Read the last dimension
    ReadAtomicVar(*iter,dimSize(numDim-1));
    iter++;
    ReadAbort(iter->Str != ")", iter->LineNumber, "Expected ) not found\n");
    iter++;
    ReadAbort(iter->Str!="=",iter->LineNumber,"Expected = not found\n");
    iter++;
  
    IOVarBase *newVar = NewASCIIVar (myName, myType, numDim, dimSize);
  
    if (numDim == 1) {
      if (myType == DOUBLE_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<double,1> *)newVar)->ArrayValue);
      else if (myType == INT_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<int,1> *)   newVar)->ArrayValue);
      else if (myType == STRING_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<string,1> *)newVar)->ArrayValue);
      else if (myType == BOOL_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<bool,1> *)  newVar)->ArrayValue);
      else if (myType == COMPLEX_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<complex<double>,1> *)  newVar)->ArrayValue);
    }
    if (numDim == 2) {
      if (myType == DOUBLE_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<double,2> *)newVar)->ArrayValue);
      else if (myType == INT_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<int,2> *)   newVar)->ArrayValue);
      else if (myType == STRING_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<string,2> *)newVar)->ArrayValue);
      else if (myType == BOOL_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<bool,2> *)  newVar)->ArrayValue);
      else if (myType == COMPLEX_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<complex<double>,2> *) newVar)->ArrayValue);
    }
    if (numDim == 3) {
      if (myType == DOUBLE_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<double,3> *)newVar)->ArrayValue);
      else if (myType == INT_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<int,3> *)   newVar)->ArrayValue);
      else if (myType == STRING_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<string,3> *)newVar)->ArrayValue);
      else if (myType == BOOL_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<bool,3> *)  newVar)->ArrayValue);
      else if (myType == COMPLEX_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<complex<double>,3> *)  newVar)->ArrayValue);
    }
    if (numDim == 4) {
      if (myType == DOUBLE_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<double,4> *)newVar)->ArrayValue);
      else if (myType == INT_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<int,4> *)   newVar)->ArrayValue);
      else if (myType == STRING_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<string,4> *)newVar)->ArrayValue);
      else if (myType == BOOL_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<bool,4> *)  newVar)->ArrayValue);
      else if (myType == COMPLEX_TYPE)
	ReadArrayData(iter, tokenList, ((IOVarASCII<complex<double>,4> *)  newVar)->ArrayValue);
    }
    return (newVar);
  }



  /// This function parses a variable assigment from the list of tokens,
  /// creates a new VarASCIIClass object and puts the appropriate value
  /// in that object.  It recognizes any of the atomic types or arrays
  /// of theose atomic types.
  IOVarBase* ReadASCIIVar (list<TokenClass>::iterator &iter,
			       list<TokenClass> &tokenList)
  {
    IOVarBase *newVar;
    IODataType myType=GetType(iter->Str);
    if (myType==INVALID){
      ReadAbort(iter->Str!="Array",iter->LineNumber,
		"Invalid Type: "+iter->Str+"\n");
      iter++;
      newVar = ReadArray(iter,tokenList);
    }
    else {
      iter++;
      string myName=iter->Str;

      iter++;
      ReadAbort(iter->Str!="=",iter->LineNumber,"Expected equals sign\n");
      iter++;
      TokenClass valToken=*iter;
      iter++;
      ReadAbort(iter->Str!=";",iter->LineNumber,"Expected semicolon\n");
      iter++;
      blitz::Array<int,1> dims(1);
      newVar = NewASCIIVar (myName, myType, 0, dims);
      if (myType == DOUBLE_TYPE)
	ReadAtomicVar (valToken, ((IOVarASCII<double,0> *)newVar)->Value);
      else if (myType == INT_TYPE)
	ReadAtomicVar (valToken, ((IOVarASCII<int,0> *)newVar)->Value);
      else if (myType == STRING_TYPE)
	ReadAtomicVar (valToken, ((IOVarASCII<string,0> *)newVar)->Value);
      else if (myType == BOOL_TYPE)
	ReadAtomicVar (valToken, ((IOVarASCII<bool,0> *)newVar)->Value);
      else if (myType == COMPLEX_TYPE)
	ReadAtomicVar (valToken, ((IOVarASCII<complex<double>,0> *)newVar)->Value);
    }
    return(newVar);
  }


  /// ReadSection parses a section in the input file.  It takes as
  /// arguments this sections parent, its name, a tokenlist iterator,
  /// the tokenlist, and a bool telling us if we want to look for a "}"
  /// at the end of the input.  If we don't, we keep parsing until the
  /// buffer runs out.  Calls itself recursively as necessary, builing
  /// up a tree of sections and variables.
  bool IOTreeASCIIClass::ReadSection (IOTreeClass *parent,
				      string myName,
				      list<TokenClass>::iterator &iter,
				      list<TokenClass> &tokenList,
				      bool wantEndBrace)
  {
    Parent = parent;
    Name = myName;
    while ((iter != tokenList.end()) && (iter->Str != "}")) {
      if (iter->Str == "Section") {
	IOTreeClass *newTree;
	iter++;
	ReadAbort(iter->Str != "(", iter->LineNumber, "Expected ( not found\n");
	iter++;
	string newName = iter->Str;
	iter++;
	// Check for included section
	if (iter->Str == ",") {
	  // Get filename
	  iter++;
	  string fileName = StripQuote(iter->Str);
	  iter++;
	  ReadAbort (iter->Str!=")", iter->LineNumber, "Expected ) not found\n");
	  iter++;
	  ReadAbort (iter->Str!=";", iter->LineNumber, "Expected ; not found\n");
	  iter++;	
	  newTree = ReadTree (fileName, newName, this);
	}
	else {
	  ReadAbort(iter->Str != ")", iter->LineNumber, 
		    "Expected ) not found\n");
	  iter++;
	  ReadAbort(iter->Str != "{", iter->LineNumber, 
		    "Expected { not found\n");
	  iter++;
	  newTree = new IOTreeASCIIClass();
	  ((IOTreeASCIIClass*)newTree)->ReadSection((IOTreeClass*)this,
						    newName,iter,tokenList,
						    true);         
	}
	SectionList.push_back(newTree);
      }
      else {
	IOVarBase *newVar =  ReadASCIIVar(iter, tokenList);
	VarList.push_back(newVar);
      }
    }
    if ((iter==tokenList.end()) && wantEndBrace) {
      cerr << "Unexpected end of file before } \n";
      exit (1);
    }
	    
    if (iter!=tokenList.end())  
      iter++;
    return (true);
  }

  IOTreeClass* IOTreeASCIIClass::NewSection(string name)
  {
    IOTreeClass* tempSection=new IOTreeASCIIClass();
    tempSection->Name=name;
    tempSection->Parent=this;
    tempSection->MyNumber=CurrSecNum;
    tempSection->SetUnderscores(UseUnderscores);
    CurrSecNum++;
    SectionList.push_back(tempSection);
    MarkModified();
    return tempSection;
  }

  void IOTreeASCIIClass::IncludeSection(IOTreeClass *newSection)
  {
    newSection->MyNumber=CurrSecNum++;
    SectionList.push_back(newSection);
    MarkModified();
  }



  bool IOTreeASCIIClass::NewFile (string fileName,
				  string mySectionName,
				  IOTreeClass *parent)
  {
    FileName=fileName;
    Parent=parent;
    Name=mySectionName;
    return true;
  }



  /// OpenFile takes a filename to open, the name of this section and
  /// the parent of this section.  It reads the file into a buffer,
  /// converts it to a list of tokens, then parses the tokens,
  /// constructing a tree of sections containing variables lists.  
  bool IOTreeASCIIClass::OpenFile(string fileName, string myName, 
				  IOTreeClass *parent)
  {
    blitz::Array<char,1> buffer;
    bool success = ReadWithoutComments(fileName,buffer);
    if (!success)
      return false;
    list<TokenClass> tokenList;
    Tokenize(buffer,tokenList);
    list<TokenClass>::iterator iter=tokenList.begin();
    ReadSection(parent,myName,iter,tokenList, false);
    return true;
  }



  void IOTreeASCIIClass::WriteSection(ofstream &outFile,int indentNum)
  {
    list<IOVarBase*>::iterator varIter=VarList.begin();
    while (varIter!=VarList.end()){
      ASCIIPrintIndent(indentNum,outFile);
      (*varIter)->Print(outFile);    
      varIter++;
    }
    list<IOTreeClass*>::iterator secIter=SectionList.begin();
    while (secIter!=SectionList.end()){
      if ((*secIter)->FileName==""){
	ASCIIPrintIndent(indentNum,outFile);
	outFile<<"Section ("<<(*secIter)->Name<<")\n";
	ASCIIPrintIndent(indentNum,outFile);
	outFile<<"{\n";
	((IOTreeASCIIClass*)(*secIter))->WriteSection(outFile,indentNum+1);
	ASCIIPrintIndent(indentNum,outFile);
	outFile<<"}\n\n";
      }
      else {
	ASCIIPrintIndent(indentNum,outFile);
	outFile<<"Section ("<<(*secIter)->Name<<", \"";
	outFile<<(*secIter)->FileName<<"\");"<<endl;
	(*secIter)->FlushFile();
      }
      secIter++;
    }
  }



  void IOTreeASCIIClass::FlushFile()
  {
    ofstream outfile;
    if ((FileName!="") && IsModified){
      outfile.open(FileName.c_str());
      WriteSection(outfile,0);
    }

    list<IOTreeClass*>::iterator iter = SectionList.begin();
    while (iter != SectionList.end()) {
      (*iter)->FlushFile();
      iter++;
    }
  }

  /// CloseFile recursively destroys the tree of data we constructed.
  void IOTreeASCIIClass::CloseFile()
  {  
    // First, free all the variables in the list
    while (!VarList.empty()) {
      delete(VarList.front());
      VarList.pop_front();
    }
  
    // Now, call all closes recursively and delete all sections
    while (!SectionList.empty())
      {
	SectionList.front()->CloseFile();
	delete SectionList.front();
	SectionList.pop_front();
      }
  }

  int 
  IOVarASCII<double,0>::GetRank()
  { return 0; }
  IODataType 
  IOVarASCII<double,0>::GetType()
  { return DOUBLE_TYPE; }
  IOFileType 
  IOVarASCII<double,0>::GetFileType()
  { return ASCII_TYPE; }
  int
  IOVarASCII<double,0>::GetExtent(int i)
  { return 1; }
  void
  IOVarASCII<double,0>::Resize(int n)
  {
    cerr << "Cannot resize atomic variable.\n";
    abort();
  }

  int 
  IOVarASCII<int,0>::GetRank()
  { return 0; }
  IODataType 
  IOVarASCII<int,0>::GetType()
  { return INT_TYPE; }
  IOFileType 
  IOVarASCII<int,0>::GetFileType()
  { return ASCII_TYPE; }
  int
  IOVarASCII<int,0>::GetExtent(int i)
  { return 1; }
  void
  IOVarASCII<int,0>::Resize(int n)
  {
    cerr << "Cannot resize atomic variable.\n";
    abort();
  }

  int 
  IOVarASCII<string,0>::GetRank()
  { return 0; }
  IODataType 
  IOVarASCII<string,0>::GetType()
  { return STRING_TYPE; }
  IOFileType 
  IOVarASCII<string,0>::GetFileType()
  { return ASCII_TYPE; }
  int
  IOVarASCII<string,0>::GetExtent(int i)
  { return 1; }
  void
  IOVarASCII<string,0>::Resize(int n)
  {
    cerr << "Cannot resize atomic variable.\n";
    abort();
  }
  bool
  IOVarASCII<string,0>::VarRead (string &val) 
  { 
    val = Value; 
    return true; 
  }

  int 
  IOVarASCII<bool,0>::GetRank()
  { return 0; }
  IODataType 
  IOVarASCII<bool,0>::GetType()
  { return BOOL_TYPE; }
  IOFileType 
  IOVarASCII<bool,0>::GetFileType()
  { return ASCII_TYPE; }
  int
  IOVarASCII<bool,0>::GetExtent(int i)
  { return 1; }
  void
  IOVarASCII<bool,0>::Resize(int n)
  {
    cerr << "Cannot resize atomic variable.\n";
    abort();
  }

  int 
  IOVarASCII<complex<double>,0>::GetRank()
  { return 0; }
  IODataType 
  IOVarASCII<complex<double>,0>::GetType()
  { return COMPLEX_TYPE; }
  IOFileType 
  IOVarASCII<complex<double>,0>::GetFileType()
  { return ASCII_TYPE; }
  int
  IOVarASCII<complex<double>,0>::GetExtent(int i)
  { return 1; }
  void
  IOVarASCII<complex<double>,0>::Resize(int n)
  {
    cerr << "Cannot resize atomic variable.\n";
    abort();
  }

}



