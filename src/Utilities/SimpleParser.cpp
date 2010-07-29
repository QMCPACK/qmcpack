//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <cstdio>
#include <cstring>


#include "Utilities/SimpleParser.h"
#include <algorithm>
using namespace std;

char* readLine(char *s, int max, istream &fp){
  char ch;
  int i = 0;
  while ( fp.get(ch) && !( ch == '\n' || ch ==';') )
  {
    if ( ch == '\\' ) // line continuation character
    {
      // check if backslash is followed by a newline
      fp.get(ch);
      if ( ch == '\n' )
      {
        // backslash followed by newline, do nothing
      }
      else
      {
        // backslash not followed by newline
        if ( i < max - 1 )
          s[i++] = '\\';
        if ( i < max - 1 )
          s[i++] = ch; //
      }
    }
    else
    {
      if (i < max - 1)
        s[i++] = ch;
    }
  }
  if (max > 0) s[i] = '\0';  // add terminating NULL

  if ( !(ch == '\n' || ch == ';' ) )
    return NULL;             // return NULL for end of file
  return s;
}


// NOTE that it only adds strings
unsigned parsewords(const char *inbuf, vector<string>& slist){
  const char* token = "=, \t\n\"";

  char *tmpstr = new char[strlen(inbuf)+1];
  strcpy(tmpstr, inbuf);
  slist.erase(slist.begin(), slist.end());

  int num = 0;
  char* tokenp = strtok(tmpstr, token);
  while(tokenp && tokenp[0] != '#') {
    num++;
    slist.push_back(string(tokenp));
    tokenp = strtok(0,token);
  }

  delete [] tmpstr;
  return num;

}

unsigned parsewords(const char *inbuf, list<string>& slist){
  const char* token = "=, \t\n";

  char *tmpstr = new char[strlen(inbuf)+1];
  strcpy(tmpstr, inbuf);
  slist.erase(slist.begin(), slist.end());

  int num = 0;
  char* tokenp = strtok(tmpstr, token);
  while(tokenp && tokenp[0] != '#') {
    num++;
    slist.push_back(string(tokenp));
    tokenp = strtok(0,token);
  }

  delete [] tmpstr;
  return num;

}

int getwords(vector<string>& slist, istream &fp, string& aline) {

  const int max = 1024;
  char s[max];

  if(readLine(s,max,fp)) {
    aline.clear();
    aline.append(s);
    return parsewords(s,slist);
  } else
    return -1;

}


int getwords(vector<string>& slist, istream &fp) {

  const int max = 1024;
  char s[max];

  if(readLine(s,max,fp)) 
    return parsewords(s,slist);
  else
    return -1;

}



void readXmol(istream& fxmol,double* data,int numvar) {

  vector<string> slist;
  int argc = getwords(slist,fxmol);
  unsigned natom = atoi(slist.front().c_str());

  argc =  getwords(slist, fxmol);

  int ii=0;
  for(int i=0; i<natom; i++) {
    argc =  getwords(slist, fxmol);
    for(int ivar=1; ivar<=numvar; ivar++) {
      data[ii++] = atof(slist[ivar].c_str());
    }
  }
}


/* \fn 
int getwords(vector<string>& slist,istream& fpos, const char* field, const char* terminate) 
* \param slist, input strings between <field> </field>
* \param fpos   istream 
* \param field  <filed> data </field>
* \param terminate string to stop searching
*/

int 
getwords(vector<string>& slist,istream& fpos, const char* field, const char* terminate) {

  slist.erase(slist.begin(), slist.end());

  vector<string> vlist;
  //char start_key[128];
  char end_key[128];
  //sprintf(start_key,"<%s>",field);
  sprintf(end_key,"</%s>",field);

//   bool start = false;
//   do {
//     int nw =getwords(vlist,fpos);
//     if(nw == 0) continue;
//     if(vlist[0] == terminate) 
//       return slist.size();
//     if(vlist[0] == start_key) {
//       slist.insert(slist.end(), vlist.begin()+1, vlist.end());
//       start = true;
//     }
//   } while(!start);

  bool start = true;
  do {
    int nw =getwords(vlist,fpos);
    if(nw == 0) continue;
    if(vlist[0] == terminate)
      return slist.size();
    if(vlist[0] == end_key) {
      start = false;
    } else {
      slist.insert(slist.end(), vlist.begin(), vlist.end());
    }
  } while(start);

  return slist.size();
}

/////////////////////////////////////////////////////////////
// insert parsed strings of istream until terminate is encountered
/////////////////////////////////////////////////////////////
int 
getwords(vector<string>& slist,istream& fpos, const char* terminate) {

  vector<string> vlist;

  // first check if the input list already contains "terminate"
  vector<string>::iterator it = find(slist.begin(), slist.end(), terminate);
  if(it != slist.end()) 
    return slist.size();

  //slist.erase(slist.begin(), slist.end()); // remove input number
  bool start = true;
  do {
    int nw =getwords(vlist,fpos);
    if(nw == 0) continue;
    if(vlist[0] == terminate)
      return slist.size();
    else 
      slist.insert(slist.end(), vlist.begin(), vlist.end());
  } while(start);

  return slist.size();
}

////////////////////////////////////////////////////////
// simple parser to get around XML parser problem
////////////////////////////////////////////////////////
unsigned parseXwords(char *inbuf, vector<string>& slist){
  const char* token = "=, <>\"\t\n";

  char *tmpstr = new char[strlen(inbuf)+1];
  strcpy(tmpstr, inbuf);
  slist.erase(slist.begin(), slist.end());

  int num = 0;
  char* tokenp = strtok(tmpstr, token);
  while(tokenp && tokenp[0] != '#') {
    num++;
    slist.push_back(string(tokenp));
    tokenp = strtok(0,token);
  }

  delete [] tmpstr;
  return num;

}

int getXwords(vector<string>& slist, istream &fp) {

  const int max = 1024;
  char s[max];

  if(readLine(s,max,fp)) 
    return parseXwords(s,slist);
  else
    return -1;

}


/////////////////////////////////////////////////////////////
// insert parsed strings of istream until terminate is encountered
/////////////////////////////////////////////////////////////
int 
getXwords(vector<string>& slist,istream& fpos, const char* terminate) {

  vector<string> vlist;

  // first check if the input list already contains "terminate"
  vector<string>::iterator it = find(slist.begin(), slist.end(), terminate);
  if(it != slist.end()) 
    return slist.size();

  //slist.erase(slist.begin(), slist.end()); // remove input number
  bool start = true;
  do {
    int nw =getXwords(vlist,fpos);
    if(nw == 0) continue;
    if(vlist[0] == terminate)
      return slist.size();
    else 
      slist.insert(slist.end(), vlist.begin(), vlist.end());
  } while(start);

  return slist.size();
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
