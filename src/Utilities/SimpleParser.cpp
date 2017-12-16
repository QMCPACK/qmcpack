//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include <cstdio>
#include <cstring>


#include "Utilities/SimpleParser.h"
#include <algorithm>

char* readLine(char *s, int max, std::istream &fp)
{
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
  if (max > 0)
    s[i] = '\0';  // add terminating NULL
  if ( !(ch == '\n' || ch == ';' ) )
    return NULL;             // return NULL for end of file
  return s;
}


// NOTE that it only adds strings
unsigned parsewords(const char *inbuf, std::vector<std::string>& slist, const std::string &extra_tokens /* = "" */)
{
  std::string token = "=, \t\n\"";
  token.append(extra_tokens);

  char *tmpstr = new char[strlen(inbuf)+1];
  strcpy(tmpstr, inbuf);
  slist.erase(slist.begin(), slist.end());
  int num = 0;
  char* tokenp = strtok(tmpstr,token.c_str());
  while(tokenp && tokenp[0] != '#')
  {
    num++;
    slist.push_back( std::string(tokenp));
    tokenp = strtok(0,token.c_str());
  }
  delete [] tmpstr;
  return num;
}

unsigned parsewords(const char *inbuf, std::list<std::string>& slist)
{
  const char* token = "=, \t\n";
  char *tmpstr = new char[strlen(inbuf)+1];
  strcpy(tmpstr, inbuf);
  slist.erase(slist.begin(), slist.end());
  int num = 0;
  char* tokenp = strtok(tmpstr, token);
  while(tokenp && tokenp[0] != '#')
  {
    num++;
    slist.push_back( std::string(tokenp));
    tokenp = strtok(0,token);
  }
  delete [] tmpstr;
  return num;
}

int getwords(std::vector<std::string>& slist, std::istream &fp, std::string& aline)
{
  const int max = 1024;
  char s[max];
  if(readLine(s,max,fp))
  {
    aline.clear();
    aline.append(s);
    return parsewords(s,slist);
  }
  else
    return -1;
}

/* dummy argument present so function's type signature can be distinguished
   from previous function */

int getwords(std::vector<std::string>& slist, std::istream &fp, int dummy /* = 0*/,
             const std::string &extra_tokens /* ="" */)
{
  const int max = 1024;
  char s[max];
  if(readLine(s,max,fp))
    return parsewords(s,slist,extra_tokens);
  else
    return -1;
}



void readXmol( std::istream& fxmol,double* data,int numvar)
{
  std::vector<std::string> slist;
  int argc = getwords(slist,fxmol);
  unsigned natom = atoi(slist.front().c_str());
  argc =  getwords(slist, fxmol);
  int ii=0;
  for(int i=0; i<natom; i++)
  {
    argc =  getwords(slist, fxmol);
    for(int ivar=1; ivar<=numvar; ivar++)
    {
      data[ii++] = atof(slist[ivar].c_str());
    }
  }
}


/* \fn
int getwords(std::vector<std::string>& slist,std::istream& fpos, const char* field, const char* terminate)
* \param slist, input strings between <field> </field>
* \param fpos   std::istream
* \param field  <filed> data </field>
* \param terminate std::string to stop searching
*/

int
getwords(std::vector<std::string>& slist,std::istream& fpos, const char* field, const char* terminate)
{
  slist.erase(slist.begin(), slist.end());
  std::vector<std::string> vlist;
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
  do
  {
    int nw =getwords(vlist,fpos);
    if(nw == 0)
      continue;
    if(vlist[0] == terminate)
      return slist.size();
    if(vlist[0] == end_key)
    {
      start = false;
    }
    else
    {
      slist.insert(slist.end(), vlist.begin(), vlist.end());
    }
  }
  while(start);
  return slist.size();
}

/////////////////////////////////////////////////////////////
// insert parsed strings of std::istream until terminate is encountered
/////////////////////////////////////////////////////////////
int
getwords(std::vector<std::string>& slist,std::istream& fpos, const char* terminate)
{
  std::vector<std::string> vlist;
  // first check if the input list already contains "terminate"
  std::vector<std::string>::iterator it = find(slist.begin(), slist.end(), terminate);
  if(it != slist.end())
    return slist.size();
  //slist.erase(slist.begin(), slist.end()); // remove input number
  bool start = true;
  do
  {
    int nw =getwords(vlist,fpos);
    if(nw == 0)
      continue;
    if(vlist[0] == terminate)
      return slist.size();
    else
      slist.insert(slist.end(), vlist.begin(), vlist.end());
  }
  while(start);
  return slist.size();
}

////////////////////////////////////////////////////////
// simple parser to get around XML parser problem
////////////////////////////////////////////////////////
unsigned parseXwords(char *inbuf, std::vector<std::string>& slist)
{
  const char* token = "=, <>\"\t\n";
  char *tmpstr = new char[strlen(inbuf)+1];
  strcpy(tmpstr, inbuf);
  slist.erase(slist.begin(), slist.end());
  int num = 0;
  char* tokenp = strtok(tmpstr, token);
  while(tokenp && tokenp[0] != '#')
  {
    num++;
    slist.push_back( std::string(tokenp));
    tokenp = strtok(0,token);
  }
  delete [] tmpstr;
  return num;
}

int getXwords(std::vector<std::string>& slist, std::istream &fp)
{
  const int max = 1024;
  char s[max];
  if(readLine(s,max,fp))
    return parseXwords(s,slist);
  else
    return -1;
}


/////////////////////////////////////////////////////////////
// insert parsed strings of std::istream until terminate is encountered
/////////////////////////////////////////////////////////////
int
getXwords(std::vector<std::string>& slist,std::istream& fpos, const char* terminate)
{
  std::vector<std::string> vlist;
  // first check if the input list already contains "terminate"
  std::vector<std::string>::iterator it = find(slist.begin(), slist.end(), terminate);
  if(it != slist.end())
    return slist.size();
  //slist.erase(slist.begin(), slist.end()); // remove input number
  bool start = true;
  do
  {
    int nw =getXwords(vlist,fpos);
    if(nw == 0)
      continue;
    if(vlist[0] == terminate)
      return slist.size();
    else
      slist.insert(slist.end(), vlist.begin(), vlist.end());
  }
  while(start);
  return slist.size();
}

