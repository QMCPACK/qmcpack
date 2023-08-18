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


#include <array>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <regex>
#include <string>

#include "SimpleParser.h"
#include <algorithm>

char* readLine(char* s, int max, std::istream& fp)
{
  char ch;
  int i = 0;
  while (fp.get(ch) && !(ch == '\n' || ch == ';'))
  {
    if (ch == '\\') // line continuation character
    {
      // check if backslash is followed by a newline
      fp.get(ch);
      if (ch == '\n')
      {
        // backslash followed by newline, do nothing
      }
      else
      {
        // backslash not followed by newline
        if (i < max - 1)
          s[i++] = '\\';
        if (i < max - 1)
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
    s[i] = '\0'; // add terminating NULL
  if (!(ch == '\n' || ch == ';'))
    return nullptr; // return NULL for end of file
  return s;
}


// NOTE that it only adds strings
unsigned parsewords(const char* inbuf, std::vector<std::string>& slist, const std::string& extra_tokens /* = "" */)
{
  std::string token = "=, \t\n\"";
  token.append(extra_tokens);

  std::string tmpstr(inbuf);
  slist.clear();
  int num      = 0;
  char* tokenp = strtok(tmpstr.data(), token.c_str());
  while (tokenp && tokenp[0] != '#')
  {
    num++;
    slist.push_back(std::string(tokenp));
    tokenp = strtok(nullptr, token.c_str());
  }
  return num;
}

unsigned parsewords(const char* inbuf, std::list<std::string>& slist)
{
  const char* token = "=, \t\n";
  std::string tmpstr(inbuf);
  slist.clear();
  unsigned num = 0;
  char* tokenp = strtok(tmpstr.data(), token);
  while (tokenp && tokenp[0] != '#')
  {
    num++;
    slist.push_back(std::string(tokenp));
    tokenp = strtok(nullptr, token);
  }
  return num;
}

int getwords(std::vector<std::string>& slist, std::istream& fp, std::string& aline)
{
  const int max = 1024;
  char s[max];
  if (readLine(s, max, fp))
  {
    aline.clear();
    aline.append(s);
    return parsewords(s, slist);
  }
  else
    return -1;
}

/* dummy argument present so function's type signature can be distinguished
   from previous function */

int getwords(std::vector<std::string>& slist,
             std::istream& fp,
             int dummy /* = 0*/,
             const std::string& extra_tokens /* ="" */)
{
  const int max = 1024;
  char s[max];
  if (readLine(s, max, fp))
    return parsewords(s, slist, extra_tokens);
  else
    return -1;
}

// Variation of getwords that splits merged numbers due to fortran fixed format
// Handles unambiguous cases with minus signs only "123-456" -> "123 -456"
int getwordsWithMergedNumbers(std::vector<std::string>& slist,
                              std::istream& fp,
                              int dummy /* = 0*/,
                              const std::string& extra_tokens /* ="" */)
{
  const int max = 1024;
  char s[max];
  if (readLine(s, max, fp))
  {
    std::regex dash("-");
    const std::string space_dash(" -");
    std::string merged(s);
    std::string unmerged = std::regex_replace(merged, dash, space_dash, std::regex_constants::format_default);
    return parsewords(unmerged.c_str(), slist, extra_tokens);
  }
  else
  {
    return -1;
  }
}

void readXmol(std::istream& fxmol, double* data, int numvar)
{
  std::vector<std::string> slist;
  getwords(slist, fxmol);
  unsigned natom = atoi(slist.front().c_str());
  getwords(slist, fxmol);
  int ii         = 0;
  for (int i = 0; i < natom; i++)
  {
    getwords(slist, fxmol);
    for (int ivar = 1; ivar <= numvar; ivar++)
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

int getwords(std::vector<std::string>& slist, std::istream& fpos, const char* field, const char* terminate)
{
  slist.clear();
  std::array<char, 128> end_key;
  if (std::snprintf(end_key.data(), end_key.size(), "</%s>", field) < 0)
    throw std::runtime_error("Error extract end_key from field.");

  std::vector<std::string> vlist;
  while (true)
  {
    if (getwords(vlist, fpos) == 0)
      continue;
    if (vlist[0] == terminate || vlist[0] == end_key.data())
      break;
    slist.insert(slist.end(), std::make_move_iterator(vlist.begin()), std::make_move_iterator(vlist.end()));
  };
  return slist.size();
}

////////////////////////////////////////////////////////
// simple parser to get around XML parser problem
////////////////////////////////////////////////////////
unsigned parseXwords(const char* inbuf, std::vector<std::string>& slist)
{
  slist.clear();

  const char* token = "=, <>\"\t\n";
  std::string tmpstr(inbuf);
  unsigned num = 0;
  char* tokenp = strtok(tmpstr.data(), token);
  while (tokenp && tokenp[0] != '#')
  {
    num++;
    slist.push_back(std::string(tokenp));
    tokenp = strtok(nullptr, token);
  }
  return num;
}

int getXwords(std::vector<std::string>& slist, std::istream& fp)
{
  const int max = 1024;
  char s[max];
  if (readLine(s, max, fp))
    return parseXwords(s, slist);
  else
    return -1;
}
