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
#ifndef OHMMS_SIMPLEPARSER_H
#define OHMMS_SIMPLEPARSER_H

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>
#if (__GNUC__ == 2)
#include <strstream.h>
#else
#include <sstream>
#endif
using namespace std;

char* readLine(char *s, int max, istream &fp);

int getwords(vector<string>& slist, istream &fp);
int getwords(vector<string>& slist,istream& fpos, const char* field, const char* terminate);
int getwords(vector<string>& slist,istream& fpos, const char* terminate);
int getXwords(vector<string>& slist, istream &fp);
int getXwords(vector<string>& slist,istream& fpos, const char* terminate);

unsigned parsewords(char *inbuf, vector<string>& slist);
unsigned parsewords(char *inbuf, list<string>& slist);

void readXmol(istream& ,double*,int);

#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
