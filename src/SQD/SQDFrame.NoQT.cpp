//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/****************************************************************************
** Form implementation generated from reading ui file 'sqd.ui'
**
** Created: Mon May 31 19:07:56 2004
**      by: The User Interface Compiler ($Id$)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#include <set>
#include "SQD/SQDFrame.h"

SQDFrame::SQDFrame()
  : HFSolver(NULL)
{
}

/*
 *  Destroys the object and frees any allocated resources
 */
SQDFrame::~SQDFrame()
{
  // no need to delete child widgets, Qt does it all for us
  if(HFSolver)
    delete HFSolver;
}

