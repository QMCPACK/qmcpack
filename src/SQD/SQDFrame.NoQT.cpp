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
  if(HFSolver) delete HFSolver;
}

