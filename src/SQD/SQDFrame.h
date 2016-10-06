//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/****************************************************************************
** Form interface generated from reading ui file 'sqd.ui'
**
** Created: Mon May 31 19:07:43 2004
**      by: The User Interface Compiler ($Id$)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#ifndef SQDFRAME_H
#define SQDFRAME_H

#include "SQD/HartreeFock.h"
#ifdef HAVE_QT
#include <qvariant.h>
#include <qwidget.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QLabel;
class QPushButton;
class NxyCurvePlot;
class QTextEdit;

/** A main application to solve Hartree-Fock equation in a spherical potential.
 */
class SQDFrame : public QWidget
{
  Q_OBJECT

public:

  SQDFrame( QWidget* parent = 0, const char* name = 0, WFlags fl = 0 );
  ~SQDFrame();

  /**
     @param fname the name of an input file
     @brief solve Hartree-Fock eqution with an input xml file fname
  */
  bool solve(const char* fname);

public slots:

  void slotSave();

protected:

  /**
    @brief Show the potential and wave functions
  */
  void showEvent ( QShowEvent * );

protected slots:

  virtual void languageChange();

private:
  QLabel* potLabel;
  QLabel* wfsLabel;
  NxyCurvePlot* potFrame;
  NxyCurvePlot* wfsFrame;
  QPushButton* saveButton;
  QPushButton* quitButton;
  QTextEdit *d_log;

  std::string elementType;
  std::string potType;
  std::string gridType;

  ohmmshf::RadialPotentialSet Pot;
  ohmmshf::SphericalOrbitalTraits::BasisSetType Psi;
  ohmmshf::HartreeFock *HFSolver;

  xmlDocPtr m_doc;

  //xmlpp::DomParser myParser;
  //xmlpp::Node *docRoot;

};
#else

/** A main application to solve Hartree-Fock equation in a spherical
    potential.
 */
class SQDFrame
{

public:

  SQDFrame();
  ~SQDFrame();

  /**
     @param fname the name of an input file
     @brief solve Hartree-Fock eqution with an input xml file fname
  */
  bool solve(const char* fname);

  void show();

  ohmmshf::RadialPotentialSet Pot;
  ohmmshf::SphericalOrbitalTraits::BasisSetType Psi;
  ohmmshf::HartreeFock *HFSolver;

  xmlDocPtr m_doc;
};
#endif
#endif // SQDFRAME_H
