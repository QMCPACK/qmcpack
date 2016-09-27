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
** Form implementation generated from reading ui file 'sqd.ui'
**
** Created: Mon May 31 19:07:56 2004
**      by: The User Interface Compiler ($Id$)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#include "PseudoGen/SQDFrame.h"

#ifdef HAVE_QT
#include <set>
#include <qvariant.h>
#include <qlabel.h>
#include <qframe.h>
#include <qpushbutton.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>
#include <qapplication.h>
#include <qsplitter.h>
#include <qvbox.h>
#include <qtextedit.h>
#include <qfile.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qlayout.h>
#include <qwt_plot_canvas.h>
#include "SQD/NxyCurvePlot.h"

/*
 *  Constructs a SQDFrame as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 */
SQDFrame::SQDFrame( QWidget* parent, const char* name, WFlags fl )
  : QWidget( parent, name, fl ), PPSolver(NULL)
{
  if ( !name )
    setName( "SQDFrame" );
  potLabel = new QLabel( this, "poteLabel" );
  potLabel->setGeometry( QRect( 240, 20, 260, 24 ) );
  potFrame = new NxyCurvePlot(this,"potPlot");
  potFrame->setGeometry( QRect( 20, 60, 540, 200 ) );
  potFrame->setFrameShape( QFrame::StyledPanel );
  potFrame->setFrameShadow( QFrame::Raised );
  wfsLabel = new QLabel( this, "wfsLabel" );
  wfsLabel->setGeometry( QRect( 160, 280, 400, 24 ) );
  wfsFrame = new NxyCurvePlot(this,"wfsPlot");
  wfsFrame->setGeometry( QRect( 20, 320, 540, 270 ) );
  d_log = new QTextEdit(this);
  d_log->setGeometry( QRect( 40, 600, 510, 200 ));
  d_log->setTextFormat(Qt::LogText);
  QFrame *hb = new QFrame(this,"buttonbox");
  hb->setGeometry(QRect(40,805,540,60));
  QHBoxLayout *layout = new QHBoxLayout(hb);
  quitButton = new QPushButton(hb, "quitButton" );
  layout->addItem(new QSpacerItem( 50, 20, QSizePolicy::Expanding, QSizePolicy::Minimum ));
  layout->addWidget(quitButton);
  layout->addItem(new QSpacerItem( 50, 20, QSizePolicy::Fixed, QSizePolicy::Minimum ));
  saveButton = new QPushButton( hb, "saveButton" );
  layout->addWidget(saveButton);
  layout->addItem(new QSpacerItem( 50, 20, QSizePolicy::Expanding, QSizePolicy::Minimum ));
  connect( quitButton, SIGNAL(clicked()), qApp, SLOT(quit()) );
  connect( saveButton, SIGNAL(clicked()), this, SLOT(slotSave()) );
  resize( QSize(600, 860).expandedTo(minimumSizeHint()) );
  languageChange();
  clearWState( WState_Polished );
}

/*
 *  Destroys the object and frees any allocated resources
 */
SQDFrame::~SQDFrame()
{
  // no need to delete child widgets, Qt does it all for us
  if(PPSolver)
    delete PPSolver;
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void SQDFrame::languageChange()
{
  setCaption( tr( "Effective-mass simulation of a spherical quantum dot" ) );
  setIconText( tr( "sqd" ) );
  potLabel->setText( tr( "<font size=\"+1\">Radial potential</font>" ) );
  wfsLabel->setText( tr( "<font size=\"+1\">Solutions of a radial Schrodinger Eq. <fon size=\"+2\"><i>u<sub>l</sub></i></font></font>" ) );
  quitButton->setText( tr( "Quit" ) );
  saveButton->setText( tr( "Save" ) );
  QToolTip::add(quitButton, "Leave the application");
  QToolTip::add(saveButton, "Save the results of the application");
}


/**
   @overload void SQDFrame::showEvent(QShowEvent *event)
   @brief Fill the frame with the results from Hartree-Fock calculations
   *
   *The summary is given in the log frame and potential and the solution of the orbitals
   \f$u_{nl} = r R_{nl}\f$ are plotted with respect to the radius.
*/
void SQDFrame::showEvent( QShowEvent *event)
{
  int max_rad_all = PPSolver->report();
  QString msg
  = QString("Potential type: %1\nGrid type: %2\nEffective mass: %3").arg(potType).arg(gridType).arg(Pot.getMass());
  d_log->append(msg);
  d_log->append("Eigen energies of a spherical quantum dot");
  for(int i=0; i<Psi.size(); i++)
  {
    msg
    = QString("\tOrbital n=%1.l=%2.m=%3.s=%4  %5").arg(Psi.N[i]).arg(Psi.L[i]).arg(Psi.M[i]).arg(Psi.S[i]).arg(PPSolver->getE(i));
    d_log->append(msg);
  }
  long y1 = potFrame->insertCurve("Potential Energy");
  potFrame->appendCurveData(y1,
                            Psi.m_grid->data(),
                            Pot.Vext->data(),
                            max_rad_all);
  potFrame->replot();
  int orbindex = 0;
  wfsFrame->setAxisTitle(QwtPlot::xBottom, "Radius in AU");
  for(int orb=0; orb<Psi.NumUniqueOrb; orb++)
  {
    int n = Psi.N[orbindex];
    int l = Psi.L[orbindex];
    QString cname= QString("n=%1.l=%2").arg(n,l);
    y1 = wfsFrame->insertCurve(cname);
    wfsFrame->appendCurveData(y1,Psi.m_grid->data(),Psi(orbindex).data(),max_rad_all);
    orbindex += Psi.IDcount[orb];
  }
  wfsFrame->replot();
//   double e = HFSolver->getE(0);
//   int orb = 0;
//   int ng = Psi.m_grid->size();
//   int n = Psi.N[orb];
//   int l = Psi.L[orb];
//   QString cname= QString("n=%1.l=%2").arg(n,l);
//   y1 = wfsFrame->insertCurve(cname);
//   wfsFrame->appendCurveData(y1,Psi.m_grid->data()+1,Psi(orb).data()+1,ng-1);
//   for(orb=1; orb<Psi.NumUniqueOrb; orb++){
//     int n_cur = Psi.N[orb];
//     int l_cur = Psi.L[orb];
//     //if(std::abs(e-HFSolver->getE(orb))>1e-6 || n_cur != n || l_cur != l ) {
//     if(n_cur != n || l_cur != l) {
//       n=n_cur; l=l_cur;
//       cname= QString("n=%1.l=%2").arg(n,l);
//       y1 = wfsFrame->insertCurve(cname);
//       wfsFrame->appendCurveData(y1,Psi.m_grid->data()+1,Psi(orb).data()+1,ng-1);
//       e = HFSolver->getE(orb);
//     }
//   }
//   wfsFrame->replot();
}

/**
 @brief save results
*/
void SQDFrame::slotSave()
{
  QString fname = QFileDialog::getSaveFileName(".",QString::null,this);
  bool save2file=false;
  if(!fname.isEmpty())
  {
    QFileInfo fi(fname);
    if(fi.exists())
    {
      QMessageBox mb("Save as",
                     "File exists. Do you really want to overwrite the file?",
                     QMessageBox::Information,
                     QMessageBox::Yes | QMessageBox::Default,
                     QMessageBox::No,
                     QMessageBox::Cancel | QMessageBox::Escape);
      mb.setButtonText(QMessageBox::Yes,"Save");
      mb.setButtonText(QMessageBox::No,"Don't Save");
      save2file = (mb.exec() == QMessageBox::Yes);
    }
    else
    {
      save2file = true;
    }
  }
  if(save2file)
  {
    QFile file(fname); // Write the text to a file
    if ( file.open( IO_WriteOnly ) )
    {
      QTextStream stream( &file );
      stream << d_log->text();
      d_log->setModified( FALSE );
    }
  }
}
#else
SQDFrame::SQDFrame() : PPSolver(NULL)
{
}

/*
 *  Destroys the object and frees any allocated resources
 */
SQDFrame::~SQDFrame()
{
  // no need to delete child widgets, Qt does it all for us
  if(PPSolver)
    delete PPSolver;
}


void SQDFrame::show()
{
  PPSolver->report();
}
#endif
