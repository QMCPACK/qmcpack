//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#ifndef OHMMS_XMLTREEBUILDER_H
#define OHMMS_XMLTREEBUILDER_H

/*!\class OhmmsQt
 * \brief A generic application. 
 */
#include <map>
#include <string>
#include <iostream>
#include <qmainwindow.h>
#include "OhmmsData/libxmldefs.h"

class QTextEdit;
class QListView;
class QScrollView;


class XmlTreeBuilder: public QMainWindow {
    
  Q_OBJECT
  
public:

  ///constructor
  XmlTreeBuilder();
  
  ///destructor
  ~XmlTreeBuilder();
  
  ///parse a file
  bool parse(const char* fname);
  
protected:
  
  void closeEvent( QCloseEvent* );

protected slots:

void resizeEvent(QResizeEvent *);
  
private slots:

  void load( const QString &fileName );
  void save();
  void choose();
  void saveAs();
  void about();

private:

    xmlDocPtr m_doc;
    xmlXPathContextPtr m_context;

    QTextEdit *e;
    QString d_filename;
    QListView* nodeView;
    QWidget* d_view;
  };

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
