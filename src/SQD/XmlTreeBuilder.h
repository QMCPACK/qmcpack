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


class XmlTreeBuilder: public QMainWindow
{

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
