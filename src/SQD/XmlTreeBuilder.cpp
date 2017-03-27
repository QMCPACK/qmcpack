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
    
    
#ifdef HAVE_CONFIG_H
#include "ohmms-config.h"
#endif
#include "SQD/XmlTreeBuilder.h"
#include <qimage.h>
#include <qpixmap.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtextedit.h>
#include <qfile.h>
#include <qfiledialog.h>
#include <qstatusbar.h>
#include <qmessagebox.h>
#include <qprinter.h>
#include <qapplication.h>
#include <qaccel.h>
#include <qtextstream.h>
#include <qpainter.h>
#include <qpaintdevicemetrics.h>
#include <qwhatsthis.h>
#include <qsimplerichtext.h>
#include <qsplitter.h>
#include <qlistview.h>
#include <qlabel.h>
#include <qvbox.h>
#include "SQD/filesave.xpm"
#include "SQD/fileopen.xpm"
#include "SQD/fileprint.xpm"
#include "SQD/QtUtility.h"
#include "OhmmsData/libxmldefs.h"
/*! \brief Constructor
 * Register OhmmsElementBase objects corresponding to the major tags of input xml files.
 * OhmmsElementSet contains default objects \ref ohmmsdefaults and the constructor
 * adds OhmmsElementBase subclasses that depend on the default objects.
 *
 */
XmlTreeBuilder::XmlTreeBuilder()
  : QMainWindow( 0, "ohmms application main window", WDestructiveClose | WGroupLeader ),
    m_doc(NULL), m_context(NULL)
{
  QPixmap openIcon, saveIcon, printIcon;
  QToolBar * fileTools = new QToolBar( this, "file operations" );
  fileTools->setLabel( "File Operations" );
  openIcon = QPixmap( fileopen_xpm );
  QToolButton * fileOpen
  = new QToolButton( openIcon, "Open File", QString::null,
                     this, SLOT(choose()), fileTools, "open file" );
  saveIcon = QPixmap( filesave );
  QToolButton * fileSave
  = new QToolButton( saveIcon, "Save File", QString::null,
                     this, SLOT(save()), fileTools, "save file" );
  (void)QWhatsThis::whatsThisButton( fileTools );
  const char * fileOpenText = "<p><img source=\"fileopen\"> "
                              "Click this button to open a <em>new file</em>.<br>"
                              "You can also select the <b>Open</b> command "
                              "from the <b>File</b> menu.</p>";
  QWhatsThis::add( fileOpen, fileOpenText );
  QMimeSourceFactory::defaultFactory()->setPixmap( "fileopen", openIcon );
  const char * fileSaveText = "<p>Click this button to save the file you "
                              "are editing. You will be prompted for a file name.\n"
                              "You can also select the <b>Save</b> command "
                              "from the <b>File</b> menu.</p>";
  QWhatsThis::add( fileSave, fileSaveText );
  const char * filePrintText = "Click this button to print the file you "
                               "are editing.\n"
                               "You can also select the Print command "
                               "from the File menu.";
  QPopupMenu * file = new QPopupMenu( this );
  menuBar()->insertItem( "&File", file );
  int id;
  id = file->insertItem( openIcon, "&Open...",
                         this, SLOT(choose()), CTRL+Key_O );
  file->setWhatsThis( id, fileOpenText );
  id = file->insertItem( saveIcon, "&Save",
                         this, SLOT(save()), CTRL+Key_S );
  file->setWhatsThis( id, fileSaveText );
  id = file->insertItem( "Save &As...", this, SLOT(saveAs()) );
  file->setWhatsThis( id, fileSaveText );
  file->insertSeparator();
  file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_W );
  file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );
  menuBar()->insertSeparator();
  QPopupMenu * help = new QPopupMenu( this );
  menuBar()->insertItem( "&Help", help );
  help->insertItem( "&About", this, SLOT(about()), Key_F1 );
  help->insertSeparator();
  help->insertItem( "What's &This", this, SLOT(whatsThis()), SHIFT+Key_F1 );
  //d_view = new QSplitter( Qt::Vertical, this);
  //d_view = new QScrollView(this);
  d_view = new QVBox(this);
  setCentralWidget(d_view);
  nodeView = new QListView(d_view,"whatever");
  nodeView->addColumn("Node Name");
  nodeView->addColumn("Value");
  nodeView->addColumn("Unit");
  nodeView->setSorting(-1);
  nodeView->setRootIsDecorated(TRUE);
  resize( 450, 600 );
}

XmlTreeBuilder::~XmlTreeBuilder()
{
  if(m_doc)
    xmlFreeDoc(m_doc);
}

void XmlTreeBuilder::closeEvent( QCloseEvent* ce )
{
  switch( QMessageBox::information( this, "Ohmms Session",
                                    "Do you want to save the changes"
                                    " to the document?",
                                    "Yes", "No", "Cancel",
                                    0, 1 ) )
  {
  case 0:
    save();
    ce->accept();
    break;
  case 1:
    ce->accept();
    break;
  case 2:
  default: // just for sanity
    ce->ignore();
    break;
  }
}

void XmlTreeBuilder::resizeEvent(QResizeEvent *e)
{
  d_view->resize(e->size());
}

void XmlTreeBuilder::choose()
{
  QString fn = QFileDialog::getOpenFileName( QString::null, QString::null,
               this);
  if ( !fn.isEmpty() )
    load( fn );
  else
    statusBar()->message( "Loading aborted", 2000 );
}

void XmlTreeBuilder::load( const QString &fileName )
{
  parse(fileName.ascii());
  setCaption( fileName );
  statusBar()->message( "Loaded document " + fileName, 2000 );
}


void XmlTreeBuilder::save()
{
  if(m_doc)
    xmlSaveFormatFile(d_filename.ascii(),m_doc,1);
}

void XmlTreeBuilder::saveAs()
{
  QString fn = QFileDialog::getSaveFileName( QString::null, QString::null,
               this );
  if ( !fn.isEmpty() )
  {
    d_filename = fn;
    save();
  }
  else
  {
    statusBar()->message( "Saving aborted", 2000 );
  }
}

void XmlTreeBuilder::about()
{
  QMessageBox::about( this, "Ohmms Session",
                      "This example demonstrates simple use of "
                      "QMainWindow,\nQMenuBar and QToolBar.");
}

/*!\fn bool XmlTreeBuilder::parse(const char* fname) {
 *\param fname an input xml file to be parsed.
 *\return true if successful.
 */
bool XmlTreeBuilder::parse(const char* fname)
{
  d_filename = fname; //QString("%1.%2").arg(fname).arg(1);
  QFileInfo fi(fname);
  //free an existing document
  if(m_doc)
    xmlFreeDoc(m_doc);
  // build an XML tree from a the file;
  m_doc = xmlParseFile(fname);
  if (m_doc == NULL)
    return false;
  // Check the document is of the right kind
  xmlNodePtr cur = xmlDocGetRootElement(m_doc);
  if (cur == NULL)
  {
    fprintf(stderr,"empty document\n");
    xmlFreeDoc(m_doc);
    m_doc = NULL;
    return false;
  }
  QListViewItem *a=new QListViewItem(nodeView,fi.fileName());
  SimpleTreeBuilder(a,cur);
  nodeView->setOpen(a,true);
  nodeView->show();
  return true;
}


