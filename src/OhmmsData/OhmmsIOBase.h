//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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

#ifndef OHMMS_OHMMSIOBASE_H
#define OHMMS_OHMMSIOBASE_H

class OhmmsIOBase {

public:

  OhmmsIOBase():stride(-1), OwnBuffer(false) {}

  virtual void open(const char*,ios_base::openmode iomode) = 0;
  virtual void close() = 0;
  virtual void flush() = 0;

  virtual bool read() = 0;
  virtual bool write(int iter) = 0;

protected:

  ios_base::openmode Mode;
  int stride;

private:
  ///no copy construction is available
  OhmmsIOBase(OhmmsIOBase&) { }
  
};

template<class T>
class StreamIO: public OhmmsIOBase {

  T& ref_;
  iostream* m_buffer;
  bool OwnBuffer;  

public:

  StreamIO(T& a, ios_base::openmode mode): ref_(T),
					   m_buffer(NULL), 
					   Mode(mode){}
  ~StreamIO() { 
    if(OwnBuffer && m_buffer) delete m_buffer;
  }

  void open(const char* fname, ios_mode::openmode mode) {
    if(m_buffer) {
      if(Mode != mode) m_buffer->close(); ///mode change for an existing buffer
    } else {
      m_buffer = new fstream;
      OwnBuffer = true;
    }
    Mode = mode;
    m_buffer->open(fname,mode);
  }

  inlin void setBuffer(iostream* ebuffer) {
    OwnBuffer = false;
    m_buffer = ebuffer;
  }
  inline void flush() {
    if(m_buffer) m_buffer->flush();
  }

  inline void close() {
    if(m_buffer) { m_buffer->close();}
  }

  inline bool read() {
    *m_buffer >> ref_;
    return true;
  }  

  inline bool write(int iter) {
    if(iter%stride == 0)  *m_buffer << ref_;
    return true;
  }  
};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
