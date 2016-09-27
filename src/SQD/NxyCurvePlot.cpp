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
    
    


#include <iostream>
#include <stdlib.h>
#include <qfile.h>
#include <qstringlist.h>
#include <qvaluevector.h>
#include <qptrvector.h>
#include <qwt_legend.h>
#include <qwt_plot_canvas.h>
#include "SQD/NxyCurvePlot.h"

std::map<long,QColor> NxyCurvePlot::ColorList;

NxyCurvePlot::NxyCurvePlot(QWidget *parent, const char* aname):
  QwtPlot(parent,aname)
{
  ColorList[0] = black;
  ColorList[1] = black;
  ColorList[2] = red;
  ColorList[3] = green;
  ColorList[4] = blue;
  ColorList[5] = cyan;
  ColorList[6] = magenta;
  setFrameStyle(QFrame::Box);
  setLineWidth(0);
  setCanvasLineWidth(2);
  setCanvasBackground(white);
  enableGridX(TRUE);
  enableGridY(TRUE);
  setGridMajPen(QPen(gray, 0, DotLine));
  //setAutoLegend(TRUE);
  //setLegendPos(Qwt::Right);
  initCurve();
  connect(this, SIGNAL(plotMousePressed(const QMouseEvent&)),
          this, SLOT(modifyCurves(const QMouseEvent&)));
}

QSize NxyCurvePlot::sizeHint() const
{
  return QSize(540,270);
}

void NxyCurvePlot::initCurve()
{
  removeCurves();
}

void NxyCurvePlot::clear()
{
  initCurve();
  replot();
}

/*!
 *\brief append a dataset to a QwtplotCurve curveId
 */
void NxyCurvePlot::appendCurveData(long curveId,
                                   double *x,
                                   double *y,
                                   int size)
{
  QwtPlotCurve *curve = NxyCurvePlot::curve(curveId);
  //curve->setPen(darkRed);
  curve->setPen(QPen(ColorList[curveId],2));
  curve->setStyle(QwtCurve::Spline);
  ///attach a curve to the specifies arrays without copyings
  //curve->setRawData(x,y,size);
  setCurveData(curveId,x,y,size);
  //d_legend->setPen(curveId,ColorList[curveId]);
}

void NxyCurvePlot::modifyCurves(const QMouseEvent &e)
{
}

