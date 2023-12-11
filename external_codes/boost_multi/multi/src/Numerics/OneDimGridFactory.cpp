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


#include "OneDimGridFactory.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{
std::unique_ptr<OneDimGridFactory::GridType> OneDimGridFactory::createGrid(xmlNodePtr cur)
{
  std::unique_ptr<GridType> agrid;
  RealType ri     = 1e-5;
  RealType rf     = 100.0;
  RealType ascale = -1.0e0;
  RealType astep  = 1.25e-2;
  IndexType npts  = 1001;
  std::string gridType("log");
  std::string gridID("invalid");
  OhmmsAttributeSet radAttrib;
  radAttrib.add(gridType, "type");
  radAttrib.add(npts, "npts");
  radAttrib.add(ri, "ri");
  radAttrib.add(rf, "rf");
  radAttrib.add(ascale, "ascale");
  radAttrib.add(astep, "astep");
  radAttrib.add(ascale, "scale");
  radAttrib.add(astep, "step");
  radAttrib.add(gridID, "id");
  radAttrib.add(gridID, "name");
  radAttrib.add(gridID, "ref");
  if (cur != NULL)
    radAttrib.put(cur);
  if (gridType == "log")
  {
    if (ascale > 0.0)
    {
      LOGMSG("Using log grid with default values: scale = " << ascale << " step = " << astep << " npts = " << npts)
      agrid = std::make_unique<LogGridZero<RealType>>();
      agrid->set(astep, ascale, npts);
    }
    else
    {
      LOGMSG("Using log grid with default values: ri = " << ri << " rf = " << rf << " npts = " << npts)
      if (ri < std::numeric_limits<RealType>::epsilon())
      {
        ri = std::numeric_limits<RealType>::epsilon();
        app_error() << "   LogGrid cannot accept r=0 for the initial point. Using ri=" << ri << std::endl;
      }
      agrid = std::make_unique<LogGrid<RealType>>();
      agrid->set(ri, rf, npts);
    }
  }
  else if (gridType == "linear")
  {
    LOGMSG("Using linear grid with default values: ri = " << ri << " rf = " << rf << " npts = " << npts)
    agrid = std::make_unique<LinearGrid<RealType>>();
    agrid->set(ri, rf, npts);
  }
  else
  {
    throw UniformCommunicateError("Unknown gridtype. Valid settings are \"log\" and \"linear\"\n");
  }
  return agrid;
}
} // namespace qmcplusplus
