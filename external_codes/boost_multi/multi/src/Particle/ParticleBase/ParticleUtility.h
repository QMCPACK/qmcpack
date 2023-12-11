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


#ifndef OHMMS_PARTICLEUTILITY_H
#define OHMMS_PARTICLEUTILITY_H

namespace qmcplusplus
{

////////////////////////////////////////////////////////////////
// Iterator is exposed. Parallel Implementation requires special care
////////////////////////////////////////////////////////////////
template<class PL, class PV>
void convert(const PL& lat, const PV& pin, PV& pout)
{
  if (pin.InUnit == pout.InUnit)
  {
    pout = pin;
    return;
  }
  if (pin.InUnit)
  {
    for (int i = 0; i < pin.size(); i++)
      pout[i] = lat.toCart(pin[i]);
    return;
  }
  else
  {
    for (int i = 0; i < pin.size(); i++)
      pout[i] = lat.toUnit(pin[i]);
    return;
  }
}

////////////////////////////////////////////////////////////////
// Iterator is exposed. Parallel Implementation requires special care
////////////////////////////////////////////////////////////////
template<class PL, class PV>
void convert2Cart(const PL& lat, PV& pin)
{
  if (pin.InUnit)
  {
    PV tmp(pin.size());
    tmp        = pin;
    pin.InUnit = false;
    for (int i = 0; i < pin.size(); i++)
      pin[i] = lat.toCart(pin[i]);
  }
}

template<class PL, class PV>
void convert2Unit(const PL& lat, PV& pin)
{
  if (!pin.InUnit)
  {
    PV tmp(pin.size());
    tmp        = pin;
    pin.InUnit = true;
    for (int i = 0; i < pin.size(); i++)
      pin[i] = lat.toUnit(pin[i]);
  }
}

////////////////////////////////////////////////////////////////
// Apply BC conditions to put the position type in lattice box [0,1)
////////////////////////////////////////////////////////////////
template<class PL, class PV>
void wrapAroundBox(const PL& lat, const PV& pin, PV& pout)
{
  if (pin.InUnit)
  {
    if (pout.InUnit)
    {
      for (int i = 0; i < pin.size(); i++)
      {
        pout[i] = lat.BConds.wrap(pin[i]); //unit -> unit
      }
    }
    else
    {
      for (int i = 0; i < pin.size(); i++)
        pout[i] = lat.toCart(lat.BConds.wrap(pin[i])); //unit -> cart
    }
  }
  else
  {
    if (pout.InUnit)
    {
      for (int i = 0; i < pin.size(); i++)
        pout[i] = lat.BConds.wrap(lat.toUnit(pin[i])); //cart -> unit
    }
    else
    {
      for (int i = 0; i < pin.size(); i++)
        pout[i] = lat.toCart(lat.BConds.wrap(lat.toUnit(pin[i]))); //cart -> cart
    }
  }
}

} // namespace qmcplusplus
#endif
