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
    
    




#ifndef OHMMS_TAG_MAKER_H
#define OHMMS_TAG_MAKER_H
/*!\class TagMaker
 * \brief Assign a unique tag whenver TagMaker::TagMaker() is called.
 */
class TagMaker
{
public:

  TagMaker()
  {
    MyTag = (++CurrentTag);
  }
  ~TagMaker() {}
  int operator()()const
  {
    return MyTag;
  }
private:
  int MyTag;
  static int CurrentTag;
};
#endif // OHMMS_TAG_MAKER_H

