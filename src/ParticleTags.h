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
    
    


#ifndef OHMMS_PARTICLETAGS_H
#define OHMMS_PARTICLETAGS_H

/**@defgroup XMLTags Reserved XML tags
 * @brief Reserved names for xml elements and attributes.
 */

#include <string>
/** Define the tags for the data types and attributes of particle classes
 *
 *The tags are grouped into
 *- the data type of an attribute (type name)
 *- the name of an attribute (object name)
 *- xml node and attribute names
 */
struct ParticleTags
{

  /** enumeraton for datatypes for particle attributes
   *
   *Each attribute of type T is stored in a container ParticleAttrib<T>.
   *ParticleBase::AttribTypeMap is used to map the std::string and the enumeration.
   */
  enum {PA_IndexType=0, /*!< index (int) type */
        PA_ScalarType,  /*!< scalar (float|double) type */
        PA_StringType,  /*!< std::string type */
        PA_PositionType,/*!< position (D-dim vector) type */
        PA_TensorType, /*!< tensor type */
        PA_UnkownType /*!< not known type */
       };

  /** @ingroup XMLTags
   *@{
   * @brief type names of particle attributes, e.g., posArray
   */
  ///the tag for uninitialized datatype and name of attributes
  static std::string null_tag;

  ///the datatype tag for the index-type attribute
  static std::string indextype_tag;

  ///the datatype tag for the scalar-type attribute
  static std::string scalartype_tag;

  ///the datatype tag for the string-type attribute
  static std::string stringtype_tag;

  ///the datatype tag for the pos-type attribute
  static std::string postype_tag;

  ///the datatype tag for the gradient-type attribute
  static std::string gradtype_tag;

  ///the datatype tag for the laplacian-type attribute
  static std::string laptype_tag;

  ///the datatype tag for the tensor-type attribute
  static std::string tensortype_tag;

  ///the datatype tag for the xmol-type attribute
  static std::string xmoltype_tag;
  /*@}*/

  /** @ingroup XMLTags
   *@{
   * @brief names of particle attributes, e.g., position
   */
  ///the name for position attribute
  static std::string position_tag;

  ///the name for id attribute
  static std::string id_tag;

  ///the name for ionid attribute
  static std::string ionid_tag;

  ///the name for trajectory attribute
  static std::string trajectory_tag;

  ///the name for force attribute
  static std::string force_tag;

  ///the name for velocity attribute
  static std::string velocity_tag;

  ///the name for energy attribute
  static std::string energy_tag;

  ///the name for sumbc (accumulative boundary crossings) attribute
  static std::string sumbc_tag;
  /*@}*/

  /** @ingroup XMLTags
   *@{
   * @brief xml element names
   */
  ///the name of particle node
  static std::string root_tag;
  ///the name of particle attribute node
  static std::string attrib_tag;
  ///the attribute name of name
  static std::string name_tag;
  ///the attribute name of datatype
  static std::string datatype_tag;
  ///the attribute name of condition
  static std::string condition_tag;
  ///the attribute name of size
  static std::string size_tag;
  ///the attribute name of format
  static std::string format_tag;
  ///the attribute name of role
  static std::string role_tag;
  /*@}*/
};
#endif
