///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/fqmc/engine/var_dependencies.cpp
///
/// \brief   implementation file for the variable dependencies class
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include<vector>
#include<map>
#include<utility>
#include<algorithm>
#include<complex>

#include<formic/utils/exception.h>
#include<formic/utils/archive.h>
#include<formic/utils/lmyengine/var_dependencies.h>

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Constructs the object from a list of doubles, which should contain triplets of numbers
///        with the first number being the index of the dependent variable, the second number
///        being the scalar factor relating that variable to the independend variable, and the
///        third number being the index of the independend variable.
///
///        That is to say, for the triplet of numbers d0, d1, d2, we have
///
///          variable[d0] = d1 * variable[d2]
///
/// \param[in]      ntot       the total number of variables
/// \param[in]      input_vec  the vector of numbers telling how variables depend on each other
///
///////////////////////////////////////////////////////////////////////////////////////////////////
formic::VarDeps::VarDeps(const int ntot, const std::vector<double> & input_vec) {

  // ensure that the total number of variables is positive
  if ( ntot < 0 )
    throw formic::Exception("total number of variables (%i) must be positive in formic::VarDeps::VarDeps") % ntot;

  // ensure that the vector length is a multiple of three and is not too long
  if ( input_vec.size() % 3 != 0 )
    throw formic::Exception("length of input_vec (%i) was not a multiple of three in formic::VarDeps::VarDeps") % input_vec.size();
  if ( input_vec.size() / 3 >= ntot )
    throw formic::Exception("number of dependencies in input_vec (%i) is equal to or greater than the number of variables (%i) in formic::VarDeps::VarDeps") % (input_vec.size()/3) % ntot;

  // get the number of independend variables (easy since each dependency in the input vec reduces it by 1)
  const int nind = ntot - input_vec.size() / 3;

  // size the data members appropriately
  m_deps_per_ind.assign(nind, 0);
  m_indices.assign(ntot, -1);
  m_scalars.assign(ntot, 0.0);

  // process the input data to get a map from dependent to independent variables and a map from each independent variable to those that depend on it and the corresponding scalars
  std::map<int, int> dep_to_ind_map;
  std::map<int, std::vector<std::pair<double,int> > > ind_to_dep_map;
  for (int i = 0; i < input_vec.size(); i += 3) {

    // get and check the dependent variable index
    const int dep = formic::round(input_vec.at(i+0));
    if ( dep < 0 || dep >= ntot )
      throw formic::Exception("dependent variable index %i is not in the allowable range [0,%i) in formic::VarDeps::VarDeps") % dep % ntot;
    if ( dep_to_ind_map.count(dep) > 0 )
      throw formic::Exception("variable %i was listed as dependend multiple times in formic::VarDeps::VarDeps") % dep;
    if ( ind_to_dep_map.count(dep) > 0 )
      throw formic::Exception("variable %i was listed as both dependend and independend in formic::VarDeps::VarDeps") % dep;

    // get and check the scalar multiplier for the dependency
    double scl = input_vec.at(i+1);
    //if ( scl == 0.0 )
    //  throw formic::Exception("the scalar by which one variable depends on another cannot be zero in formic::VarDeps::VarDeps");

    // get and check the independent variable index
    int ind = formic::round(input_vec.at(i+2));
    if ( ind < 0 || ind >= ntot )
      throw formic::Exception("independent variable index %i is not in the allowable range [0,%i) in formic::VarDeps::VarDeps") % ind % ntot;

    // if the "indepenedent" variable is actually dependent on another variable, update the independent index and the scalar appropriately
    if ( dep_to_ind_map.count(ind) > 0 ) {
      const int actual_ind = dep_to_ind_map[ind];
      if ( ind_to_dep_map.count(actual_ind) < 1 )
        throw formic::Exception("expected actual_ind (%i) to be in the independent-to-dependent map in formic::VarDeps::VarDeps") % actual_ind;
      bool found = false;
      for (int j = 0; j < ind_to_dep_map[actual_ind].size(); j++) {
        if ( ind_to_dep_map[actual_ind].at(j).second == ind ) {
          if ( found )
            throw formic::Exception("expected to find ind (%i) only once in the independent-to-dependent map for actual_ind = %i in formic::VarDeps::VarDeps") % ind % actual_ind;
          scl *= ind_to_dep_map[actual_ind].at(j).first;
          found = true;
        }
      }
      if ( !found )
        throw formic::Exception("failed to find ind (%i) in the independent-to-dependent map for actual_ind = %i in formic::VarDeps::VarDeps") % ind % actual_ind;
      ind = actual_ind;
    }

    // add to the dependent-to-independent mapping
    dep_to_ind_map[dep] = ind;

    // add to the independent-to-dependent mapping
    ind_to_dep_map[ind].push_back(std::pair<double,int>(scl, dep));

    // if this was a newly found independent variable, make its first dependency itself
    if ( ind_to_dep_map[ind].size() == 1 )
      ind_to_dep_map[ind].insert(ind_to_dep_map[ind].begin(), std::pair<double,int>(1.0, ind));

  }

  // make dependent on themselves any variables that did not appear in the mapping
  for (int i = 0; i < ntot; i++)
    if ( dep_to_ind_map.count(i) == 0 && ind_to_dep_map.count(i) == 0 )
      ind_to_dep_map[i].push_back(std::pair<double,int>(1.0, i));

  // check that we found all the independent variables
  if ( ind_to_dep_map.size() != nind )
    throw formic::Exception("size of ind_to_dep_map (%i) did not equal nind (%i) in formic::VarDeps::VarDeps") % ind_to_dep_map.size() % nind;

  // populate the vectors that will be used to convert back and forth
  std::map<int, std::vector<std::pair<double,int> > >::const_iterator it = ind_to_dep_map.begin();
  for (int i = 0, j = 0; i < nind; i++, it++) {
    m_deps_per_ind.at(i) = (it->second).size();
    for (std::vector<std::pair<double,int> >::const_iterator jt = (it->second).begin(); jt != (it->second).end(); j++, jt++) {
      if ( j >= ntot )
        throw formic::Exception("unexpectedly reached the end of m_indices and m_scalars in formic::VarDeps::VarDeps");
      m_indices.at(j) = jt->second;
      m_scalars.at(j) = jt->first;
      //std::cout << boost::format("var %4i     =     %14.6e    *    var %4i") % jt->second % jt->first % it->first << std::endl;
    }
  }

  // ensure we covered all the indices
  if ( std::find(m_indices.begin(), m_indices.end(), -1) != m_indices.end() )
    throw formic::Exception("not all variables were accounted for in formic::VarDeps::VarDeps");

  //// ensure all the indices are distinct
  //for (int i = 0; i < m_indices.size(); i++) 
  //  std::cout << m_indices.at(i) << std::endl;
  //std::set<int> chris_check(m_indices.begin(), m_indices.end());
  //if ( chris_check.size() != m_indices.size() )
  //  throw formic::Exception("repeated index in m_dices in formic::VarDeps::VarDeps");
  //if ( std::set<int>(m_indices.begin(), m_indices.end()).size() != m_indices.size() )
  //  throw formic::Exception("repeated index in m_indices in formic::VarDeps::VarDeps");

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Constructs the object from two separate VarDeps objects.  The indices within the
///        second object are all incremented by the number of variables in the first object
///        to produce a unified index set in which the index sets from the two objects are
///        entirely distinct.
///
/// \param[in]      vd0        the first VarDeps object to combine
/// \param[in]      vd1        the sedond VarDeps object to combine
///
///////////////////////////////////////////////////////////////////////////////////////////////////
formic::VarDeps::VarDeps(const formic::VarDeps & vd0, const formic::VarDeps & vd1) {

  // get the new total number of variables
  const int ntot = vd0.n_tot() + vd1.n_tot();

  // get the new number of independent variables
  const int nind = vd0.n_ind() + vd1.n_ind();

  // size the data members appropriately
  m_deps_per_ind.assign(nind, 0);
  m_indices.assign(ntot, -1);
  m_scalars.assign(ntot, 0.0);

  // add the data for the first object
  int j = 0;
  for (int i = 0; i < vd0.n_ind(); i++) {
    m_deps_per_ind.at(i) = vd0.m_deps_per_ind.at(i);
    for (int k = 0; k < m_deps_per_ind.at(i); k++, j++) {
      m_indices.at(j) = vd0.m_indices.at(j);
      m_scalars.at(j) = vd0.m_scalars.at(j);
    }
  }

  // check that the correct number of indices and scalars were added
  if ( j != vd0.n_tot() )
    throw formic::Exception("added %i indices and scalars from vd0 instead of the expected %i in formic::VarDeps::VarDeps") % j % vd0.n_tot();

  // add the data for the second object
  for (int i = 0; i < vd1.n_ind(); i++) {
    m_deps_per_ind.at(vd0.n_ind()+i) = vd1.m_deps_per_ind.at(i);
    for (int k = 0; k < vd1.m_deps_per_ind.at(i); k++, j++) {
      m_indices.at(j) = vd1.m_indices.at(j-vd0.n_tot()) + vd0.n_tot();
      m_scalars.at(j) = vd1.m_scalars.at(j-vd0.n_tot());
    }
  }

  // check that the correct number of indices and scalars were added
  if ( j != ntot )
    throw formic::Exception("added %i indices and scalars from vd0 and vd1 instead of the expected %i in formic::VarDeps::VarDeps") % j % ntot;

  //for (int i = 0; i < m_deps_per_ind.size(); i++) 
  //  std::cout << boost::format("%d ") % m_deps_per_ind[i];
  //std::cout << std::endl;
  
  //for (int i = 0; i < m_indices.size(); i++) 
  //  std::cout << boost::format("%d ") % m_indices[i];
  //std::cout << std::endl;

  //for ( int i = 0; i < m_scalars.size(); i++) 
  //  std::cout << boost::format("%10.8e ") % m_scalars[i];
  //std::cout << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Replaces this variable dependency scheme with a combination of this one and
///         another in which this one comes first.
///
/// \param[in]      other      the other dependency scheme to combine with
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::VarDeps::append(const formic::VarDeps & other) {
  *this = formic::VarDeps(*this, other);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Fills in the values of all the variables given the values of the independent variables.
///
/// \param[in]      ind_vars   a length this->n_ind() array holding the ind var values
/// \param[out]     all_vars   a length this->n_tot() array to be filled with the values for all
///                            the variables
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void formic::VarDeps::expand_ind_to_all(const S * const ind_vars, S * const all_vars) const {
  for (int i = 0, j = 0; i < this->n_ind(); i++)
    for (int k = j + m_deps_per_ind[i]; j < k; j++)
      all_vars[m_indices[j]] = m_scalars[j] * ind_vars[i];
}

template void formic::VarDeps::expand_ind_to_all(const double * const ind_vars, double * const all_vars) const;
template void formic::VarDeps::expand_ind_to_all(const std::complex<double> * const ind_vars, std::complex<double> * const all_vars) const;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Given the partial derivatives of some function w.r.t. each of all the variables,
///        computes the total derivatives w.r.t. each of the independent variables
///
/// \param[in]      all_ders   a length this->n_tot() array holding the partial derivatives
/// \param[out]     ind_ders   a length this->n_ind() array to be filled with the total derivs
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void formic::VarDeps::compute_ind_derivs(const S * const all_ders, S * const ind_ders) const {
  for (int i = 0, j = 0; i < this->n_ind(); i++) {
    ind_ders[i] = formic::zero(S());
    for (int k = j + m_deps_per_ind[i]; j < k; j++)
      ind_ders[i] += m_scalars[j] * all_ders[m_indices[j]];
  }
}

template void formic::VarDeps::compute_ind_derivs(const double * const all_ders, double * const ind_ders) const;
template void formic::VarDeps::compute_ind_derivs(const std::complex<double> * const all_ders, std::complex<double> * const ind_ders) const;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Given the values of all the variables, extracts the values of the independent variables
///
/// \param[in]   all_vals   a length this->n_tot() array holding all the variables
/// \param[out]  ind_vals   a length this->n_ind() array to be filled with the independent vars
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void formic::VarDeps::extract_ind_vals(const S * const all_vals, S * const ind_vals) const { 
  for (int i = 0, j = 0; i < this->n_ind(); i++) {
    ind_vals[i] = all_vals[m_indices[j]];
    j += m_deps_per_ind[i];
  }
}

template void formic::VarDeps::extract_ind_vals(const double * const, double * const) const;
template void formic::VarDeps::extract_ind_vals(const std::complex<double> * const, std::complex<double> * const) const;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Sets the dependent variables equal to their multipliers times their independent variable
///
/// \param[in,out]  vars       a length this->n_tot() array holding the variables
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void formic::VarDeps::enforce_dependencies(S * const vars) const {
  for (int i = 0, j = 0; i < this->n_ind(); i++) {
    const S ind_val = vars[m_indices[j++]];
    for (int k = j + m_deps_per_ind[i] - 1; j < k; j++)
      vars[m_indices[j]] = m_scalars[j] * ind_val;
  }
}

template void formic::VarDeps::enforce_dependencies(double * const vars) const;
template void formic::VarDeps::enforce_dependencies(std::complex<double> * const vars) const;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Adopts the sign structure of the provided data for the dependency scalars
///
/// \param[in]      vars       a length this->n_tot() array holding the data in question
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class S> void formic::VarDeps::adopt_signs_from(const S * const vars) {
  for (int i = 0, j = 0; i < this->n_ind(); i++) {
    const S ind_val = vars[m_indices[j++]];
    for (int k = j + m_deps_per_ind[i] - 1; j < k; j++) {
      if ( ind_val != formic::zero(S()) ) { // only adopt a sign if the independent variable is nonzero
        const S dep_val = vars[m_indices[j]];
        const S s = ( dep_val / ind_val ) / ( dep_val == formic::zero(S()) ? 1.0 : std::abs( dep_val / ind_val ) );
        m_scalars[j] = formic::real(s) * std::abs(m_scalars[j]);
      }
    }
  }
}

template void formic::VarDeps::adopt_signs_from(const double * const vars);
template void formic::VarDeps::adopt_signs_from(const std::complex<double> * const vars);

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Writes to the supplied archive the command used to create this object's dependencies
///        in the named data object
///
/// \param[in,out]  arch       the archive to write the command to
/// \param[in]      data_name  the name of the data object
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::VarDeps::to_command(formic::Archive & arch, const std::string & data_name) const {

  // write keyword
  arch << "SetDataDependencies" << std::endl;

  // write data object name
  arch << data_name << std::endl;

  // write instructions not to absorb sign structures
  arch << false << std::endl;

  // prepare input vector
  std::vector<double> input_vec;
  for (int i = 0, j = 0; i < this->n_ind(); i++) {
    if ( m_deps_per_ind[i] > 1 ) {
      const int ind = m_indices[j++];
      for (int k = j + m_deps_per_ind[i] - 1; j < k; j++) {
        input_vec.push_back(double(m_indices[j]));
        input_vec.push_back(double(m_scalars[j]));
        input_vec.push_back(double(ind));
      }
    } else {
      j++;
    }
  }

  // write input vector
  arch << input_vec << std::endl;

}
