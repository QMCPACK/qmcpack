///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/fqmc/engine/var_dependencies.h
///
/// \brief   header file for the variable dependencies class
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FORMIC_VAR_DEPENDENCIES_HEADER
#define FORMIC_VAR_DEPENDENCIES_HEADER

#include<vector>
#include<string>

namespace formic {

  class Archive;

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   Class to keep track of which variables are independent and which variables depend
  ///          on others.
  ///
  /// This class keep track of simple variable dependencies of the form
  /// variable_2 = scalar * variable_1.  It is intended to be used in two ways.
  ///
  /// First, it can take a vector of the independent variables and fill in a vector of all the
  /// variables using the dependencies.
  ///
  /// Second, it can take a vector of derivatives of a function w.r.t. all the variables, in which
  /// the derivatives were taken assuming all variables were independent, and compute from it
  /// the vector of derivatives w.r.t. only the independent variables.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  class VarDeps {

    private:

      /// \brief  for each independent variable, the number of variables that depend on it
      std::vector<int> m_deps_per_ind;

      /// \brief  lists the indices that tell which variables depend on each independent variable
      std::vector<int> m_indices;

      /// \brief  the scalar multiples telling how variables depend on each other, e.g. dep_var_3 = scalar * ind_var_1
      std::vector<double> m_scalars;

    public:

      /// \brief  Default constructor
      VarDeps() {}

      VarDeps(const int ntot, const std::vector<double> & input_vec);

      VarDeps(const formic::VarDeps & vd0, const formic::VarDeps & vd1);

      /// \brief  Returns the number of total variables
      int n_tot() const { return m_indices.size(); }

      /// \brief  Returns the number of independent variables
      int n_ind() const { return m_deps_per_ind.size(); }

      void append(const formic::VarDeps & other);

      template<class S> void expand_ind_to_all(const S * const ind_vars, S * const all_vars) const;

      template<class S> void compute_ind_derivs(const S * const all_ders, S * const ind_ders) const;

      template<class S> void extract_ind_vals(const S * const all_vals, S * const ind_vals) const; 

      template<class S> void enforce_dependencies(S * const vars) const;

      template<class S> void adopt_signs_from(const S * const vars);

      void to_command(formic::Archive & arch, const std::string & data_name) const;

  };

}

#endif
