/////////////////////////////////////////////////////////////////////////////////////////////////
// \file formic/utils/lmyengine/UpdaterBase.h
//
//
// \brief header file for harmonic davidson method updater
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ENGINE_UPDATER_BASE_HEADER
#define ENGINE_UPDATER_BASE_HEADER

#include<vector>

#include<boost/shared_ptr.hpp>

#include"formic/utils/zero_one.h"
#include"formic/utils/mpi_interface.h"
#include"formic/utils/lmyengine/matrix_builder.h"
#include"formic/utils/lmyengine/eigen_solver.h"
#include"formic/utils/lmyengine/davidson_solver.h"
#include"formic/utils/lmyengine/spam_solver.h"
#include"formic/utils/lmyengine/energy_target_accu.h"
#include"formic/utils/numeric.h"
#include"formic/utils/matrix.h"
#include"formic/utils/lmyengine/var_dependencies.h"

namespace cqmc {

  namespace engine {

    ////////////////////////////////////////////////////////////////////////////////////////////
    // \brief A class to perform harmonic linear method update of the current wave function 
    //
    //
    ////////////////////////////////////////////////////////////////////////////////////////////
    template<typename S> class UpdaterBase {

      private:
        
      /// \brief harmonic Davidson shift(omega)
      double _omega;

      /// \brief |value/guiding|^2 list
      const std::vector<double> & _vgs;

      /// \brief weight list
      const std::vector<double> & _weight;

      /// \brief vector to store linear method shift scale
      const std::vector<double> & _shift_scale;

      /// \brief a flag to tell whether we are doing ground state calculation
      bool _ground_state;

      public:
        
      // constructor
      UpdaterBase(const std::vector<double> & vgs,
                  const std::vector<double> & weight,
                  const std::vector<double> & shift_scale,
                  const double omega,
                  const bool ground_state)
      :_vgs(vgs),
      _weight(weight),
      _shift_scale(shift_scale),
      _omega(omega),
      _ground_state(ground_state),
      {}

      // default destructor
      ~UpdaterBase() {}

      virtual void perform_update(const formic::VarDeps * dep_ptr, 
                                  std::vector<S> & vf_var, 
                                  std::vector<bool> & good_solve, 
                                  std::vector<int> & shift_solved, 
                                  std::ostream & output) = 0;


   }; // end of class definition
 } // end of engine
} // end of cqmc

#endif
