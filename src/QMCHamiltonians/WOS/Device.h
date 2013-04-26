#ifndef OHMMS_QMC_DEVICE_H
#define OHMMS_QMC_DEVICE_H
#include <vector>
#include "Numerics/Spline3D/Config.h"
#include "QMCHamiltonians/WOS/Domain.h"
#include "QMCHamiltonians/WOS/WOSParticles.h"

using namespace qmcplusplus;

class Device
{

public:

  /// width of skin region
  double m_skin;

  /// dielectric of the device
  double m_eps;

  /// factor for continuous charge density
  double q_fac;

  /// continuous charge density
  double rho0;

  /// boundaries of the device : capacitor
  posvec_t r_min;
  posvec_t r_max;

  /// boundary conditions: the applied Gate voltages
  std::vector<double> Vapp;


  /// initialise and construct the Device
  Device(double delta,
         double eps,
         double qdensity,
         const std::vector<double>& appV,
         const posvec_t& min,
         const posvec_t& max)
  {
    m_skin = delta;
    m_eps = eps;
    rho0 = qdensity;
    r_min = min;
    r_max = max;
    q_fac = rho0 / ( 6.0 * m_eps );
    Vapp.resize(appV.size());
    Vapp = appV;
    //    flush();
  }

  /// make the maximum sphere, i.e. find the nearest distance to boundary
  void MaximumSphere( Domain& );

  /// check for first passage out of the boundary
  double passage( Domain& );
  double OC_passage( const double,
                     Domain&,
                     WOSParticles*);

  double contrib0( int, const Domain&, WOSParticles* );
  double contribk( const Domain&, const WOSParticles* );
  double OC_contrib0( double,
                      const posvec_t&,
                      WOSParticles*);

  void flush()
  {
    cout << m_skin << '\t' << m_eps << '\t' << rho0 << endl;
    cout << r_min << endl;
    cout << r_max << endl;
    cout << Vapp[4] << '\t' << Vapp[5] << endl;
  }


};
#endif
