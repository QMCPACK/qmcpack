#include "QMCHamiltonians/WOS/WOSPotential.h"

inline ValueType WOSPotential::evaluate(ParticleBase& P) {
  
  RealType PE = 0.0; RealType PEsq = 0.0;
  
  RealType vHartree = 0.0;
  
  /// create a domain
  Domain domain;
  
  int itermax = 1;
  //      Random.init(0,1,1356823);
  
  for(int iat = 0; iat < P.getTotalNum(); iat++){
    
    for(int iter = 0; iter < itermax ; iter++){
      
      RealType vbare = 0.0;
      RealType Gee = 0.0;
      
      /// start runner from particle position
      domain.runner = P.R[iat];
      domain.in_device = true;
      //cout << "runner: " << domain.runner << endl;
      
      /// dielectric at runner position
      QDot->epsilon(domain);
      //cout << "eps " << domain.eps_d << endl;
      
      /// initialise weight_bc
      QDot->weight_init();
      //cout << "wt: " << QDot->weight_bc << endl;
      
      /// create spherical domain within layer
      QDot->MakeLimitedSphere( domain );
      //cout << "radius: " << domain.radius << endl;
      
      vbare = QDot->Image_Contribution( iat, domain ,P );
      //	  vbare += QDot->Domain_Contribution( domain, P );
      cout << "vbare: " << vbare << endl;
      
      QDot->sample_point( domain );
      //cout << "r: " << domain.runner << endl;
      
      QDot->MakeMaximumSphere( domain );
      //cout << "radius " << domain.radius << endl;
      
      vbare += QDot->Domain_Contribution( domain, P );
      cout << "vbare: " << vbare << endl;
	  
      vbare += QDot->passage( domain );
      cout << "vbare: " << vbare << '\t' << domain.in_device << endl;
      
      while( domain.in_device ) {
	
	QDot->epsilon( domain );
	//cout << "eps" << domain.eps_d << endl;
	
	QDot->sample_point( domain );
	//cout << "r: " << domain.runner << endl;
	
	QDot->MakeMaximumSphere( domain );
	//cout << "d " << domain.radius << endl;
	
	vbare += QDot->Domain_Contribution( domain, P );
	cout << "vbare: " << vbare << endl;
	
	vbare += QDot->passage( domain );
	cout << "vbare: " << vbare << '\t' << domain.in_device << endl;
	
      }
      cout <<  endl;
      
      PE += vbare;
      PEsq += vbare*vbare;
      //	  cout << iter << '\t' << vbare << endl;
      //cout << "poten:: " << vbare << '\t' << QDot->weight_bc << endl;
    }
    
    PE /= 1.0*itermax; PEsq /= 1.0*itermax;
	
    cout <<  P.R[0][0] << '\t' << -PE << '\t' << sqrt(fabs(PEsq-PE*PE)/(1.0*itermax)) << endl;
    //      //cout << " I am here:: very good " << P.R[0] << endl;
    exit(-1);
    //cout << endl << endl << endl << endl ;

  }
  
  
  return PE;
}

