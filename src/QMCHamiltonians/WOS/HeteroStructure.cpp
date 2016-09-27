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
    
    
#include "QMCHamiltonians/WOS/HeteroStructure.h"

#include <fstream>
#include <iomanip>

void HeteroStructure::construct()
{
  /// add bottom interface
  add_Interface( new Interface( intervals_d[0], 0.0, 0.0 ) );
  /// add bottom layer
  posvec_t bot_corner = r_min;
  bot_corner[2] = intervals_d[0];
  posvec_t top_corner = r_max;
  top_corner[2] = intervals_d[1];
  add_Layer( new Layer(eps_d[0], offset_d[0],bot_corner,top_corner ) );
  for( int ilayer = 1; ilayer < nlayers_d; ilayer++ )
  {
    /// add a layer
    bot_corner = r_min;
    bot_corner[2] = intervals_d[ilayer];
    top_corner = r_max;
    top_corner[2] = intervals_d[ilayer + 1];
    add_Layer( new Layer(eps_d[ ilayer ],
                         offset_d[ ilayer ],bot_corner,top_corner ) );
    /// add the interface above the layer
    double diff_eps = layers[ilayer - 1]->eps_d - layers[ilayer]->eps_d;
    double av_eps = 0.5 * ( layers[ilayer-1]->eps_d + layers[ilayer]->eps_d );
    add_Interface( new Interface(intervals_d[ilayer], diff_eps, av_eps) );
  }
  add_Interface( new Interface(intervals_d[nlayers_d], 0.0, 0.0 ) );
  return;
}


void HeteroStructure::MakeLimitedSphere( Domain& domain )
{
  int id = 0;
  int nlayers = layers.size();
  for( int ilayer = 0; ilayer < nlayers; ilayer++ )
  {
    if( !layers[ id ]->find( domain ) )
    {
      id++;
    }
    else
    {
      layers[ id ]->makeSphere( domain );
      return;
    }
  }
  id = 0;
  int nifc = interfaces.size();
  for(int ifc = 0; ifc < nifc; ifc++)
  {
    if(!interfaces[ id ]->find( domain ) )
    {
      id++;
    }
    else
    {
      domain.inter_face = true;
      MakeIfcSphere( domain, id );
      return;
    }
  }
}


void HeteroStructure::MakeIfcSphere(Domain& domain, int ifc)
{
  double d[6];
  d[0] = domain.runner[0] - r_min[0];
  d[1] = r_max[0] - domain.runner[0];
  d[2] = domain.runner[1] - r_min[1];
  d[3] = r_max[1]-domain.runner[1];
  d[4] = domain.runner[2] - interfaces[ifc-1]->z_val;
  d[5] = interfaces[ifc+1]->z_val - domain.runner[2];
  double* it = std::min_element(d,d+6);
  domain.radius = (*it);
  return;
}

void HeteroStructure::MakeMaximumSphere(Domain& domain)
{
  double d[6];
  d[0] = domain.runner[0]-r_min[0];
  d[1] = r_max[0]-domain.runner[0];
  d[2] = domain.runner[1]-r_min[1];
  d[3] = r_max[1]-domain.runner[1];
  d[4] = domain.runner[2]-r_min[2];
  d[5] = r_max[2]-domain.runner[2];
  double* it = std::min_element(d,d+6);
  domain.radius = (*it);
  return;
}

void HeteroStructure::add_Layer(Layer* aLayer)
{
  int id = layers.size();
  layers.push_back( aLayer );
  layers[ id ]->setid(id);
  r_min[0] = std::min( layers[id]->r_min[0], r_min[0] );
  r_max[0] = std::max( layers[id]->r_max[0], r_max[0] );
}

void HeteroStructure::add_Interface(Interface* anInterface)
{
  int id = interfaces.size();
  interfaces.push_back(anInterface);
  interfaces[id]->setid(id);
}


void HeteroStructure::calc_dfrac(Domain& domain)
{
  for( int i = 0; i < interfaces.size(); i++)
    interfaces[i]->calc_dfrac(domain.runner[2],domain.radius);
}

void HeteroStructure::sample_prob(const Domain& domain)
{
  double Stotal = 0.0;
  int nlayers = layers.size();
  /// probability ratios
  for(int ilayer = 0; ilayer < nlayers; ilayer++)
  {
    double frac1 = interfaces[ilayer+1]->d_frac;
    double frac2 = interfaces[ilayer]->d_frac;
    layers[ilayer]->sample_prob(domain,frac1,frac2);
    interfaces[ilayer]->sample_prob(domain);
    Stotal += layers[ilayer]->prob_d + interfaces[ilayer]->prob_d;
  }
  /// Normalise them
  double cum_prob = 0.0;
  for(int ilayer = 0; ilayer < nlayers; ilayer++)
  {
    layers[ilayer]->prob_d /= Stotal;
    cum_prob += layers[ilayer]->prob_d;
    layers[ilayer]->prob_d = cum_prob;
  }
  for(int ilayer = 0; ilayer < nlayers; ilayer++)
  {
    interfaces[ilayer]->prob_d /= Stotal;
    cum_prob += interfaces[ilayer]->prob_d;
    interfaces[ilayer]->prob_d = cum_prob;
  }
  weight_bc *= Stotal;
  return;
}


double HeteroStructure::epsilon(const posvec_t& r)
{
  int nlayers = layers.size();
  for( int ilayer = 0; ilayer < nlayers; ilayer++ )
  {
    if(layers[ilayer]->find(r))
    {
      return layers[ilayer]->eps_d;
    }
  }
  for( int ilayer = 1; ilayer < nlayers; ilayer++ )
  {
    if(interfaces[ilayer]->find(r))
    {
      return interfaces[ilayer]->eps_d;
    }
  }
}

void HeteroStructure::epsilon(Domain& domain)
{
  domain.eps_d = epsilon(domain.runner);
  return;
}

void HeteroStructure::sample_point(Domain& domain)
{
  calc_dfrac(domain);
  sample_prob(domain);
  int nlayers = layers.size();
  double phi = 2.0 * M_PI * Random();
  double rnd = Random();
  double hop = Random();
  for( int ilayer = 0; ilayer < nlayers; ilayer++)
  {
    if ( hop <= layers[ilayer]->prob_d)
    {
      double theta = acos( ( 1.0 - rnd ) * interfaces[ ilayer + 1 ]->d_frac
                           + rnd * interfaces[ ilayer ]->d_frac );
      domain.WalkOnSphere( theta, phi );
      return;
    }
  }
  for(int ilayer = 0; ilayer < nlayers; ilayer++)
  {
    if(hop <= interfaces[ilayer]->prob_d)
    {
      double d = std::abs( interfaces[ilayer]->d_frac );
      double denom = 1.0 - ( 1.0 - d ) * rnd;
      denom = denom * denom;
      double rho = domain.radius * d * sqrt( 1.0 / denom - 1.0 );
      weight_bc *= interfaces[ ilayer ]->xsign;
      domain.WalkOnDisk( rho, phi, interfaces[ ilayer ]->z_val );
      return;
    }
  }
}

double HeteroStructure::Image_Contribution(int ipt,
    const Domain& domain,
    const qmcplusplus::ParticleBase& P)
{
  double Gself = -Gfac/(domain.eps_d*domain.radius + 1e-30);
  double Gee = 0.0;
  /// calculate distance of domain from particle
  for(int iat = 0; iat < P.getTotalNum(); iat++)
  {
    if( iat != ipt )
    {
      posvec_t rs = P.R[iat] - domain.runner;
      double dist = sqrt(dot(rs,rs));
      if( dist < domain.radius )
      {
        Gee += ( domain.radius - dist ) /
               ( domain.eps_d * domain.radius * dist + 1e-30);
      }
    }
  }
  //  std::cout << "Gself & Gee: " << Gself << '\t' << Gee << std::endl;
  Gself += Gee;
  return 0.5*Gself;
}

double HeteroStructure::Domain_Contribution(const Domain& domain,
    const qmcplusplus::ParticleBase& P)
{
  double d = domain.radius;
  //  rho * asq / 6 eps0
  //  rho = n e = n * 1.60219e-19, asq = (.529e-10)^2; eps0 = 8.8543373e-12
  //  1/4pie0 = 8.9874e9, and aB = 0.52917725e-10m
  // gives energy in eV, so / by 27.11
  // 1 Ha = 27.211396eV
  //  double units = 3.1122e-31;
  //double units = 3.0158289e-27;  // correct
  //double units = 3.06158289e-27;   // for eV and nm
  double units = 3.04158289e-27;   // for eV and aB/nm rad scaled
  //double units = 3.1506304e-31;      // Ha and aB
  /// sample point according to Green's function
  double phi = 2.0 * M_PI * Random();
  double cosq = 2.0 * Random() - 1.0;
  //double cosq = 2.0 * 0.701 - 1.0;
  /// solve for d by Newton's method
  double dr = 0.5;
  double xi = Random();
  for( int i = 0; i <= 10; i++ )
  {
    double numer = dr * dr * ( 2.0 * dr - 3.0 ) + xi;
    double denom = 6.0 * dr * ( dr - 1.0 );
    dr -= numer / ( denom + 1.e-10 );
  }
  dr *= d;
  /// sample the point r
  posvec_t r;
  double drsinq = dr * sqrt ( 1.0 - cosq * cosq );
  r[0] = domain.runner[0] + drsinq * cos ( phi );
  r[1] = domain.runner[1] + drsinq * sin ( phi );
  r[2] = domain.runner[2] + dr * cosq;
  /// find the index on Grid3D
  gridvec_t ir = Full_Grid->ptn(r);
  int index = Full_Grid->index(Full_Grid->ptn(r));
  double rhov = rhoq[index];
  //double rhov = rho_Spline->evaluate(r);
  double radius = domain.radius*0.052917725;
  //  std::cout << domain.radius << '\t' << radius << std::endl;
  rhov = units * rhov * radius * radius / domain.eps_d;
  //cout << "rhoContrib: " << r << '\t' << temp << '\t' << rhov << std::endl;
  double Gee = 0.0;
  /// calculate distance of domain from particle
  for(int iat = 0; iat < P.getTotalNum(); iat++)
  {
    posvec_t rs = P.R[iat] - domain.runner;
    double dist = sqrt(dot(rs,rs));
    if( dist < domain.radius )
    {
      Gee += ( domain.radius - dist ) /
             ( domain.eps_d * domain.radius * dist + 1e-30);
    }
  }
  Gee *= Gfac;
  //  std::cout << "rhov & Gee " << rhov << '\t' << Gee << '\t' << ir << std::endl;
  /// only half of the e-e interaction to prevent double counting
  double contribution = rhov + 0.5*Gee;
  return contribution;
}

void HeteroStructure::read_file(int sign,
                                const char* file,
                                std::vector<double>& vec)
{
  std::ifstream infile(file, ios_base::in);
  double val;
  double eVtoHa = 0.036749309;
  double units = 1.0;//eVtoHa;
  for(int i = 0; i < Full_Grid->size; i++)
  {
    infile >> val;
    vec[i] = sign*units*val;
  }
  infile.close();
  return;
}

double HeteroStructure::passage(Domain& domain)
{
  /// sampled boundary Voltage
  double vapp = 0.0;
  if( domain.radius <= skin_d )
  {
    /// nearest Grid point and its index
    gridvec_t ir = Full_Grid->ptn(domain.runner);
    int index = Full_Grid->index(ir);
    /// wieghted (different eps) and sampled potential
    vapp = weight_bc * bare_potential[index];
    /// runner now on surface
    domain.in_device = false;
    //cout << "witihin skin : " << domain.runner << std::endl;
    //    std::cout << "gridpt: " << ir << " index: " << index << '\t' << bare_potential[index] << std::endl;
  }
  return vapp;
}
