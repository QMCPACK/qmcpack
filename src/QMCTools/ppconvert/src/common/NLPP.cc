//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "NLPP.h"
#include "GKIntegration.h"
#include "MatrixOps.h"

bool
NLPPClass::IsNonlocal()
{
  return true;
}

void 
ChannelPotential::Read(IOSectionClass &in, Grid *grid)
{
  assert (in.ReadVar("l", l));
  assert (in.ReadVar("n_principal", n_principal));
  assert (in.ReadVar("Cutoff", rc));
  R0 = 1.75 * rc;
  assert (in.ReadVar("Occupation", Occupation));
  assert (in.ReadVar("Eigenvalue", Eigenvalue));
  Array<double,1> V_data, u_data;
  assert (in.ReadVar("Vl", V_data));
  assert (in.ReadVar("ul", u_data));
  V.Init (grid, V_data);
  u.Init (grid, u_data);
}

void
ChannelPotential::Write(IOSectionClass &out)
{
  out.WriteVar("l", l);
  out.WriteVar("n_principal", n_principal);
  out.WriteVar("Cutoff", rc);
  out.WriteVar("Occupation", Occupation);
  out.WriteVar("Eigenvalue", Eigenvalue);
  out.WriteVar("Vl", V.Data());
  out.WriteVar("ul", u.Data());
}

void
NLPPClass::Read(IOSectionClass &in)
{
  assert(in.ReadVar("AtomicNumber", AtomicNumber));
  assert(in.ReadVar("LocalChannel", lLocal));
  assert(in.ReadVar("ValenceCharge", Zion));
  assert(in.ReadVar("Symbol", Symbol));
  
  int numChannels = in.CountSections("lChannel");
  Vl.resize(numChannels);
//   vector<Array<double,1> > vl(numChannels), ul(numChannels);
//   vector<double> rc(numChannels);
  assert (in.OpenSection("PotentialGrid"));
  PotentialGrid = ReadGrid (in);
  in.CloseSection(); // "PotentialGrid"
  for (int i=0; i<numChannels; i++) {
    assert (in.OpenSection("lChannel", i));
    int l;
    assert (in.ReadVar("l", l));
    if (l > numChannels) {
      cerr << "Skipped channels in NLPPClass read.\n";
      abort();
    }
    Vl[l].Read(in, PotentialGrid);
    in.CloseSection (); // "lChannel"
  }
  for (int l=0; l<numChannels; l++) {
    Array<double,1> deltav(PotentialGrid->NumPoints);
    deltav = Vl[l].V.Data()- Vl[lLocal].V.Data();
    Vl[l].DeltaV.Init (PotentialGrid, deltav);
  }
}


void
NLPPClass::Write(IOSectionClass &out)
{
  out.WriteVar ("Type", "NLPP");
  out.WriteVar ("LocalChannel", lLocal);
  out.WriteVar ("AtomicNumber", AtomicNumber);
  out.WriteVar ("Symbol", Symbol);
  out.WriteVar ("ValenceCharge", Zion);
  out.NewSection("PotentialGrid");
  PotentialGrid->Write(out);
  out.CloseSection();  // "PotentialGrid"
  for (int i=0; i<Vl.size(); i++) {
    out.NewSection("lChannel");
    Vl[i].Write (out);
    out.CloseSection(); // "lChannel"
  }
}



// This function computes the Kleinman-Bylander projector for the
// nonlocal channels of the pseudopotential.  It just calls
// SetupProjector on each of the nonlocal channels.
void
NLPPClass::SetupProjectors(double G_max, double G_FFT)
{
  for (int i=0; i<Vl.size(); i++) 
    if (Vl[i].l != lLocal) 
      Vl[i].SetupProjector(G_max, G_FFT);

  FILE *fout = fopen ("zeta_r.dat", "w");
  for (double r=0.0; r<50.0; r+=0.001) {
    fprintf (fout, "%12.16e ", r);
      for (int l=0; l<Vl.size(); l++)
	if (l != lLocal)
	  fprintf (fout, "%12.16e ", Vl[l].zeta_r(r));
    fprintf (fout, "\n");
  }
  fclose (fout);
  
  fout = fopen ("zeta_q.dat", "w");
  for (double q=0.0; q<G_FFT; q+=0.001) {
    fprintf (fout, "%12.16e ", q);
      for (int l=0; l<Vl.size(); l++)
	if (l != lLocal)
	  fprintf (fout, "%12.16e ", Vl[l].zeta_q(q));
    fprintf (fout, "\n");
  }
  fclose (fout);

  fout = fopen ("chi_q.dat", "w");
  for (double q=0.0; q<G_FFT; q+=0.001) {
    fprintf (fout, "%12.16e ", q);
      for (int l=0; l<Vl.size(); l++)
	if (l != lLocal)
	  fprintf (fout, "%12.16e ", Vl[l].chi_q(q));
    fprintf (fout, "\n");
  }
  fclose (fout);

  fout = fopen ("chi_r.dat", "w");
  for (double r=0.0; r<50.0; r+=0.001) {
    fprintf (fout, "%12.16e ", r);
      for (int l=0; l<Vl.size(); l++)
	if (l != lLocal)
	  fprintf (fout, "%12.16e ", Vl[l].chi_r(r));
    fprintf (fout, "\n");
  }
  fclose (fout);

}


// A(q,q') = q^2 q'^2 \int_0^R0 j_l(qr) r^2 j_l(q'r) dr
double
ChannelPotential::A(double q, double qp)
{
  if (l == 0) {
    if (q == qp) 
      return -0.25*q*(-2.0*q*R0 + sin(2.0*q*R0));
    else
      return -1.0/(q*q - qp*qp) * 
	q * qp *(q*cos(q*R0)*sin(qp*R0) - qp*cos(qp*R0)*sin(q*R0));
  }
  else if (l == 1){
    if (q == qp) 
      return (-2.0 + 2.0*q*q*R0*R0 + 2.0*cos(2.0*q*R0) + q*R0*sin(2.0*q*R0))/(4.0*R0);
    else
      return q*qp/(q*q - qp*qp) *
	(q*cos(qp*R0)*sin(q*R0) - qp*cos(q*R0)*sin(qp*R0))
	- 1.0/R0 * sin(q*R0)*sin(qp*R0);
  }
  else if (l == 2) {
    if (qp == 0.0 || q == 0.0)
      return 0.0;
    else if (q == qp) 
      return 1.0/(4*q*q*R0*R0*R0)*
	(-6.0 -6.0*q*q*R0*R0 + 2.0*q*q*q*q*R0*R0*R0*R0 + 
	 (6.0-6.0*q*q*R0*R0)*cos(2.0*q*R0) + 
	 q*R0*(12.0-q*q*R0*R0)*sin(2.0*q*R0));
    else 
      return 1.0/(q*qp*(q*q-qp*qp)*R0*R0*R0) *
	(sin(q*R0)*(qp*R0*(-3.0*qp*qp + q*q*(3.0+qp*qp*R0*R0))*cos(qp*R0) +
		    3.0*(qp*qp-q*q)*sin(qp*R0)) - q*R0*cos(q*R0) *
	 (3.0*qp*(q*q - qp*qp)*R0*cos(qp*R0) + 
	  (3.0*qp*qp + q*q*(-3.0 + qp*qp * R0*R0))*sin(qp*R0)));
  }
  else
    return sqrt(-1.0);
}

// double
// ChannelPotential::A(double q, double qp)
// {
//   double sum = 0.0;
//   double delta = 0.01;
//   for (double r=0.0; r<R0; r+= delta)
//     sum += jl(l,q*r)*jl(l,qp*r)*r*r;
//   sum *= delta*q*q*qp*qp;
//   return sum;
// }


// This function computes the Kleinman-Bylander projector from the
// DeltaV and u splines.  It then applies the methods of King-Smith et
// al to filter out the high Fourier components, so that the
// projector can be applied without error in real-space.
void
ChannelPotential::SetupProjector (double G_max, double G_FFT)
{
  Grid &grid = *u.grid;
  // First, compute zeta_r, normalization, and E_KB
  Job = NORM;
  GKIntegration<ChannelPotential> integrator(*this);
  double norm = integrator.Integrate(0.0, grid.End, 1.0e-12);
  ProjectorNorm = 1.0/sqrt(norm);
  Job = EKB;
  E_KB = norm/integrator.Integrate(0.0, grid.End, 1.0e-12);
  cerr << "l = " << l << "  Norm is " << norm 
       << "  E_KB is " << E_KB << "  R0 = " << R0 << endl;
  
  // Compute zeta(r)
  Array<double,1> zeta(grid.NumPoints);
  zeta(0) = ProjectorNorm * DeltaV(0)*u(1.0e-8)*1.0e8;
  for (int i=1; i<grid.NumPoints; i++) 
    zeta(i) = ProjectorNorm * DeltaV(i)*u(i)/grid(i);
  zeta_r.Init (&grid, zeta);

  // Compute zeta(q)
  Job = ZETA_Q;
  qGrid.Init (0.0, G_FFT, 1000);
  zeta.resize(qGrid.NumPoints);
  for (int i=0; i<qGrid.NumPoints; i++) {
    qCurrent = qGrid(i);
    zeta(i) = integrator.Integrate(0.0, grid.End, 1.0e-12);
  }
  zeta_q.Init (&qGrid, zeta);


  double gamma = G_FFT - G_max;
  // Zero out zeta_q above gamma;
  Array<double,1> chi_q_data(qGrid.NumPoints);
  for (int i=0; i<qGrid.NumPoints; i++) 
    chi_q_data(i) = (qGrid(i) >= gamma) ? 0.0 : zeta_q(i);
  chi_q.Init  (&qGrid, chi_q_data);  

  // Now for the magic:  We adjust chi_q between G_max and gamma so
  // that the real-space oscillations outside R0 are damped out
  // See King-Smith et al, PRB 44 13063
  // Find index of gamma
  int gammaIndex = qGrid.ReverseMap(gamma);
  int G_maxIndex  = qGrid.ReverseMap(G_max)+1;
  double delta = qGrid(1)-qGrid(0);
  int nb = gammaIndex - G_maxIndex;
  Array<double,1> b(nb), x(nb);
  // First, create the b vector
  b = 0.0;
  for (int i=0; i<nb; i++) {
    double q = qGrid(G_maxIndex+i);
    for (int j=0; j<G_maxIndex; j++) {
      double qp = qGrid (j);
      b(i) -= delta * A(q, qp)*zeta_q(j);
    }
  }
  // Now, create the M matrix
  Array<double,2> M(nb, nb);
  for (int i=0; i<nb; i++) {
    double q = qGrid(G_maxIndex+i);
    for (int j=0; j<nb; j++) {
      double qp = qGrid(G_maxIndex+j);
      M(i,j) = delta*A(q, qp);
    }
    M(i,i) -= 0.5*M_PI*q*q;
  }
//   if (l==2) {
//     FILE *fout = fopen ("M.dat", "w");
//     for (int i=0; i<nb; i++) {
//       for (int j=0; j<nb; j++) 
// 	fprintf (fout, "%24.16e ", M(i,j));
//       fprintf (fout, "\n");
//     }
//     fclose (fout);
//     fout = fopen ("b.dat", "w");
//     for (int i=0; i<nb; i++)
//       fprintf (fout, "%24.16e\n", b(i));
//     fclose (fout);
//   }
  // Now solve Mx = b
  Array<int,1> perm;
  double sign;
  LUdecomp (M, perm, sign);
  LUsolve (M, perm, b);
  x = b;
  for (int i=0; i<nb; i++) 
    chi_q(G_maxIndex+i) = x(i);
  chi_q(G_maxIndex+nb) = 0.0;
  
//   FILE *fout = fopen("chicheck.dat","w");
//   for (double q=0.0; q<qGrid.End; q+=0.001)
//     fprintf (fout, "%24.16e %24.16e\n", q, chi_q(q));
//   fclose (fout);

  // Now transform back to real-space, computing chi(r)
  Job = CHI_R;
  Array<double,1> chi_r_data(grid.NumPoints);
  for (int i=1; i<grid.NumPoints; i++) {
    rCurrent = grid(i);
    chi_r_data(i) = integrator.Integrate(0.0, gamma, 1.0e-10);
  }
  chi_r.Init (&grid, chi_r_data);

  // Finally, check to see if chi_r is small outside R0
  Job = CHECK_CHI_R;
  double norm2 = integrator.Integrate(0.0, grid.End, 1.0e-8);
  double error = integrator.Integrate( R0, grid.End, 1.0e-8);
  if (error > 1.0e-16)
    cerr << "Fractional error in real-space projection = "
	 << (error / norm2) << endl;
}


double
NLPPClass::V(int l, double r)
{
  if (l < Vl.size())
    return Vl[l].V(r);
  else
    return Vl[lLocal].V(r);
}

double
NLPPClass::dVdr(int l, double r)
{
  if (l < Vl.size())
    return Vl[l].V.Deriv(r);
  else
    return Vl[lLocal].V.Deriv(r);
}

double
NLPPClass::d2Vdr2(int l, double r)
{
  if (l < Vl.size())
    return Vl[l].V.Deriv2(r);
  else
    return Vl[lLocal].V.Deriv2(r);
}

double
NLPPClass::V(double r)
{
  return Vl[lLocal].V(r);
}

double
NLPPClass::dVdr(double r)
{
  return Vl[lLocal].V.Deriv(r);
}

double
NLPPClass::d2Vdr2(double r)
{
  return Vl[lLocal].V.Deriv2(r);
}
