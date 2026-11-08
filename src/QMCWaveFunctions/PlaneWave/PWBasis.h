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


/** @file PWBasis.h
 * @brief Plane-wave basis values and Cartesian derivative work buffers.
 */
#ifndef QMCPLUSPLUS_PLANEWAVEBASIS_BLAS_H
#define QMCPLUSPLUS_PLANEWAVEBASIS_BLAS_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Message/Communicate.h"
#include "CPU/e2iphi.h"
#include "hdf/hdf_archive.h"

/** If defined, use recursive method to build the basis set for each position
 *
 * performance improvement is questionable: load vs sin/cos
 */
//#define PWBASIS_USE_RECURSIVE

namespace qmcplusplus
{
/** Evaluate primitive-cell plane waves exp(i (k+G).r) at Cartesian positions.
 *
 * The reduced twist and integer G vectors are converted using the primitive
 * lattice, even when the target ParticleSet is a tiled supercell. Evaluation
 * updates reusable buffers for contraction with orbital coefficients; concurrent
 * evaluations require separate basis instances.
 */
class PWBasis : public QMCTraits
{
public:
  using ParticleLayout = qmcplusplus::Lattice;
  using GIndex_t       = TinyVector<IndexType, 3>;

private:
  /// Initialization marker; stores the largest recursion extent in recursive mode.
  int maxmaxg;
  /// Largest absolute G index in each dimension for recursive phase generation.
  GIndex_t maxg;
  /// Legacy cutoff value; the current implementation retains every input G vector.
  RealType ecut;
  /// Primitive-cell twist in reduced reciprocal coordinates.
  PosType twist;
  /// Primitive-cell twist converted to a Cartesian wave vector.
  PosType twist_cart;

  /// Integer reciprocal-lattice coefficients of the input G vectors.
  std::vector<GIndex_t> gvecs;
  /// Nonnegative recursion indices, gvecs[ig][idim]+maxg[idim].
  std::vector<GIndex_t> gvecs_shifted;

  /// Negative squared Cartesian lengths of k+G, used for Laplacians.
  std::vector<RealType> minusModKplusG2;
  /// Cartesian k+G vectors used for phases and derivatives.
  std::vector<PosType> kplusgvecs_cart;

  /// Per-dimension phase factors for the optional recursive evaluation path.
  Matrix<ComplexType> C;
  /// Legacy real-basis sign flags; current real projections belong to PWRealOrbitalSet.
  std::vector<int> negative;

public:
  /// Column layout of Z: value, Laplacian, and Cartesian gradient components.
  enum
  {
    PW_VALUE,
    PW_LAP,
    PW_GRADX,
    PW_GRADY,
    PW_GRADZ,
    PW_MAXINDEX
  };

  /// Basis VGL data indexed by (G-vector, PW_* column), filled by evaluateAll().
  Matrix<ComplexType> Z;

  /// Basis values indexed by G-vector, filled by evaluate().
  Vector<ComplexType> Zv;

  /// Hessians indexed by (G-vector, row*OHMMS_DIM+column), filled by evaluateHessian().
  Matrix<ComplexType> Zh;

  /// Cartesian phases (k+G).r used by the direct exponential evaluator.
  Vector<RealType> phi;

  /** Map an input coefficient index to its retained basis row.
   *
   * Orbital coefficient loaders honor negative entries as excluded G vectors.
   * The current trimforecut() retains all vectors, so this map is an identity.
   */
  std::vector<int> inputmap;

  ///total number of basis functions
  int NumPlaneWaves;

  /// Primitive lattice used to convert reduced reciprocal coordinates.
  ParticleLayout Lattice;

  ///default constructor
  PWBasis() : maxmaxg(0), NumPlaneWaves(0) {}

  /// Select an initial reduced twist; readbasis() supplies the lattice and G vectors.
  PWBasis(const PosType& twistangle) : maxmaxg(0), twist(twistangle), NumPlaneWaves(0) {}

  ~PWBasis() {}

  /// Set the reduced primitive twist and rebuild the basis work arrays if it changes.
  void setTwistAngle(const PosType& tang);

  /// Prepare shifted reciprocal vectors and resize the evaluation work arrays.
  void reset();

  /** Load G vectors from /electrons/kpoint_0/gvectors in an orbital archive.
   *
   * The file's /supercell/primitive_vectors overrides lat when present. All input
   * G vectors are retained; ecutoff does not currently truncate this list.
   *
   * @param h5basisgroup open orbital HDF archive
   * @param ecutoff legacy cutoff value to store
   * @param lat fallback lattice if primitive vectors are absent from the file
   * @param pwname legacy group-name argument; the ESHDF G-vector path is fixed
   * @param pwmultname legacy multiplier tag used in diagnostic output
   * @param resizeContainer whether to call reset() after loading the vectors
   * @return number of input plane waves
   * @throws std::runtime_error on an HDF read failure or an empty G-vector list
   */
  int readbasis(hdf_archive& h5basisgroup,
                RealType ecutoff,
                const ParticleLayout& lat,
                const std::string& pwname     = "planewaves",
                const std::string& pwmultname = "multipliers",
                bool resizeContainer          = true);

  /** Prepare Cartesian k+G vectors, Laplacian factors, and the coefficient map.
   *
   * Despite its legacy name, this routine currently retains every G vector;
   * the kinetic-energy cutoff branch is disabled.
   */
  void trimforecut();

#if defined(PWBASIS_USE_RECURSIVE)
  /** Fill the recursion coefficients matrix.
   *
   * @todo Generalize to non-orthorohmbic cells
   */
  inline void BuildRecursionCoefs(const PosType& pos)
  {
    PosType tau_red(Lattice.toUnit(pos));
//      RealType phi=TWOPI*tau_red[0];
//      RealType nphi=maxg0*phi;
//      ComplexType ct0(std::cos(phi),std::sin(phi));
//      ComplexType t(std::cos(nphi),-std::sin(nphi));
//      C0[0]=t;
//      for(int n=1; n<=2*maxg0; n++) C0[n] = (t *= ct0);
//
//      phi=TWOPI*tau_red[1];
//      nphi=maxg1*phi;
//      ct0=ComplexType(std::cos(phi),std::sin(phi));
//      t=ComplexType(std::cos(nphi),-std::sin(nphi));
//      C1[0]=t;
//      for(int n=1; n<=2*maxg1; n++) C1[n] = (t *= ct0);
//
//      phi=TWOPI*tau_red[2];
//      nphi=maxg2*phi;
//      ct0=ComplexType(std::cos(phi),std::sin(phi));
//      t=ComplexType(std::cos(nphi),-std::sin(nphi));
//      C2[0]=t;
//      for(int n=1; n<=2*maxg2; n++) C2[n] = (t *= ct0);
#pragma ivdep
    for (int idim = 0; idim < 3; idim++)
    {
      int ng        = maxg[idim];
      RealType phi  = TWOPI * tau_red[idim];
      RealType nphi = ng * phi;
      ComplexType Ctemp(std::cos(phi), std::sin(phi));
      ComplexType t(std::cos(nphi), -std::sin(nphi));
      ComplexType* restrict cp_ptr = C[idim];
      *cp_ptr++                    = t;
      for (int n = 1; n <= 2 * ng; n++)
      {
        *cp_ptr++ = (t *= Ctemp);
      }
    }
    //Base version
    //#pragma ivdep
    //      for(int idim=0; idim<3; idim++){
    //        RealType phi=TWOPI*tau_red[idim];
    //        ComplexType Ctemp(std::cos(phi),std::sin(phi));
    //        int ng=maxg[idim];
    //        ComplexType* restrict cp_ptr=C[idim]+ng;
    //        ComplexType* restrict cn_ptr=C[idim]+ng-1;
    //        *cp_ptr=1.0;
    //        for(int n=1; n<=ng; n++,cn_ptr--){
    //          ComplexType t(Ctemp*(*cp_ptr++));
    //          *cp_ptr = t;
    //          *cn_ptr = conj(t);
    //        }
    //      }
    //Not valid for general supercell
    //      // Cartesian of twist for 1,1,1 (reduced coordinates)
    //      PosType G111(1.0,1.0,1.0);
    //      G111 = Lattice.k_cart(G111);
    //
    //      //Precompute a small number of complex factors (PWs along b1,b2,b3 lines)
    //      //using a fast recursion algorithm
    //#pragma ivdep
    //      for(int idim=0; idim<3; idim++){
    //        //start the recursion with the 111 vector.
    //        RealType phi = pos[idim] * G111[idim];
    //        register ComplexType Ctemp(std::cos(phi), std::sin(phi));
    //        int ng=maxg[idim];
    //        ComplexType* restrict cp_ptr=C[idim]+ng;
    //        ComplexType* restrict cn_ptr=C[idim]+ng-1;
    //        *cp_ptr=1.0;
    //        for(int n=1; n<=ng; n++,cn_ptr--){
    //          ComplexType t(Ctemp*(*cp_ptr++));
    //          *cp_ptr = t;
    //          *cn_ptr = conj(t);
    //        }
    //      }
  }

  /// Fill Zv at a Cartesian position using recursive phase generation.
  inline void evaluate(const PosType& pos)
  {
    BuildRecursionCoefs(pos);
    RealType twistdotr = dot(twist_cart, pos);
    ComplexType pw0(std::cos(twistdotr), std::sin(twistdotr));
    //Evaluate the planewaves for particle iat.
    for (int ig = 0; ig < NumPlaneWaves; ig++)
    {
      //PW is initialized as exp(i*twist.r) so that the final basis evaluations are for (twist+G).r
      ComplexType pw(pw0); //std::cos(twistdotr),std::sin(twistdotr));
      for (int idim = 0; idim < 3; idim++)
        pw *= C(idim, gvecs_shifted[ig][idim]);
      //pw *= C0[gvecs_shifted[ig][0]];
      //pw *= C1[gvecs_shifted[ig][1]];
      //pw *= C2[gvecs_shifted[ig][2]];
      Zv[ig] = pw;
    }
  }
  /** Fill Z with basis values, gradients, and Laplacians at P.activeR(iat).
   * This recursive path does not populate Zv or the Hessian buffer.
   */
  inline void evaluateAll(const ParticleSet& P, int iat)
  {
    const PosType& r(P.activeR(iat));
    BuildRecursionCoefs(r);
    RealType twistdotr = dot(twist_cart, r);
    ComplexType pw0(std::cos(twistdotr), std::sin(twistdotr));
    //Evaluate the planewaves and derivatives.
    ComplexType* restrict zptr = Z.data();
    for (int ig = 0; ig < NumPlaneWaves; ig++, zptr += 5)
    {
      //PW is initialized as exp(i*twist.r) so that the final basis evaluations
      //are for (twist+G).r
      ComplexType pw(pw0);
      // THE INDEX ORDER OF C DOESN'T LOOK TOO GOOD: this could be fixed
      for (int idim = 0; idim < 3; idim++)
        pw *= C(idim, gvecs_shifted[ig][idim]);
      //pw *= C0[gvecs_shifted[ig][0]];
      //pw *= C1[gvecs_shifted[ig][1]];
      //pw *= C2[gvecs_shifted[ig][2]];
      zptr[0] = pw;
      zptr[1] = minusModKplusG2[ig] * pw;
      zptr[2] = kplusgvecs_cart[ig][0] * ComplexType(-pw.imag(), pw.real());
      zptr[3] = kplusgvecs_cart[ig][1] * ComplexType(-pw.imag(), pw.real());
      zptr[4] = kplusgvecs_cart[ig][2] * ComplexType(-pw.imag(), pw.real());
    }
  }
#else
  /// Fill Zv at a Cartesian position using direct complex exponentials.
  inline void evaluate(const PosType& pos)
  {
    //Evaluate the planewaves for particle iat.
    for (int ig = 0; ig < NumPlaneWaves; ig++)
      phi[ig] = dot(kplusgvecs_cart[ig], pos);
    eval_e2iphi(NumPlaneWaves, phi.data(), Zv.data());
  }
  /// Fill Zv and Z with values and VGL data at P.activeR(iat); Hessian Zh is unchanged.
  inline void evaluateAll(const ParticleSet& P, int iat)
  {
    const PosType& r(P.activeR(iat));
    evaluate(r);
    ComplexType* restrict zptr = Z.data();
    for (int ig = 0; ig < NumPlaneWaves; ig++, zptr += 5)
    {
      //PW is initialized as exp(i*twist.r) so that the final basis evaluations
      //are for (twist+G).r
      ComplexType& pw = Zv[ig];
      zptr[0]         = pw;
      zptr[1]         = minusModKplusG2[ig] * pw;
      zptr[2]         = kplusgvecs_cart[ig][0] * ComplexType(-pw.imag(), pw.real());
      zptr[3]         = kplusgvecs_cart[ig][1] * ComplexType(-pw.imag(), pw.real());
      zptr[4]         = kplusgvecs_cart[ig][2] * ComplexType(-pw.imag(), pw.real());
    }
  }
#endif

  /** Fill Zh with Cartesian Hessians -K_row*K_column*Z(ig, PW_VALUE), K=k+G.
   *
   * @pre evaluateAll() has populated Z for the desired particle position.
   * This routine reuses those values and does not update Zv or Z.
   */
  void evaluateHessian();

  //    /** Fill the recursion coefficients matrix.
  //     *
  //     * @todo Generalize to non-orthorohmbic cells
  //     */
  //    void BuildRecursionCoefsByAdd(const PosType& pos)
  //    {
  //      // Cartesian of twist for 1,1,1 (reduced coordinates)
  //      PosType G111(1.0,1.0,1.0);
  //      G111 = Lattice.k_cart(G111);
  //      //PosType redP=P.Lattice.toUnit(P.R[iat]);
  //      //Precompute a small number of complex factors (PWs along b1,b2,b3 lines)
  //      for(int idim=0; idim<3; idim++){
  //        //start the recursion with the 111 vector.
  //        RealType phi = pos[idim] * G111[idim];
  //        int ng(maxg[idim]);
  //        RealType* restrict cp_ptr=logC[idim]+ng;
  //        RealType* restrict cn_ptr=logC[idim]+ng-1;
  //        *cp_ptr=0.0;
  //        //add INTEL vectorization
  //        for(int n=1; n<=ng; n++,cn_ptr--){
  //          RealType t(phi+*cp_ptr++);
  //          *cp_ptr = t;
  //          *cn_ptr = -t;
  //        }
  //      }
  //    }
};
} // namespace qmcplusplus
#endif
