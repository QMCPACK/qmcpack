//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
N = 256;
k = 13;

A = rand(N,N);
Ainv = A^(-1);

AinvT = Ainv';

delta = rand(1,N);

AinvTp = AinvT;

rowk = AinvT(k,:);

ratio = 1.0 + (delta*Ainv)(k)
ratio2 = dot(Ainv(:,k),A(k,:)+delta)
rInv = 1.0/ratio

for row=[1:N]
  AinvTp(row,:) -= rInv*dot(AinvT(row,:),delta)*rowk;
end

Anew = A;
Anew(k,:) += delta;

AnewInvT = (Anew^(-1))';

max(max(abs(AnewInvT - AinvTp)))
