************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      Subroutine NucAtt(mGrid,Rho,nRho,P2_ontop,
     &                  nP2_ontop,iSpin,F_xc,dF_dRho,
     &                  ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_XX)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      use nq_Grid, only: Grid
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
#include "nq_info.fh"
      Real*8  Rho(nRho,mGrid), dF_dRho(ndF_dRho,mGrid),
     &        dF_dP2ontop(ndF_dP2ontop,mGrid),
     &        P2_ontop(nP2_ontop,mGrid), F_xc(mGrid)
      Real*8, Allocatable:: RA(:,:), ZA(:), Eff(:)
      Integer, Allocatable:: nStab(:)
*
      Call Get_nAtoms_All(mCenter)
      Call mma_allocate(RA,3,mCenter,Label='RA')
      Call Get_Coord_All(RA,mCenter)
*
      Call mma_allocate(ZA,mCenter,Label='ZA')
      Call Get_iScalar('Unique atoms',nCenter)
*
      Call mma_allocate(nStab,nCenter,Label='nStab')
      Call Get_iArray('nStab',nStab,nCenter)
      Call mma_allocate(Eff,nCenter,Label='Eff')
      Call Get_dArray('Effective Nuclear Charge',Eff,nCenter)
      Call Get_iScalar('nSym',nSym)
*
      iOff = 1
      Do i = 1, nCenter
         n=nSym/nStab(i)
         call dcopy_(n,[Eff(i)],0,ZA(iOff),1)
         iOff = iOff + n
      End Do
*
      Call mma_deallocate(Eff)
      Call mma_deallocate(nStab)
*
      Call Do_NucAtt_(mGrid,Rho,nRho,P2_ontop,nP2_ontop,
     &                iSpin,F_xc,dF_dRho,ndF_dRho,
     &                dF_dP2ontop,ndF_dP2ontop,Grid,
     &                RA,ZA,mCenter)
*
      Call mma_deallocate(ZA)
      Call mma_deallocate(RA)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(T_XX)
      End
      Subroutine Do_NucAtt_(mGrid,Rho,nRho,P2_ontop,nP2_ontop,
     &                      iSpin,F_xc,dF_dRho,ndF_dRho,
     &                      dF_dP2ontop,ndF_dP2ontop,Grid,RA,ZA,mCenter)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid), dF_dRho(ndF_dRho,mGrid),
     &        dF_dP2ontop(ndF_dP2ontop,mGrid),
     &        P2_ontop(nP2_ontop,mGrid), F_xc(mGrid)
      Real*8 Grid(3,mGrid),RA(3,mCenter),ZA(mCenter)
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=1
*
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
*
         d_alpha=Rho(1,iGrid)
         DTot=Two*d_alpha
*
*------- Accumulate contributions to the nuclear attraction energy
*
         Attr=Zero
         Do i = 1, mCenter
            x=Grid(1,iGrid)-RA(1,i)
            y=Grid(2,iGrid)-RA(2,i)
            z=Grid(3,iGrid)-RA(3,i)
            Fact=ZA(i)/Sqrt(x**2+y**2+z**2)
            Attr=Attr+Fact
         End Do
         F_xc(iGrid)=F_xc(iGrid)-Attr*DTot
*
         dF_dRho(ipR,iGrid)=-Attr
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=/=1
*
      Else
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
*
         d_alpha=Rho(1,iGrid)
         d_beta =Rho(2,iGrid)
         DTot=d_alpha+d_beta
*
*------- Accumulate contributions to the nuclear attraction energy
*
         Attr=Zero
         Do i = 1, mCenter
            x=Grid(1,iGrid)-RA(1,i)
            y=Grid(2,iGrid)-RA(2,i)
            z=Grid(3,iGrid)-RA(3,i)
            Fact=ZA(i)/Sqrt(x**2+y**2+z**2)
            Attr=Attr+Fact
         End Do
         F_xc(iGrid)=F_xc(iGrid)-Attr*DTot
*
         dF_dRho(ipRa,iGrid)=-Attr
         dF_dRho(ipRb,iGrid)=-Attr
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End
