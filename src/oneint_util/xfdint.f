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
* Copyright (C) 1995, Roland Lindh                                     *
************************************************************************
      SubRoutine XFdInt(
#define _CALLING_
#include "int_interface.fh"
     &                 )
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of nuclear attraction     *
*         integrals.                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, April '95                               *
************************************************************************
      use external_centers
      use Phase_Info
      Implicit Real*8 (A-H,O-Z)
      External TNAI, Fake, XCff2D, XRys2D
#include "itmax.fh"
#include "real.fh"
#include "print.fh"

#include "int_interface.fh"

*-----Local variables
      Real*8 C(3), TC(3), Coori(3,4), CoorAC(3,2),
     &       ZFd((iTabMx+1)*(iTabMx+2)/2), ZRFd((iTabMx+1)*(iTabMx+2)/2)
      Logical EQ, NoLoop, NoSpecial
      Integer iAnga(4), iDCRT(0:7), iStb(0:7), jCoSet(8,8)
      Character ChOper(0:7)*3
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
C     nElem(ixyz) = 2*ixyz+1
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      iRout = 151
      iPrint = nPrint(iRout)
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*
*---- Loop over charges and dipole moments in the external field
*
      nData=3
      Do iOrdOp = 0, nOrdOp
*
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = iOrdOp
      iAnga(4) = 0
      call dcopy_(3, A,1,Coori(1,1),1)
      call dcopy_(3,RB,1,Coori(1,2),1)
      mabMin = nabSz(Max(la,lb)-1)+1
      mabMax = nabSz(la+lb)
      If (EQ(A,RB)) mabMin=nabSz(la+lb-1)+1
      mcdMin = nabSz(iOrdOp-1)+1
      mcdMax = nabSz(iOrdOp)
      lab=(mabMax-mabMin+1)
      kab=nElem(la)*nElem(lb)
      lcd=(mcdMax-mcdMin+1)
      labcd=lab*lcd
*
*     Compute FLOP's and size of work array which Hrr will use.
*
      Call mHrr(la,lb,nFLOP,nMem)
*
*---- Distribute the work array
*
      ip2 = 1
      ip1 = ip2 + nZeta*Max(labcd,lcd*nMem)
      mArr = nArr - Max(labcd,lcd*nMem)
*
*     Find center to accumulate angular momentum on. (HRR)
*
      If (la.ge.lb) Then
       call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
       call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
*
*     Loop over centers of the external field.
*
      iDum=0
      Do iFd = 1, nXF
*
         NoLoop=.True.
         Do jElem = 1, nElem(iOrdOp)
            ZFd(jElem)=XF(nData+jElem,iFd)
*        Divide quadrupole diagonal by 2 due to different normalisation
            if((iOrdOp.eq.2).and.
     &           (jElem.eq.1.or.jElem.eq.4.or.jElem.eq.6))
     &           ZFd(jElem)=ZFd(jElem)*0.5D0
            NoLoop = NoLoop .and.  ZFd(jElem).eq.Zero
         End Do
*
         If (NoLoop) Go To 111
*------- Pick up the center coordinates
         C(1:3)=XF(1:3,iFd)

         If (iPrint.ge.99) Call RecPrt('C',' ',C,1,3)
*
*------- Generate stabilizor of C
*
         iChxyz=iChAtm(C)
         Call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)
*
*--------Find the DCR for M and S
*
         Call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
         Fact = DBLE(nStabM) / DBLE(LmbdT)
*
         If (iPrint.ge.99) Then
            Write (6,*) ' m      =',nStabM
            Write (6,'(9A)') '(M)=',(ChOper(iStabM(ii)),
     &            ii = 0, nStabM-1)
            Write (6,*) ' s      =',nStb
            Write (6,'(9A)') '(S)=',(ChOper(iStb(ii)),
     &            ii = 0, nStb-1)
            Write (6,*) ' LambdaT=',LmbdT
            Write (6,*) ' t      =',nDCRT
            Write (6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),
     &            ii = 0, nDCRT-1)
         End If

*
         Do lDCRT = 0, nDCRT-1
            Call OA(iDCRT(lDCRT),C,TC)
*
            jElem=0
            Do ix = iOrdOp, 0, -1
               If (Mod(ix,2).eq.0) Then
                  Factx=One
               Else
                  Factx=DBLE(iPhase(1,iDCRT(lDCRT)))
               End If
               Do iy = iOrdOp-ix, 0, -1
                  If (Mod(iy,2).eq.0) Then
                     Facty=One
                  Else
                     Facty=DBLE(iPhase(2,iDCRT(lDCRT)))
                  End If
                  iz = iOrdOp-ix-iy
                  If (Mod(iz,2).eq.0) Then
                     Factz=One
                  Else
                     Factz=DBLE(iPhase(3,iDCRT(lDCRT)))
                  End If
*
                  jElem = jElem + 1
                  ZRFd(jElem)=Factx*Facty*Factz*ZFd(jElem)
               End Do
            End Do
*
            call dcopy_(3,TC,1,CoorAC(1,2),1)
            call dcopy_(3,TC,1,Coori(1,3),1)
            call dcopy_(3,TC,1,Coori(1,4),1)
*
*           Compute integrals with the Rys quadrature.
*
            nT = nZeta
            NoSpecial=.True.
            Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &               [One],[One],1,P,nZeta,
     &               TC,1,rKappa,[One],Coori,Coori,CoorAC,
     &               mabmin,mabmax,mcdMin,mcdMax,
     &               Array(ip1),mArr*nZeta,
     &               TNAI,Fake,XCff2D,XRys2D,NoSpecial)
*
*---------- The integrals are now ordered as ijkl,e,f
*
*           a) Change the order to f,ijkl,e
*           b) Unfold e to ab, f,ijkl,ab
*           c) Change the order back to ijkl,ab,f
*
*a)--------
*
            Call DGeTMO(Array(ip1),nZeta*lab,nZeta*lab,lcd,
     &                  Array(ip2),lcd)
*
*b)---------Use the HRR to unfold e to ab
*
            Call HRR(la,lb,A,RB,Array(ip2),lcd*nZeta,nMem,ipIn)
            ip3=ip2-1+ipIn
*
*c)--------
*
            Call DGeTMO(Array(ip3),lcd,lcd,nZeta*kab,Array(ip1),
     &                  nZeta*kab)
*
*-----------Accumulate contributions to the symmetry adapted operator
*
            nOp = NrOpr(iDCRT(lDCRT))
            ipI=ip1
*
            Do i = 1, nElem(iOrdOp)
               If (ZRFd(i).ne.Zero)
     &            Call SymAdO(Array(ipI),nZeta,la,lb,nComp,Final,nIC,
     &                        nOp         ,lOper,iChO,-Fact*ZRFd(i))
               ipI=ipI+nZeta*nElem(la)*nElem(lb)
            End Do
*
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
            Write (6,*) (Fact*ZFd(i),i = 1, nElem(iOrdOp))
            Call RecPrt('Array(ip1)',' ',Array(ip1),nZeta,
     &              (la+1)*(la+2)/2*(lb+1)*(lb+2)/2*nElem(iOrdOp))
            Call RecPrt('Final',' ',Final,
     &              nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2*nIC)
#endif

*
         End Do  ! End loop over DCRs
*
111      Continue
      End Do     ! iFd
*
          nData = nData + nElem(iOrdOp)
      End Do     !  iOrdOp
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_integer(nHer)
         Call Unused_real_array(CCoor)
         Call Unused_real_array(PtChrg)
         Call Unused_integer(iAddPot)
      End If
      End
