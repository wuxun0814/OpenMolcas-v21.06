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
* Copyright (C) 1995,2001,2008, Roland Lindh                           *
************************************************************************
      SubRoutine PCMHss(
#define _CALLING_
#include "hss_interface.fh"
     &                 )
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of nuclear attraction     *
*         integrals.                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, May 1995                                *
*                                                                      *
*             Modified to PCM gradients September 2001, Lund, by       *
*             R. Lindh.                                                *
*             Modified to PCM Hessian February 2008, Lund by           *
*             R. Lindh.                                                *
************************************************************************
      use PCM_arrays, only: PCM_SQ, PCMTess
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
      External TNAI1, Fake, XCff2D
#include "real.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
#include "disp2.fh"
#include "rctfld.fh"

#include "hss_interface.fh"

*     Local variables
      Integer iDCRT(0:7), Index(3,4)
      Real*8 Coori(3,4), CoorAC(3,2), C(3), TC(3)
      Logical NoLoop, JfGrd(0:2,0:3),
     &        JfHss(0:3,0:2,0:3,0:2),
     &        IfG(0:3), Tr(0:3)
      Integer iAnga(4), iStb(0:7), JndGrd(0:2,0:3,0:7),
     &        JndHss(0:3,0:2,0:3,0:2,0:7),
     &         mOp(4), iuvwx(4)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*                                                                      *
************************************************************************
*                                                                      *
*    We will have five terms here!
*
*    1) Sum(i) q_i V_ie^xy
*    2) Sum(ij) V_in^y Q_ij V_je^x
*    3) Sum(ij) V_ie^y Q_ij V_je^x
*    Maurizio to add comments for the last two terms!
*
      iRout = 151
      iPrint = nPrint(iRout)
*
      nip = 1
      ipA = nip
      nip = nip + nAlpha*nBeta
      ipB = nip
      nip = nip + nAlpha*nBeta
      ipDAO = nip
      nip = nip + nAlpha*nBeta*nElem(la)*nElem(lb)*nElem(nOrdOp)
      If (nip-1.gt.nZeta*nArr) Then
         Write (6,*) 'nip-1.gt.nZeta*nArr'
         Call ErrTra
         Call Abend()
      End If
      nArray = nZeta*nArr - nip + 1
*
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = nOrdOp
      iAnga(4) = 0
      call dcopy_(3, A,1,Coori(1,1),1)
      call dcopy_(3,RB,1,Coori(1,2),1)
      If (la.ge.lb) Then
       call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
       call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
      iuvwx(1) = dc(mdc)%nStab
      iuvwx(2) = dc(ndc)%nStab
      mOp(1) = nOp(1)
      mOp(2) = nOp(2)
*
      ipAOff = ipA
      Do iBeta = 1, nBeta
         call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
         ipAOff = ipAOff + nAlpha
      End Do
*
      ipBOff = ipB
      Do iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
         ipBOff = ipBOff + 1
      End Do
*
*---- Modify the density matrix with the prefactor
*
      nDAO = nElem(la) * nElem(lb)
      Do iDAO = 1, nDAO
         Do iZeta = 1, nZeta
            Fact = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
            DAO(iZeta,iDAO) = Fact * DAO(iZeta,iDAO)
         End Do
      End Do
      If (iPrint.ge.99) Call RecPrt('DAO',' ',DAO,nZeta,nDAO)
*
*---- Generate stabilizor of C, i.e. a center of a tile.
*
      nStb=1
      iStb(0)=0
*
*     Loop over the tiles
*
      Do iTs = 1, nTs
         q_i=PCM_SQ(1,iTs)+PCM_SQ(2,iTs)
         NoLoop = q_i.eq.Zero
         If (NoLoop) Go To 111
*------- Pick up the tile coordinates
         C(1:3)=PCMTess(1:3,iTs)

         If (iPrint.ge.99) Call RecPrt('C',' ',C,1,3)
         Call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
         Fact = -q_i*DBLE(nStabM) / DBLE(LmbdT)
*
         Call DYaX(nZeta*nDAO,Fact,DAO,1,Array(ipDAO),1)
*
         iuvwx(3) = nStb
         iuvwx(4) = nStb
*
         Do lDCRT = 0, nDCRT-1
            mOp(3) = NrOpr(iDCRT(lDCRT))
            mOp(4) = mOp(3)
            Call OA(iDCRT(lDCRT),C,TC)
            call dcopy_(3,TC,1,CoorAC(1,2),1)
            call dcopy_(3,TC,1,Coori(1,3),1)
            call dcopy_(3,TC,1,Coori(1,4),1)
*
*           Initialize JfGrd, JndGrd, JfHss, and JndHss.
*
            Call LCopy(12,[.False.],0,JfGrd,1)
            Call ICopy(nSym*4*3,[0],0,JndGrd,1)
            Call LCopy(144,[.False.],0,JfHss,1)
            Call ICopy(nSym*16*9,[0],0,JndHss,1)
*
*           Overwrite with information in IfGrd, IndGrd, IfHss,
*           and IndHss. This sets up the info for the first two
*           centers! Make sure that no translational invariance
*           is used.

            Do iAtom = 0, 1
               Do iCar  = 0, 2
                  JfGrd(iCar,iAtom) = Ifgrd(iCar,iAtom)
                  Do iIrrep=0,nSym-1
                     JndGrd(iCar,iAtom,iIrrep)=
     &                  Abs(IndGrd(iCar,iAtom,iIrrep))
                  End Do
                  Do jAtom = 0, 1
                     Do jCar = 0, 2
                        JfHss(iAtom,iCar,jAtom,jCar) =
     &                    IfHss(iAtom,iCar,jAtom,jCar)
                        Do iIrrep=0,nSym-1
                           JndHss(iAtom,iCar,jAtom,jCar,iIrrep) =
     &                       Abs(IndHss(iAtom,iCar,jAtom,jCar,iIrrep))
                        End Do
                     End Do
                  End Do
               End Do
            End Do
*
*-----------Derivatives with respect to the operator is computed via
*           the translational invariance.
*           Note: We want no such thing!
*
*
*           The third center is calculated by translational invariance.
*           This requires the 2nd derivatives on the other centers.
*           Note: We want no such thing!
*
            Call LCopy(4,[.False.],0,Tr,1)
*
            IfG(0)=.True.
            IfG(1)=.True.
            IfG(2)=.False.
            IfG(3)=.False.
            Call LCopy(12,[.False.],0,JfGrd,1)
*
*           Compute integrals with the Rys quadrature.
*
            nDiff=2
            mRys=(la+lb+2+nDiff+nOrdOp)/2
            Eta=One
            EInv=One
            nFinal=0
            Call Rysg2(iAnga,mRys,nZeta,
     &                 Array(ipA),Array(ipB),[One],[One],
     &                 Zeta,ZInv,nZeta,[Eta],[EInv],1,
     &                 P,nZeta,TC,1,Coori,Coori,CoorAC,
     &                 Array(nip),nArray,
     &                 TNAI1,Fake,XCff2D,
     &                 Array(ipDAO),nDAO*nElem(nOrdOp),
     &                 Hess,nHess,
     &                 JfGrd,JndGrd,
     &                 JfHss,JndHss,mOp,iuvwx,IfG,
     &                 nFinal,Index,.False.,.True.,Tr)
*
*           Call RecPrt(' In PCMHss:Hess',' ',Hess,nHess,1)
         End Do  ! End loop over DCRs
*
111      Continue
      End Do     ! End loop over centers in the external field
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(Final)
        Call Unused_integer(nHer)
        Call Unused_real_array(Ccoor)
        Call Unused_integer_array(lOper)
      End If
      End
