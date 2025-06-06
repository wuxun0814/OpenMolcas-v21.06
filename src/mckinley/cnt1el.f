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
* Copyright (C) 1990, Roland Lindh                                     *
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Cnt1El(Kernel,KrnlMm,Label,
     &                 iDCnt,iDCar,loper,rHrmt,DiffOp,dens,
     &                 Lab_Dsk,iadd)
************************************************************************
*                                                                      *
* Object: to compute the one-electron integrals. The method employed at*
*         this point is not necessarily the fastest. However, the total*
*         time for the computation of integrals will depend on the time*
*         spent in computing the two-electron integrals.               *
*         The memory at this point is assumed to be large enough to do *
*         the computation in core.                                     *
*         The data is structured with respect to four indices, two (my *
*         ny or i j) refer to primitives or basis functions and two (a *
*         b) refer to the components of the cartesian or spherical     *
*         harmonic gaussians.                                          *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
*             Rewritten for gradients needed in hessian calculations   *
*             and general operator treatment                           *
*             May '95 By:                                              *
*             Anders Bernhardsson , Dept. of Theoretical Chemistry,    *
*             University  of Lund, SWEDEN.                             *
************************************************************************
      use Real_Spherical
      use iSD_data
      use Basis_Info
      use Center_Info
      use Sizes_of_Seward, only:S
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
*     External Kernel, KrnlMm
      External KrnlMm
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "disp2.fh"
#include "nsd.fh"
#include "setup.fh"
* log trans   integer dcent
      Real*8 A(3), B(3), RB(3),CCoor(3),dens(*)
      Character Label*8
      Integer nOp(2), ip(8), iDCRR(0:7), iDCRT(0:7),
     &        iStabM(0:7), iStabO(0:7), IndGrd(0:7)
      Logical IfGrd(3,2),EQ,DiffOP,DiffCnt,Trans(2)
      Integer, Parameter:: iTwoj(0:7)=[1,2,4,8,16,32,64,128]
      Character(LEN=8) Lab_dsk
      Real*8, Allocatable:: Integrals(:), Zeta(:), ZI(:), Kappa(:),
     &                      PCoor(:,:), Fnl(:), Kern(:), ScrSph(:),
     &                      SO(:), Scr(:)
      Logical, External :: TF
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      Subroutine Kernel(
#define _CALLING_
#include "grd_mck_interface.fh"
     &                 )
#include "grd_mck_interface.fh"
      End Subroutine Kernel
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the number of blocks from each component of the operator
*     and the irreps it will span.
*
*
C Differentiated symmetry-unique center IDCNT
C Derivative wrt component IDCAR=1,2,3 (d/dx,d/dy,d/dz)
C INDDSP(IDCNT,IIRREP) is the number of displacements in
C earlier center/irrep. Thus it is an offset.
*
      nOrdOp=0
      Call iCopy(nIrrep,[0],0,IndGrd,1)
      loper=0
#ifdef _DEBUGPRINT_
      iprint=99
#endif
      nnIrrep=nIrrep
      If (sIrrep) nnIrrep=1
      Do iIrrep=0,nnIrrep-1
         nDisp = IndDsp(iDcnt,iIrrep)
C First set NDISP=ordering number of this displacement.
C Then loop over directions d/dx,d/dy,d/dz
         Do iCar=1,3
            iComp = 2**(iCar-1)
            If ( TF(iDCnt,iIrrep,iComp)) Then
               ndisp=ndisp+1
C NDISP is now the ordering number of this displacement.
               If (iDCar.eq.icar) Then
                  loper=loper+2**iIrrep
                  IndGrd(iIrrep) = nDisp
               End If
            End If
         End Do
      End Do
      nIC=0
      If (loper.eq.0) Return
C For the displacement represented by this symmetry-unique
C center IDCNT and this component IDCAR, the differentiation
C operator has components with irreps that have been marked
C with '1' in LOPER, regarded as a flag array.
C INDGRD(IIRREP) will be zero, except for those irreps, and
C will then contain the ordering number of the displacement.

C Allocate one integral array for each of these irreps.
C The address is kept in array IP().
      nIC=0
      Call ICopy(nIrrep,[0],0,ip,1)

      iStart=1
      Do iIrrep =0,nIrrep-1
         If (iAnd(2**iIrrep,loper).ne.0) Then
            LenInt=nFck(iIrrep)
            nIc=nIC+1
            ip(NIC)=iStart
            iStart=iStart+LenInt
         End If
      End Do
      LenInt_Tot=iStart-1
      Call mma_allocate(Integrals,LenInt_Tot,Label='Integrals')
      Integrals(:)=Zero

C Obtain ISTABO, the stabilizer of the totally symmetric irrep(!)
C Note: 3rd parameter is bit-packed set of irreps
C so '1' contains only irrep nr 0.
C But then ISTABO will be the whole group!? and NSTABO=NIRREP?!
      Call SOS(iStabO,nStabO,1)
*
*-----Auxiliary memory allocation.
*
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Basis_Mode('Valence')
      Call Nr_Shells(nSkal)
      Call Setup_iSD()
*                                                                      *
************************************************************************
*                                                                      *
*-----Double loop over shells.
*
      Do iS = 1, nSkal
         iShll  = iSD( 0,iS)
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         iCnttp = iSD(13,iS)
         iCnt   = iSD(14,iS)
         A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
*
         Do jS = 1, iS
            jShll  = iSD( 0,jS)
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jPrim  = iSD( 5,jS)
            jAO    = iSD( 7,jS)
            mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
            jCnttp = iSD(13,jS)
            jCnt   = iSD(14,jS)
            B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
*-------Call kernel routine to get memory requirement. Observe, however
*       that kernels which will use the HRR will allocate that
*       memory internally.
*
        maxi=S%maxPrm(iAng)*S%maxprm(jang)
        Call mma_allocate(Zeta,maxi,Label='Zeta')
        Call mma_allocate(ZI,maxi,Label='ZI')
        Call mma_allocate(Kappa,maxi,Label='Kappa')
        Call mma_allocate(PCoor,maxi,3,Label='PCoor')
        Call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)
*
*       Memory requirements for contraction and Symmetry
*       adoption of derivatives.
*
        lFinal = S%MaxPrm(iAng) * S%MaxPrm(jAng) *
     &           nElem(iAng)*nElem(jAng)*nIrrep
*
        MemKrn=Max(MemKer*Maxi,lFinal)
        Call mma_Allocate(Kern,MemKrn,Label='Kern')
*
*            Save some memory and use Scrt area for
*            transformation
*
*       Allocate memory for the final integrals all in the
*       primitive basis.
*
        Call mma_allocate(Fnl,lFinal,Label='Fnl')
*
*       Scratch area for the transformation to spherical gaussians
*
        nScr1=S%MaxBas(iAng)*S%MaxBas(jAng)*nElem(iAng)*nElem(jAng)*nIC
        Call mma_allocate(ScrSph,nScr1,Label='ScrSph')
*
*         At this point we can compute Zeta.
*         This is now computed in the ij or ji order.
*
          Call ZXia(Zeta,ZI,
     &              iPrim,jPrim,Shells(iShll)%Exp,
     &                          Shells(jShll)%Exp)
*
*
            DiffCnt=((mdci.eq.iDCnt).or.(mdcj.eq.iDCnt))
            If ((.not.DiffCnt).and.(.not.DiffOp)) Goto 131
            Call lCopy(6,[.false.],0,IfGrd,1)
C Logical trans(2)
C trans(iCnt) is true means there will be a sign shift in the SYMADO
C routine for the contribution to the integral from the
C differentiation wrt center iCnt
            Call lCopy(2,[.false.],0,trans,1)
            If (mdci.eq.iDCnt) Then
                IfGrd(idCar,1)=.true.
            End If
            If (mdcj.eq.iDCnt) Then
                IfGrd(idCar,2)=.true.
            End If
*
            If (IfGrd(iDCar,1).and.IfGrd(iDCar,2).and.
     &          (.not.DiffOp)) Then
              IfGrd(iDCar,2)=.false.
              Trans(2)=.true.
            End If
            If (Label.eq.'CONNECTI') Trans(2)=.false.
            If (Label.eq.'OVRGRDA') Trans(2)=.false.
*
*           Allocate memory for SO integrals that will be generated by
*           this batch of AO integrals.
*
            nSO=0
            Do iIrrep=0,nIrrep-1
                If (iAnd(loper,2**iIrrep).ne.0) Then
                 iSmLbl=2**iIrrep
                 nSO=nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
               End If
            End Do
#ifdef _DEBUGPRINT_
            If (iPrint.ge.29) Write (6,*) ' nSO=',nSO
#endif
            If (nSO.eq.0) Go To 131
            Call mma_Allocate(SO,nSO*iBas*jBas,Label='SO')
            SO(:)=Zero
*
*           Find the DCR for A and B
*
            Call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,
     &                     dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
*
*           Find the stabilizer for A and B
*
            Call Inter(dc(mdci)%iStab,dc(mdci)%nStab,
     &                 dc(mdcj)%iStab,dc(mdcj)%nStab,
     &                 iStabM,nStabM)
*
            Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
*
*           Compute normalization factor
*
            iuv = dc(mdci)%nStab*dc(mdcj)%nStab
            Fact = DBLE(iuv*nStabO) / DBLE(nIrrep**2 * LmbdT)
            If (MolWgh.eq.1) Then
               Fact = Fact * DBLE(nIrrep)**2 / DBLE(iuv)
            Else If (MolWgh.eq.2) Then
               Fact = sqrt(DBLE(iuv))*DBLE(nStabO)/DBLE(nIrrep*LmbdT)
            End If
*
*           Loops over symmetry operations acting on the basis.
*
            nOp(1) = NrOpr(0)
            if(jBas.lt.-999999) write(6,*) 'gcc overoptimization',nDCRR
            Do 140 lDCRR = 0, nDCRR-1
             Call OA(iDCRR(lDCRR),B,RB)
             nOp(2) = NrOpr(iDCRR(lDCRR))
             If (Label.ne.'CONNECTI'
     &           .and.EQ(A,RB).and. (.Not.DiffOp)) Go To 140
*
*            Compute kappa and P.
*
             Call Setup1(Shells(iShll)%Exp,iPrim,
     &                   Shells(jShll)%Exp,jPrim,
     &                   A,RB,Kappa,PCoor,ZI)
*
*            Compute AO integrals.
*            for easy implementation of NA integrals.
*
             call dcopy_(lFinal,[0.0d0],0,Fnl,1)
             Call Kernel(Shells(iShll)%Exp,iPrim,
     &                   Shells(jShll)%Exp,jPrim,
     &                   Zeta,ZI,
     &                   Kappa,PCoor,
     &                   Fnl,iPrim*jPrim,
     &                   iAng,jAng,A,RB,nOrder,Kern,
     &                   MemKrn,Ccoor,nOrdOp,IfGrd,IndGrd,nop,
     &                   loper,dc(mdci)%nStab,
     &                   dc(mdcj)%nStab,nic,idcar,idcnt,
     &                   iStabM,nStabM,trans,nIrrep)
*
*
*        Transform from primitive to contracted basis functions.
*        Order of transformation is fixed. It has been shown through
*        testing that the index order ij,ab will give a performance
*        that is up to 20% faster than the ab,ij index order.
*
*
*            Transform i,jabx to jabx,I
             kk=nElem(iAng)*nElem(jAng)*nIC
             Call DGEMM_('T','N',
     &                   jPrim*kk,iBas,iPrim,
     &                   1.0d0,Fnl,iPrim,
     &                         Shells(iShll)%pCff,iPrim,
     &                   0.0d0,Kern,jPrim*kk)
*
*            Transform j,abxI to abxI,J
*
             Call DGEMM_('T','N',
     &                   kk*iBas,jBas,jPrim,
     &                   1.0d0,Kern,jPrim,
     &                         Shells(jShll)%pCff,jPrim,
     &                   0.0d0,Fnl,kk*iBas)
*
*            Transform to spherical gaussians if needed.
*
                 kk=nElem(iAng)*nElem(jAng)
*
                 If (Shells(iShll)%Transf.or.Shells(jShll)%Transf) Then
*
*             Result comes back as IJAB or IJAb
*
                   Call CarSph(Fnl,kk,iBas*jBas*nIC,
     &                         Kern,nScr1,
     &                    RSph(ipSph(iAng)),iAng,
     &                    Shells(iShll)%Transf,
     &                    Shells(iShll)%Prjct,
     &                    RSph(ipSph(jAng)),jAng,
     &                    Shells(jShll)%Transf,
     &                    Shells(jShll)%Prjct,
     &                    ScrSph,iCmp*jCmp)
*
                  Call DGeTmO(ScrSph,nIC,nIC,
     &                    iBas*jBas*iCmp*jCmp,
     &                    Kern,iBas*jBas*iCmp*jCmp)
*
                Else
*
*             Transpose abx,IJ back to IJ,abx
*
                    Call DGeTmO(Fnl,kk*nIC,kk*nIC,
     &                   iBas*jBas,Kern,iBas*jBas)
                End If
*
*            At this point accumulate the batch of integrals onto the
*            final symmetry adapted integrals.
*
#ifdef _DEBUGPRINT_
                If (iPrint.ge.99) Then
                  Call RecPrt (' Accumulated SO integrals, so far...',
     &                               ' ',SO,iBas*jBas,nSO)
                End If
#endif
*
*------------Symmetry adapt component by component
*
             iSOBlk = 1
             iIC=1
             Do iIrrep = 0, nIrrep-1
                iSmLbl=iAnd(lOper,iTwoj(iIrrep))
                mSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
                If (mSO.eq.0) Then
                   Do jIrrep = 0, nIrrep-1
                      If (iAnd(iSmLbl,iTwoj(jIrrep)).ne.0) iIC = iIC + 1
                   End Do
                Else
                   Call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,
     &                         iShell,jShell,iShll,jShll,iAO,jAO,
     &                         Kern,iBas,jBas,nIC,iIC,
     &                         SO(iSOBlk),mSO,nOp)
                   iSOBlk = iSOBlk + mSO*iBas*jBas
                End If
             End Do
*
 140        Continue
*
*           Multiply with factors due to projection operators
*
            If (Fact.ne.One)
     &       Call DScal_(nSO*iBas*jBas,Fact,SO,1)
*
*           Scatter the SO's on to the non-zero blocks of the
*           lower triangle.
*
             iSOBlk=1
             iiC=0
             Do  iIrrep = 0, nIrrep-1
               If (iAnd(lOper,2**iIrrep).ne.0) Then
                 iSmlbl=2**iIrrep
                 iiC=iiC+1
                 mSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
                 If (nfck(iIrrep).ne.0.and.mSO.ne.0)
     &            Call SOSctt(SO(iSOBlk),iBas,jBas,mSO,
     &                    Integrals(ip(iIC)),nFck(iIrrep),iSmLbl,
     &                    iCmp,jCmp,iShell,jShell,
     &                    iAO,jAO,
     &                    nIC,Label,2**iIrrep,rHrmt)
                 iSOBlk = iSOBlk + mSO*iBas*jBas
               End If
             End Do
*
            Call mma_deallocate(SO)
 131        Continue
         Call mma_deallocate(ScrSph)
         Call mma_deallocate(Kern)
         Call mma_deallocate(Fnl)
         Call mma_deallocate(PCoor)
         Call mma_deallocate(Kappa)
         Call mma_deallocate(ZI)
         Call mma_deallocate(Zeta)
         End Do
      End Do
      Call Free_iSD()
*
*     Compute properties or write integrals to disc and
*     deallocate core.
*
      nDens=0
      nDenssq=0
      Do iI=0,nIrrep-1
         nDenssq = nDenssq + nBas(ii)**2+nBas(ii)
         nDens   = nDens   + nBas(iI)*(nBas(iI)+1)/2
      End Do
      nrOp=0

      Call mma_allocate(Scr,2*nDenssq,Label='Scr')
      Do 16 iIrrep = 0, nIrrep-1
         iSmLbl = 2**iIrrep
         If (iAnd(2**iIrrep,loper).ne.0) Then
            nrOp=nrOp+1
            jdisp=indgrd(iIrrep)
            kOper=2**iIrrep
            If (show.and.iIrrep.eq.0) Then
               Write(6,*) Label,': ',
     &               ddot_(nDens,Dens,1,Integrals(ip(nrop)),1)
               Write(6,*) 'oper: ',
     &               ddot_(nDens,Integrals(ip(nrop)),1,
     &                           Integrals(ip(nrop)),1)
               Write(6,*) 'Dens: ',ddot_(nDens,Dens,1,Dens,1)
            Else If (show) Then
               mDens=nFck(iIrrep)
               Write(6,*) Label
               Write(6,'(A,G20.10)') 'oper: ',
     &               ddot_(mDens,Integrals(ip(nrop)),1,
     &                           Integrals(ip(nrop)),1)
            End if
*
            If (iadd.ne.0) Then
               irc=-1
               iopt=0
               call dRdMck(irc,iOpt,Lab_dsk,jdisp,Scr,koper)
               If (irc.ne.0) Call SysAbendMsg('cnt1el',
     &                            'error during read in rdmck',' ')
               Call DaXpY_(nfck(iIrrep),one,Scr,1,
     &                     Integrals(ip(nrop)),1)
            End If
            irc=-1
            iopt=0
#ifdef _DEBUGPRINT_
            Write(6,'(2A,2I8)')'Lab_dsk,jdisp,koper',Lab_dsk,jdisp,koper
#endif
            Call dWrMck(irc,iOpt,Lab_dsk,jdisp,Integrals(ip(nrop)),
     &                  koper)
            If (irc.ne.0) Call SysAbendMsg('cnt1el',
     &                            'error during write in dwrmck',' ')
         End If
 16   Continue
      Call mma_deallocate(Scr)
      Call mma_deallocate(Integrals)
*
      Return
      End
