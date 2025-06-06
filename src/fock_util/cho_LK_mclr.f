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
* Copyright (C) Mickael G. Delcey                                      *
************************************************************************
      SUBROUTINE CHO_LK_MCLR(DLT,DI,DA,G2,Kappa,JI,
     &                      KI,JA,KA,FkI,FkA,
     &                      MO_Int,QVec,Ash,CMO,CMO_inv,
     &                      nOrb,nAsh,nIsh,doAct,Fake_CMO2,
     &                      LuAChoVec,LuIChoVec,iAChoVec)

**********************************************************************
*  Author : M. G. Delcey based on cho_LK_rassi_x
*
*  Note:  this routine differs from CHO_LK_RASSI_X because it can
*         handle inactive or active matrix and 1-index transformed
*         densities
*
*
C *************** INACTIVE AO-BASIS FOCK MATRIX **********************
C
C   FI(ab) = 2 * sum_J  Lab,J * U(J)  -  sum_Jk  Yka,J * Xkb,J
C
C      U(J) = sum_gd  Lgd,J * DI(gd)
C
C      a,b,g,d:  AO-index
C      k:        MO-index   belonging to (Inactive)
C      v,w,x,y:  MO-indeces belonging to (Active)
C
**********************************************************************
      use ChoArr, only: nBasSh, nDimRS
      use ChoSwp, only: nnBstRSh, InfVec, IndRed
      use Data_Structures, only: DSBA_Type, Allocate_DSBA,
     &                           Deallocate_DSBA
      use Data_Structures, only: SBA_Type
      use Data_Structures, only: Allocate_SBA, Deallocate_SBA
      use Data_Structures, only: NDSBA_Type, Allocate_NDSBA,
     &                           Deallocate_NDSBA
      use Data_Structures, only: G2_Type, Allocate_G2, Deallocate_G2
      use Data_Structures, only: Allocate_L_Full, Deallocate_L_Full,
     &                           L_Full_Type
      use Data_Structures, only: Allocate_Lab, Deallocate_Lab,
     &                           Lab_Type

#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit Real*8 (a-h,o-z)
#include "warnings.fh"

      Real*8 G2(*), MO_Int(*)
      Integer   kOff(8),kaOff(8)
      Integer   iASQ(8,8,8)
      Integer   LuAChoVec(8),LuIChoVec(8)
      Real*8    tread(2),tcoul(2),texch(2),tintg(2),tact(2)
      Real*8    tint1(2),tint2(2),tint3(2),tQmat(2)
      Real*8    tmotr(2),tscrn(2)
      Integer   nChMo(8)

      Type (DSBA_Type) DLT, DI, DA, Kappa, JI, KI, JA, KA, FkI, FkA,
     &                 QVec, Ash(2), CMO, CMO_Inv, Tmp(2), QTmp(2),
     &                 CM(2)
      Type (DSBA_Type) JALT
      Type (SBA_Type) Lpq(3)
      Type (NDSBA_Type) DiaH
      Type (G2_Type) MOScr
      Type (L_Full_Type) L_Full
      Type (Lab_Type) Lab

#ifdef _DEBUGPRINT_
      Logical   Debug
#endif
      Logical   timings,DoScreen,add
      Real*8    thrv(2),xtau(2),norm
      Character*50 CFmt
      Character(LEN=14), Parameter :: SECNAM = 'CHO_LK_MCLR'
#include "chomclr.fh"
#include "real.fh"
      Logical Fake_CMO2,DoAct,ReadInter
      Save nVec_
*
#include "cholesky.fh"
      Integer nOrb(8),nAsh(8),nIsh(8)

      Logical, Parameter ::  DoRead = .false.
      Real*8, Parameter :: FactCI = -Two, FactXI = Half, xone=-One
#include "choorb.fh"
#include "stdalloc.fh"

      Real*8 LKThr

      Character*6 mode
      Integer, External :: Cho_LK_MaxVecPerBatch
      Real*8,  External :: Cho_LK_ScreeningThreshold

      Real*8, Allocatable :: Lrs(:,:), VJ(:), Drs(:), Frs(:,:)

      Integer, Allocatable:: nnBfShp(:,:), kOffSh(:,:),
     &                       iShp_rs(:), Indx(:,:,:)
      Real*8, Allocatable:: SvShp(:,:),Diag(:),AbsC(:), SumAClk(:,:,:),
     &                      Ylk(:,:,:), MLk(:,:,:), Faa(:), Fia(:)
#if defined (_MOLCAS_MPP_)
      Real*8, Allocatable:: DiagJ(:)
#endif

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
************************************************************************

#ifdef _DEBUGPRINT_
      Debug=.false.! to avoid double printing in CASSCF-debug
#endif

      ! Allow LT-format access to JA although it is in SQ-format
      Call Allocate_DSBA(JALT,nBas,nBas,nSym,Case='TRI',Ref=JA%A0)
      timings=.false.
*
*
      IREDC = -1  ! unknown reduced set in core

      nDen = 2  ! the two bi-orthonormal sets of orbitals
      If (Fake_CMO2) nDen = 1  ! MO1 = MO2
      kDen=nDen

      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      ! 1 --> CPU   2 --> Wall
      tread(:) = zero  !time read/transform vectors
      tcoul(:) = zero  !time for computing Coulomb
      texch(:) = zero  !time for computing Exchange
      tintg(:) = zero  !time for computing (tw|xy) integrals
      tmotr(:) = zero  !time for the half-transf of vectors
      tscrn(:) = zero  !time for screening overhead
      tint1(:) = zero
      tint2(:) = zero
      tint3(:) = zero
      tQmat(:) = zero
      tact(:)  = zero

C ==================================================================

c --- Various offsets
c --------------------
      nnO=0
      nnA=0
      nA2=0
      kOff(1)=0
      kAOff(1)=0
      MaxB=nBas(1)
      MaxAct=nAsh(1)
      DO ISYM=2,NSYM
        MaxB=Max(MaxB,nBas(iSym))
        MaxAct=Max(MaxAct,nAsh(iSym))
        nnO = nnO + nOrb(iSym-1)
        nnA = nnA + nAsh(iSym-1)
        nA2 = nA2 + nAsh(iSym-1)**2
        kOff(iSym)=nnO
        kAOff(iSym)=nnA
      END DO
      nnO = nnO + nOrb(nSym)
      nnA = nnA + nAsh(nSym)
      nA2 = nA2 + nAsh(nSym)**2

      nChMO(:) = nOrb(:)

      If (DoAct) Then
        ioff=0
        Do ijsym=1,nsym
          Do isym=1,nsym
            jsym=iEOR(isym-1,ijsym-1)+1
            Do kSym=1,nSym
               lSym=iEOR(kSym-1,ijSym-1)+1
               iASQ(isym,jsym,kSym)=ioff
               ioff=ioff+nASh(iSym)*nAsh(jSym)*
     &                   nASh(kSym)*nASh(lSym)
            End Do
          End Do
        End Do

C *** memory for the Q matrices --- temporary array
        Call Allocate_DSBA(QTmp(1),nBas,nAsh,nSym)
        Call Allocate_DSBA(QTmp(2),nBas,nAsh,nSym)
        QTmp(1)%A0(:)=Zero
        QTmp(2)%A0(:)=Zero

        iCase=1
        Call Allocate_G2(MOScr,nAsh,nSym,iCase)
        MOScr%A0(:)=Zero
      End If
**************************************************
      If (Deco) Then
         Call Allocate_DSBA(CM(1),nBas,nBas,nSym)
         Call Allocate_DSBA(CM(2),nBas,nBas,nSym)

         Call Allocate_DSBA(Tmp(1),nBas,nBas,nSym)
         Call Allocate_DSBA(Tmp(2),nBas,nBas,nSym)

         Do iS=1,nSym
*
**       Create Cholesky orbitals from DI
*
           Call CD_InCore(DI%SB(iS)%A2,nBas(iS),CM(1)%SB(iS)%A2,
     &                    nBas(iS),nChMO(iS),1.0d-12,irc)
           If (.not.Fake_CMO2) Then
*
**         MO transform
*

             Call DGEMM_('T','T',nChMO(iS),nBas(iS),nBas(iS),
     &                  1.0D0,CM(1)%SB(iS)%A2,nBas(iS),
     &                        CMO_inv%SB(iS)%A2,nBas(iS),
     &                  0.0d0,Tmp(2)%SB(iS)%A2,nChMO(iS))
*
**       Create one-index transformed Cholesky orbitals
*
             Call DGEMM_('N','N',nChMO(iS),nBas(iS),nBas(iS),
     &                  1.0D0,Tmp(2)%SB(iS)%A2,nChMO(iS),
     &                        Kappa%SB(is)%A2,nBas(iS),
     &                  0.0d0,Tmp(1)%SB(iS)%A2,nChMO(iS))
*
**         AO transform
*
             Call DGEMM_('N','T',nBas(iS),nChMO(iS),nBas(iS),
     &                    1.0D0,CMO%SB(iS)%A2,nBas(iS),
     &                          Tmp(1)%SB(iS)%A2,nChMO(iS),
     &                    0.0d0,CM(2)%SB(iS)%A2,nBas(iS))
           EndIf
         End Do
         Call Deallocate_DSBA(Tmp(2))
         Call Deallocate_DSBA(Tmp(1))
      Else

        Write (6,*) 'Cho_LK_MCLR: this will not work'
        Call Abend()

      EndIf

**************************************************


C --- Define the max number of vectors to be treated in core at once

      MaxVecPerBatch=Cho_LK_MaxVecPerBatch()

C --- Define the screening threshold

      LKThr=Cho_LK_ScreeningThreshold(-1.0d0)
      dmpk=1.0d-2
*     dmpk=0.0d0

C --- Vector MO transformation screening thresholds
      NumVT=NumChT
      Call GAIGOP_SCAL(NumVT,'+')
      thrv(1) = (LKThr/(Max(1,nnO)*NumVT))*dmpk**2
      xtau(1) = Sqrt((LKThr/Max(1,nnO))*dmpk)
      xtau(2) = xtau(1) ! dummy init

      If (.not.Fake_CMO2) Then
        norm=sqrt(ddot_(SIZE(Kappa%A0),Kappa%A0,1,Kappa%A0,1))
        xtau(2)=Sqrt((LKThr/Max(1,nnO))*dmpk)*norm
        dmpk=min(norm,1.0d-2)
        thrv(2)=(LKThr/(Max(1,nnO)*NumVT))*dmpk**2
      EndIf
      tau = (LKThr/Max(1,nnO))*dmpk

      MaxRedT=MaxRed
      Call GAIGOP_SCAL(MaxRedT,'+')

      If (Estimate) Then
         xtau(1)=xtau(1)/Sqrt(1.0d0*MaxRedT)
         If (.not.Fake_CMO2) xtau(2)=xtau(2)/Sqrt(1.0d0*MaxRedT)
         tau    =tau    /MaxRedT
      EndIf

      CALL mma_allocate(DIAG,NNBSTRT(1),Label='DIAG')

#if defined (_MOLCAS_MPP_)
      If (nProcs.gt.1 .and. Update .and. Is_Real_Par()) Then
         NNBSTMX=0
         Do i=1,nSym
            NNBSTMX = Max(NNBSTMX,NNBSTR(i,1))
         End Do
         CALL mma_allocate(diagJ,NNBSTMX,Label='diagJ')
         diagJ(:)=Zero
      EndIf
#endif
C *************** Read the diagonal integrals (stored as 1st red set)
      If (Update) CALL CHO_IODIAG(DIAG,2) ! 2 means "read"

c --- allocate memory for sqrt(D(a,b)) stored in full (squared) dim
      Call Allocate_NDSBA(DiaH,nBas,nBas,nSym)
      DiaH%A0(:)=Zero

c --- allocate memory for the abs(C(l)[k])
      Call mma_allocate(AbsC,MaxB,Label='AbsC')

c --- allocate memory for the Y(l)[k] vectors
      Call mma_allocate(Ylk,MaxB,nnO,nDen,Label='Ylk')

c --- allocate memory for the ML[k] lists of largest elements
c --- in significant shells
      Call mma_allocate(MLk,nShell,nnO,nDen,Label='MLk')

c --- allocate memory for the lists of  S:= sum_l abs(C(l)[k])
c --- for each shell
      Call mma_allocate(SumAClk,nShell,nnO,nDen)

c --- allocate memory for the Index arrays
      Call mma_allocate(Indx,[0,nShell],[1,nnO],[1,nDen],Label='Indx')

c --- allocate memory for kOffSh
      Call mma_allocate(kOffSh,nShell,nSym,Label='kOffh')

c --- allocate memory for nnBfShp
      nnShl_2=nShell**2
      Call mma_allocate(nnBfShp,nnShl_2,nSym,Label='nnBfShp')

c --- allocate memory for iShp_rs
      Call mma_allocate(iShp_rs,nnShl_tot,Label='iShp_rs')

c --- allocate memory for the shell-pair Frobenius norm of the vectors
      Call mma_allocate(SvShp,nnShl,2,Label='SvShp')


C *** Compute Shell Offsets ( MOs and transformed vectors)

      MxBasSh = 0

      Do iSyma=1,nSym

         LKsh=0

         Do iaSh=1,nShell    ! kOffSh(iSh,iSym)

            kOffSh(iaSh,iSyma) = LKsh

            LKsh = LKsh + nBasSh(iSyma,iaSh)

            MxBasSh = Max(MxBasSh,nBasSh(iSyma,iaSh))

         End Do

      End Do


C --- allocate memory for the diagonal elements of the Fock matrix
      Call mma_allocate(Fia,MxBasSh,Label='Fia')
      Call mma_allocate(Faa,nShell,Label='Faa')
      Fia(:)=Zero
      Faa(:)=Zero

C *** Determine S:= sum_l C(l)[k]^2  in each shell of C(a,k)
      Do jDen=1,nDen
         Do kSym=1,nSym

            Do jK=1,nOrb(kSym)
               jK_a = jK + kOff(kSym)

               Do iaSh=1,nShell

                  SKsh=zero
                  iS = kOffSh(iaSh,kSym) + 1
                  iE = kOffSh(iaSh,kSym) + nBasSh(kSym,iaSh)
                  Do ik=iS,iE
                     SKsh = SKsh
     &                    + CM(jDen)%SB(kSym)%A2(ik,jK)**2
                  End Do

                  SumAClk(iaSh,jk_a,jDen) = Sksh

               End Do

            End Do
         End Do
      End Do

C *** Compute Shell-pair Offsets in the K-matrix

      Do iSyma=1,nSym

         LKshp=0

         Do iaSh=1,nShell

          Do ibSh=1,nShell

            iShp = nShell*(iaSh-1) + ibSh

            nnBfShp(iShp,iSyma) = LKShp

            LKShp = LKShp + nBasSh(iSyma,iaSh)*nBasSh(iSyma,ibSh)

          End Do

         End Do

      End Do

C *** Mapping shell pairs from the full to the reduced set

      Call Mk_iShp_rs(iShp_rs,nShell)

      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1


C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym

        iAdr=0

        NumCV=NumCho(jSym)
        Call GAIGOP_SCAL(NumCV,'max')
        If (NumCV .lt. 1) Cycle

        JNUM=1
        Call Allocate_L_Full(L_Full,nShell,iShp_rs,JNUM,JSYM,nSym,
     &                       Memory=LFULL)

        iLoc = 3 ! use scratch location in reduced index arrays

C ****************     MEMORY MANAGEMENT SECTION    *****************
C ------------------------------------------------------------------
C --- compute memory needed to store at least 1 vector of JSYM
C --- and do all the subsequent calculations
C ------------------------------------------------------------------
         mTvec1= 0
         MxB=0
         do l=1,nSym
            k=Muld2h(l,JSYM)
            Mmax = Max(0,nOrb(k))
            If (Mmax.gt.0) MxB = Max(MxB,nBas(l))
            if (DoAct) Then
              mTvec1= mTvec1+ nAsh(k)*nBas(l)
            EndIf
         end do

         LFMAX = Max(3*mTvec1,LFULL) ! re-use memory for the active vec
         mTvec = nDen*Max(MxB,1) ! mem for storing half-transformed vec

C ------------------------------------------------------------------
C ------------------------------------------------------------------

         JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
         JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec
#if defined (_MOLCAS_MPP_)
         myJRED1=JRED1 ! first red set present on this node
         ntv0=0
#endif

c --- entire red sets range for parallel run
         Call GAIGOP_SCAL(JRED1,'min')
         Call GAIGOP_SCAL(JRED2,'max')

         kscreen=1
         DoScreen=.true.

         Do JRED=JRED1,JRED2

            CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

            If (nVrs.eq.0) GOTO 999  ! no vectors in that (jred,jsym)

            if (nVrs.lt.0) then
               Write(6,*)SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
     &                            ,nVrs
               call Abend
            endif

            Call Cho_X_SetRed(irc,iLoc,JRED)
c           !set index arrays at iLoc
            if(irc.ne.0)then
              Write(6,*)SECNAM//'cho_X_setred non-zero return code.'//
     &                          ' rc= ',irc
              call Abend
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            If(JSYM.eq.1)Then

               Call mma_allocate(Drs,nRS,Label='Drs')
               If (DoAct) Then
                 Call mma_allocate(Frs,nRS,2,Label='Frs')
               Else
                 Call mma_allocate(Frs,nRS,1,Label='Frs')
               End If
               Drs(:)=Zero
               Frs(:,:)=Zero

            EndIf

            Call mma_maxDBLE(LWORK)

            nVec = min(LWORK/(nRS+mTvec+LFMAX),min(nVrs,MaxVecPerBatch))

*Store nVec to make sure the routine always uses the same
            If (iAChoVec.eq.1) nVec_=nVec
            ReadInter=(iAChoVec.eq.2).and.(nVec.eq.nVec_)
!           nVec.ne.nVec_ should happen only if lack of memory
*           ReadInter=.false.

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) 'LWORK= ',LWORK
               WRITE(6,*) 'min. mem. need= ',nRS+mTvec+LFMAX
               WRITE(6,*) 'nRS= ',nRS
               WRITE(6,*) 'mTvec= ',mTvec
               WRITE(6,*) 'LFMAX= ',LFMAX
               WRITE(6,*) 'jsym= ',jsym
               CALL Quit(_RC_MEMORY_ERROR_)
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

            If(JSYM.eq.1)Then
C --- Transform the density to reduced storage
               mode = 'toreds'
               add = .False.
               nMat=1
               Call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,
     &                           [DLT],Drs,mode,add)
            EndIf

C --- BATCH over the vectors ----------------------------

            nBatch = (nVrs-1)/nVec + 1

            DO iBatch=1,nBatch

               If (iBatch.eq.nBatch) Then
                  JNUM = nVrs - nVec*(nBatch-1)
               else
                  JNUM = nVec
               endif


               JVEC = nVec*(iBatch-1) + iVrs
               IVEC2 = JVEC - 1 + JNUM

               CALL CWTIME(TCR1,TWR1)

               CALL CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,
     &                        NUMV,IREDC,MUSED)

               If (NUMV.le.0 .or.NUMV.ne.JNUM ) then
                  RETURN
               End If

               CALL CWTIME(TCR2,TWR2)
               tread(1) = tread(1) + (TCR2 - TCR1)
               tread(2) = tread(2) + (TWR2 - TWR1)

               If(JSYM.eq.1)Then
C ************ (alpha+beta) COULOMB CONTRIBUTION  ****************
C
C --- Contraction with the density matrix
C ---------------------------------------
C --- V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * DI(rs)
C==========================================================
C
                  CALL CWTIME(TCC1,TWC1)

                  Call mma_allocate(VJ,JNUM,Label='VJ')

                  CALL DGEMV_('T',nRS,JNUM,
     &                 ONE,Lrs,nRS,
     &                 Drs,1,ZERO,VJ,1)

C --- FI(rs){#J} <- FI(rs){#J} + FactCI * sum_J L(rs,{#J})*V{#J}
C===============================================================

                  Fact = dble(min(jVec-iVrs,1))

                  CALL DGEMV_('N',nRS,JNUM,
     &                 FactCI,Lrs,nRS,
     &                 VJ,1,Fact,Frs(:,1),1)


                  CALL CWTIME(TCC2,TWC2)
                  tcoul(1) = tcoul(1) + (TCC2 - TCC1)
                  tcoul(2) = tcoul(2) + (TWC2 - TWC1)

                  Call mma_deallocate(VJ)

               EndIf  ! Coulomb contribution

C *************** EXCHANGE CONTRIBUTIONS  ***********************
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
               Call Allocate_L_Full(L_Full,nShell,iShp_rs,JNUM,JSYM,
     &                              nSym)
               Call Allocate_Lab(Lab,JNUM,nBasSh,nBas,nShell,nSym,nDen)

               CALL CWTIME(TCS1,TWS1)
C ---------------------------------------------------------------------
C --- Estimate the diagonals :   D(a,b) = sum_J (Lab,J)^2
C
               If (Estimate) Then

                  Call Fzero(DIAG(1+iiBstR(jSym,1)),NNBSTR(jSym,1))

                  Do krs=1,nRS

                     mrs = iiBstR(JSYM,iLoc) + krs
                     jrs = IndRed(mrs,iLoc) ! address in 1st red set

                     Do jvc=1,JNUM

                        Diag(jrs) = Diag(jrs) + Lrs(krs,jvc)**2

                     End Do

                  End Do

               EndIf

               CALL CWTIME(TCS2,TWS2)
               tscrn(1) = tscrn(1) + (TCS2 - TCS1)
               tscrn(2) = tscrn(2) + (TWS2 - TWS1)

               CALL CWTIME(TCX1,TWX1)

C *** Reorder vectors to Full-dimensions
C ***
C *** Vectors are returned in the storage LaJ,b with the restriction:
C ***
C ***    Sym(a).ge.Sym(b)
C ***
C *** and blocked in shell pairs

               CALL CHO_getShFull(Lrs,lread,JNUM,JSYM,IREDC,L_Full,
     &                            SvShp,nnShl,iShp_rs,nnShl_tot)


               CALL CWTIME(TCX2,TWX2)
               texch(1) = texch(1) + (TCX2 - TCX1)
               texch(2) = texch(2) + (TWX2 - TWX1)


               IF (DoScreen) THEN

                   CALL CWTIME(TCS1,TWS1)

c --- Compute DH(a,b)=sqrt(D(a,b)) from the updated diagonals.
c ---                              Only the symmetry blocks with
c ---                              compound symmetry JSYM are computed
c --------------------------------------------------------------------
                   ired1 = 1 ! location of the 1st red set
                   Call swap_tosqrt(irc,ired1,NNBSTRT(1),JSYM,
     &                               DIAH,DIAG)

                   CALL CWTIME(TCS2,TWS2)
                   tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                   tscrn(2) = tscrn(2) + (TWS2 - TWS1)

               ENDIF


               Do kSym=1,nSym

                  lSym=MulD2h(JSYM,kSym)

                  Do jK=1,nOrb(kSym)

                    jk_a = kOff(kSym)+jK

                    Lab%A0(1:nDen*nBas(lSym)*JNUM)=Zero

                   IF (DoScreen) THEN

                     CALL CWTIME(TCS1,TWS1)
C------------------------------------------------------------------
C --- Setup the screening
C------------------------------------------------------------------

                     Do jDen=1,nDen

                        Do ik=1,nBas(kSym)
                           AbsC(ik) = abs(CM(jDen)%SB(kSym)%A2(ik,jK))
                        End Do

                        If (lSym.ge.kSym) Then
c --------------------------------------------------------------
C --- Y(l)[k] = sum_n  DH(l,n) * |C(n)[k]|
C===============================================================
                           Mode(1:1)='N'
                           n1 = nBas(lSym)
                           n2 = nBas(kSym)

                        Else
c --------------------------------------------------------------
C --- Y(l)[k] = sum_n  DH(n,l) * |C(n)[k]|
C===============================================================
                           Mode(1:1)='T'
                           n1 = nBas(kSym)
                           n2 = nBas(lSym)

                        EndIf

                        If (n1>0)
     &                  CALL DGEMV_(Mode(1:1),n1,n2,
     &                             ONE,DiaH%SB(lSym,kSym)%A2,n1,
     &                                 AbsC,1,
     &                            ZERO,Ylk(1,jK_a,jDen),1)

                     End Do

C --- List the shells present in Y(l)[k] by the largest element
                     Do jDen=1,nDen
                        Do ish=1,nShell
                           YshMax=zero
                           Do ibs=1,nBasSh(lSym,ish)
                              YshMax = Max(YshMax,
     &                          Ylk(koffSh(ish,lSym)+ibs,jK_a,jDen))
                           End Do
                           MLk(ish,jK_a,jDen) = YshMax
                        End Do
                     End Do


C --- Sort the lists ML[k]
                     Do jDen=1,nDen
                        Do ish=1,nShell
                           Indx(iSh,jk_a,jDen) = ish
                        End Do
                     End Do

C ****  The Max in the MO set 1 is used as reference
                     numSh1=0  ! # of significant shells in MO set 1
                     YMax=MLk(1,jK_a,1)
                     jmlmax=1
                     Do iml=2,nShell  ! get the max in the MO set 1
                        If (MLk(iml,jK_a,1).gt.YMax) then
                           YMax = MLk(iml,jK_a,1)
                           jmlmax = iml
                        Endif
                     End Do
                     If (jmlmax.ne.1) then  ! swap positions
                        xTmp = MLk(1,jK_a,1)
                        iTmp = Indx(1,jk_a,1)
                        MLk(1,jK_a,1) = YMax
                        Indx(1,jk_a,1) = Indx(jmlmax,jk_a,1)
                        MLk(jmlmax,jK_a,1) = xTmp
                        Indx(jmlmax,jk_a,1) = iTmp
                     Endif

C **** Sort the list for the MO set 2   iff  MOs1.ne.MOs2
                     If (.not.Fake_CMO2) Then
                       numSh2=0  ! # of significant shells in MO set 2
                       jml=1
                       Do while (jml.le.nShell)

                         YMax=MLk(jml,jK_a,2)
                         jmlmax=jml

                         Do iml=jml+1,nShell  ! get the max
                           If (MLk(iml,jK_a,2).gt.YMax) then
                              YMax = MLk(iml,jK_a,2)
                              jmlmax = iml
                           Endif
                         End Do

                         If(jmlmax.ne.jml) then  ! swap positions
                          xTmp = MLk(jml,jK_a,2)
                          iTmp = Indx(jml,jk_a,2)
                          MLk(jml,jK_a,2) = YMax
                          Indx(jml,jk_a,2) = Indx(jmlmax,jk_a,2)
                          MLk(jmlmax,jK_a,2) = xTmp
                          Indx(jmlmax,jk_a,2) = iTmp
                         Endif

c --- Exact bounds (quadratic scaling of the MO transformation)
c --- Note that in true RASSI the exchange matrix is not
c --- positive definite.
c
                         If(MLk(jml,jK_a,2).ge.xtau(2))then
                           numSh2 = numSh2 + 1
                         else
                           jml=nShell  ! exit the loop
                         endif

                         jml=jml+1

                       End Do

                       Indx(0,jk_a,2) = numSh2
                       numSh1 = 1

                     Else ! fake biorthonormal basis

                       numSh2 = 6669666 ! dummy assignement

                       If (MLk(1,jK_a,1) .ge. xtau(1))  numSh1 = 1

                     EndIf

C **** Sort the list for the MO set 1 only if needed
                     If(numSh2.gt.0) then
                       jml=2 ! the 1st element has already been treated
                       Do while (jml.le.nShell)

                          YMax=MLk(jml,jK_a,1)
                          jmlmax=jml
                          Do iml=jml+1,nShell  ! get the max
                             If (MLk(iml,jK_a,1).gt.YMax) then
                                YMax = MLk(iml,jK_a,1)
                                jmlmax = iml
                             Endif
                          End Do

                          If(jmlmax.ne.jml) then  ! swap positions
                            xTmp = MLk(jml,jK_a,1)
                            iTmp = Indx(jml,jk_a,1)
                            MLk(jml,jK_a,1) = YMax
                            Indx(jml,jk_a,1) = Indx(jmlmax,jk_a,1)
                            MLk(jmlmax,jK_a,1) = xTmp
                            Indx(jmlmax,jk_a,1) = iTmp
                          Endif

                          If( .not.Fake_CMO2  .and.
     &                       MLk(jml,jK_a,1).ge.xtau(1))then
                             numSh1 = numSh1 + 1

c --- Here we use a non-exact bound for the exchange matrix because a
c     fake rassi (MOs1=MOs2) has a positive definite exchange
                          ElseIf ( Fake_CMO2  .and.
     &                            MLk(jml,jK_a,1).ge.xtau(1) ) then
                             numSh1 = numSh1 + 1
                          Else
                             jml=nShell  ! exit the loop
                          Endif

                          jml=jml+1

                       End Do
                     Else
                       numSh1 = 0
                     EndIf

                     Indx(0,jk_a,1) = numSh1

                     CALL CWTIME(TCS2,TWS2)
                     tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                     tscrn(2) = tscrn(2) + (TWS2 - TWS1)
C------------------------------------------------------------------
                   ENDIF    ! Screening setup

C --- Transform vectors for shells in the lists ML[k]
C
C --- Screening based on the Frobenius norm: sqrt(sum_ij A(i,j)^2)
C
C ---  || La,J[k] ||  .le.  || Lab,J || * || Cb[k] ||

                      CALL CWTIME(TCT1,TWT1)

                      Do jDen=1,nDen

*
*MGD maybe check that NumSh1 matches the one on the record?
*
                         Do iSh=1,Indx(0,jk_a,jDen)

                           iaSh = Indx(iSh,jk_a,jDen)

                           Lab%Keep(iaSh,jDen)=.True.

                           ibcount=0

                           Do ibSh=1,nShell

                            iOffShb = kOffSh(ibSh,kSym)

                            iShp = iTri(iaSh,ibSh)

                            If (iShp_rs(iShp)<=0) Cycle

                            IF (lSym.ge.kSym) Then

                              If(nnBstRSh(JSym,iShp_rs(iShp),iLoc)*
     &                           nBasSh(lSym,iaSh)*
     &                           nBasSh(kSym,ibSh) .gt. 0
     &                        .and. abs(SumAClk(ibSh,jK_a,jDen)*
     &                        SvShp(iShp_rs(iShp),1) ).ge. thrv(jDen))
     &                           Then

                              If (ReadInter.and.(jDen.eq.1)) Then
*
**   Read vectors
*

                                 ibcount=1

                                 Exit
                              Else
*
**   Or compute them
*
*
                                 ibcount = ibcount + 1

                                 l1=1
                                 If (iaSh<ibSh) l1=2


C ---  LaJ,[k] = sum_b  L(aJ,b) * C(b)[k]
C ---------------------------------------

                                 CALL DGEMV_('N',nBasSh(lSym,iaSh)*JNUM,
     &                                        nBasSh(kSym,ibSh),
     &                        ONE,L_Full%SPB(lSym,iShp_rs(iShp),l1)%A21,
     &                                        nBasSh(lSym,iaSh)*JNUM,
     &                            CM(jDen)%SB(kSym)%A2(iOffShb+1:,jK),1,
     &                                  ONE,Lab%SB(iaSh,lSym,jDen)%A,1)

                              End If


                              EndIf

                            Else   ! lSym < kSym
                              If(nnBstRSh(JSym,iShp_rs(iShp),iLoc)*
     &                           nBasSh(lSym,iaSh)*
     &                           nBasSh(kSym,ibSh) .gt. 0
     &                        .and. abs(SumAClk(ibSh,jK_a,jDen)*
     &                        SvShp(iShp_rs(iShp),1) ).ge. thrv(jDen))
     &                           Then
                                 ibcount = ibcount + 1

                                 l1 = 1
                                 If (ibSh<iaSh) l1 = 2


C ---  LJa,[k] = sum_b  L(b,Ja) * C(b)[k]
C ---------------------------------------

                                  CALL DGEMV_('T',nBasSh(kSym,ibSh),
     &                                     JNUM*nBasSh(lSym,iaSh),
     &                      One,L_Full%SPB(kSym,iShp_rs(iShp),l1)%A12,
     &                                      nBasSh(kSym,ibSh),
     &                            CM(jDen)%SB(kSym)%A2(iOffShb+1:,jK),1,
     &                                  ONE,Lab%SB(iaSh,lSym,jDen)%A,1)


                            EndIf
                            EndIf

                           End Do ! ibSh

                           If (lSym>=kSym) Then

                           If (ReadInter.and.(jDen.eq.1)) Then
**
                             If (ibcount.gt.0) Then
                               lvec=nBasSh(lSym,iaSh)*JNUM
                               call DDAFILE(LuIChoVec(Jsym),2,
     &                           Lab%SB(iaSh,lSym,jDen)%A,
     &                         lvec,iAdr)
                             EndIf

                           Else
*
**   Store vectors
*
                             If ((iAChoVec.eq.1).and.(jDen.eq.1)
     &                           .and.(ibcount.gt.0)) Then
                               lvec=nBasSh(lSym,iaSh)*JNUM
                               call DDAFILE(LuIChoVec(Jsym),1,
     &                           Lab%SB(iaSh,lSym,jDen)%A,
     &                         lvec,iAdr)
*
                             EndIf
                           End If
                           End If

c --- The following re-assignement is used later on to check if the
c --- iaSh vector LaJ[k] can be neglected because identically zero

                             If (ibcount==0)
     &                          Lab%Keep(iaSh,jDen) = .False.

                         End Do ! iSh

                         Do iSh=Indx(0,jk_a,jDen)+1,nshell
                            iaSh = Indx(iSh,jk_a,jDen)
                            Lab%Keep(iaSh,jDen) = .False.
                         End Do

                      End Do  ! jDen

                      CALL CWTIME(TCT2,TWT2)
                      tmotr(1) = tmotr(1) + (TCT2 - TCT1)
                      tmotr(2) = tmotr(2) + (TWT2 - TWT1)

C --- Prepare the J-screening

                      CALL CWTIME(TCS1,TWS1)

                      Do iSh=1,Indx(0,jk_a,1)

                         iaSh = Indx(iSh,jk_a,1)

                         iaSkip=Merge(1,0,Lab%Keep(iaSh,1))

                         jaSkip=Merge(1,0,Lab%Keep(iaSh,kDen))

                         If (iaSkip*jaSkip==0) Cycle

                         IF (lSym.ge.kSym) Then

C ---  Faa,[k] = sum_J  LaJ[k2]*LaJ[k1]
C -------------------------------------
                            Inc=nBasSh(lSym,iaSh)
                            n1=1

                         Else   ! lSym < kSym

C ---  Faa,[k] = sum_J  LJa[k2]*LJa[k1]
C -------------------------------------
                            Inc=1
                            n1=JNUM

                         End If

                         Temp=Zero
                         Do ia=1,nBasSh(lSym,iaSh)
                            Fia(ia)=DDot_(JNUM,
     &                       Lab%SB(iaSh,lSym,kDen)%A(1+n1*(ia-1):),Inc,
     &                       Lab%SB(iaSh,lSym,   1)%A(1+n1*(ia-1):),Inc)
                            Temp=Max(Abs(Fia(ia)),Temp)
                         End Do

                         Faa(iaSh)=Temp

                      End Do

                      CALL CWTIME(TCS2,TWS2)
                      tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                      tscrn(2) = tscrn(2) + (TWS2 - TWS1)


C------------------------------------------------------------
C --- Compute exchange matrix for the interacting shell pairs
C------------------------------------------------------------

                      CALL CWTIME(TCX1,TWX1)

                      Do lSh=1,Indx(0,jk_a,1)

                         iaSh = Indx(lSh,jk_a,1)

                         iaSkip=Merge(1,0,Lab%Keep(iaSh,kDen))

                         mSh = 1

                         Do while (mSh.le.Indx(0,jk_a,kDen))

                            ibSh = Indx(mSh,jk_a,kDen)

                            ibSkip=Merge(1,0,Lab%Keep(ibSh,   1))

                            iShp = nShell*(iaSh-1) + ibSh

                            iOffShb = kOffSh(ibSh,lSym)

                            iOffAB = nnBfShp(iShp,lSym)

                            xFab = sqrt(abs(Faa(iaSh)*Faa(ibSh)))

                            If (MLk(lSh,jK_a,1)*MLk(mSh,jK_a,kDen)
     &                          .lt.tau) Then

                                mSh = Indx(0,jk_a,kDen) !skip rest

                            ElseIf ( xFab.ge.tau/MaxRedT
     &                              .and. iaSkip*ibSkip.eq.1) Then

                               nBsa = nBasSh(lSym,iaSh)
                               IF (lSym.ge.kSym) Then

C ---  F(a,b)[k] = F(a,b)[k] + FactXI * sum_J  X2(a,J)[k] * X1(b,J)[k]
C --------------------------------------------------------------------
                                  n1 = nBasSh(lSym,iaSh)
                                  n2 = nBasSh(lSym,ibSh)
                                  Mode(1:1)='N'
                                  Mode(2:2)='T'

                               ELSE   ! lSym < kSym

C ---  F(a,b)[k] = F(a,b)[k] + FactXI * sum_J  X2(J,a)[k] * X1(J,b)[k]
C --------------------------------------------------------------------
                                  n1 = JNUM
                                  n2 = JNUM
                                  Mode(1:1)='T'
                                  Mode(2:2)='N'

                               EndIf

                               CALL DGEMM_(Mode(1:1),Mode(2:2),
     &                         nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,
     &                            FactXI,Lab%SB(iaSh,lSym,kDen)%A,n1,
     &                                   Lab%SB(ibSh,lSym,   1)%A,n2,
     &                               ONE,KI%SB(lSym)%A1(1+iOffAB:),nBsa)

                            EndIf


                            mSh = mSh + 1  ! update shell counter


                         End Do

                      End Do

                      CALL CWTIME(TCX2,TWX2)
                      texch(1) = texch(1) + (TCX2 - TCX1)
                      texch(2) = texch(2) + (TWX2 - TWX1)


                  End Do  ! loop over k MOs


               End Do   ! loop over MOs symmetry

               Call Deallocate_Lab(Lab)
               Call Deallocate_L_Full(L_Full)
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
               DoScreen=.false. ! avoid redo screening inside batch loop

C --- Diagonals updating. It only makes sense if Nscreen > 0

               If (Update .and. Nscreen .gt. 0) Then

                  CALL CWTIME(TCS1,TWS1)
C ---------------------------------------------------------------------
C --- update the diagonals :   D(a,b) = D(a,b) - sum_J (Lab,J)^2
C
C --- subtraction is done in the 1st reduced set
#if defined (_MOLCAS_MPP_)
                  If (nProcs .gt. 1 .and. Is_Real_Par()) then

                   Do krs=1,nRS

                     mrs = iiBstR(JSYM,iLoc) + krs
                     jrs = IndRed(mrs,iLoc) - iiBstR(JSYM,1)

                     Do jvc=1,JNUM

                        DiagJ(jrs) = DiagJ(jrs) + Lrs(krs,jvc)**2
                     End Do

                   End Do

                  Else

                   Do krs=1,nRS

                     mrs = iiBstR(JSYM,iLoc) + krs
                     jrs = IndRed(mrs,iLoc) ! address in 1st red set

                     Do jvc=1,JNUM

                        Diag(jrs) = Diag(jrs) - Lrs(krs,jvc)**2
                     End Do

                   End Do

                  EndIf

#else
                  Do krs=1,nRS

                     mrs = iiBstR(JSYM,iLoc) + krs
                     jrs = IndRed(mrs,iLoc) ! address in 1st red set

                     Do jvc=1,JNUM

                        Diag(jrs) = Diag(jrs) - Lrs(krs,jvc)**2
                     End Do

                  End Do
#endif

                  CALL CWTIME(TCS2,TWS2)
                  tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                  tscrn(2) = tscrn(2) + (TWS2 - TWS1)

               EndIf

C ************  END EXCHANGE CONTRIBUTION  ****************

C --------------------------------------------------------------------
C --- First half Active transformation  Lvb,J = sum_a  C1(v,a) * Lab,J
C --------------------------------------------------------------------

               CALL CWTIME(TCINT1,TWINT1)

C --- Set up the skipping flags
C --- The memory used before for the full-dimension AO-vectors
C ---     is now re-used to store half and full transformed
C ---     vectors in the active space
C -------------------------------------------------------------
               If (DoAct) Then

               iSwap = 0  ! Lvb,J are returned
               ! Lvb,J
               ! Lvi,J i general MO index
               ! L~vi,J ~ transformed index
               Call Allocate_SBA(Lpq(1),nAsh,nBas,nVec,JSYM,nSym,iSwap)
               Call Allocate_SBA(Lpq(2),nAsh,nAsh,nVec,JSYM,nSym,iSwap)
               Call Allocate_SBA(Lpq(3),nAsh,nAsh,nVec,JSYM,nSym,iSwap)

*MGD should we compute only if there are active orbitals in this sym?
*
**             Read vectors
*
               kMOs = 1  !
               nMOs = 1  ! Active MOs (1st set)
               If (iAChoVec.eq.2) Then
                 ioff=0
                 Do i=1,nSym
                    k = Muld2h(i,JSYM)
                    lvec=nAsh(k)*nBas(i)*JNUM
                    iAdr2=(JVEC-1)*nAsh(k)*nBas(i)+ioff
                    call DDAFILE(LuAChoVec(Jsym),2,Lpq(1)%SB(k)%A3,
     &                           lvec,iAdr2)
                    ioff=ioff+nAsh(k)*nBas(i)*NumCho(jSym)
                 End Do
               Else
* Lrs * MO
                 CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,
     &                           JSYM,iSwap,IREDC,nMOs,kMOs,ASh(1),
     &                           Lpq,DoRead)
               EndIf

               if (irc.ne.0) then
                  RETURN
               endif

*
**             Store vectors
*
               ioff=0
               If (iAChoVec.eq.1) Then
                 Do i=1,nSym
                    k = Muld2h(i,JSYM)
                    lvec=nAsh(k)*nBas(i)*JNUM
                    iAdr2=(JVEC-1)*nAsh(k)*nBas(i)+ioff
                    call DDAFILE(LuAChoVec(Jsym),1,Lpq(1)%SB(k)%A3,
     &                           lvec,iAdr2)
                    ioff=ioff+nAsh(k)*nBas(i)*NumCho(jSym)
                 End Do
               EndIf
               CALL CWTIME(TCINT2,TWINT2)
               tint1(1) = tint1(1) + (TCINT2 - TCINT1)
               tint1(2) = tint1(2) + (TWINT2 - TWINT1)

C --------------------------------------------------------------------
C --- Active-Active transformation  Lvw,J = sum_b  Lvb,J * C2(w,b)
C --------------------------------------------------------------------
               Do iSymb=1,nSym

                  iSymv = MulD2h(JSYM,iSymb)
                  NAv = nAsh(iSymv)
                  NAw = nAsh(iSymb)

                  If(NAv*Naw.ne.0)Then

                   Do JVC=1,JNUM

                    CALL CWTIME(TCINT2,TWINT2)
*Lv~w
                    If (.not.Fake_CMO2) Then
                     CALL DGEMM_('N','T',NAv,Naw,NBAS(iSymb),
     &                          One,Lpq(1)%SB(iSymv)%A3(:,:,JVC),NAv,
     &                              Ash(2)%SB(iSymb)%A2,NAw,
     &                         Zero,Lpq(3)%SB(iSymv)%A3(:,:,JVC),NAv)
                    CALL CWTIME(TCINT4,TWINT4)
               tint1(1) = tint1(1) + (TCINT4 - TCINT2)
               tint1(2) = tint1(2) + (TWINT4 - TWINT2)
C --------------------------------------------------------------------
C --- Formation of the Q matrix Qpx = Lpy Lv~w Gxyvw
C --------------------------------------------------------------------
*~Lxy=Lv~w Gxyvw
*MGD probably additional nSym loop
                     CALL DGEMV_('N',NAv*Naw,NAv*Naw,
     &                  ONE,G2,NAv*Naw,
     &                      Lpq(3)%SB(iSymv)%A3(:,:,JVC),1,
     &                  ZERO,Lpq(2)%SB(iSymv)%A3(:,:,JVC),1)
*Qpx=Lpy ~Lxy
                     Call DGEMM_('T','N',NBAS(iSymb),NAw,Nav,
     &                          Two,Lpq(1)%SB(iSymv)%A3(:,:,JVC),NAv,
     &                              Lpq(2)%SB(iSymv)%A3(:,:,JVC),Naw,
     &                          ONE,QTmp(2)%SB(iSymb)%A2,NBAS(iSymb))
                      CALL CWTIME(TCINT3,TWINT3)
                      tQmat(1) = tQmat(1) + (TCINT3 - TCINT4)
                      tQmat(2) = tQmat(2) + (TWINT3 - TWINT4)
                    EndIf
                   CALL CWTIME(TCINT3,TWINT3)
C --------------------------------------------------------------------
C --- Active-Active transformation  Lvw,J = sum_b  Lvb,J * C2(w,b)
C --------------------------------------------------------------------
*Lvw
                    CALL DGEMM_('N','T',NAv,Naw,NBAS(iSymb),
     &                         One,Lpq(1)%SB(iSymv)%A3(:,:,JVC),NAv,
     &                             Ash(1)%SB(iSymb)%A2,NAw,
     &                        Zero,Lpq(2)%SB(iSymv)%A3(:,:,JVC),NAv)

                   CALL CWTIME(TCINT2,TWINT2)
                   tint1(1) = tint1(1) + (TCINT2 - TCINT3)
                   tint1(2) = tint1(2) + (TWINT2 - TWINT3)
                   End Do
*
C
C
C *************** EVALUATION OF THE (TW|XY) INTEGRALS ***********
                   If (iSymv.gt.iSymb) Go to 20
                   Do isymx=1,iSymb
                     iSymy=MulD2h(JSYM,iSymx)
                     If (iSymy.gt.iSymx.or.
     &                  (iSymb.eq.isymx.and.iSymy.gt.iSymv)) Cycle
                     Nax=nAsh(iSymx)
                     Nay=nAsh(iSymy)
                     If(NAx*Nay.ne.0)Then

                      i = 2
                      If (.not.Fake_CMO2) i = 3

* (tu|v~w) = Ltu*Lv~w
                        Call DGEMM_('N','T',Nav*Naw,NAx*Nay,JNUM,
     &                       One, Lpq(2)%SB(iSymv)%A3,NAv*Naw,
     &                            Lpq(i)%SB(iSymx)%A3,NAx*Nay,
     &                       One,MOScr%SB(iSymb,iSymv,iSymx)%A2,NAv*Naw)
                     ENdIf
                   EndDo
                   CALL CWTIME(TCINT3,TWINT3)
                   tint2(1) = tint2(1) + (TCINT3 - TCINT2)
                   tint2(2) = tint2(2) + (TWINT3 - TWINT2)
 20                Continue
                  EndIf
               EndDo
C
C
C ************ EVALUATION OF THE ACTIVE FOCK MATRIX *************
*Coulomb term
               If (JSYM.eq.1) Then

                 Call mma_allocate(VJ,JNUM,Label='VJ')
                 VJ(:)=Zero

                 Do iSymb=1,nSym

                    iSymv = MulD2h(JSYM,iSymb)
                    NAv = nAsh(iSymv)
                    NAw = nAsh(iSymb)

                    If(NAv*Naw.ne.0)Then
                     i=3
                     If (Fake_CMO2) i = 2

                     CALL DGEMV_('T',Nav*Naw,JNUM,
     &                          ONE,Lpq(i)%SB(iSymv)%A3,Nav*Naw,
     &                         DA%SB(iSymB)%A2,1,ONE,VJ,1)
                    EndIf
                 EndDo
*
                 CALL DGEMV_('N',nRS,JNUM,
     &                     -FactCI,Lrs,nRS,
     &                     VJ,1,1.0d0,Frs(:,2),1)

                 Call mma_deallocate(VJ)

               EndIf
                   CALL CWTIME(TCINT2,TWINT2)
                   tact(1) = tact(1) + (TCINT2 - TCINT3)
                   tact(2) = tact(2) + (TWINT2 - TWINT3)
C --------------------------------------------------------------------
C --- Formation of the Q matrix Qpx = L~py Lvw Gxyvw
C --------------------------------------------------------------------
               Do iSymb=1,nSym

                  iSymv = MulD2h(JSYM,iSymb)
                  NAv = nAsh(iSymv)
                  NAw = nAsh(iSymb)

                  If(NAv*Naw.ne.0)Then

                    Lpq(3)%SB(iSymv)%A3(:,:,:)=Zero

                    Do iSymx=1,nSym
                      iSymy=MulD2h(JSYM,iSymx)
                      Nax=nAsh(iSymx)
                      Nay=nAsh(iSymy)

                      ipG    = 1 +iASQ(isymb,iSymv,iSymx)

                      If(NAx*Nay.ne.0)Then
                       Call DGEMM_('N','N',Nav*Naw,JNUM,NAx*Nay,
     &                              One,G2(ipG),NAv*Naw,
     &                                  Lpq(2)%SB(iSymy)%A3,NAx*Nay,
     &                              ONE,Lpq(3)%SB(iSymv)%A3,Nav*Naw)
                      EndIf

                    End Do

                    Do JVC=1,JNUM
                      Call DGEMM_('T','N',NBAS(iSymb),NAw,Nav,
     &                           One,Lpq(1)%SB(iSymv)%A3(:,:,JVC),NAv,
     &                               Lpq(3)%SB(iSymv)%A3(:,:,JVC),Nav,
     &                           ONE,QTmp(1)%SB(iSymb)%A2,NBAS(iSymb))
                    End Do

*MGD check timing loops are correct (not always timed from
*previous reference)
                 EndIf
               ENd Do
               CALL CWTIME(TCINT3,TWINT3)
               tQmat(1) = tQmat(1) + (TCINT3 - TCINT2)
               tQmat(2) = tQmat(2) + (TWINT3 - TWINT2)

               Call Deallocate_SBA(Lpq(2))

               Call Allocate_SBA(Lpq(2),nAsh,nBas,nVec,JSYM,nSym,iSwap)
C
C
C ************ EVALUATION OF THE ACTIVE FOCK MATRIX *************
*Exchange term
               Do iSymb=1,nSym

                  iSymv = MulD2h(JSYM,iSymb)
                  NAv = nAsh(iSymv)
                  NAw = nAsh(iSymb)

                  If(NAv*nBas(iSymb).ne.0)Then

                   Do JVC=1,JNUM
                     Call DGEMM_('T','N',NBAS(iSymb),Nav,Nav,
     &                           ONE,Lpq(1)%SB(iSymv)%A3(:,:,JVC),Nav,
     &                               DA%SB(iSymb)%A2,Nav,
     &                    ZERO,Lpq(2)%SB(iSymv)%A3(:,:,JVC),NBAS(iSymb))
                   End Do
                   CALL CWTIME(TCINT2,TWINT2)
                   tact(1) = tact(1) + (TCINT2 - TCINT3)
                   tact(2) = tact(2) + (TWINT2 - TWINT3)

C --------------------------------------------------------------------
C --- First half Active transformation  L~vb,J = sum_a  ~C(v,a) * Lab,J
C --------------------------------------------------------------------
                   If (.not.Fake_CMO2) Then
                     CALL CWTIME(TCINT2,TWINT2)

                     CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,
     &                           JSYM,iSwap,IREDC,nMOs,kMOs,Ash(2),
     &                           Lpq,DoRead)
                     CALL CWTIME(TCINT3,TWINT3)
                     tint1(1) = tint1(1) + (TCINT3 - TCINT2)
                     tint1(2) = tint1(2) + (TWINT3 - TWINT2)
C --------------------------------------------------------------------
C --- Formation of the Q matrix Qpx = Lp~y Lvw Gxyvw
C --------------------------------------------------------------------
                     Do JVC=1,JNUM
                      Call DGEMM_('T','N',NBAS(iSymb),NAw,Nav,
     &                           One,Lpq(1)%SB(iSymb)%A3(:,:,JVC),NAv,
     &                               Lpq(3)%SB(iSymv)%A3(:,:,JVC),Naw,
     &                           ONE,QTmp(2)%SB(iSymb)%A2,NBAS(iSymb))
                     End Do
                     CALL CWTIME(TCINT2,TWINT2)
                     tQmat(1) = tQmat(1) + (TCINT2 - TCINT3)
                     tQmat(2) = tQmat(2) + (TWINT2 - TWINT3)
                   EndIf
                 EndIf
               EndDO ! iSymb
C
C
C ************ EVALUATION OF THE ACTIVE FOCK MATRIX *************
*Exchange term
*MGD what if Naw=0 but not NBas(b)?
               Do iSymb=1,nSym

                  iSymv = MulD2h(JSYM,iSymb)
                  NAv = nAsh(iSymv)
                  NAw = nAsh(iSymb)

                  If(NAv*Naw.ne.0)Then

                   Do JVC=1,JNUM
                     Do is=1,NBAS(iSymb)
                      CALL DGEMV_('N',NBAS(iSymb),Nav,
     &                 -FactXI,Lpq(2)%SB(iSymb)%A3(:,:,JVC),NBAS(iSymb),
     &                         Lpq(1)%SB(iSymv)%A3(:,is,JVC),1,
     &                     ONE,KA%SB(iSymB)%A2(:,is),1)

                    EndDo
                   End Do
                   CALL CWTIME(TCINT3,TWINT3)
                   tact(1) = tact(1) + (TCINT3 - TCINT2)
                   tact(2) = tact(2) + (TWINT3 - TWINT2)
*
                  EndIf

               End Do ! iSymb
*
               Call Deallocate_SBA(Lpq(3))
               Call Deallocate_SBA(Lpq(2))
               Call Deallocate_SBA(Lpq(1))

               EndIf ! If (DoAct)

               CALL CWTIME(TCINT2,TWINT2)
               tintg(1) = tintg(1) + (TCINT2 - TCINT1)
               tintg(2) = tintg(2) + (TWINT2 - TWINT1)

               if (irc.ne.0) then
                  RETURN
               endif

C ---------------- END (TW|XY) EVALUATION -----------------------

            END DO  ! end batch loop


            If(JSYM.eq.1)Then
c --- backtransform fock matrix to full storage
               mode = 'tofull'
               add = .True.
               nMat = 1
               Call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,
     &                           [JI],Frs(:,1),mode,add)
               If (DoAct) Then
                  Call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,
     &                              [JALT],Frs(:,2),mode,add)
               End If
            EndIf

C --- free memory
            Call mma_deallocate(Lrs)

            If(JSYM.eq.1)Then
              Call mma_deallocate(Frs)
              Call mma_deallocate(Drs)
            EndIf


999         Continue

C --- Screening control section
            DoScreen = kscreen.eq.Nscreen

            if (.not.DoScreen) then
                kscreen = kscreen + 1
            else
                kscreen = 1
            endif

#if defined (_MOLCAS_MPP_)
            If (nProcs.gt.1 .and. Update .and. DoScreen
     &          .and. Is_Real_Par()) Then
               Call GaDsum(DiagJ,nnBSTR(JSYM,1))
               Call Daxpy_(nnBSTR(JSYM,1),xone,DiagJ,1,
     &                    Diag(1+iiBstR(JSYM,1)),1)
               Call Fzero(DiagJ,nnBSTR(JSYM,1))
            EndIf
C--- Need to activate the screening to setup the contributing shell
C--- indeces the first time the loop is entered .OR. whenever other nodes
C--- have performed screening in the meanwhile
            If (nProcs.gt.1 .and. .not.DoScreen .and. nVrs.eq.0
     &          .and. Is_Real_Par()) Then
               ntv0=ntv0+1
               DoScreen = (JRED.lt.myJRED1 .or. ntv0.ge.Nscreen)
               if (DoScreen) ntv0=0
            EndIf
#endif

         END DO   ! loop over red sets

      END DO  ! loop over JSYM

* --- Accumulate Coulomb and Exchange contributions
      Do iSym=1,nSym

         Do iaSh=1,nShell

            ioffa = kOffSh(iaSh,iSym)

            Do ibSh=1,nShell

               iShp = nShell*(iaSh-1) + ibSh
               iOffAB = nnBfShp(iShp,iSym)

               iShp = nShell*(ibSh-1) + iaSh
               jOffAB = nnBfShp(iShp,iSym)

               ioffb = kOffSh(ibSh,iSym)

               Do ib=1,nBasSh(iSym,ibSh)

                Do ia=1,nBasSh(iSym,iaSh)

                  iab = nBasSh(iSym,iaSh)*(ib-1) + ia
*MGD warning with sym
                  jab = nBasSh(iSym,ibSh)*(ia-1) + ib

                  iag = ioffa + ia
                  ibg = ioffb + ib

                  iabg = iTri(iag,ibg)

                  FkI%SB(iSym)%A2(iag,ibg)
     &                     = JI%SB(iSym)%A1(iabg)
     &                     + KI%SB(iSym)%A1(iOffAB+iab)
     &                     + KI%SB(iSym)%A1(jOffAB+jab)

                  FkA%SB(iSym)%A2(iag,ibg)
     &                     = JALT%SB(iSym)%A1(iabg)
     &                     + KA%SB(iSym)%A2(iag,ibg)
     &                     + KA%SB(iSym)%A2(ibg,iag)

                End Do

               End Do

            End Do

         End Do

      End Do
      CALL CWTIME(TCINT1,TWINT1)

      Call Deallocate_DSBA(JALT)

      If (DoAct) Then
*
**Compute MO integrals
*
        Do iS=1,nSym
          Do jS=1,nsym
            ijS=iEOR(is-1,js-1)+1
            Do kS=1,nSym
              ls=iEOr(ijs-1,ks-1)+1

              Do iAsh=1,nAsh(is)
                Do jAsh=1,nAsh(js)
                  iij=itri(iAsh+kAOff(is),jAsh+kAOff(jS))
                  Fac1=1.0d0
                  If (iAsh+kAOff(is).eq.jAsh+kAoff(jS)) Fac1=2.0d0

                  Do kAsh=1,nAsh(ks)
                    Do lAsh=1,nAsh(ls)
                      ikl=itri(lAsh+kAOff(lS),kAsh+kAOff(kS))
                      Fac2=1.0d0
                      If (lAsh+kAOff(lS).eq.kAsh+kAOff(kS)) Fac2=2.0d0
                      If (iij.ne.ikl) Fac2=Fac2*0.5d0

                      ipG=itri(iij,ikl)

                      MO_Int(ipG)=MO_Int(ipG)+Fac1*Fac2
     &                       *MOScr%SB(iS,jS,kS)%A4(iAsh,jAsh,kAsh,lAsh)
                    End Do
                  End Do
                End Do
              End Do

            End Do
          End Do
        End Do
      EndIf
*
**Transform Fock and Q matrix to MO basis
*
      Do iS=1,nSym
        jS=iS
        If (nBas(iS).ne.0) Then
          If (DoAct) Then
            Call DGEMM_('T','N',nBas(jS),nBas(iS),nBas(iS),
     &                  1.0d0,FkA%SB(iS)%A2,nBas(iS),
     &                        CMO%SB(iS)%A2,nBas(iS),
     &                  0.0D0,JA%SB(iS)%A2,nBas(jS))
            Call DGEMM_('T','N',nBas(jS),nBas(jS),nBas(iS),
     &                  1.0d0,JA%SB(iS)%A2,nBas(iS),
     &                        CMO%SB(jS)%A2,nBas(jS),
     &                  0.0d0,FkA%SB(iS)%A2,nBas(jS))
          EndIf
          Call DGEMM_('T','N',nBas(jS),nBas(iS),nBas(iS),
     &                1.0d0,FkI%SB(iS)%A2,nBas(iS),
     &                      CMO%SB(iS)%A2,nBas(iS),
     &                0.0d0,JA%SB(iS)%A2,nBas(jS))
          Call DGEMM_('T','N',nBas(jS),nBas(jS),nBas(iS),
     &                1.0d0,JA%SB(iS)%A2,nBas(iS),
     &                      CMO%SB(jS)%A2,nBas(jS),
     &                0.0d0,FkI%SB(iS)%A2,nBas(jS))
          If (DoAct) Then
            If (Fake_CMO2) Then
              Call DGEMM_('T','N',nBas(jS),nAsh(iS),nBas(jS),
     &                     1.0d0,CMO%SB(iS)%A2,nBas(jS),
     &                           QTmp(1)%SB(js)%A2,nBas(jS),
     &                     0.0d0,QVec%SB(js)%A2,nBas(jS))
            Else
              Call DGEMM_('T','N',nBas(jS),nAsh(iS),nBas(jS),
     &                     1.0d0,CMO%SB(iS)%A2,nBas(jS),
     &                           QTmp(2)%SB(jS)%A2,nBas(jS),
     &                     0.0d0,QVec%SB(jS)%A2,nBas(jS))
              Call DGEMM_('T','N',nBas(jS),nAsh(iS),nBas(jS),
     &                     1.0d0,CMO%SB(iS)%A2,nBas(jS),
     &                           QTmp(1)%SB(jS)%A2,nBas(jS),
     &                     0.0d0,QTmp(2)%SB(jS)%A2,nBas(jS))
              Call DGEMM_('N','N',nBas(jS),nAsh(iS),nBas(jS),
     &                    -1.0d0,Kappa%SB(iS)%A2,nBas(jS),
     &                           QTmp(2)%SB(jS)%A2,nBas(jS),
     &                     1.0d0,QVec%SB(jS)%A2,nBas(jS))
            EndIf
          EndIf
        EndIf
      End Do

      CALL CWTIME(TCINT2,TWINT2)
      tint3(1) = tint3(1) + (TCINT2 - TCINT1)
      tint3(2) = tint3(2) + (TWINT2 - TWINT1)


      If (DoAct) Then
        Call Deallocate_DSBA(QTmp(2))
        Call Deallocate_DSBA(QTmp(1))
        Call Deallocate_G2(MOScr)
      EndIf

      Call mma_deallocate(Faa)
      Call mma_deallocate(Fia)
      Call mma_deallocate(SvShp)
      Call mma_deallocate(iShp_rs)
      Call mma_deallocate(nnBfShp)
      Call mma_deallocate(kOffSh)
      Call mma_deallocate(Indx)
      Call mma_deallocate(SumAClk)
      Call mma_deallocate(MLk)
      Call mma_deallocate(Ylk)
      Call mma_deallocate(AbsC)
      Call Deallocate_NDSBA(DiaH)
#if defined (_MOLCAS_MPP_)
      If (nProcs.gt.1 .and. Update .and. Is_Real_Par())
     &    CALL mma_deallocate(DiagJ)
#endif
      Call mma_deallocate(Diag)

      If (Deco) Then
         Call Deallocate_DSBA(CM(2))
         Call Deallocate_DSBA(CM(1))
      End If

      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1

*
*---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky MCLR timing from '//SECNAM
      Write(6,CFmt)'----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Fock matrix construction        CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ VECTORS                     '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COULOMB                          '
     &                           //'         ',tcoul(1),tcoul(2)
         Write(6,'(2x,A26,2f10.2)')'EXCHANGE                         '
     &                           //'         ',
     &                             tscrn(1)+tmotr(1)+texch(1),
     &                             tscrn(2)+tmotr(2)+texch(2)
         Write(6,'(2x,A26,2f10.2)')'  SCREENING                      '
     &                           //'         ',tscrn(1),tscrn(2)
         Write(6,'(2x,A26,2f10.2)')'  MO TRANSFORM                   '
     &                           //'         ',tmotr(1),tmotr(2)
         Write(6,'(2x,A26,2f10.2)')'  FORMATION                      '
     &                           //'         ',texch(1),texch(2)
         Write(6,'(2x,A26,2f10.2)')'ACTIVE INT, Q AND FOCK MATRIX    '
     &                           //'         ',
     &                     tint1(1)+tint2(1)+tint3(1)+tQmat(1)+tact(1),
     &                     tint1(2)+tint2(2)+tint3(2)+tQmat(2)+tact(2)
         Write(6,'(2x,A26,2f10.2)')'  MO TRANSFORM                   '
     &                           //'         ',tint1(1),tint1(2)
         Write(6,'(2x,A26,2f10.2)')'  INTEGRAL                       '
     &                           //'         ',tint2(1),tint2(2)
         Write(6,'(2x,A26,2f10.2)')'  Q MATRIX                       '
     &                           //'         ',tQmat(1),tQmat(2)
         Write(6,'(2x,A26,2f10.2)')'  ACTIVE FOCK MATRIX             '
     &                           //'         ',tact(1),tact(2)
         Write(6,'(2x,A26,2f10.2)')'  MO BACK TRANSFORM              '
     &                           //'         ',tint3(1),tint3(2)
         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif


c Print the Fock-matrix
#ifdef _DEBUGPRINT_
      if(Debug) then !to avoid double printing in RASSI-debug

      WRITE(6,'(6X,A)')'TEST PRINT FROM '//SECNAM
      WRITE(6,'(6X,A)')
      WRITE(6,'(6X,A)')'***** INACTIVE FOCK MATRIX ***** '
      DO ISYM=1,NSYM
        IF( NBAS(ISYM).GT.0 ) THEN
          WRITE(6,'(6X,A)')
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          call CHO_OUTPUT(FkI%SB(ISYM)%A2,1,NBAS(ISYM),1,NBAS(ISYM),
     &                    NBAS(ISYM),NBAS(ISYM),1,6)
        ENDIF
      END DO

      endif

#endif

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(nIsh)
      END

**************************************************************
