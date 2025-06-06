************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine DrvDFT(h1,TwoHam,D,RepNuc,nh1,First,Dff,
     &                  lRF,KSDFT,ExFac,Do_Grad,Grad,nGrad,iSpin,
     &                  D1I,D1A,nD1,DFTFOCK)
      use KSDFT_Info, only: KSDFA, funcaa, funcbb, funccc
      Implicit Real*8 (a-h,o-z)
      External LSDA, Overlap, BLYP, BPBE, B3LYP, HFS, HFB,
     &         XAlpha, LSDA5, B3LYP5, B2PLYP, TLYP, NLYP,
     &         NucAtt, NEWF, NEWF1, OLYP, O3LYP, OPBE,
     &         PBE, PBE0, PBEsol, M06L, M06, M062X, HFO,
     &         M06HF, Checker, SSBSW, SSBD, HFG, GLYP, GPBE,
     &         HFB86, B86LYP, B86PBE, BWIG, KT3,
     &         O2PLYP,  KT2,  RGE2, REVPBE,
     &         PTCA,S12G, S12H
#include "real.fh"
#include "stdalloc.fh"
#include "nq_info.fh"
#include "debug.fh"
#include "pamint.fh"
#include "ksdft.fh"
      Real*8 h1(nh1), TwoHam(nh1), D(nh1,2), Grad(nGrad), Vxc_ref(2)
      Real*8 D1I(nD1),D1A(nD1)
      Logical First, Dff, lRF,  Do_Grad
      Logical Do_MO,Do_TwoEl, Found
      Character*(*) KSDFT
      Character*4 DFTFOCK
      Real*8, Allocatable:: D_DS(:,:), F_DFT(:,:)
*
      KSDFA = KSDFT
      lKSDFT=LEN(KSDFT)
      Debug=.False.
*                                                                      *
************************************************************************
*                                                                      *
c     Call SetQue('Trace=on')
*                                                                      *
************************************************************************
*                                                                      *
      Call Put_iScalar('Multiplicity',iSpin)
      Call Get_iScalar('nSym',mIrrep)
      Call Get_iArray('nBas',mBas,mIrrep)
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*                                                                      *
      Call Get_dScalar('DFT exch coeff',CoefX)
      Call Get_dScalar('DFT corr coeff',CoefR)
*
************************************************************************
*                                                                      *
      If (Do_Grad) Call FZero(Grad,nGrad)
*                                                                      *
************************************************************************
*                                                                      *
      If (iSpin.eq.1) Then
         nD=1
      Else
         nD=2
      End If
*
*     What is this?
*
      If (DFTFOCK.eq.'DIFF') nD=2
      If (DFTFOCK.eq.'ROKS') nD=2
      Call mma_allocate(D_DS,nh1,nD,Label='D_DS')
*
*---- Get the total density
*
      Call Get_D1ao(D_DS,nh1)
*     Call RecPrt('D1ao',' ',D_DS(:,1),nh1,1)
*
*
*---- Get the spin density
*
      If (nD.ne.1) Then
         Call Get_D1Sao(D_DS(:,2),nh1)
*        Call RecPrt('D1Sao',' ',D_DS(:,2),nh1,1)
      End If
*
*---- Compute alpha and beta densities
*
*     Call RecPrt('DTot',' ',D_DS(:,1),nh1,1)
*     Call RecPrt('DSpn',' ',D_DS(:,2),nh1,1)
      If (nD.eq.1) Then
        D_DS(:,1)=Half*D_DS(:,1)
      Else
         Do i = 1, nh1
            DTot=D_DS(i,1)
            DSpn=D_DS(i,2)
            d_Alpha=Half*(DTot+DSpn)
            d_Beta =Half*(DTot-DSpn)
            D_DS(i,1)=d_Alpha
            D_DS(i,2)=d_Beta
         End Do
      End If
*     Call RecPrt('Da',' ',D_DS(:,1),nh1,1)
*     Call RecPrt('Db',' ',D_DS(:,2),nh1,1)
*
      If(KSDFT(1:3).ne.'SCF') Then
        Call Get_iArray('nIsh',nIsh,mIrrep)
        Call Get_iArray('nFro',nFro,mIrrep)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     DFT functionals, compute integrals over the potential
*
      Func            =Zero
      Funcaa          =Zero
      Funcbb          =Zero
      Funccc          =Zero
      Dens_I          =Zero
      Dens_a1         =Zero
      Dens_b1         =Zero
      Dens_a2         =Zero
      Dens_b2         =Zero
      Dens_t1         =Zero
      Dens_t2         =Zero
      Grad_I          =Zero
      Tau_I           =Zero
      Do_MO           =.False.
      Do_TwoEl        =.False.
*
      If (nD.eq.2.and.DFTFOCK.eq.'DIFF') Then
         numAO=0
         If(KSDFT(1:3).ne.'SCF') Then
           Do iIrrep=0,mIrrep-1
             nAsh(iIrrep)=0
           End Do
           Call qpg_iArray('nAsh',Found,nData)
           If(Found .and. nData.eq.mIrrep) Then
             Call Get_iArray('nAsh',nAsh,mIrrep)
           End If
           Do iIrrep=0,mIrrep-1
             numAO=numAO+nAsh(iIrrep)
           End Do
        End If
        If (numAO.ne.0)
     &  Do_TwoEl        =.True.
        Do_MO           =.True.

      End If
*
*     nFckDim: number of different types of Fock matrices. Normally for
*     conventional functionals we have one Fock matrix for closed shell
*     calculations and two (F_alpha and F_beta) for open shell systems.
*     For CASDFT we have always two (F_inactive and F_active)
*
************************************************************************
*                                                                      *
*      LSDA LDA SVWN  GLM stuff                                        *
*                                                                      *
       If (KSDFT.eq.'LSDA ' .or.
     &     KSDFT.eq.'LDA '  .or.
     &     KSDFT.eq.'TLSDA'  .or. !GLM
     &     KSDFT.eq.'FTLSDA'  .or. !AMS
     &     KSDFT.eq.'SVWN ') Then
         If(KSDFT.eq.'TLSDA'
     &     .or.KSDFT.eq.'FTLSDA') Do_MO=.true. !GLM
         If(KSDFT.eq.'TLSDA'
     &     .or.KSDFT.eq.'FTLSDA') Do_TwoEl=.true. !GLM
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(LSDA   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
c          write(6,*) 'Func in drvdft :', Func
*                                                                      *
************************************************************************
*                                                                      *
*      LSDA5 LDA5 SVWN5                                                *
*                                                                      *
       Else If (KSDFT.eq.'LSDA5' .or.
     &          KSDFT.eq.'LDA5'  .or.
     &          KSDFT.eq.'TLSDA5 '.or.
     &          KSDFT.eq.'SVWN5') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(LSDA5  ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
c         write(6,*) 'Func in drvdft :', Func
*                                                                      *
************************************************************************
*                                                                      *
*     HFB                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFB') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(HFB    ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     HFO                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFO') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(HFO    ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     HFG                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFG') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(HFG    ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     HFB86                                                            *
*                                                                      *
       Else If (KSDFT.eq.'HFB86') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(HFB86   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      HFS                                                             *
*                                                                      *
       Else If (KSDFT.eq.'HFS') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(HFS    ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      XALPHA                                                          *
*                                                                      *
       Else If (KSDFT.eq.'XALPHA') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(XAlpha ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     Overlap                                                          *
*                                                                      *
      Else If (KSDFT.eq.'Overlap') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(Overlap,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     NucAtt                                                           *
*                                                                      *
      Else If (KSDFT.eq.'NucAtt') Then
         ExFac=One
         Functional_type=LDA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(NucAtt,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     BWIG                                                             *
*                                                                      *
      Else If (KSDFT.eq.'BWIG') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(BWIG   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     BLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'BLYP'
     &       .or.  KSDFT.eq.'TBLYP' !GLM
     &       .or.  KSDFT.eq.'FTBLYP' !AMS
     &        ) Then
       If(KSDFT.eq.'TBLYP'
     &       .or. KSDFT.eq.'FTBLYP') Do_MO=.true.
       If(KSDFT.eq.'TBLYP'
     &       .or. KSDFT.eq.'FTBLYP') Do_TwoEl=.true.
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(BLYP   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     OLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'OLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(OLYP   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     KT3                                                              *
*                                                                      *
      Else If (KSDFT.eq.'KT3') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(KT3   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     KT2                                                              *
*                                                                      *
      Else If (KSDFT.eq.'KT2') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(KT2   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     GLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'GLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(GLYP   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B86LYP                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B86LYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(B86LYP   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     BPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'BPBE') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(BPBE   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     OPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'OPBE'
     &     .or.KSDFT.eq.'TOPBE' !GLM
     &     .or.KSDFT.eq.'FTOPBE'!AMS
     &         ) then
         If(KSDFT.eq.'TOPBE'
     &     .or.KSDFT.eq.'FTOPBE') Do_MO=.true.
         If(KSDFT.eq.'TOPBE'
     &     .or.KSDFT.eq.'FTOPBE') Do_TwoEl=.true.
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(OPBE   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     GPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'GPBE') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(GPBE   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B86PBE                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B86PBE') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(B86PBE  ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     TLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'TLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(TLYP   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     NLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'NLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(NLYP   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP                                                            *
*                                                                      *
      Else If (KSDFT.eq.'B3LYP ') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(B3LYP  ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     O3LYP                                                            *
*                                                                      *
      Else If (KSDFT.eq.'O3LYP ') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(O3LYP  ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B2PLYP                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B2PLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(B2PLYP  ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     O2PLYP                                                           *
*                                                                      *
      Else If (KSDFT.eq.'O2PLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(O2PLYP  ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP5                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B3LYP5') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(B3LYP5 ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     PBE                                                              *
*                                                                      *
      Else If (KSDFT.eq.'PBE'
     &     .or.KSDFT.eq.'TPBE' !GLM
     &     .or.KSDFT.eq.'FTPBE'!AMS
     &         ) then
         If(KSDFT.eq.'TPBE'
     &     .or.KSDFT.eq.'FTPBE') Do_MO=.true.
         If(KSDFT.eq.'TPBE'
     &     .or.KSDFT.eq.'FTPBE') Do_TwoEl=.true.
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(PBE   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     revPBE                                                           *
*                                                                      *
      Else If (KSDFT.eq.'REVPBE'
     &     .or.KSDFT.eq.'TREVPBE' !GLM
     &     .or.KSDFT.eq.'FTREVPBE'!AMS
     &         ) then
         If(KSDFT.eq.'TREVPBE'
     &     .or.KSDFT.eq.'FTREVBPE') Do_MO=.true.
         If(KSDFT.eq.'TREVPBE'
     &     .or.KSDFT.eq.'FTREVPBE') Do_TwoEl=.true.
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(REVPBE,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     SSBSW                                                              *
*                                                                      *
      Else If (KSDFT.eq.'SSBSW'
     &     .or.KSDFT.eq.'TSSBSW') Then !GLM
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(SSBSW   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     SSBD                                                             *
*                                                                      *
      Else If (KSDFT.eq.'SSBD'
     &     .or.KSDFT.eq.'TSSBD') Then !GLM
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(SSBD ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     S12H                                                             *
*                                                                      *
      Else If (KSDFT.eq.'S12H') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(S12H  ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     S12G                                                             *
*                                                                      *
      Else If (KSDFT.eq.'S12G'
     &     .or.KSDFT.eq.'TS12G') Then !GLM
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(S12G  ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     PBEsol                                                           *
*                                                                      *
      Else If (KSDFT.eq.'PBESOL') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(PBEsol   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     RGE2                                                             *
*                                                                      *
      Else If (KSDFT.eq.'RGE2') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(RGE2   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     PTCA                                                             *
*                                                                      *
      Else If (KSDFT.eq.'PTCA') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(PTCA   ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     PBE0                                                             *
*                                                                      *
      Else If (KSDFT.eq.'PBE0') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(PBE0  ,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06-L                                                            *
*                                                                      *
      Else If (KSDFT.eq.'M06L') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(M06L,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06                                                              *
*                                                                      *
      Else If (KSDFT.eq.'M06 ') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(M06,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06-2X                                                           *
*                                                                      *
      Else If (KSDFT.eq.'M062X') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(M062X,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06-HF                                                           *
*                                                                      *
      Else If (KSDFT.eq.'M06HF') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(M06HF,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     Checker                                                          *
*                                                                      *
      Else If (KSDFT.eq.'CHECKER') Then
         ExFac=Zero
         Functional_type=meta_GGA_type2
         nFckDim = nD
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         Call DrvNQ(Checker,F_DFT,nFckDim,Func,
     &              D_DS,nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
*     CASDFT functionals:                                              *
*                                                                      *
      Else If (KSDFT(1:4).eq.'NEWF') Then
*                                                                      *
*        These functionals are still under construction.               *
*        The code, written by S.G. & C., will be now modified by       *
*                      Giovanni Ghigo (CGG)                            *
*                                                                      *
         If (DFTFOCK.ne.'DIFF') Then
            Call WarningMessage(2,
     &                 ' This is CASDFT type functional !!!;'
     &               //' You cannot use it in this calculation.')
            Call Quit_OnUserError()
         End If
         Do_Grad  = .False.
         Do_MO    = .True.
         Do_TwoEl = .True.
*
         ExFac=Get_ExFac(KSDFT)
         Functional_type=CASDFT_type
         nFckDim = 2
         Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
         F_DFT(:,:)=Zero
         If ( KSDFT(5:5).eq.'0' )
     &      Call DrvNQ(NEWF ,F_DFT,nFckDim,Func,
     &                 D_DS,nh1,nD,Do_Grad,
     &                 Grad,nGrad,
     &                 Do_MO,Do_TwoEl,DFTFOCK)
         If ( KSDFT(5:5).eq.'1' )
     &      Call DrvNQ(NEWF1 ,F_DFT,nFckDim,Func,
     &                 D_DS,nh1,nD,Do_Grad,
     &                 Grad,nGrad,
     &                 Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
      Else
         Call WarningMessage(2,
     &               ' DrvDFT: Undefined functional type!')
         Write (6,*) '         Functional=',KSDFT(1:lKSDFT)
         Call Quit_OnUserError()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Energy_integrated=Func
*                                                                      *
************************************************************************
*                                                                      *
      If (KSDFT.eq.'Overlap'.or.KSDFT.eq.'NucAtt') Then
         call dcopy_(nh1,F_DFT,1,h1,1)
         If (KSDFT.eq.'NucAtt') Energy_integrated=Func
      Else
*
*        Put out the integrated DFT energy and the DFT Fock matrices
*        on the RUNFILE
*
*        Call Put_DFT_Energy(Energy_integrated)
         Call Poke_dScalar('KSDFT energy',Energy_integrated)
         Call Put_dScalar('CASDFT energy',Energy_integrated)
         Call Put_dExcdRa(F_DFT,nFckDim*nh1)
*         Write(6,'(a,f22.16)') " Energy in drvdft ",Energy_integrated
#ifdef _DEBUGPRINT_
         Write(6,'(a,f22.16)') " Energy ",Energy_integrated
         If (nFckDim.eq.1) Then
            Do i=1,nh1
               Write(6,'(i4,f22.16)') i,F_DFT(i,1)
            End Do
         Else
            Do i=1,nh1
              Write(6,'(i4,3f22.16)') i,F_DFT(i,1),
     &                                  F_DFT(i,2),
     &        F_DFT(i,1)+F_DFT(i,2)/2.0d0
            End Do
         End If
#endif

*
*        In the SCF program (traclc.f) the program computes the trace
*        of the one-electron hamiltonian over a set of densities. The
*        DFT contribution is not linear with respect to variations of
*        the density. However, with the following term we can include
*        the linear component in that code.
*
         Fact = Two
         If (nD.ne.1) Fact=One
         Vxc_ref(1)=Fact*DDot_(nh1,F_DFT(:,1),1,D_DS,1)
         If (nD.ne.1) Then
           Vxc_ref(2)=DDot_(nh1,F_DFT(:,2),1,D_DS(:,2),1)
         Else
            Vxc_ref(2)=Zero
         End If
         Call Put_Temp('Vxc_ref ',Vxc_ref,2)
      End If
*
      Call mma_deallocate(F_DFT)
      Call mma_deallocate(D_DS)
      Call Free_iSD()
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(TwoHam)
         Call Unused_real_array(D)
         Call Unused_real(RepNuc)
         Call Unused_logical(First)
         Call Unused_logical(Dff)
         Call Unused_logical(lRF)
         Call Unused_real_array(D1I)
         Call Unused_real_array(D1A)
      End If
      End
