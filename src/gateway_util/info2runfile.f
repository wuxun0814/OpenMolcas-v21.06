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
* Copyright (C) 2006, Roland Lindh                                     *
*               2019, Ignacio Fdez. Galvan                             *
************************************************************************
      SubRoutine Info2Runfile()
************************************************************************
*                                                                      *
*     Object: dump misc. information to the runfile.                   *
*                                                                      *
*     Author: Roland Lindh, Dept Chem. Phys., Lund University, Sweden  *
*             October 2006                                             *
************************************************************************
      use Period
      use Basis_Info
      use Center_Info
      use external_centers, only: iXPolType, XF
      use Temporary_Parameters, only: Expert, DirInt
      use Sizes_of_Seward, only: S
      use RICD_Info, only: Do_RI, Cholesky, Cho_OneCenter, LocalDF
      use Real_Info, only: CoC, CoM
      use Logical_Info, only: DoFMM
      use Symmetry_Info, only: nIrrep, VarR, VarT
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cholesky.fh"
#include "real.fh"
#include "rctfld.fh"
#include "stdalloc.fh"
#include "print.fh"
      Character xLblCnt(MxAtom)*(LENIN)
      Logical Pseudo
      Integer nDel(8)
      Logical Found
*
#include "embpcharg.fh"
      Real*8, Dimension(:,:), Allocatable :: DCo
      Real*8, Dimension(:), Allocatable :: DCh, DCh_Eff
      Integer, Allocatable :: NTC(:), ICh(:), IsMM(:), nStab(:)
************************************************************************
*                                                                      *
      Call ICopy(8,[0],0,nDel,1)
      Call Put_iArray('nFro',nDel,nIrrep) ! put to 0
      Call qpg_iArray('nDel',Found,nData)
      If (.Not.Found) Then
         Call Put_iArray('nDel',nDel,nIrrep)
      End If
      Do i=1,nIrrep ! note that in info.fh is nBas(0:7)
         nDel(i)=nBas(i-1)-nDel(i)
      End Do
      Call Put_iArray('nOrb',nDel,nIrrep) ! nDel is corrupted here!
*
*     VarR and VarT
*
      If (lRF.and..not.PCM) VarT=.True.
      Pseudo=.False.
      Do iCnttp = 1, nCnttp
         Pseudo = Pseudo .or. (dbsc(iCnttp)%pChrg .and.
     &                         dbsc(iCnttp)%Fixed)
      End Do
      If (.not.DoEMPC) Then
         If (Allocated(XF).or.Pseudo) Then
            VarR=.True.
            VarT=.True.
         End If
      End If
*
*     Manipulate the option flag
*
      iOption=0
      If (DirInt) iOption=iOr(iOption,1)
      If (Expert) iOption=iOr(iOption,2)
      If (lRF) iOption=iOr(iOption,4)
      If (lLangevin.or.iXPolType.gt.0) iOption=iOr(iOption,8)
      If (PCM) Then
         iOption=iOr(iOption,16)
         nPCM_Info=0
         Call Put_iScalar('PCM info length',nPCM_Info)
      End If
      iOption=iOr(iOption,32)
*  2el-integrals from the Cholesky vectors
      If (Cholesky.or.Do_RI) iOption=iOr(iOption,2**9)
*  RI-Option
      If (Do_RI) Then
         iOption=iOr(iOption,2**10)
*  Local or non-local
         If (LocalDF) Then
            Call Put_LDFAccuracy()
            Call Put_LDFConstraint()
            iLocalDF=1
         Else
            iLocalDF=0
         End If
         Call Put_iScalar('DF Mode',iLocalDF)
      End If
*  1C-CD
      If (Cholesky.and.Cho_1Center) iOption=iOr(iOption,2**12)
      Cho_OneCenter=Cho_1Center
      Call Put_iScalar('System BitSwitch',iOption)

      Call Put_iScalar('Highest Mltpl',S%nMltpl)
      iGO = 0
      Call Put_iScalar('Grad ready',iGO)
      iFMM = 0
      If (DoFMM) iFMM = 1
      Call Put_iScalar('FMM',iFMM)
*
      Call qpg_iScalar('Saddle Iter',Found)
      If (.NOT.Found) Then
         iter_S=0
         Call Put_iScalar('Saddle Iter',iter_S)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Generate list of unique centers (atoms+pseudo)
*
      nNuc = 0
      Do iCnttp = 1, nCnttp
         If (.Not.dbsc(iCnttp)%Frag.and.
     &       .Not.dbsc(iCnttp)%Aux) nNuc = nNuc + dbsc(iCnttp)%nCntr
      End Do
*
      Call mma_allocate(DCo,3,nNuc,label='DCo')
      Call mma_allocate(ICh,nNuc,label='ICh')
      Call mma_allocate(DCh_Eff,nNuc,label='DCh_Eff')
      mdc = 0
      iNuc = 0
      Do iCnttp = 1, nCnttp
         If (.Not.dbsc(iCnttp)%Frag.and.
     &       .Not.dbsc(iCnttp)%Aux) Then
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               mdc = mdc + 1
               iNuc = iNuc+ 1
               DCo(1:3,iNuc)=dbsc(iCnttp)%Coor(1:3,iCnt)
               DCh_Eff(iNuc)=dbsc(iCnttp)%Charge
               ICh(iNuc)    =dbsc(iCnttp)%AtmNr
               xLblCnt(iNuc)=dc(mdc)%LblCnt(1:LENIN)
            End Do
         Else
            mdc  = mdc + dbsc(iCnttp)%nCntr
         End If
      End Do
*
      Call Put_iScalar('Unique centers',nNuc)
      Call Put_dArray('Un_cen Coordinates',DCo,3*nNuc)
      Call Put_iArray('Un_cen charge',ICh,nNuc)
      Call Put_dArray('Un_cen effective charge',DCh_Eff,nNuc)
      Call Put_cArray('Un_cen Names',xLblCnt(1),(LENIN)*nNuc)
*
      Call mma_deallocate(DCh_Eff)
      Call mma_deallocate(ICh)
      Call mma_deallocate(DCo)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate list of unique atoms
*
      nNuc = 0
      Do iCnttp = 1, nCnttp
         If (.Not.dbsc(iCnttp)%pChrg.and.
     &       .Not.dbsc(iCnttp)%Frag.and.
     &       .Not.dbsc(iCnttp)%Aux) nNuc = nNuc + dbsc(iCnttp)%nCntr
      End Do
*
      Call mma_allocate(DCo,3,nNuc,label='DCo')
      Call mma_allocate(DCh,nNuc,label='DCh')
      Call mma_allocate(DCh_Eff,nNuc,label='DCh_Eff')
      Call mma_allocate(nStab,nNuc,label='nStab')
      mdc = 0
      iNuc = 0
      Do iCnttp = 1, nCnttp
         If (.Not.dbsc(iCnttp)%pChrg.and.
     &       .Not.dbsc(iCnttp)%Frag.and.
     &       .Not.dbsc(iCnttp)%Aux) Then
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               mdc = mdc + 1
               iNuc = iNuc+ 1
               DCo(1:3,iNuc)=dbsc(iCnttp)%Coor(1:3,iCnt)
               DCh_Eff(iNuc)=dbsc(iCnttp)%Charge
               DCh(iNuc)    =DBLE(dbsc(iCnttp)%AtmNr)
               xLblCnt(iNuc)=dc(mdc)%LblCnt(1:LENIN)
               nStab(iNuc)=dc(mdc)%nStab
            End Do
         Else
            mdc  = mdc + dbsc(iCnttp)%nCntr
         End If
      End Do
*
      Call Put_iScalar('Unique atoms',nNuc)
*
      Call mma_allocate(IsMM,SIZE(dbsc),Label='IsMM')
      Do i = 1, SIZE(dbsc)
         IsMM(i)=dbsc(i)%IsMM
      End Do
      Call Put_iArray('IsMM',IsMM,SIZE(dbsc))
      Call mma_deallocate(IsMM)
*
      Call Put_dArray('Unique Coordinates',DCo,3*nNuc)
      Call Put_dArray('Center of Mass',CoM,3)
      Call Put_dArray('Center of Charge',CoC,3)
      Call Put_dArray('Nuclear charge',DCh,nNuc)
      Call Put_dArray('Effective nuclear Charge',DCh_Eff,nNuc)
      Call Put_cArray('Unique Atom Names',xLblCnt(1),(LENIN)*nNuc)
      Call Put_iArray('nStab',nStab,nNuc)
      If (Allocated(XF)) Call Put_iScalar('nXF',nXF)
      If (Cell_l) Then
         Call Put_dArray('Unit Cell Vector',VCell,9)
         Call Put_iArray('Spread of Coord.',ispread,3)
         Call Put_iScalar('Unit Cell NAtoms',lthCell)
         Call Put_iArray('Unit Cell Atoms',AdCell,lthCell)
      End IF
*
*     Initiate entry to zero.
      call dcopy_(nNuc,[0.0D0],0,DCh_Eff,1)
      Call Put_dArray('Mulliken Charge',DCh_Eff,nNuc)
*
      Call mma_deallocate(nStab)
      Call mma_deallocate(DCh_Eff)
      Call mma_deallocate(DCh)
      Call mma_deallocate(DCo)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate a translation array iNuc -> iCnttp (for espf)
*
      Call mma_allocate(NTC,nNuc,label='NTC')
*
      iNTC = 0
      Do iCnttp = 1, nCnttp
         If (.Not.dbsc(iCnttp)%pChrg.and.
     &       .Not.dbsc(iCnttp)%Frag.and.
     &       .Not.dbsc(iCnttp)%Aux) Then
            Do iNuc = 1, dbsc(iCnttp)%nCntr
               NTC(iNTC+iNuc) = iCnttp
            End Do
            iNTC = iNTC + dbsc(iCnttp)%nCntr
         End If
      End Do
*
      Call Put_iArray('Atom -> Basis',NTC,nNuc)
      Call mma_deallocate(NTC)
*                                                                      *
************************************************************************
*                                                                      *
*     Coordinate list as above but for all centers with proper basis
*     functions.
*
      nNuc = 0
      Do iCnttp = 1, nCnttp
         If (.Not.dbsc(iCnttp)%Frag.and.
     &       .Not.dbsc(iCnttp)%Aux) nNuc = nNuc + dbsc(iCnttp)%nCntr
      End Do
*
      Call mma_allocate(DCo,3,nNuc,label='DCo')
      iNuc = 0
      Do iCnttp = 1, nCnttp
         If (.Not.dbsc(iCnttp)%Frag.and.
     &       .Not.dbsc(iCnttp)%Aux) Then
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               iNuc = iNuc+ 1
               DCo(1:3,iNuc)=dbsc(iCnttp)%Coor(1:3,iCnt)
            End Do
         End If
      End Do
*
      Call Put_iScalar('Bfn atoms',nNuc)
      Call Put_dArray('Bfn Coordinates',DCo,3*nNuc)
*
      Call mma_deallocate(DCo)
*                                                                      *
************************************************************************
*                                                                      *
*     Coordinate list as above but for only pseudo centers
*
      nNuc = 0
      Do iCnttp = 1, nCnttp
         If (dbsc(iCnttp)%pChrg) nNuc = nNuc + dbsc(iCnttp)%nCntr
      End Do
*
      Call mma_allocate(DCo,3,nNuc,label='DCo')
      Call mma_allocate(DCh,nNuc,label='DCh')
      iNuc = 0
      Do iCnttp = 1, nCnttp
         If (dbsc(iCnttp)%pChrg)Then
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               iNuc = iNuc+ 1
               DCo(1:3,iNuc)=dbsc(iCnttp)%Coor(1:3,iCnt)
               DCh(iNuc)   = DBLE(dbsc(iCnttp)%AtmNr)
            End Do
         Else
         End If
      End Do
*
      Call Put_iScalar('Pseudo atoms',nNuc)
      Call Put_dArray('Pseudo Coordinates',DCo,3*nNuc)
      Call Put_dArray('Pseudo charge',DCh,nNuc)
*
      Call mma_deallocate(DCh)
      Call mma_deallocate(DCo)
*                                                                      *
************************************************************************
*                                                                      *
      Call Mk_ChDisp()
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
      Subroutine Put_LDFAccuracy()
      Implicit None
#include "localdf.fh"
      Call Put_dScalar('LDF Accuracy',Thr_Accuracy)
      End
      Subroutine Put_LDFConstraint()
      Implicit None
#include "localdf.fh"
      Call Put_iScalar('LDF Constraint',LDF_Constraint)
      End
