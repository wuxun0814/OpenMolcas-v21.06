!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

Module Interfaces_SCF

Implicit None
Private
Public :: dOne_SCF, MinDns, OccDef, PMat_SCF, TraClc_i, vOO2OV

Interface

  SubRoutine dOne_SCF(nSym,nBas,nOrb,nFro,CMO,nCMO,Occ,Dlt,alpha_density)
    Integer :: nSym,nCMO,nBas(nSym),nOrb(nSym),nFro(nSym)
    Real*8, Target:: CMO(nCMO), Occ(*), Dlt(*)
    Logical :: alpha_density
  End SubRoutine dOne_SCF

  SubRoutine MinDns(Dens,mBT,NumD,XCff,ltXCff,nD)
    Integer :: mBT,NumD,ltXCff,nD
    Real*8 :: XCff(ltXCff,nD)
    Real*8, Target :: Dens(mBT,nD,NumD)
  End SubRoutine MinDns

  SubRoutine OccDef(Occ,mmB,nD,CMO,mBB)
    Integer :: mmB,nD,mBB
    Real*8 :: Occ(mmB,nD)
    Real*8, Target :: CMO(mBB,nD)
  End SubRoutine OccDef

  SubRoutine PMat_SCF(Dens,OneHam,TwoHam,nDT,NumDT,nXCf,FstItr,XCf,nD,E_DFT,nE_DFT,Vxc,Fock)
    Integer :: nDT,NumDT,nXCf,nD,nE_DFT
    Real*8, Target :: Dens(nDT,nD,NumDT)
    Real*8 :: OneHam(nDT),XCf(nXCf,nD),E_DFT(nE_DFT),Vxc(nDT,nD,NumDT),Fock(nDT,nD)
    Real*8, Target :: TwoHam(nDT,nD,NumDT)
    Logical :: FstItr
  End SubRoutine PMat_SCF

  SubRoutine TraClc_i(OneHam,Dens,TwoHam,Vxc,nDT,NumDT,iterLw,TrDh,TrDP,TrDD,nTr,nD)
    Integer :: nDT,NumDT,iterLw,nTr,nD
    Real*8 :: OneHam(nDT),TrDh(nTr,nTr,nD),TrDP(nTr,nTr,nD),TrDD(nTr,nTr,nD)
    Real*8, Target :: Dens(nDT,nD,NumDT),TwoHam(nDT,nD,NumDT),Vxc(nDT,nD,NumDT)
  End SubRoutine TraClc_i

  SubRoutine vOO2OV(v1,n1,v2,n2,nD)
    Integer :: n1,n2,nD
    Real*8, Target :: v1(n1,nd),v2(n2,nD)
  End SubRoutine vOO2OV

End Interface

End Module Interfaces_SCF
