!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990,1991,1993,1998, Roland Lindh                      *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drv2El(Integral_WrOut,ThrAO)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals.                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Modified for k2 loop. August '91                         *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             Modified driver. Jan. '98                                *
!***********************************************************************

use iSD_data, only: iSD
use Basis_Info, only: dbsc
use Real_Info, only: CutInt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Eight
use Definitions, only: wp, iwp

implicit none
external :: Integral_WrOut
real(kind=wp), intent(in) :: ThrAO
real(kind=wp) :: A_int, Disc, Dix_Mx, ExFac(1), P_Eff, PP_Count, PP_Eff, PP_Eff_delta, S_Eff, ST_Eff, T_Eff, TCpu1, TCpu2, Thize, &
                 TMax_all, TskCount, TskHi, TskLw, TWall1, Twall2
integer(kind=iwp) :: iCnttp, ijS, iOpt, iS, iTOffs(8,8,8), jCnttp, jS, kCnttp, klS, kS, lCnttp, lS, nij, Nr_Dens, nSkal
logical(kind=iwp) :: Verbose, Indexation, FreeK2, W2Disc, PreSch, DoIntegrals, DoFock, DoGrad, FckNoClmb(1), FckNoExch(1), &
                     Triangular
character(len=72) :: SLine
real(kind=wp), allocatable :: Dens(:), Fock(:), TInt(:), TMax(:,:)
integer(kind=iwp), parameter :: nTInt = 1, mDens = 1
integer(kind=iwp), allocatable :: Pair_Index(:,:)
logical(kind=iwp), external :: Rsv_GTList

!                                                                      *
!***********************************************************************
!                                                                      *
SLine = 'Computing 2-electron integrals'
call StatusLine(' Seward:',SLine)
!                                                                      *
!***********************************************************************
!                                                                      *
ExFac = One
Nr_Dens = 1
DoIntegrals = .true.
DoFock = .false.
DoGrad = .false.
FckNoClmb = .false.
FckNoExch = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Basis_Mode('Valence')
call Setup_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize for 2-electron integral evaluation. Do not generate
! tables for indexation.

Indexation = .false.
call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
!                                                                      *
!***********************************************************************
!                                                                      *
Thize = Zero               ! Not used for conventional integrals
PreSch = .true.            ! Not used for conventional integrals

Disc = Zero
Dix_Mx = Zero
TskHi = Zero
TskLw = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute entities for prescreening at shell level

call mma_allocate(TMax,nSkal,nSkal)
call Shell_MxSchwz(nSkal,TMax)
TMax_all = Zero
do iS=1,nSkal
  do jS=1,iS
    TMax_all = max(TMax_all,TMax(iS,jS))
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of non-vanishing pairs

call mma_allocate(Pair_Index,2,nSkal*(nSkal+1)/2)
nij = 0
do iS=1,nSkal
  do jS=1,iS
    if (TMax_All*TMax(iS,jS) >= CutInt) then
      nij = nij+1
      Pair_Index(1,nij) = iS
      Pair_Index(2,nij) = jS
    end if
  end do
end do
P_Eff = real(nij,kind=wp)
!                                                                      *
!***********************************************************************
!                                                                      *
Triangular = .true.
call Init_TList(Triangular,P_Eff)
call Init_PPList()
call Init_GTList()
iOpt = 0

PP_Eff = P_Eff**2
PP_Eff_delta = 0.1_wp*PP_Eff
PP_Count = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu1,TWall1)

! big loop over individual tasks distributed over individual nodes

do
  ! make reservations of a tesk in global task list and get task range
  ! in return. Function will be false if no more tasks to execute.
  if (.not. Rsv_GTlist(TskLw,TskHi,iOpt,W2Disc)) exit
  W2Disc = .false.

  ! Now do a quadruple loop over shells

  ijS = int((One+sqrt(Eight*TskLw-Three))/Two)
  iS = Pair_Index(1,ijS)
  jS = Pair_Index(2,ijS)
  klS = int(TskLw-real(ijS,kind=wp)*(real(ijS,kind=wp)-One)/Two)
  kS = Pair_Index(1,klS)
  lS = Pair_Index(2,klS)
  TskCount = TskLw

  do while (TskCount-TskHi <= 1.0e-10_wp)

    ! Logic to avoid computing integrals in a mixed muonic and
    ! electronic basis.

    iCnttp = iSD(13,iS)
    jCnttp = iSD(13,jS)
    if (dbsc(iCnttp)%fMass == dbsc(jCnttp)%fMass) then
      kCnttp = iSD(13,kS)
      lCnttp = iSD(13,lS)
      if (dbsc(kCnttp)%fMass == dbsc(lCnttp)%fMass) then

        S_Eff = real(ijS,kind=wp)
        T_Eff = real(klS,kind=wp)
        ST_Eff = S_Eff*(S_Eff-One)/Two+T_Eff
        if (ST_Eff >= PP_Count) then
          write(SLine,'(A,F5.2,A)') 'Computing 2-electron integrals,',ST_Eff/PP_Eff*100.0_wp,'% done so far.'
          call StatusLine(' Seward:',SLine)
          PP_Count = PP_Count+PP_Eff_delta
        end if

        A_int = TMax(iS,jS)*TMax(kS,lS)
        if (A_Int >= CutInt) then
          ! from Dens are dummy arguments
          call mma_allocate(TInt,nTInt,label='TInt')
          call mma_allocate(Dens,mDens,label='Dens')
          call mma_allocate(Fock,mDens,label='Fock')
          call Eval_Ints_New_Inner(iS,jS,kS,lS,TInt,nTInt,iTOffs,Integral_WrOut,Dens,Fock,mDens,ExFac,Nr_Dens,FckNoClmb,FckNoExch, &
                                   Thize,W2Disc,PreSch,Dix_Mx,Disc,TskCount,DoIntegrals,DoFock)
          call mma_deallocate(TInt)
          call mma_deallocate(Dens)
          call mma_deallocate(Fock)
        end if
      end if
    end if
    TskCount = TskCount+One
    if (TskCount-TskHi > 1.0e-10_wp) exit
    klS = klS+1
    if (klS > ijS) then
      ijS = ijS+1
      klS = 1
    end if
    iS = Pair_Index(1,ijS)
    jS = Pair_Index(2,ijS)
    kS = Pair_Index(1,klS)
    lS = Pair_Index(2,klS)
  end do

  ! Use a time slot to save the number of tasks and shell
  ! quadrupltes process by an individual node
  call SavStat(1,One,'+')
  call SavStat(2,TskHi-TskLw+One,'+')
end do
! End of big task loop
call CWTime(TCpu2,TWall2)
call SavTim(1,TCpu2-TCpu1,TWall2-TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_GTList()
call Free_PPList()
call Free_TList()

call mma_deallocate(Pair_Index)
call mma_deallocate(TMax)
!                                                                      *
!***********************************************************************
!                                                                      *
! Terminate integral environment.

Verbose = .false.
FreeK2 = .true.
call Term_Ints(Verbose,FreeK2)
call Free_iSD()

return

end subroutine Drv2El
