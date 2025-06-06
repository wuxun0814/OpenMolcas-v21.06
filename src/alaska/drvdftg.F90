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
! Copyright (C) 2002, Roland Lindh                                     *
!***********************************************************************

subroutine DrvDFTg(Grad,Temp,nGrad)
!***********************************************************************
!                                                                      *
! Object: driver for computation of gradient with respect to the DFT   *
!         energy.                                                      *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chem. Phys.                       *
!             University of Lund, SWEDEN                               *
!             August 2002                                              *
!***********************************************************************

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Para_Info, only: King
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
#include "Molcas.fh"
#include "print.fh"
#include "rctfld.fh"
#include "disp.fh"
#include "nq_info.fh"
integer(kind=iwp) :: iDFT, iDumm, iEnd, iIrrep, iOpt, iPrint, iRout, iSpin, jPrint, LuWr, nDens
real(kind=wp) :: Dummy1(1), Dummy2(1), Dummy3(1), Dummy4, Dumm0(1), Dumm1(1), ExFac, TCpu1, TCpu2, TWall1, TWall2
logical(kind=iwp) :: First, Dff, Do_Grad, l_casdft
character(len=80) :: Label
character(len=16) :: KSDFT
character(len=4) :: DFTFOCK

!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu1,TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
!...  Prologue
DFTFOCK = 'SCF '
iRout = 131
iPrint = nPrint(iRout)
LuWr = u6

nDens = 0
do iIrrep=0,nIrrep-1
  nDens = nDens+nBas(iIrrep)*(nBas(iIrrep)+1)/2
end do

!     D F T - g r a d i e n t                                          *
!***********************************************************************
!8)                                                                    *
!     D F T - g r a d i e n t

!call Get_iOption(iDFT)

call Get_cArray('DFT functional',KSDFT,16)
l_casdft = (KSDFT(1:5) == 'TLSDA') .or. (KSDFT(1:6) == 'TLSDA5') .or. (KSDFT(1:5) == 'TBLYP') .or. (KSDFT(1:5) == 'TOPBE') .or. &
           (KSDFT(1:6) == 'TSSBSW') .or. (KSDFT(1:5) == 'TSSBD') .or. (KSDFT(1:5) == 'TS12G') .or. (KSDFT(1:4) == 'TPBE') .or. &
           (KSDFT(1:5) == 'FTPBE') .or. (KSDFT(1:7) == 'TREVPBE') .or. (KSDFT(1:8) == 'FTREVPBE') .or. &
           (KSDFT(1:6) == 'FTLSDA') .or. (KSDFT(1:6) == 'FTOPBE') .or. (KSDFT(1:6) == 'FTBLYP)')

if (l_casdft) then
  DFTFOCK = 'ROKS'
  call Get_iScalar('System BitSwitch',iOpt)
  iOpt = ibset(iOpt,6)
  call Put_iScalar('System BitSwitch',iOpt)
end if

call Get_iScalar('System BitSwitch',iDFT)
if (btest(iDFT,6)) then

  call StatusLine(' Alaska:',' Computing DFT gradients')

  First = .true.
  Dff = .false.
  call Get_cArray('DFT functional',KSDFT,16)
  ExFac = Zero ! Set to proper value at retrun!
  Do_Grad = .true.
  call Get_iScalar('Multiplicity',iSpin)
  !write(LuWr,*) 'DrvDFTg: KSDFT=',KSDFT
  !write(LuWr,*) 'DrvDFTg: ExFac=',ExFac
  iDumm = 1
  call DrvDFT(Dummy1,Dummy2,Dummy3,Dummy4,nDens,First,Dff,lRF,KSDFT,ExFac,Do_Grad,Temp,nGrad,iSpin,Dumm0,Dumm1,iDumm,DFTFOCK)

  iEnd = 1
  do
    if (KSDFT(iEnd:iEnd) == ' ') then
      iEnd = iEnd-1
      exit
    else
      iEnd = iEnd+1
    end if
  end do
  Label = 'The DFT('//KSDFT(1:iEnd)//') contribution'
  jPrint = nPrint(112)
  !AMS
  !jPrint = 15
  if (jPrint >= 15) call PrGrad(Label,Temp,nGrad,ChDisp,5)
  if (king()) call DaXpY_(nGrad,One,Temp,1,Grad,1)
  if (iPrint >= 6) then
    write(LuWr,*)
    if (Grid_Type == Moving_Grid) then
      write(LuWr,*) 'DFT contribution computed for a moving grid.'
    else
      write(LuWr,*) 'DFT contribution computed for a fixed grid.'
    end if
    write(LuWr,*)
  end if

end if
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu2,TWall2)
call SavTim(5,TCpu2-TCpu1,TWall2-TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine DrvDFTg
