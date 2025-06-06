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

subroutine TRAONE_MOTRA(PAO,PMO,TEMP,CMO)
! Transformation program: one-electron section
!
! Objective: transforms a one-electron matrix PAO in AO-basis
!            to a molecular orbital matrix PMO.

#include "intent.fh"

use motra_global, only: nBas, nDel, nFro, nOrb, nSym
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: PAO(*), CMO(*)
real(kind=wp), intent(_OUT_) :: PMO(*), TEMP(*)
integer(kind=iwp) :: IAO, ICMO, IMO, IOFF, ISYM

ICMO = 1
IAO = 1
IMO = 1
do ISYM=1,NSYM
  ICMO = ICMO+NBAS(ISYM)*NFRO(ISYM)
  IOFF = 1+NBAS(ISYM)*NBAS(ISYM)
  call SQUARE(PAO(IAO),TEMP(1),1,NBAS(ISYM),NBAS(ISYM))
  !call MXMA(CMO(ICMO),NBAS(ISYM),1,TEMP,1,NBAS(ISYM),TEMP(IOFF),1,NORB(ISYM),NORB(ISYM),NBAS(ISYM),NBAS(ISYM))
  call DGEMM_('T','N',NORB(ISYM),NBAS(ISYM),NBAS(ISYM),One,CMO(ICMO),NBAS(ISYM),TEMP,NBAS(ISYM),Zero,TEMP(IOFF),max(1,NORB(ISYM)))
  call MXMT(TEMP(IOFF),1,NORB(ISYM),CMO(ICMO),1,NBAS(ISYM),PMO(IMO),NORB(ISYM),NBAS(ISYM))
  ICMO = ICMO+NBAS(ISYM)*(NORB(ISYM)+NDEL(ISYM))
  IAO = IAO+NBAS(ISYM)*(NBAS(ISYM)+1)/2
  IMO = IMO+NORB(ISYM)*(NORB(ISYM)+1)/2
end do

return

end subroutine TRAONE_MOTRA
