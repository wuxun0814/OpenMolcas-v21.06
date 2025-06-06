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

module Gateway_global

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: G_Mode = 1, &
                                S_Mode = 2, &
                                GS_Mode = 3
integer(kind=iwp) :: Run_Mode

public :: Run_Mode, G_Mode, S_Mode, GS_Mode

end module Gateway_global
