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

subroutine embPotMem( &
#define _CALLING_
#include "mem_interface.fh"
)

#include "macros.fh"

use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
unused_var(la)
unused_var(lb)
unused_var(lr)

nHer = 0
Mem = 0

return

end subroutine embPotMem
