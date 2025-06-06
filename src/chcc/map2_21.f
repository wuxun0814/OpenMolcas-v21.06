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
        subroutine Map2_21 (A,B,d1,d2)
c
c       this routine do:
c       map B(21) <- A(12)
c
        implicit none
        integer d1,d2
        real*8 A(1:d1,1:d2)
        real*8 B(1:d2,1:d1)
c
c       help variables
        integer i1,i2
c
        do i1=1,d1
        do i2=1,d2
        b(i2,i1)=a(i1,i2)
        end do
        end do
c
        return
        end
