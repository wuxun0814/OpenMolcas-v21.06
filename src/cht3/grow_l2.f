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
        subroutine grow_l2(A,B,nc,nv,dima,dimb,lasta,lastb)
c
c this routine do :
c
c grow A(A,B,m) from the blocked cholesky vectors
c B(a',b',m)
c
        implicit none
        integer i1,i2,i3,dima,dimb,nc,nv
        integer lasta,lastb
        real*8 A(nv,nv,nc),B(dima,dimb,nc)
c
        do i3=1,nc
        do i1=1,dima
        do i2=1,dimb
        A(lasta+i1,lastb+i2,i3)=B(i1,i2,i3)
        A(lastb+i2,lasta+i1,i3)=B(i1,i2,i3)
        end do
        end do
        end do
c
        return
        end
