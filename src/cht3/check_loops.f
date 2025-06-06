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
        subroutine check_loops(nv,vblock,nla,nlb)
c
        integer nv,vblock,nla,nlb
        integer nuga,nga,ngb,ngc

        nuga=nv/vblock
cmp! pridavok
        if((nuga*vblock).lt.nv)nuga=nuga+1
c
        nla=0
        do nga=1,nuga
        do ngb=1,nga
        do ngc=1,ngb
        nla=nla+1
        end do
        end do
        end do
c
        nlb=0
        do nga=1,nuga
        do ngb=1,nga
        do ngc=1,nuga
        nlb=nlb+1
        end do
        end do
        end do
c
        return
        end
