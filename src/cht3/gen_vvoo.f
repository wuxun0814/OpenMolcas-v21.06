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
        subroutine gen_vvoo (w,l1,tmp,l2)
c
c this routine do
c
c regenerate (ab,ij) integrals from blocked
c MO cholesky vectors
c
c <vv|oo> = (vo|vo)
c
c --------
c
c       L1(m,I,A')
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        real*8 tmp(*),l1(*),w(*),l2(*)
        integer a,b,dima,dimb,length,lasta,lastb
        integer a_tmp,b_tmp
c
        do a=1,NvGrp
c
c1        read tmp(m,I,A')
c
        dima=DimGrpaR(a)
         length=nc*no*dima
         call GetX_t3 (tmp,length,LunAux,L1Name(a),1,1)
c
c5        map l1 (A',I,m) <- tmp (m,I,A')
c
         call Map3_321_t3 (tmp,l1,nc,no,dima)
c
c  ----- read tmp(m,I,B')
c
        do b=1,a
c
        dimb=DimGrpaR(b)
         length=nc*no*dimb
         call GetX_t3 (tmp,length,LunAux,L1Name(b),1,1)
c
c4        map l2 (m,B',I) <- tmp (m,I,B')
c
        call Map3_132_t3 (tmp,l2,nc,no,dimb)
c
c        zero tmp
c
        call zeroma (tmp,1,dima*no*dimb*no)
c
c7      mult tmp(A',I,B',J) <- l1 (A',I,m) . l2(m,B',J)
c
        call mc0c1a3b (
     & dima*no,nc,nc,dimb*no,
     & dima*no,dimb*no,
     & dima*no,nc,dimb*no,l1,l2,tmp)
c
          lasta=0
        if (a.gt.1) then
          do a_tmp=1,a-1
          lasta=lasta+DimGrpaR(a_tmp)
          end do
        end if
c
          lastb=0
        if (b.gt.1) then
          do b_tmp=1,b-1
          lastb=lastb+DimGrpaR(b_tmp)
          end do
        end if
c
c        write (6,'(A,4(i4,x))') 'BB1 dima, dimb, lasta, lastb ',
c     & dima,dimb,lasta,lastb
        call grow_vvoo(w,tmp,no,nv,dima,dimb,lasta,lastb)
c
        if (a.ne.b) then
c
        call Map4_3412_t3 (tmp,l2,dima,no,dimb,no)
c
c        write (6,'(A,4(i4,x))') 'BB2 dima, dimb, lasta, lastb ',
c     & dimb,dima,lastb,lasta
        call grow_vvoo(w,l2,no,nv,dimb,dima,lastb,lasta)
c
        end if
c
c3        end loop over B'
c
        end do
c
c3        end loop over A'
c
        end do
c
        return
        end
