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
       subroutine exth4 (a,b,dimp,dimpq,dimr,p)
c
c     this routine extract A(pq,r) -> B_p(q,r)
c
c     a     - matrxi a (Input)
c     b     - matrix b (Output)
c     dimp  - dimension of p (and also q) (Input)
c     dimpq - dimension of pq (Input)
c     dimr  - dimension of r (Input)
c     p     - value of index p (Input)
c
#include "t31.fh"
       integer dimp,dimpq,dimr,p
       real*8 a(1:dimpq,1:dimr)
       real*8 b(1:dimp,dimr)
c
c     help variables
c
       integer q,r,qp,pq0
c
       if (p.eq.0) then
       return
       end if
c
c     q>p part
       if (p.gt.1) then
       pq0=nshf(p)
       do 20 r=1,dimr
       do 21 q=1,p-1
       b(q,r)=a(pq0+q,r)
 21     continue
 20     continue
       end if
c
c     q=p part
       do 40 r=1,dimr
       b(p,r)=0.0d0
 40     continue
c
c     r<p part
       if (p.lt.dimp) then
       do 60 r=1,dimr
       do 61 q=p+1,dimp
       qp=nshf(q)+p
       b(q,r)=-a(qp,r)
 61     continue
 60     continue
       end if
c
       return
       end
