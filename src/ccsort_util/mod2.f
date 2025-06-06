************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1994, Per Ake Malmqvist                                *
*               1995,1996, Pavel Neogrady                              *
************************************************************************
       subroutine mod2 (nsym,nish,nash,nssh,norb,fi,eps)
c
c     this routine define fi(p,q) = delta(p,q).eps(p)
c     and redefine nish=nish+nash, nash=0
c
c     this is suitable for closed shell case
c
       integer nsym
       integer nish(1:8)
       integer nash(1:8)
       integer nssh(1:8)
       integer norb(1:8)
       real*8 fi(*)
       real*8 eps(*)
c
c     help variables
c
       integer p,q,isym,pq,padd
c
c1    redefine foki
c
       pq=0
       padd=0
       do 200 isym=1,nsym
c
       do 100 p=1,norb(isym)
       do 101 q=1,p
       pq=pq+1
       if (p.eq.q) then
       fi(pq)=eps(padd+p)
       else
       fi(pq)=0.0d0
       end if
 101    continue
 100    continue
c
       padd=padd+norb(isym)
 200    continue
c
c2    redefine n's
c
       do 300 isym=1,nsym
       nish(isym)=nish(isym)+nash(isym)
       nash(isym)=0
 300    continue
c
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(nssh)
       end
