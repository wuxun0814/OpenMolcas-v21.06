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
       subroutine cct3_expand (wrk,wrksize,
     & nind,exptyp,mapda,mapia,ssa,possb0,mapdb,
     & mapib,rc)
c
c     this routine realize expansion
c
c     A(pqrs) -> B(pqrs)
c
c     nind   - # of indexex in matrix A  (Input)
c     exptyp - type of expansion :  (Input)
c     1 - pq,r,s -> p,q,r,s  pq,r -> p,q,r  pq -> p,q
c     2 - p,qr,s -> p,q,r,s  p,qr -> p,q,r
c     3 - p,q,rs -> p,q,r,s
c     4 - pq,rs  -> p,q,r,s
c     5 - pq,rs  -> p,q,rs
c     6 - pq,rs  -> pq,r,s
c     mapda  - direct map matrix corresponding to A  (Input)
c     mapia  - inverse map matrix corresponding to A  (Input)
c     ssa    - overall symetry state  of matrix A  (Input)
c     possb0 - initial possition of matrix B in WRK  (Input)
c     mapdb  - direct map matrix corresponding to B  (Output)
c     mapib  - inverse map matrix corresponding to B  (Output)
c     rc     - return (error) code  (Output)
c
c     Table of expansions
c
c     nind  exptyp          Operation             Implementation
c     4       0     A(p,q,r,s) -> B(p,q,r,s)     Realized in map
c     4       1     A(pq,r,s)  -> B(p,q,r,s)           Yes
c     4       2     A(p,qr,s)  -> B(p,q,r,s)           Yes
c     4       3     A(p,q,rs)  -> B(p,q,r,s)           Yes
c     4       4     A(pq,rs)   -> B(p,q,r,s)           Yes
c     4       5     A(pq,rs)   -> B(p,q,rs)            Yes
c     4       6     A(pq,rs)   -> B(pq,r,s)            Yes
c
c     3       0     A(p,q,r)   -> B(p,q,r)       Realized in map
c     3       1     A(pq,r)    -> B(p,q,r)             Yes
c     3       2     A(p,qr)    -> B(p,q,r)             Yes
c
c     2       0     A(p,q)     -> B(p,q)         Realized in map
c     2       1     A(pq)      -> B(p,q)               Yes
c
c     1       0     A(p)       -> B(p)           Realized in map
c
c
#include "t31.fh"
#include "wrk.fh"
c
       integer nind,exptyp,ssa,possb0,rc
       integer mapda(0:512,1:6),mapdb(0:512,1:6)
       integer mapia(1:8,1:8,1:8),mapib(1:8,1:8,1:8)
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6
       integer sa1,sa2,sa3,sa4
       integer na,ia,ib1,ib2,ib3,ib4
       integer typa,posst
c
       rc=0
       na=mapda(0,5)
       typa=mapda(0,6)
c
c     general tests
c
       if (exptyp.eq.0) then
c     RC=1  : exptyp=0 (for exptyp=0, there is no sopystical expansion, NCI)
       rc=1
       return
       end if
c
c     get mapdb,mapib
c
       if ((nind.eq.4).and.(exptyp.eq.5)) then
       nhelp1=3
       else if ((nind.eq.4).and.(exptyp.eq.6)) then
       nhelp1=1
       else
       nhelp1=0
       end if
c
       call cct3_grc0 (nind,nhelp1,mapda(0,1),mapda(0,2),mapda(0,3),
     & mapda(0,4),ssa,
     & possb0,posst,mapdb,mapib)
c
c
       if (nind.lt.2) then
c     RC=2 - number of indexes < 2 (NCI)
       rc=2
       end if
c
       if (nind.eq.2) then
c
c     ********** 2 index *********
c
       if (exptyp.eq.1) then
c
c2.1  expand A(pq) -> B(p,q)
c
c     tests
c
       if (typa.ne.1) then
c     RC=3 : nind=2, exptyp=1 (typA is not 1, Stup)
       rc=3
       return
       end if
c
       do 100 ia=1,na
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
c
       ib1=mapib(sa1,1,1)
       ib2=mapib(sa2,1,1)
c
       if (sa1.gt.sa2) then
c
c     map A(p,q) -> B(p,q)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_map11 (wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)
c
c     map A(p,q) -> - B(q,p)
c
c     possB2
       nhelp2=mapdb(ib2,1)
c     dimp,dimq
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp4=dimm(mapda(0,2),sa2)
c
       call cct3_map21 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,2,1,-1)
c
       else
c     sa1=sa2
c
c     expand A(pq) -> B(p,q)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB
       nhelp2=mapdb(ib1,1)
c     dimp
       nhelp3=dimm(mapda(0,1),sa1)
c
       call cct3_expand0 (wrk(nhelp1),wrk(nhelp2),mapda(ia,2),nhelp3)
c
       end if

 100    continue
c
       else
c     RC=4 : nind=2, exptyp=@ (Inproper exptyp, Stup)
       rc=4
       return
       end if
c
       else if (nind.eq.3) then
c
c     ********** 3 index *********
c
       if (exptyp.eq.1) then
c
c3.1  expand A(pq,r) -> B(p,q,r)
c
c     tests
c
       if (typa.ne.1) then
c     RC=5 : nind=3, exptyp=1 (typA is not 1, Stup)
       rc=5
       return
       end if
c
       do 300 ia=1,na
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
c
       ib1=mapib(sa1,sa2,1)
       ib2=mapib(sa2,sa1,1)
c
       if (sa1.gt.sa2) then
c
c     map A(p,q,r) -> B(p,q,r)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_map11 (wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)
c
c     map A(p,q,r) -> - B(q,p,r)
c
c     possB2
       nhelp2=mapdb(ib2,1)
c     dimp,dimq,dimr
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp4=dimm(mapda(0,2),sa2)
       nhelp5=dimm(mapda(0,3),sa3)
c
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,
     & nhelp5,2,1,3,-1)
c
       else
c     sa1=sa2
c
c     expand A(pq,r) -> B(p,q,r)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB
       nhelp2=mapdb(ib1,1)
c     dimpq
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp3=nhelp3*(nhelp3-1)/2
c     dimp,dimr
       nhelp4=dimm(mapda(0,1),sa1)
       nhelp5=dimm(mapda(0,3),sa3)
c
       call cct3_expand1 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5)
c
       end if

 300    continue
c
       else if (exptyp.eq.2) then
c
c3.2  expand A(p,qr) -> B(p,q,r)
c
c     tests
c
       if (typa.ne.2) then
c     RC=6 : nind=3, exptyp=2 (typA is not 2, Stup)
       rc=6
       return
       end if
c
       do 400 ia=1,na
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
c
       ib1=mapib(sa1,sa2,1)
       ib2=mapib(sa1,sa3,1)
c
       if (sa2.gt.sa3) then
c
c     map A(p,q,r) -> B(p,q,r)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_map11 (wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)
c
c     map A(p,q,r) -> - B(q,p,r)
c
c     possB2
       nhelp2=mapdb(ib2,1)
c     dimp,dimq,dimr
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp4=dimm(mapda(0,2),sa2)
       nhelp5=dimm(mapda(0,3),sa3)
c
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,
     & nhelp5,1,3,2,-1)
c
       else
c     sa2=sa3
c
c     expand A(p,qr) -> B(p,q,r)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB
       nhelp2=mapdb(ib1,1)
c     dimqr
       nhelp3=dimm(mapda(0,2),sa1)
       nhelp3=nhelp3*(nhelp3-1)/2
c     dimp,dimq
       nhelp4=dimm(mapda(0,1),sa1)
       nhelp5=dimm(mapda(0,2),sa2)
c
       call cct3_expand3 (wrk(nhelp1),wrk(nhelp2),nhelp4,nhelp3,nhelp5)
c
       end if

 400    continue
c
       else
c     RC=7 : nind=3, exptyp=@ (Inproper exptyp, Stup)
       rc=7
       return
       end if
c
       else if (nind.eq.4) then
c
c     ********** 4 index *********
c
       if (exptyp.eq.1) then
c
c4.1  expand A(pq,r,s) -> B(p,q,r,s)
c
c     tests
c
       if (typa.ne.1) then
c     RC=8 : nind=4, exptyp=1 (typA is not 1, Stup)
       rc=8
       return
       end if
c
       do 500 ia=1,na
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
       sa4=mapda(ia,6)
c
       ib1=mapib(sa1,sa2,sa3)
       ib2=mapib(sa2,sa1,sa3)
c
       if (sa1.gt.sa2) then
c
c     map A(p,q,r,s) -> B(p,q,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_map11 (wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)
c
c     map A(p,q,r,s) -> - B(q,p,r,s)
c
c     possB2
       nhelp2=mapdb(ib2,1)
c     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp4=dimm(mapda(0,2),sa2)
       nhelp5=dimm(mapda(0,3),sa3)
       nhelp6=dimm(mapda(0,4),sa4)
c
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,
     & nhelp6,2,1,3,4,-1)
c
       else
c     sa1=sa2
c
c     expand A(pq,r_s) -> B(p,q,r_s)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB
       nhelp2=mapdb(ib1,1)
c     dimpq
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp3=nhelp3*(nhelp3-1)/2
c     dimp,dimr_s
       nhelp4=dimm(mapda(0,1),sa1)
       nhelp5=dimm(mapda(0,3),sa3)*dimm(mapda(0,4),sa4)
c
       call cct3_expand1 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5)
c
       end if

 500    continue
c
       else if (exptyp.eq.2) then
c
c4.2  expand A(p,qr,s) -> B(p,q,r,s)
c
c     tests
c
       if (typa.ne.2) then
c     RC=9 : nind=4, exptyp=2 (typA is not 2, Stup)
       rc=9
       return
       end if
c
       do 600 ia=1,na
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
       sa4=mapda(ia,6)
c
       ib1=mapib(sa1,sa2,sa3)
       ib2=mapib(sa1,sa3,sa2)
c
       if (sa2.gt.sa3) then
c
c     map A(p,q,r,s) -> B(p,q,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_map11 (wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)
c
c     map A(p,q,r,s) -> - B(p,r,q,s)
c
c     possB2
       nhelp2=mapdb(ib2,1)
c     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp4=dimm(mapda(0,2),sa2)
       nhelp5=dimm(mapda(0,3),sa3)
       nhelp6=dimm(mapda(0,4),sa4)
c
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,
     &  nhelp6,1,3,2,4,-1)
c
       else
c     sa2=sa3
c
c     expand A(p,qr,s) -> B(p,q,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB
       nhelp2=mapdb(ib1,1)
c     dimqr
       nhelp3=dimm(mapda(0,2),sa1)
       nhelp3=nhelp3*(nhelp3-1)/2
c     dimp,dims,dimq
       nhelp4=dimm(mapda(0,1),sa1)
       nhelp5=dimm(mapda(0,4),sa4)
       nhelp5=dimm(mapda(0,2),sa2)
c
       call cct3_expand2 (wrk(nhelp1),wrk(nhelp2),nhelp4,nhelp3,nhelp5,
     &               nhelp6)
c
       end if

 600    continue
c
       else if (exptyp.eq.3) then
c
c4.3  expand A(p,q,rs) -> B(p,q,r,s)
c
c     tests
c
       if (typa.ne.3) then
c     RC=10: nind=4, exptyp=3 (typA is not 3, Stup)
       rc=10
       return
       end if
c
       do 700 ia=1,na
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
       sa4=mapda(ia,6)
c
       ib1=mapib(sa1,sa2,sa3)
       ib2=mapib(sa1,sa2,sa4)
c
       if (sa3.gt.sa4) then
c
c     map A(p,q,r,s) -> B(p,q,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_map11 (wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)
c
c     map A(p,q,r,s) -> - B(q,p,r,s)
c
c     possB2
       nhelp2=mapdb(ib2,1)
c     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp4=dimm(mapda(0,2),sa2)
       nhelp5=dimm(mapda(0,3),sa3)
       nhelp6=dimm(mapda(0,4),sa4)
c
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,
     & nhelp6,1,2,4,3,-1)
c
       else
c     sa1=sa2
c
c     expand A(p_q,rs) -> B(p_q,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB
       nhelp2=mapdb(ib1,1)
c     dimrs
       nhelp3=dimm(mapda(0,3),sa3)
       nhelp3=nhelp3*(nhelp3-1)/2
c     dimr,dimp_q
       nhelp4=dimm(mapda(0,3),sa3)
       nhelp5=dimm(mapda(0,1),sa1)*dimm(mapda(0,2),sa2)
c
       call cct3_expand3 (wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp3,nhelp4)
c
       end if

 700    continue
c
       else if (exptyp.eq.4) then
c
c4.4  expand A(pq,rs) -> B(p,q,r,s)
c
c     tests
c
       if (typa.ne.4) then
c     RC=11: nind=4, exptyp=4 (typA is not 4, Stup)
       rc=11
       return
       end if
c
       do 800 ia=1,na
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
       sa4=mapda(ia,6)
c
       ib1=mapib(sa1,sa2,sa3)
       ib2=mapib(sa2,sa1,sa3)
       ib3=mapib(sa1,sa2,sa4)
       ib4=mapib(sa2,sa1,sa4)
c
       if ((sa1.gt.sa2).and.(sa3.gt.sa4)) then
c
c     map A(p,q,r,s) -> B(p,q,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_map11 (wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)
c
c     map A(p,q,r,s) -> - B(q,p,r,s)
c     map A(p,q,r,s) -> - B(p,q,s,r)
c     map A(p,q,r,s) -> + B(q,p,s,r)
c
c     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp4=dimm(mapda(0,2),sa2)
       nhelp5=dimm(mapda(0,3),sa3)
       nhelp6=dimm(mapda(0,4),sa4)
c
c     possB2
       nhelp2=mapdb(ib2,1)
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,
     & nhelp5,nhelp6,2,1,3,4,-1)
c
c     possB3
       nhelp2=mapdb(ib3,1)
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,
     &  nhelp5,nhelp6,1,2,4,3,-1)
c
c     possB4
       nhelp2=mapdb(ib4,1)
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,
     &  nhelp6,2,1,4,3,1)
c
       else if ((sa1.eq.sa2).and.(sa3.eq.sa4)) then
c
c     expand A(pq,rs) -> B(p,q,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB
       nhelp2=mapdb(ib1,1)
c     dimp,dimq
       nhelp5=dimm(mapda(0,1),sa1)
       nhelp6=dimm(mapda(0,3),sa3)
c     dimpq,dimrs
       nhelp3=nhelp5*(nhelp5-1)/2
       nhelp4=nhelp6*(nhelp6-1)/2
c
       call cct3_expand40 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,
     &                nhelp6)
c
       else if (sa1.eq.sa2) then
c
c     expand A(pq,r_s) -> B(p,q,r_s)
c
c     possA
       nhelp1=mapda(ia,1)
c     dimr,dims
       nhelp5=dimm(mapda(0,3),sa3)
       nhelp6=dimm(mapda(0,4),sa4)
c     dimpq
       nhelp4=dimm(mapda(0,1),sa1)
       nhelp3=nhelp4*(nhelp4-1)/2
c
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_expand1 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5*nhelp6,
     &               nhelp4)
c
c     expand A(pq,r,s) -> - B(p,q,s,r)
c
c     possB3
       nhelp2=mapdb(ib3,1)
       call cct3_expand41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5,nhelp6,
     &                nhelp4)
c
       else if (sa3.eq.sa4) then
c
c     expand A(p_q,rs) -> B(p_q,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     dimp,dimq
       nhelp5=dimm(mapda(0,1),sa1)
       nhelp6=dimm(mapda(0,2),sa2)
c     dimrs
       nhelp4=dimm(mapda(0,3),sa3)
       nhelp3=nhelp4*(nhelp4-1)/2
c
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_expand3 (wrk(nhelp1),wrk(nhelp2),nhelp5*nhelp6,nhelp3,
     &               nhelp4)
c
c     expand A(p,q,rs) -> - B(q,p,r,s)
c
c     possB4
       nhelp2=mapdb(ib2,1)
       call cct3_expand41 (wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp6,nhelp3,
     &                nhelp4)
c
       end if
c
 800    continue
c
       else if (exptyp.eq.5) then
c
c4.5  expand A(pq,rs) -> B(p,q,rs)
c
c     tests
c
       if (typa.ne.4) then
c     RC=12: nind=4, exptyp=5 (typA is not 4, Stup)
       rc=12
       return
       end if
c
       do 900 ia=1,na
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
       sa4=mapda(ia,6)
c
       ib1=mapib(sa1,sa2,sa3)
       ib2=mapib(sa2,sa1,sa3)
c
       if (sa3.gt.sa4) then
c
c     map A(p,q,r,s) -> B(p,q,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_map11 (wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)
c
c     map A(p,q,r,s) -> - B(q,p,r,s)
c
c     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp4=dimm(mapda(0,2),sa2)
       nhelp5=dimm(mapda(0,3),sa3)
       nhelp6=dimm(mapda(0,4),sa4)
c
c     possB2
       nhelp2=mapdb(ib2,1)
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,
     & nhelp6,2,1,3,4,-1)
c
       else if ((sa1.eq.sa2).and.(sa3.eq.sa4)) then
c
c     expand A(pq,rs) -> B(p,q,rs)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB
       nhelp2=mapdb(ib1,1)
c     dimp,dimr
       nhelp5=dimm(mapda(0,1),sa1)
       nhelp6=dimm(mapda(0,3),sa3)
c     dimpq,dimrs
       nhelp3=nhelp5*(nhelp5-1)/2
       nhelp4=nhelp6*(nhelp6-1)/2
c
       call cct3_expand1 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5)
c
       else if (sa1.eq.sa2) then
c
c     expand A(pq,r_s) -> B(p,q,r_s)
c
c     possA
       nhelp1=mapda(ia,1)
c     dimr,dims
       nhelp5=dimm(mapda(0,3),sa3)
       nhelp6=dimm(mapda(0,4),sa4)
c     dimpq
       nhelp4=dimm(mapda(0,1),sa1)
       nhelp3=nhelp4*(nhelp4-1)/2
c
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_expand1 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5*nhelp6,
     &               nhelp4)
c
       else if (sa3.eq.sa4) then
c
c     map A(p,q,rs) -> B(p,q,rs)
c
c     possA
       nhelp1=mapda(ia,1)
c     dimp,dimq
       nhelp5=dimm(mapda(0,1),sa1)
       nhelp6=dimm(mapda(0,2),sa2)
c     dimrs
       nhelp4=dimm(mapda(0,3),sa3)
       nhelp3=nhelp4*(nhelp4-1)/2
c
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp6,
     &  nhelp3,1,2,3,1)
c
c     map A(p,q,rs) -> - B(q,p,rs)
c
c     possB4
       nhelp2=mapdb(ib2,1)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp5,nhelp6,
     & nhelp3,2,1,3,-1)
c
       end if
c
 900    continue
c
       else if (exptyp.eq.6) then
c
c4.6  expand A(pq,rs) -> B(pq,r,s)
c
c     tests
c
       if (typa.ne.4) then
c     RC=13: nind=4, exptyp=6 (typA is not 4, Stup)
       rc=13
       return
       end if
c
       do 1000 ia=1,na
c
       sa1=mapda(ia,3)
       sa2=mapda(ia,4)
       sa3=mapda(ia,5)
       sa4=mapda(ia,6)
c
       ib1=mapib(sa1,sa2,sa3)
       ib3=mapib(sa1,sa2,sa4)
c
       if ((sa1.gt.sa2).and.(sa3.gt.sa4)) then
c
c     map A(p,q,r,s) -> B(p,q,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_map11 (wrk(nhelp1),wrk(nhelp2),mapda(ia,2),1)
c
c     map A(p,q,r,s) -> - B(q,p,r,s)
c     map A(p,q,r,s) -> - B(p,q,s,r)
c     map A(p,q,r,s) -> + B(q,p,s,r)
c
c     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp4=dimm(mapda(0,2),sa2)
       nhelp5=dimm(mapda(0,3),sa3)
       nhelp6=dimm(mapda(0,4),sa4)
c
c     possB3
       nhelp2=mapdb(ib3,1)
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,
     & nhelp6,1,2,4,3,-1)
c
       else if ((sa1.eq.sa2).and.(sa3.eq.sa4)) then
c
c     expand A(pq,rs) -> B(pq,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     possB
       nhelp2=mapdb(ib1,1)
c     dimp,dimq
       nhelp5=dimm(mapda(0,1),sa1)
       nhelp6=dimm(mapda(0,3),sa3)
c     dimpq,dimrs
       nhelp3=nhelp5*(nhelp5-1)/2
       nhelp4=nhelp6*(nhelp6-1)/2
c
c@@?  call cct3_expand4 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,nhelp6)
       call cct3_expand3 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp6)
c
       else if (sa1.eq.sa2) then
c
c     map A(pq,r,s) -> B(pq,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     dimp,dimr,dims
       nhelp3=dimm(mapda(0,1),sa1)
       nhelp4=dimm(mapda(0,3),sa3)
       nhelp5=dimm(mapda(0,4),sa4)
c     dimpq
       nhelp6=nhelp3*(nhelp3-1)/2
c
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp6,nhelp4,
     &  nhelp5,1,2,3,1)
c
c     map A(pq,r,s) -> - B(pq,s,r)
c
c     possB3
       nhelp2=mapdb(ib3,1)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp6,
     & nhelp4,nhelp5,1,3,2,-1)
c
       else if (sa3.eq.sa4) then
c
c     expand A(p,q,rs) -> B(p,q,r,s)
c
c     possA
       nhelp1=mapda(ia,1)
c     dimp,dimq
       nhelp5=dimm(mapda(0,1),sa1)
       nhelp6=dimm(mapda(0,2),sa2)
c     dimrs
       nhelp4=dimm(mapda(0,3),sa3)
       nhelp3=nhelp4*(nhelp4-1)/2
c
c     possB1
       nhelp2=mapdb(ib1,1)
       call cct3_expand3 (wrk(nhelp1),wrk(nhelp2),nhelp5*nhelp6,nhelp3,
     &               nhelp4)
c
       end if
c
 1000   continue
c
       else
c     RC=14: nind=4, exptyp=@ (Inproper exptyp, Stup)
       rc=14
       end if
c
       else
c     RC=15- nind=@ (number of indexes >4, NCI)
       rc=15
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(mapia)
       end
