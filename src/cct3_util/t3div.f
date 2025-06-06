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
       subroutine t3div (wrk,wrksize,
     & mapdw,mapdv,ssw,mapdd1,mapid1,mapdd2,mapid2,
     & typdiv,i,j,k,symi,symj,symk,ec,rc)
c
c     mapdw  - direct map matrix of W (Input)
c     mapdv  - direct map matrix of V (Input)
c     ssw    - overall symmetry state of matrix W (V) (Input)
c     mapdd1 - direct map matrix of 1.st diagonal el. (Input)
c     mapid1 - inverse map matrix of 1.st diagonal el. (Input)
c     mapdd2 - direct map matrix of 2.st diagonal el. (Input)
c     mapid2 - inverse map matrix of 2.st diagonal el. (Input)
c     (if there is only one spin, use any map's for 2)
c     typdiv - typ of operation (see Table) (Input)
c     i      - value of occupied index i (Inlut)
c     j      - value of occupied index j (Inlut)
c     k      - value of occupied index k (Inlut)
c     symi   - symmetry of index i (Input)
c     symj   - symmetry of index j (Input)
c     symk   - symmetry of index k (Input)
c     ec     - energy contribution from this part (Output)
c     rc     - return (error) code (Output)
c
c     this routine realize division by denominators and
c     calculate corresponding energy contribution
c
c     ec = sum (abc) [ W(abc).V(abc) / Dijkabc ]
c
c     for following types of W,V
c
c     typdiv         Operation                Implemented
c     1     W(abc)  . V(abc) /Dijkabc           Yes
c     2     W(ab,c) . W(ab,c)/Dijkabc           Yes
c     3     W(a,bc) . V(a,bc)/Dijkabc           Yes
c
c     N.B. spin combinations aaa,bbb for 1; aab for 2; and abb for 3
c     are authomatically assumed
c
#include "t31.fh"
#include "wrk.fh"
c
       integer ssw,typdiv,i,j,k,symi,symj,symk,rc
       integer mapdw(0:512,1:6)
       integer mapdv(0:512,1:6)
       integer mapdd1(0:512,1:6)
       integer mapdd2(0:512,1:6)
c       integer mapiw(1:8,1:8,1:8)
c       integer mapiv(1:8,1:8,1:8)
       integer mapid1(1:8,1:8,1:8)
       integer mapid2(1:8,1:8,1:8)
       real*8 ec
c
c     help variables
c
       integer iw,possw,possv
       integer id1,id2,id3,possd1,possd2,possd3
       integer syma,symb,symc,dima,dimb,dimc
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6
       integer index
       real*8 eco,denijk
c
c
c0.*  some tests
c
       if (mapdw(0,6).ne.mapdv(0,6)) then
c     RC=1 : typW is not equal to typV (Stup)
       rc=1
       return
       end if
c
       if ((typdiv.eq.1).and.(mapdw(0,6).ne.5)) then
c     RC=2 : typdiv=1, typW is not 5 (Stup)
       rc=2
       return
       end if
c
       if ((typdiv.eq.2).and.(mapdw(0,6).ne.1)) then
c     RC=3 : typdiv=2, typW is not 1 (Stup)
       rc=3
       return
       end if
c
       if ((typdiv.eq.3).and.(mapdw(0,6).ne.2)) then
c     RC=4 : typdiv=3, typW is not 2 (Stup)
       rc=4
       return
       end if
c
c0.*  vanish ec
c
       ec=0.0d0
c
c0.*  def denijk
c
       if (typdiv.eq.1) then
c     cases aaa,bbb
c
c     diagonal part i
       id1=mapid1(symi,1,1)
       possd1=mapdd1(id1,1)
       index=possd1+i-1
       denijk=wrk(index)
c
c     diagonal part j
       id1=mapid1(symj,1,1)
       possd1=mapdd1(id1,1)
       index=possd1+j-1
       denijk=denijk+wrk(index)
c
c     diagonal part k
       id1=mapid1(symk,1,1)
       possd1=mapdd1(id1,1)
       index=possd1+k-1
       denijk=denijk+wrk(index)
c
       else if (typdiv.eq.2) then
c     case aab
c
c     diagonal part i
       id1=mapid1(symi,1,1)
       possd1=mapdd1(id1,1)
       index=possd1+i-1
       denijk=wrk(index)
c
c     diagonal part j
       id1=mapid1(symj,1,1)
       possd1=mapdd1(id1,1)
       index=possd1+j-1
       denijk=denijk+wrk(index)
c
c     diagonal part k
       id2=mapid2(symk,1,1)
       possd2=mapdd2(id2,1)
       index=possd2+k-1
       denijk=denijk+wrk(index)
c
       else if (typdiv.eq.3) then
c     case abb
c
c     diagonal part i
       id1=mapid1(symi,1,1)
       possd1=mapdd1(id1,1)
       index=possd1+i-1
       denijk=wrk(index)
c
c     diagonal part j
       id2=mapid2(symj,1,1)
       possd2=mapdd2(id2,1)
       index=possd2+j-1
       denijk=denijk+wrk(index)
c
c     diagonal part k
       id2=mapid2(symk,1,1)
       possd2=mapdd2(id2,1)
       index=possd2+k-1
       denijk=denijk+wrk(index)
c
       end if
c
c
       if (typdiv.eq.1) then
c1    case W(pqr). V(pqr)
c
       do 100 iw=1,mapdw(0,5)
c
c1.*  def possition of W,V
       possw=mapdw(iw,1)
       possv=mapdv(iw,1)
c
c1.*  def symmetry status
       syma=mapdw(iw,3)
       symb=mapdw(iw,4)
       symc=mapdw(iw,5)
c
c1.*  def dimensions
       dima=dimm(mapdw(0,1),syma)
       dimb=dimm(mapdw(0,2),symb)
       dimc=dimm(mapdw(0,3),symc)
c
c1.*  realize packing
c
       if (syma.eq.symc) then
c1.a  case syma=symb=symc
c
c1.a.*find address for d1,2,3 (the same one)
       id1=mapid1(syma,1,1)
c
c1.a.*def possition of d1,2,3 (the same one)
       possd1=mapdd1(id1,1)
c
c1.a.*def additional dimensions
       nhelp1=dima*(dima-1)*(dima-2)/6
       nhelp2=dimm(5,syma)
       if (mapdw(0,1).eq.3) then
c     alpha cese
       nhelp3=noa(syma)
       else
c     beta case
       nhelp3=nob(syma)
       end if
c
c1.a.*do packing
       call t3dhlp4 (wrk(possw),wrk(possv),dima,nhelp1,
     & denijk,eco,
     & wrk(possd1),
     & nhelp2,nhelp3)
       ec=ec+eco
c
       else if (syma.eq.symb) then
c1.b  case syma=symb.ne.symc
c
c1.b.*find address for d1,3
       id1=mapid1(syma,1,1)
       id3=mapid1(symc,1,1)
c
c1.b.*def possition of d1,3
       possd1=mapdd1(id1,1)
       possd3=mapdd1(id3,1)
c
c1.b.*def additional dimensions
       nhelp1=dima*(dima-1)/2
       nhelp2=dimm(5,syma)
       nhelp3=dimm(5,symc)
       if (mapdw(0,1).eq.3) then
c     alpha cese
       nhelp4=noa(syma)
       nhelp5=noa(symc)
       else
c     beta case
       nhelp4=nob(syma)
       nhelp5=nob(symc)
       end if
c
c1.b.*do packing
       call t3dhlp2 (wrk(possw),wrk(possv),dima,nhelp1,dimc,
     & denijk,eco,
     & wrk(possd1),wrk(possd3),
     & nhelp2,nhelp3,nhelp4,nhelp5)
       ec=ec+eco
c
       else if (symb.eq.symc) then
c1.c  case syma.ne.symb=symc
c
c1.c.*find address for d1,2
       id1=mapid1(syma,1,1)
       id2=mapid1(symb,1,1)
c
c1.c.*def possition of d1,2
       possd1=mapdd1(id1,1)
       possd2=mapdd1(id2,1)
c
c1.c.*def additional dimensions
       nhelp1=dimb*(dimb-1)/2
       nhelp2=dimm(5,syma)
       nhelp3=dimm(5,symb)
       if (mapdw(0,1).eq.3) then
c     alpha cese
       nhelp4=noa(syma)
       nhelp5=noa(symb)
       else
c     beta case
       nhelp4=nob(syma)
       nhelp5=nob(symb)
       end if
c
c1.c.*do packing
       call t3dhlp3 (wrk(possw),wrk(possv),dima,dimb,nhelp1,
     & denijk,eco,
     & wrk(possd1),wrk(possd2),
     & nhelp2,nhelp3,nhelp4,nhelp5)
       ec=ec+eco
c
       else
c1.d  case syma.ne.symb.ne.symc
c
c1.d.*find address for d1,2,3
       id1=mapid1(syma,1,1)
       id2=mapid1(symb,1,1)
       id3=mapid1(symc,1,1)
c
c1.d.*def possition of d1,2,3
       possd1=mapdd1(id1,1)
       possd2=mapdd1(id2,1)
       possd3=mapdd1(id3,1)
c
c1.d.*def additional dimensions
       nhelp1=dimm(5,syma)
       nhelp2=dimm(5,symb)
       nhelp3=dimm(5,symc)
       if (mapdw(0,1).eq.3) then
c     alpha cese
       nhelp4=noa(syma)
       nhelp5=noa(symb)
       nhelp6=noa(symc)
       else
c     beta case
       nhelp4=nob(syma)
       nhelp5=nob(symb)
       nhelp6=nob(symc)
       end if
c
c1.d.*do packing
       call t3dhlp1 (wrk(possw),wrk(possv),dima,dimb,dimc,
     & denijk,eco,
     & wrk(possd1),wrk(possd2),wrk(possd3),
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6)
       ec=ec+eco
c
       end if
c
c
 100    continue

       else if (typdiv.eq.2) then
c2    case W(pq,r) . V(pq,r)
c
       do 200 iw=1,mapdw(0,5)
c
c2.*  def possition of W,W
       possw=mapdw(iw,1)
       possv=mapdv(iw,1)
c
c2.*  def symmetry status
       syma=mapdw(iw,3)
       symb=mapdw(iw,4)
       symc=mapdw(iw,5)
c
c2.*  def dimensions
       dima=dimm(mapdw(0,1),syma)
       dimb=dimm(mapdw(0,2),symb)
       dimc=dimm(mapdw(0,3),symc)
c
c2.*  realize packing
c
       if (syma.eq.symb) then
c2.a  case syma=symb,symc
c
c2.a.*find address for d1,3
       id1=mapid1(syma,1,1)
       id3=mapid2(symc,1,1)
c
c2.a.*def possition of d1,3
       possd1=mapdd1(id1,1)
       possd3=mapdd2(id3,1)
c
c2.a.*def additional dimensions
       nhelp1=dima*(dima-1)/2
       nhelp2=dimm(5,syma)
       nhelp3=dimm(5,symc)
       nhelp4=noa(syma)
       nhelp5=nob(symc)
c
c2.a.*do packing
       call t3dhlp2 (wrk(possw),wrk(possv),dima,nhelp1,dimc,
     & denijk,eco,
     & wrk(possd1),wrk(possd3),
     & nhelp2,nhelp3,nhelp4,nhelp5)
       ec=ec+eco
c
       else
c2.b  case syma.ne.symb,symc
c
c2.b.*find address for d1,2,3
       id1=mapid1(syma,1,1)
       id2=mapid1(symb,1,1)
       id3=mapid2(symc,1,1)
c
c2.b.*def possition of d1,2,3
       possd1=mapdd1(id1,1)
       possd2=mapdd1(id2,1)
       possd3=mapdd2(id3,1)
c
c2.b.*def additional dimensions
       nhelp1=dimm(5,syma)
       nhelp2=dimm(5,symb)
       nhelp3=dimm(5,symb)
       nhelp4=noa(syma)
       nhelp5=noa(symb)
       nhelp6=nob(symc)
c
c2.b.*do packing
       call t3dhlp1 (wrk(possw),wrk(possv),dima,dimb,dimc,
     & denijk,eco,
     & wrk(possd1),wrk(possd2),wrk(possd3),
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6)
       ec=ec+eco
c
       end if
c
 200    continue
c
c
       else if (typdiv.eq.3) then
c3    case B(p,qr) = B(p,qr) + ns*(A(p,r,q)-A(p,q,r))
c
       do 300 iw=1,mapdw(0,5)
c
c3.*  def possition of W,V
       possw=mapdw(iw,1)
       possv=mapdv(iw,1)
c
c3.*  def symmetry status
       syma=mapdw(iw,3)
       symb=mapdw(iw,4)
       symc=mapdw(iw,5)
c
c3.*  def dimensions
       dima=dimm(mapdw(0,1),syma)
       dimb=dimm(mapdw(0,2),symb)
       dimc=dimm(mapdw(0,3),symc)
c
c3.*  realize packing
c
       if (symb.eq.symc) then
c3.a  case syma,symb=symc
c
c3.a.*find address for d1,2
       id1=mapid1(syma,1,1)
       id2=mapid2(symb,1,1)
c
c3.a.*def possition of d1,2
       possd1=mapdd1(id1,1)
       possd2=mapdd2(id2,1)
c
c3.a.*def additional dimensions
       nhelp1=dimb*(dimb-1)/2
       nhelp2=dimm(5,syma)
       nhelp3=dimm(5,symb)
       nhelp4=noa(syma)
       nhelp5=nob(symb)
c
c3.a.*do packing
       call t3dhlp3 (wrk(possw),wrk(possv),dima,dimb,nhelp1,
     & denijk,eco,
     & wrk(possd1),wrk(possd2),
     & nhelp2,nhelp3,nhelp4,nhelp5)
       ec=ec+eco
c
       else
c3.b  case syma,symb.ne.symc
c
c3.b.*find address for d1,2,3
       id1=mapid1(syma,1,1)
       id2=mapid2(symb,1,1)
       id3=mapid2(symc,1,1)
c
c3.b.*def possition of d1,2,3
       possd1=mapdd1(id1,1)
       possd2=mapdd2(id2,1)
       possd3=mapdd2(id3,1)
c
c3.b.*def additional dimensions
       nhelp1=dimm(5,syma)
       nhelp2=dimm(5,symb)
       nhelp3=dimm(5,symc)
       nhelp4=noa(syma)
       nhelp5=nob(symb)
       nhelp6=nob(symc)
c
c3.b.*do packing
       call t3dhlp1 (wrk(possw),wrk(possv),dima,dimb,dimc,
     & denijk,eco,
     & wrk(possd1),wrk(possd2),wrk(possd3),
     & nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6)
       ec=ec+eco
c
       end if
c
 300    continue
c
c
       else
c     RC=5 , typdiv is not 1,2,3 (NCI)
       rc=5
       return
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(ssw)
       end
