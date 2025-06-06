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
       subroutine ireorg (wrk,wrksize,
     & symp,symq,symr,syms,typp,typq,typr,typs,
     & posspv1,possqv1,possrv1,posssv1,typpv1,typqv1,typrv1,typsv1,
     & typv2,possv10,possv20,fact)
c
c     this routine is up level routine for ireorg1 (also more detailed
c     description can be found there).
c     v1 must be of type 0, v2 can be 0,1,3 and 4
c     this routine only prepair some constants, required by ireorg1,
c     that can be deduced form input data - dimpq,dimrs,dimt-dimx
c
c     symp-s     - symetries of p-s (I)
c     typp-s     - types of indexes p-s in V2 (I)
c     possp-sv1  - possitions of p-s ind. in v1 (I)
c     typp-sv1   - types of indices, corresponding to p-s in V1 (I)
c     typv2      - type of V2 (0,1,2,4) (I)
c     possv10,20 - initial possitions of V1 and V2 in wrk (I)
c     fact       - multiplication factors (usually +-1.0d0) (I)
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer symp,symq,symr,syms,typp,typq,typr,typs
       integer posspv1,possqv1,possrv1,posssv1,typpv1,typqv1,typrv1,
     & typsv1
       integer typv2,possv10,possv20
       real*8 fact
c
c     help variables
c
       integer ind(1:4)
       integer rc,dimpq,dimrs
       integer :: nhelp=-1,mhelp=-1
c
c*    define dimensions of V1
c
       call ireorg2 (symp,typpv1,nhelp,rc)
       ind(posspv1)=nhelp
       call ireorg2 (symq,typqv1,nhelp,rc)
       ind(possqv1)=nhelp
       call ireorg2 (symr,typrv1,nhelp,rc)
       ind(possrv1)=nhelp
       call ireorg2 (syms,typsv1,nhelp,rc)
       ind(posssv1)=nhelp
c
c*    def dimpq,dimrs
c
       call ireorg2 (symp,typp,nhelp,rc)
       call ireorg2 (symq,typq,mhelp,rc)
c
       if (((typv2.eq.1).or.(typv2.eq.4)).and.(symp.eq.symq)) then
       dimpq=(nhelp*(nhelp-1))/2
       else
       dimpq=nhelp*mhelp
       end if
c
       call ireorg2 (symr,typr,nhelp,rc)
       call ireorg2 (syms,typs,mhelp,rc)
c
       if (((typv2.eq.3).or.(typv2.eq.4)).and.(symr.eq.syms)) then
       dimrs=(nhelp*(nhelp-1))/2
       else
       dimrs=nhelp*mhelp
       end if
c
c*    use ireorg1
c
       call ireorg1 (symp,symq,symr,syms,typp,typq,typr,typs,
     & posspv1,possqv1,possrv1,posssv1,typpv1,typqv1,typrv1,typsv1,
     & typv2,wrk(possv10),wrk(possv20),fact,dimpq,dimrs,
     & ind(1),ind(2),ind(3),ind(4))
c
       return
       end
