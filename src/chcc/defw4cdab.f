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
        subroutine DefW4cdab (W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,
     c                        aSGrp,bSGrp,cSGrp,dSGrp)
c
c        define W(a,b,c,d) from (cd|ab)
c
        implicit none
        integer dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp
        real*8 W(1:dima,1:dimb,1:dimc,1:dimd)
        real*8 Wx(1:cdLen,1:abLen)
c
c        help variables
        integer a,b,c,d,ab,cd
c
        if ((aSGrp.eq.bSGrp).and.(cSGrp.eq.dSGrp)) then
c        case (c=d|a=b)
          do a=2,dima
          ab=a*(a-1)/2
          do b=1,a-1
          ab=ab+1
            do c=2,dimc
            cd=c*(c-1)/2
            do d=1,c-1
            cd=cd+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(cd,ab)
              W(b,a,c,d)=W(b,a,c,d)+Wx(cd,ab)
              W(a,b,d,c)=W(a,b,d,c)+Wx(cd,ab)
              W(b,a,d,c)=W(b,a,d,c)+Wx(cd,ab)
            end do
            end do
            do c=1,dimc
            cd=c*(c+1)/2
              W(a,b,c,c)=W(a,b,c,c)+Wx(cd,ab)
              W(b,a,c,c)=W(b,a,c,c)+Wx(cd,ab)
            end do
          end do
          end do
c
          do a=1,dima
          ab=a*(a+1)/2
            do c=2,dimc
            cd=c*(c-1)/2
            do d=1,c-1
            cd=cd+1
              W(a,a,c,d)=W(a,a,c,d)+Wx(cd,ab)
              W(a,a,d,c)=W(a,a,d,c)+Wx(cd,ab)
            end do
            end do
            do c=1,dimc
            cd=c*(c+1)/2
              W(a,a,c,c)=W(a,a,c,c)+Wx(cd,ab)
            end do
          end do
c
        else if ((aSGrp.eq.bSGrp).and.(cSGrp.ne.dSGrp)) then
c        case (c,d|a=b)
          do a=2,dima
          ab=a*(a-1)/2
          do b=1,a-1
          ab=ab+1
            cd=0
            do d=1,dimd
            do c=1,dimc
            cd=cd+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(cd,ab)
              W(b,a,c,d)=W(b,a,c,d)+Wx(cd,ab)
            end do
            end do
          end do
          end do
c
          do a=1,dima
          ab=a*(a+1)/2
            cd=0
            do d=1,dimd
            do c=1,dimc
            cd=cd+1
              W(a,a,c,d)=W(a,a,c,d)+Wx(cd,ab)
            end do
            end do
          end do
c
        else if ((aSGrp.ne.bSGrp).and.(cSGrp.eq.dSGrp)) then
c        case (c=d|a,b)
          ab=0
          do b=1,dimb
          do a=1,dima
          ab=ab+1
            do c=2,dimc
            cd=c*(c-1)/2
            do d=1,c-1
            cd=cd+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(cd,ab)
              W(a,b,d,c)=W(a,b,d,c)+Wx(cd,ab)
            end do
            end do
            do c=1,dimc
            cd=c*(c+1)/2
              W(a,b,c,c)=W(a,b,c,c)+Wx(cd,ab)
            end do
          end do
          end do
c
        else
c        case (c,d|a,b)
          ab=0
          do b=1,dimb
          do a=1,dima
          ab=ab+1
            cd=0
            do d=1,dimd
            do c=1,dimc
            cd=cd+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(cd,ab)
            end do
            end do
          end do
          end do
c
        end if
c
c
        return
        end
