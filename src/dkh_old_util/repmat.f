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
* Copyright (C) 2019, Thomas J. Duignan                                *
*               2021, Rulin Feng                                       *
************************************************************************
      subroutine repmat(idbg,bInt,sInt,donorm)
      Use Basis_Info
      use Symmetry_Info, only: nIrrep
      implicit real*8(a-h,o-z)
#include "itmax.fh"
#include "Molcas.fh"
#include "rinfo.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
      integer icaddr(MxAO),numc(MxAO),ihelp(MxAtom,iTabMx),numb(MxAO)
      integer mcaddr(MxAO)
      real*8 bint(*),sint(*)
      logical New_Center,New_l,New_m, Old_Center, Old_l
      integer istart
      integer np, nc, ip, jp, kp
      logical donorm
      real*8, allocatable :: mag(:,:), u2c(:,:), scr(:,:)
      real*8, allocatable :: u2ct(:,:), fin(:,:), pa(:)
      real*8 kpp, finish
*     contracted basis, atomic basis functions
*
*     symmetry info
*
*     lant(i): number of atoms in i:th symmetry bf
*     expand the coefficient matrix into symmetry basis set
*     auxiliary
*     icaddr(i): adresses in coeff for a symmetry adapted function
*
c     do i=1,12640
c     write(67,'(d25.14)') bint(i)
c     enddo
      nSym=nIrrep
      iPrint=0
      if(iprint.ge.10.or.idbg.gt.0) then
        write(idbg,*) ' in repmat', nsym
         write(idbg,*) nSym, (nBas(i),i=0,nsym-1)
         write(idbg,*) nSym, (nrBas(i),i=1,nsym)
      endif
*
*     set up pointer
*
      k=0
      ia=0  ! center index
      ka=0  ! shell index
      Do iCnttp=1,nCnttp
         Do icnt = 1, dbsc(iCnttp)%nCntr
            ia=ia+1
            Do la=1,nAngr(ia)+1
               ka=ka+1
               ihelp(ia,la)=k
               k=k+nPrimr(ka)*nBasisr(ka)
            End Do
         End Do
      End Do
      If (iPrint.ge.10.or.idbg.gt.0) Then
         write(idbg,*) ' Help vector'
         ia=0
         Do iCnttp=1,nCnttp
            Do icnt = 1, dbsc(iCnttp)%nCntr
               ia=ia+1
               write(idbg,'(10i5)') (ihelp(ia,j),j=1,nAngr(ia)+1)
           End Do
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Loop over irreps
*
      k=0
      Do iSym = 1, nSym
         numck=1
         numcl=0
         kbias=0
*
*        Loop over basis functions in irrep
*
         Do iCont = 1, nrBas(iSym)
            k=k+1
            If (iCont.gt.1) Then
               New_Center=icent(k).ne.icent(k-1)
               New_l=lnang(k).ne.lnang(k-1)
               New_m=lmag(k).ne.lmag(k-1)
               If (New_m)  kbias=kbias-numc(k-1)
               If (New_Center.or.New_l)  kbias=0
            End If
            kbias = kbias + 1
            ka=0
            ia=0
            Do iCnttp = 1, nCnttp
               Do iCnt = 1, dbsc(iCnttp)%nCntr
                  ia=ia+1
                  Do la = 1, nAngr(ia)+1
                     ka=ka+1
                     Old_Center=icent(k).eq.ia
                     Old_l=lnang(k).eq.(la-1)
                     If (idbg.gt.0) write(idbg,*) ' at numck', k,ia,
     &                   icent(k),la-1,lnang(k),ia,New_Center,New_l
                     If (Old_Center.and.Old_l) Then
                        numc(k)=nBasisr(ka)
                        numb(k)=nPrimr(ka)
                        icaddr(k)=ihelp(ia,la)+(kbias-1)*nPrimr(ka)
                        If (k.gt.1.and.kbias.eq.1) numck=numcl+numck
                        mcaddr(k)=numck
                        numcl=nPrimr(ka)
                     End If
                  End Do  ! la
               End Do     ! iCnt
            End Do        ! iCnttp
         End Do           ! iCont
      End Do              ! iSym
      k=0
      If (iPrint.ge.10.or.idbg.gt.0) then
         ic=0
         ip=0
         Do iSym = 1, nSym
            Write (idbg,*)        ' symmetry',iSym
            Write (idbg,*)        ' numb'
            Write (idbg,'(20i4)') (numb(i+ic),  i=1,nrBas(iSym))
            Write (idbg,*)        ' numc'
            Write (idbg,'(20i4)') (numc(i+ic),  i=1,nrBas(iSym))
            Write (idbg,*)        ' Pointer to contraction vector'
            Write (idbg,'(20i4)') (icaddr(i+ic),i=1,nrBas(iSym))
            Write (idbg,*)        ' mcaddr'
            Write (idbg,'(20i4)') (mcaddr(i+ic),i=1,nrBas(iSym))
            ic=ic+nrBas(iSym)
            ip=ip+nBas(iSym-1)
         End Do
      End If
*debugdebug
c     Write (*,*) (rCof(i),i=1,4)
*debugdebug
*
*     transform
*
      If (donorm) then
        kc=0
        kcL=0
        ibasL=0
        indbL=0
        kp=0
        Do iSym = 1, nSym
*        loop over contracted
          Do iBas = 1, nrBas(iSym)
            ibasL=ibasL+1
            jbasL=kcL
            Do jbas=1,ibas
               jbasL=jbasL+1
               kc=kc+1
               sum=0.0
               ipbasL=mcaddr(ibasL) -1
*              loop over primitives in this contracted function
               Do iPrim=1,numb(ibasL)
                  ipbasL=ipbasL+1
                  jpbasL=mcaddr(jbasL) -1
                  Do jPrim=1,numb(jbasL)
                     jpbasL=jpbasL+1
                     ilarge=Max(ipbasL,jpbasL)
                     ismall=Min(ipbasL,jpbasL)
                     kp=kp+1
                     indb=indbL+(ilarge*(ilarge-1))/2 + ismall
                     sum=sum+bint(indb)*rCof(icaddr(ibasL)+iPrim)
     &                  *rCof(icaddr(jbasL)+jPrim)
                     if (idbg.gt.0) Then
                        Write (idbg,*) indb,
     &                   icaddr(ibasL)+iPrim,icaddr(jbasL)+jPrim
                        Write (idbg,*) bint(indb),
     &                                 rCof(icaddr(ibasL)+iPrim),
     &                                 rCof(icaddr(jbasL)+jPrim)
                     End If
                  End Do  ! jprim
               End Do     ! iprim
               sint(kc)=sum
c              write(66,'(d25.14)') sum
            End Do
          End Do
          kcL=kcL+nrBas(iSym)
          indbL=indbL+(nBas(iSym-1)*(nBas(iSym-1)+1))/2
c         If (idbg.gt.0) Write(idbg,*) ipbasL,jpbasL,kp
        End Do

      else

        np = nBas(0)
        nc = nrBas(1)
        ! Scratch to square the mag ints and contract
        call mma_allocate(mag,np,np)
        call mma_allocate(u2c,np,nc)
        call mma_allocate(u2ct,nc,np)
        call mma_allocate(scr,np,nc)
        call mma_allocate(fin,nc,nc)
        call mma_allocate(pa,np)

        mag(:,:)= 0.0D0
        u2c(:,:) = 0.0D0
        u2ct(:,:) = 0.0D0
        scr(:,:) = 0.0D0
        fin(:,:) = 0.0D0
        pa(:) = 0.0D0

        ! Square the uncontracted ints
        mp = 0
        do ip = 1, np
          do jp = 1, ip
            mp = mp + 1
            mag(jp, ip) = bint(mp)
            ! Suboptimal
            mag(ip, jp) = bint(mp)
          enddo
        enddo

        ! Construct contraction matrix
        do icon = 1, nc
          istart = mcaddr(icon) - 1
          do iprim = 1, numb(icon)
!            write(6,'(A,4I5,F14.10)') 'u2c loop',
!     &     icon,iprim,istart,icaddr(icon),rCof(icaddr(icon)+iprim)
            u2c(istart+iprim, icon) = rCof(icaddr(icon)+iprim)
          enddo
        enddo

        ! U2C.T * UNCON * U2C => CON
        u2ct(:,:) = transpose(u2c)

        ! Since dgemm_ is causing trouble and this is trivial
        do iprim = 1, np
          do icon = 1, nc
            kpp=0.0d0
            pa(:) = u2ct(icon,:) * mag(:,iprim)
            do ip = 1, np
              kpp = kpp + pa(ip)
            enddo
            scr(iprim, icon) = kpp
          enddo
        enddo
!
        kp = 0
        do icon = 1, nc
          do jcon = 1, nc
            pa(:) = u2c(:,icon) * scr(:,jcon)
            kpp = 0.0d0
            do ip = 1, np
              kpp = kpp + pa(ip)
            enddo
            fin(icon, jcon) = kpp
          enddo
        enddo

        kp = 0
        do icon = 1, nc
          do jcon = 1, icon
            kp = kp + 1
            sint(kp) = fin(icon, jcon)
          enddo
        enddo

        call mma_deallocate(mag)
        call mma_deallocate(u2c)
        call mma_deallocate(u2ct)
        call mma_deallocate(scr)
        call mma_deallocate(fin)
        call mma_deallocate(pa)

        call cpu_time(finish)

      endif !donorm

      Return
      End
