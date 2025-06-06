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
      subroutine t3_bta_aac(nuga,nugc,kab,kca,kac,kc,la,lxa,lxc,mi,mij,
     $adim,cdim,N,noab_a,nuab_a,noab_b,nuab_b,lu,iasblock,nga,ngc,
     $oehi,oehk,oepa,oepc,enx,vab,vca,t1aa,t1ba,t1ac,t1bc,t3a,t3b,ifvo)
      implicit none
      real*8 one,zero,den,dena,denb,denc,enx,xx,yy
      parameter (one=1.d0,zero=0.d0)
      integer nadim,adim,ncdim,cdim,i,j,k,iasblock(5),lu(6),N
      integer noab_a,nuab_a,noab_b,nuab_b,nuga,nno_a,nnoab,nugc,
     $ngab_offset,ngca_offset,ngac_offset,nuga_offset,nugc_offset
      integer ias,jk,ij,ik,kj,ki,nga,ngc,a,b,c,ab,abb,bab
      integer iasabi,iascai,iasack
      real*8 kca(adim*cdim,N,*),kac(adim*cdim,N,*)
     $,kab(adim*(adim-1)/2,N,*)
     $,kc(*),la(N*adim,*),lxa(N*adim,*),lxc(N*cdim,*)
      real*8 t3a(*),t3b(*)
      real*8 mi(cdim*adim*(adim-1)/2,*),mij(*),
     $       vab(adim*(adim-1)/2,*),vca(adim*cdim,*)
      real*8 t1aa(noab_a,*),t1ba(noab_a,*),t1ac(noab_b,*),t1bc(noab_b,*)
      real*8 oehi(*),oehk(*),oepa(*),oepc(*)
      logical ifvo
C
C iasblock(1) > ka,kb,kc   iasblock(2) > la,lb iasblock(3) > lxa,lxc,lxb
      if(adim.eq.1)return
      nno_a=noab_a*(noab_a-1)/2
      nnoab=noab_a*noab_b
      nadim=adim*(adim-1)/2
      ncdim=adim*cdim
      nuga_offset=iasblock(1)*nuga*(nuga+1)/2
      nugc_offset=iasblock(1)*nuga*nugc
      ias=iasblock(2)*(nga-1)+1
      call multi_readir(la,nno_a*adim*N,lu(2),ias)
      ias=iasblock(3)*(nga-1)+1
      call multi_readir(lxa,nnoab*adim*N,lu(5),ias)
      ias=iasblock(3)*(ngc-1)+1
      call multi_readir(lxc,nnoab*cdim*N,lu(6),ias)
C vvoo ints reading
      ngab_offset=iasblock(4)*(nga*(nga-1)/2+nga-1)+1
      ias=iasblock(2)*nuga+ngab_offset
      call multi_readir(vab,nno_a*nadim,lu(2),ias)
      ngca_offset=iasblock(5)*(nugc*(nga-1)+ngc-1)+1
      ias=iasblock(2)*nuga+iasblock(4)*nuga*(nuga+1)/2+ngca_offset
      call multi_readir(vca,nnoab*adim*cdim,lu(2),ias)
!!      ngac_offset=iasblock(5)*(nuga*(ngc-1)+nga-1)+1
C end readin vvoo ints
      ngab_offset=iasblock(1)*(nga*(nga-1)/2+nga-1)+1
      ngac_offset=iasblock(1)*(nuga*(ngc-1)+nga-1)+1
      ngca_offset=iasblock(1)*(nugc*(nga-1)+ngc-1)+1

                do i=1,noab_a
                iasabi=(i-1)*nuga_offset+ngab_offset
                call multi_readir(kab(1,1,i),N*nadim,lu(1),iasabi)
                enddo
                do i=1,noab_a
                iascai=(i-1)*nugc_offset+ngca_offset
                call multi_readir(kca(1,1,i),N*ncdim,lu(3),iascai)
                enddo
         do k=1,noab_b
         do i=1,noab_a
         ik=(k-1)*noab_a +i
         ki=(i-1)*noab_b +k
C  K_ab^ir x L_rc^ik     cba
      call DGEMM_('T','T',cdim,nadim,N,one,lxc(1,ik),N,kab(1,1,i),nadim,
     $                      zero,mi(1,i),cdim)
C
C  K_ac^ir x L_rb^ki     cab
      call DGEMM_('N','N',ncdim,adim,N,one,kca(1,1,i),ncdim,lxa(1,ki),N,
     $                    zero,t3b,ncdim)
        ab=1
      do a=2,adim
      abb=(a-1)*cdim+1
      bab=(a-1)*ncdim+1
      do b=1,a-1
      call daxpy_(cdim,-1.d0,t3b(abb),1,mi(ab,i),1)
      call daxpy_(cdim,1.d0,t3b(bab),1,mi(ab,i),1)
      ab=ab+cdim
      abb=abb+ncdim
      bab=bab+cdim
      enddo
      enddo
          enddo    ! i
! end prefactors
         iasack=(k-1)*nugc_offset+ngac_offset
         call multi_readir(kac,N*ncdim,lu(4),iasack)
         ij=0
         do i=2,noab_a
         ki=(i-1)*noab_b +k
         ik=(k-1)*noab_a +i
         kj=k-noab_b
         jk=(k-1)*noab_a
         do j=1,i-1
         ij=ij+1
         kj=kj+noab_b
         jk=jk+1
C  K_bc^kr x L_ra^ij
         call DGEMM_('N','N',ncdim,adim,N,one,kac,ncdim,la(1,ij),N,
     $                     zero,t3a,ncdim)
C transpose the first two inicesd
         ab=1
         do a=1,adim
         call transm(t3a(ab),t3b(ab),adim,cdim)
         ab=ab+ncdim
         enddo
C  K_ab^ir x L_rc^jk -K_ab^jr x L_rc^ik
         call vsub(kab(1,1,j),1,kab(1,1,i),1,kc,1,N*nadim)
         call vadd(lxc(1,jk),1,lxc(1,ik),1,mij,1,N*cdim)
         call DGEMM_('T','T',cdim,nadim,N,one,mij,N,kc,nadim,
     $                     zero,t3a,cdim)
C  K_ab^ir x L_rc^jk
!!         call DGEMM_('T','T',cdim,nadim,N,one,lxc(1,jk),N,ka,nadim,
!!     $                     zero,t3a,cdim)
C  -K_ab^jr x L_rc^ik
!!         call DGEMM_('T','T',cdim,nadim,N,-one,lxc(1,ik),N,kb,nadim,
!!     $                      one,t3a,cdim)
C
C  K_bc^ir x L_ra^kj -K_bc^jr x L_ra^ki
         call vsub(kca(1,1,j),1,kca(1,1,i),1,kc,1,N*ncdim)
         call vadd(lxa(1,kj),1,lxa(1,ki),1,mij,1,N*adim)
         call DGEMM_('N','N',ncdim,adim,N,one,kc,ncdim,mij,N,
     $                     one,t3b,ncdim)
C  K_bc^ir x L_ra^kj
!!         call DGEMM_('N','N',ncdim,adim,N,one,ka,ncdim,lxa(1,kj),N,
!!     $                     one,t3b,ncdim)
C  -K_bc^jr x L_ra^ki
!!         call DGEMM_('N','N',ncdim,adim,N,-one,kb,ncdim,lxa(1,ki),N,
!!     $                     one,t3b,ncdim)
       ab=1
      do a=2,adim
      abb=(a-1)*cdim+1
      bab=(a-1)*ncdim+1
      do b=1,a-1
      call daxpy_(cdim,-1.d0,t3b(abb),1,t3a(ab),1)
      call daxpy_(cdim,1.d0,t3b(bab),1,t3a(ab),1)
      ab=ab+cdim
      abb=abb+ncdim
      bab=bab+cdim
      enddo
      enddo
!!
      call daxpy_(nadim*cdim,-1.d0,mi(1,i),1,t3a,1)
      call daxpy_(nadim*cdim,1.d0,mi(1,j),1,t3a,1)
         den=oehi(i)+oehi(j)+oehk(k)
      ab=0
      do a=2,adim
      dena=den-oepa(a)
      do b=1,a-1
      denb=dena-oepa(b)
      do c=1,cdim
      denc=denb-oepc(c)
      ab=ab+1
      xx=t3a(ab)
      yy=xx/denc
      enx=enx+yy*xx
      t3a(ab)=yy
!!      t1aa(j,a)=t1aa(j,a)-yy*vca((a-1)*cdim+c,ki)
      enddo
      enddo
      enddo
      call expa2_uhf(t3a,cdim,adim,-1,t3b)
         call DGEMM_('N','T',1,cdim,nadim, one,vab(1,ij),1,
     $  t3a,cdim,one,t1ac(k,1),noab_b)
         call DGEMM_('N','N',1,adim,ncdim, one,vca(1,kj),1,
     $  t3b,ncdim,one,t1aa(i,1),noab_a)
         call DGEMM_('N','N',1,adim,ncdim,-one,vca(1,ki),1,
     $  t3b,ncdim,one,t1aa(j,1),noab_a)
C ccsd(T) part t2*t3
         if(ifvo) then
         call DGEMM_('N','T',1,cdim,nadim, one,kab(1,i,j),1,
     $  t3a,cdim,one,t1bc(k,1),noab_b)
         call DGEMM_('N','N',1,adim,ncdim,-one,kca(1,k,j),1,
     $  t3b,ncdim,one,t1ba(i,1),noab_a)
         call DGEMM_('N','N',1,adim,ncdim, one,kca(1,k,i),1,
     $  t3b,ncdim,one,t1ba(j,1),noab_a)
         endif
         enddo !j
         enddo !i
         enddo !k
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nuab_a)
        call Unused_integer(nuab_b)
      end if
      end
