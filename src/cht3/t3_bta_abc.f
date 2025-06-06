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
      subroutine t3_bta_abc(nuga,nugc,kab,kcb,kca,kac,kbc,kc,la,lb,lxa,
     $lxb,lxc,mi,mij,adim,bdim,cdim,N,noab_a,nuab_a,noab_b,nuab_b,lu,
     $iasblock,nga,ngb,ngc,oehi,oehk,oepa,oepb,oepc,enx,vab,vcb,vca,
     $t1aa,t1ba,t1ab,t1bb,t1ac,t1bc,t3a,t3b,ifvo)
      implicit none
cmp
c     integer imm
cmp
      real*8 one,zero,den,dena,denb,denc,enx,xx,yy
c     real*8 sumt3
      parameter (one=1.d0,zero=0.d0)
      integer nadim,adim,ncdim,cdim,bdim,nbdim,i,j,k,iasblock(5),lu(6),N
      integer noab_a,nuab_a,noab_b,nuab_b,nuga,nno_a,nnoab,nugc,
     $ngab_offset,ngca_offset,ngac_offset,nuga_offset,nugc_offset,
     $ngcb_offset,ngbc_offset
      integer ias,jk,ij,ik,kj,ki,nga,ngb,ngc,a,b,c,ab,ba
      integer iasabi,iascai,iasack,iascbi,iasbck
      real*8 kab(adim*bdim,N,*),kcb(cdim*bdim,N,*),kca(cdim*adim,N,*)
     $,kac(cdim*adim,N),kbc(cdim*bdim,N),kc(*),lxb(N*bdim,*)
      real*8 la(N*adim,*),lb(N*bdim,*),lxa(N*adim,*),lxc(N*cdim,*),
     $t3a(*),t3b(*),vab(adim*bdim,*),vca(adim*cdim,*),vcb(bdim*cdim,*)
      real*8 mi(adim*bdim*cdim,*),mij(*)
      real*8 t1aa(noab_a,*),t1ba(noab_a,*),t1ac(noab_b,*)
     $        ,t1bc(noab_b,*),t1ab(noab_a,*),t1bb(noab_a,*)
      real*8 oehi(*),oehk(*),oepa(*),oepb(*),oepc(*)
      logical ifvo
cmp
c
C iasblock(1) > ka,kb,kc   iasblock(2) > la,lb iasblock(3) > lxa,lxc,lxb
!!      sumt3=0.d0
      nno_a=noab_a*(noab_a-1)/2
      nnoab=noab_a*noab_b
      nadim=adim*bdim
      nbdim=bdim*cdim
      ncdim=adim*cdim
      nuga_offset=iasblock(1)*nuga*(nuga+1)/2
      nugc_offset=iasblock(1)*nuga*nugc
      ias=iasblock(2)*(nga-1)+1
      call multi_readir(la,nno_a*adim*N,lu(2),ias)
      ias=iasblock(2)*(ngb-1)+1
      call multi_readir(lb,nno_a*bdim*N,lu(2),ias)
      ias=iasblock(3)*(nga-1)+1
      call multi_readir(lxa,nnoab*adim*N,lu(5),ias)
      ias=iasblock(3)*(ngb-1)+1
      call multi_readir(lxb,nnoab*bdim*N,lu(5),ias)
      ias=iasblock(3)*(ngc-1)+1
      call multi_readir(lxc,nnoab*cdim*N,lu(6),ias)
C
C vvoo ints reading
      ngab_offset=iasblock(4)*(nga*(nga-1)/2+ngb-1)+1
      ias=iasblock(2)*nuga+ngab_offset
      call multi_readir(vab,nno_a*nadim,lu(2),ias)
      ngca_offset=iasblock(5)*(nugc*(nga-1)+ngc-1)+1
      ias=iasblock(2)*nuga+iasblock(4)*nuga*(nuga+1)/2+ngca_offset
      call multi_readir(vca,nnoab*ncdim,lu(2),ias)
      ngcb_offset=iasblock(5)*(nugc*(ngb-1)+ngc-1)+1
      ias=iasblock(2)*nuga+iasblock(4)*nuga*(nuga+1)/2+ngcb_offset
      call multi_readir(vcb,nnoab*nbdim,lu(2),ias)
C end readin vvoo ints
C
      ngab_offset=iasblock(1)*(nga*(nga-1)/2+ngb-1)+1
      ngac_offset=iasblock(1)*(nuga*(ngc-1)+nga-1)+1
      ngbc_offset=iasblock(1)*(nuga*(ngc-1)+ngb-1)+1
      ngca_offset=iasblock(1)*(nugc*(nga-1)+ngc-1)+1
      ngcb_offset=iasblock(1)*(nugc*(ngb-1)+ngc-1)+1
C saves reading:
      do i=1,noab_a
      iasabi=(i-1)*nuga_offset+ngab_offset
      call multi_readir(kab(1,1,i),N*nadim,lu(1),iasabi)
      enddo
      do i=1,noab_a
      iascai=(i-1)*nugc_offset+ngca_offset
      call multi_readir(kca(1,1,i),N*ncdim,lu(3),iascai)
cmp
cmp        write (*,*) 'ze tak teraz ju dam',i
cmp        call check_mat(kca(1,1,i),1,N*ncdim)
cmp
      enddo
      do i=1,noab_a
         iascbi=(i-1)*nugc_offset+ngcb_offset
         call multi_readir(kcb(1,1,i),N*nbdim,lu(3),iascbi)
      enddo
         do k=1,noab_b
         iasbck=(k-1)*nugc_offset+ngbc_offset
         call multi_readir(kbc,N*nbdim,lu(4),iasbck)
         iasack=(k-1)*nugc_offset+ngac_offset
         call multi_readir(kac,N*ncdim,lu(4),iasack)
C start calculating prefactors:
cmp
         do i=1,noab_a
         ik=(k-1)*noab_a +i
         ki=(i-1)*noab_b +k
cmp
C  K_ab^ir x L_rc^ik     cba
       call DGEMM_('T','T',cdim,nadim,N,one,lxc(1,ik),N,
     &  kab(1,1,i),nadim,
     $                      zero,mi(1,i),cdim)
C  K_bc^ir x L_ra^ki          cba
       call DGEMM_('N','N',nbdim,adim,N,one,kcb(1,1,i),nbdim,
     &  lxa(1,ki),N,
     $                     one ,mi(1,i),nbdim)
cmp
cmp        do imm=0,adim*nbdim-1
cmp        if (abs(mi(1+imm,i)).gt.10000) then
cmp          write (*,*) 'uz mi dojebane 2',imm+1,i,mi(1+imm,i)
cmp          stop
cmp        end if
cmp        end do
cmp
C
C  K_ac^ir x L_rb^ki     cab
      call DGEMM_('N','N',ncdim,bdim,N,one,kca(1,1,i),ncdim,
     &    lxb(1,ki),N,
     $                    zero,t3b,ncdim)
cmp
         ab=1
      do a=1,adim
      ba=(a-1)*cdim+1
      do b=1,bdim
      call daxpy_(cdim,-1.d0,t3b(ba),1,mi(ab,i),1)
      ab=ab+cdim
      ba=ba+ncdim
      enddo
      enddo
          enddo    ! i
cmp
! end prefactors
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
cmp
C  K_ac^kr x L_rb^ij        bac
         call DGEMM_('T','T',bdim,ncdim,N,-one,lb(1,ij),N,kac,ncdim,
     $                     zero,t3b,bdim)
C  K_bc^kr x L_ra^ij        abc
         call DGEMM_('T','T',adim,nbdim,N,one,la(1,ij),N,kbc,nbdim,
     $                     zero,t3a,adim)
C transpose the first two indices
cmp
         ab=1
         do c=1,cdim
         ba=(c-1)*nadim+1
         do b=1,bdim
         call daxpy_(adim,1.d0,t3a(ab),1,t3b(ba),bdim)
         ba=ba+1
         ab=ab+adim
         enddo
         enddo
C t3b  bac
         call transm(t3b,t3a,nadim,cdim)
C      cba in t3a
C K_ab^ir x L_rc^jk  -K_ab^jr x L_rc^ik
         call vsub(kab(1,1,j),1,kab(1,1,i),1,kc,1,N*nadim)
C         call daxpy_(N*nadim,-one,kb,1,ka,1)
         call vadd(lxc(1,jk),1,lxc(1,ik),1,mij,1,N*cdim)
         call DGEMM_('T','T',cdim,nadim,N,one,mij,N,kc,nadim,
     $                      one,t3a,cdim)
C  K_ab^ir x L_rc^jk
!!         call DGEMM_('T','T',cdim,nadim,N,one,lxc(1,jk),N,ka,nadim,
!!     $                      one,t3a,cdim)
C  -K_ab^jr x L_rc^ik
!!         call DGEMM_('T','T',cdim,nadim,N,-one,lxc(1,ik),N,kb,nadim,
!!     $                      one,t3a,cdim)
C
C  K_bc^ir x L_ra^kj -K_bc^jr x L_ra^ki         cba
         call vsub(kcb(1,1,j),1,kcb(1,1,i),1,kc,1,N*nbdim)
         call vadd(lxa(1,kj),1,lxa(1,ki),1,mij,1,N*adim)
         call DGEMM_('N','N',nbdim,adim,N,one,kc,nbdim,mij,N,
     $                     one,t3a,nbdim)
C  K_bc^ir x L_ra^kj          cba
!!         call DGEMM_('N','N',nbdim,adim,N,one,ka,nbdim,lxa(1,kj),N,
!!     $                     one,t3a,nbdim)
C  -K_bc^jr x L_ra^ki         cba
!!         call DGEMM_('N','N',nbdim,adim,N,-one,kb,nbdim,lxa(1,ki),N,
!!     $                     one,t3a,nbdim)
C  K_ac^ir x L_rb^kj  -K_ac^jr x L_rb^ki      cab
         call vsub(kca(1,1,j),1,kca(1,1,i),1,kc,1,N*ncdim)
         call vadd(lxb(1,kj),1,lxb(1,ki),1,mij,1,N*bdim)
         call DGEMM_('N','N',ncdim,bdim,N,one,kc,ncdim,mij,N,
     $                    zero,t3b,ncdim)
C  K_ac^ir x L_rb^kj      cab
!!         call DGEMM_('N','N',ncdim,bdim,N,one,ka,ncdim,lxb(1,kj),N,
!!     $                    zero,t3b,ncdim)
C  -K_ac^jr x L_rb^ki     cab
!!         call DGEMM_('N','N',ncdim,bdim,N,-one,kb,ncdim,lxb(1,ki),N,
!!     $                     one,t3b,ncdim)
cmp
         ab=1
      do a=1,adim
      ba=(a-1)*cdim+1
      do b=1,bdim
      call daxpy_(cdim,-1.d0,t3b(ba),1,t3a(ab),1)
      ab=ab+cdim
      ba=ba+ncdim
      enddo
      enddo
c
cmp        do a=1,adim
cmp        do b=1,bdim
cmp        do c=1,cdim
cmp        ab=c+(b-1)*cdim+(a-1)*bdim*cdim
cmp        ba=c+(a-1)*cdim+(b-1)*adim*cdim
cmp        t3a(ab)=t3a(ab)-t3b(ba)
cmp        end do
cmp        end do
cmp        end do
c
cmp
c
      call daxpy_(nadim*cdim,-1.d0,mi(1,i),1,t3a,1)
      call daxpy_(nadim*cdim,1.d0,mi(1,j),1,t3a,1)
         den=oehi(i)+oehi(j)+oehk(k)
      ab=0
      do a=1,adim
      ba=ab+1
      dena=den-oepa(a)
      do b=1,bdim
      denb=dena-oepb(b)
      do c=1,cdim
      denc=denb-oepc(c)
      ab=ab+1
cmp        if ((i.eq.14).and.(j.eq.1).and.(k.eq.1)) then
cmp          write (*,'(A,3(i5,x),3(f18.10,x))')
cmp     &               'a,b,c, oepa(a), oepb(b), oepc(c)',
cmp     &                a,b,c,oepa(a),oepb(b),oepc(c)
cmp          write (*,*) 'denc = ',denc
cmp          write (*,*) 'ab, t3a(ab) ',ab,t3a(ab)
cmp        end if
      xx=t3a(ab)
      yy=xx/denc
!!      sumt3=sumt3+xx
      enx=enx+yy*xx
cmp        if ((i.eq.14).and.(j.eq.1).and.(k.eq.1)) then
cmp          write (*,'(A,3(f18.10,x))') 'xx,yy,enx = ',xx,yy,enx
cmp        end if
      t3a(ab)=yy
      enddo
      enddo
      call transm(t3a(ba),t3b(ba),cdim,bdim)
      enddo
      call DGEMM_('N','T',1,cdim,nadim, one,vab(1,ij),1,
     $  t3a,cdim,one,t1ac(k,1),noab_b)
       call DGEMM_('N','N',1,adim,nbdim,-one,vcb(1,kj),1,
     $  t3a,nbdim,one,t1aa(i,1),noab_a)
       call DGEMM_('N','N',1,adim,nbdim,one,vcb(1,ki),1,
     $  t3a,nbdim,one,t1aa(j,1),noab_a)
      call DGEMM_('N','T',1,bdim,ncdim,-one,vca(1,ki),1,
     $  t3b,bdim,one,t1ab(j,1),noab_a)
      call DGEMM_('N','T',1,bdim,ncdim,one,vca(1,kj),1,
     $  t3b,bdim,one,t1ab(i,1),noab_a)
       if(ifvo) then
      call DGEMM_('N','T',1,cdim,nadim, one,kab(1,i,j),1,
     $  t3a,cdim,one,t1bc(k,1),noab_b)
       call DGEMM_('N','N',1,adim,nbdim,one,kcb(1,k,j),1,
     $  t3a,nbdim,one,t1ba(i,1),noab_a)
       call DGEMM_('N','N',1,adim,nbdim,-one,kcb(1,k,i),1,
     $  t3a,nbdim,one,t1ba(j,1),noab_a)
      call DGEMM_('N','T',1,bdim,ncdim,one,kca(1,k,i),1,
     $  t3b,bdim,one,t1bb(j,1),noab_a)
      call DGEMM_('N','T',1,bdim,ncdim,-one,kca(1,k,j),1,
     $  t3b,bdim,one,t1bb(i,1),noab_a)
      endif
cmp        write (*,'(3(i5,x),f18.10)') i,j,k,enx
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
