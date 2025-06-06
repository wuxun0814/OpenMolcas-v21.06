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
        subroutine gather_t2(t2,t2_tmp,tmp)
c
c temporary routine. In future T2 block structure will be merged
c                    with the block structure of the (T) code
c
c this routine do :
c
c cycle through T2XY files and gather them into on T2 array
c T2(nv_beta,nv_alpha,no_beta,no_alpha)
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        integer a,b,dima,dimb
cmp        integer length
        integer length
        integer lasta,lastb
c
        real*8 t2(*),tmp(*),t2_tmp(*)
        integer a_tmp,b_tmp
c
cmp        write (6,*)
cmp        write (6,*) '------ DimGrpaR ------'
cmp        write (6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
cmp        write (6,*)
c
        do a=1,NvGrp
        do b=1,a
c
cmp@@        dima=nv/NvGrp
        dima=DimGrpaR(a)
        dimb=DimGrpaR(b)
c
cmp        write (6,'(A,i3,i3,2x,A6)') 'a,b,T2Name(a,b) ',a,b,T2Name(a,b)
c
         if (a.eq.b) then  ! a=b
c open the pertinent file
c
cmp@@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
        dimb=dima
c
cmp         write (6,*) 'dima = ',dima
         length=(dima*(dima+1)*no*no)/2
cmp         write (6,*) 'length = ',length
cmp!         write (6,*) 'file size (g77) = ',16+length*8
cmp         write (6,*) 'file size (ifort) = ',8+length*8
c
        call GetX_t3 (tmp,length,LunAux,T2Name(a,b),1,1)
c
         else ! a>b
c open the pertinent file
c
cmp@@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
cmp@@         if (b.eq.NvGrp) dimb=nv-((NvGrp-1)*dimb)
c
cmp         write (6,*) 'dima = ',dima
cmp         write (6,*) 'dimb = ',dimb
         length=dima*dimb*no*no
cmp         write (6,*) 'length = ',length
cmp!         write (6,*) 'file size (g77) = ',16+length*8
cmp         write (6,*) 'file size (ifort) = ',8+length*8
c
        call GetX_t3 (tmp,length,LunAux,T2Name(a,b),1,1)
c
         end if
c
c add its contents to the t2 array
c
cmp@@        lasta=(a-1)*(nv/NvGrp)
cmp@@        lastb=(b-1)*(nv/NvGrp)
cmp@@
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
cmp@@
c
cmp        write (6,'(A,2(i5,2x),2(i3,x))') 'lasta, lastb, a, b = ',
cmp     & lasta,lastb,a,b
c
        if (a.eq.b) then ! expand and map
c expand and map l2_1 (a',b',m) <- tmp (m,ab')
        call expand4_12 (tmp,t2_tmp,dima,no,no)
        call grow_t2neq(t2,t2_tmp,dima,dimb,nv,no,
     & lasta,lastb)
        else
        call grow_t2neq(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb)
        end if
c
        end do
        end do
c
        return
        end
