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
* Copyright (C) 2016,2017, Giovanni Li Manni                           *
*               2019-2021, Oskar Weser                                 *
*               2021, Werner Dobrautz                                  *
************************************************************************

      module fciqmc_read_RDM
      use definitions, only: wp
      use para_info, only: myRank
      use general_data, only : nActEl
! Note that two_el_idx_flatten has also out parameters.
      use index_symmetry, only : two_el_idx_flatten
      use CI_solver_util, only: CleanMat
      use linalg_mod, only: abort_, verify_

      implicit none

      private
      public :: read_neci_RDM, cleanup, read_neci_GUGA_RDM
      contains

!>  @brief
!>    Read NECI RDM files
!>
!>  @author Werner Dobrautz, Oskar Weser
!>
!>  @details
!>  Read the spin-free TwoRDM file written by the GUGA-NECI
!>  implementation and transfer them to Molcas.
!>  The spin density is set to zero, because spin projection
!>  is not properly defined in the GUGA framework.

!>  @paramin[out] DMAT Average spin-free 1 body density matrix
!>  @paramin[out] DSPN spin-dependent 1-RDM (set to zero)
!>  @paramin[out] PSMAT Average spin-free 2 body density matrix
!>  @paramin[out] PAMAT 'fake' Average antisymm. 2-dens matrix
      subroutine read_neci_GUGA_RDM(DMAT, DSPN, PSMAT, PAMAT)
      real(wp), intent(out) :: DMAT(:), DSPN(:), PSMAT(:), PAMAT(:)
      integer :: file_id, isfreeunit, i
      logical :: tExist
      real(wp) :: RDMval

      if (myRank /= 0) then
          call bcast_2RDM("PSMAT")
          call bcast_2RDM("PAMAT")
          call bcast_2RDM("DMAT")
      end if

      call f_Inquire('PSMAT',tExist)
      call verify_(tExist, 'PSMAT does not exist')
      call f_Inquire('PAMAT',tExist)
      call verify_(tExist, 'PAMAT does not exist')
      call f_Inquire('DMAT',tExist)
      call verify_(tExist, 'DMAT does not exist')

      PSMAT(:) = 0.0_wp; PAMAT(:) = 0.0_wp;
      DMAT(:) = 0.0_wp; DSPN(:) = 0.0_wp;

      file_id = IsFreeUnit(11)
      call Molcas_Open(file_id, 'PSMAT')
          do while (read_line(file_id, i, RDMval))
            psmat(i) = RDMval
          end do
      close(file_id)

      file_id = IsFreeUnit(11)
      call Molcas_Open(file_id, 'PAMAT')
          do while (read_line(file_id, i, RDMval))
            pamat(i) = RDMval
          end do
      close(file_id)

      file_id = IsFreeUnit(11)
      call Molcas_Open(file_id, 'DMAT')
          do while (read_line(file_id, i, RDMval))
            dmat(i) = RDMval
          end do
      close(file_id)

      ! Clean evil non-positive semi-definite matrices,
      ! by clamping the occupation numbers between 0 and 2.
      ! DMAT is intent(inout)
      call cleanMat(DMAT)

      contains

          !> Read a line from GUGA NECI RDM file with
          !> `file_id` as unit and parse it into i and RDMval.
          !> Return true, if line was successfully read and false
          !> if EOF reached.
          !> Aborts if IO error happens.
          !> **i and RDMval are undefined, if functions returns false**!
          logical function read_line(file_id, i, RDMval)
              integer, intent(in) :: file_id
              integer, intent(out) :: i
              real(wp), intent(out) :: RDMval
              integer :: iread
              read_line = .false.
              read(file_id, "(I6,G25.17)", iostat=iread) i, RDMval
              if (iread > 0) then
                  call abort_('Error in read_next')
              else if (is_iostat_end(iread)) then
                  ! Let's try to cause an error, if code uses undefined
                  ! values for i and RDMval
                  i = huge(i)
                  RDMval = huge(RDMval)
              else
                  read_line = .true.
              end if
          end function

      end subroutine read_neci_GUGA_RDM

!>  @brief
!>    Start and control FCIQMC.
!>
!>  @author Giovanni Li Manni, Oskar Weser
!>
!>  @details
!>  Read TwoRDM files written by NECI and transfer them to Molcas.
!>  Neci can have some intermediate spin-resolved/spin-free RDMs where basically aaaa contains
!>  average of aaaa and bbbb, abab contains average of abab and baba...
!>  This is ok for CASSCF but not ok for spin-resolved properties, in which case the completely
!>  spin-resolved RDMs need to be read-in.
!>  In principle, NECI could also evaluate and store completely spin-free matrices.
!>  In that case only a reordering following Molcas convention is necessary.
!>
!>  @paramin[out] DMAT Average 1 body density matrix
!>  @paramin[out] DSPN Average spin 1-dens matrix
!>  @paramin[out] PSMAT Average symm. 2-dens matrix
!>  @paramin[out] PAMAT Average antisymm. 2-dens matrix
      subroutine read_neci_RDM(DMAT, DSPN, PSMAT, PAMAT)
      use Para_Info, only: MyRank
#include "output_ras.fh"
      real*8, intent(out) :: DMAT(:), DSPN(:), PSMAT(:), PAMAT(:)
      integer :: iUnit, isfreeunit, p, q, r, s, pq, rs, ps, rq, psrq,
     &  pqrs, iread, norb, iprlev
      logical :: tExist, switch
      real*8 :: fac, RDMval, fcnacte
      real*8 :: D_alpha(size(DMAT)), D_beta(size(DMAT))


      iprlev = iprloc(1)
      if(iprlev == debug) then
        write(6,*) 'Rank of process: ', MyRank
      end if
      switch = .true.
*     ^  This variable must become a keyword for discriminating spin-resolved from spin-free input RDMs.
*     ^  For now when .false. it is assumed that 3 files (only   aaaa, abab and abba) are fed.
*     ^  ....... when .true.  it is assumed that 6 files (adding bbbb, baba and baab) are fed.
*********************************************************************************
* Broadcasting TwoRDM generated by QMC code in master node into all processors. *
*********************************************************************************
! NOTE(Giovanni, Oskar): The suffix ".1" corresponds to the root.
!   For state averaged the .1 shall be replaced by iroot. (irdm in NECI)
      if(myRank /= 0) then
        call bcast_2RDM("TwoRDM_aaaa.1")
        call bcast_2RDM("TwoRDM_aaaa.1")
        call bcast_2RDM("TwoRDM_abab.1")
        call bcast_2RDM("TwoRDM_abba.1")
        call bcast_2RDM("TwoRDM_bbbb.1")
        call bcast_2RDM("TwoRDM_baba.1")
        call bcast_2RDM("TwoRDM_baab.1")
      end if
**********************************************************************************
******************************** existency check *********************************
**********************************************************************************
      call f_Inquire('TwoRDM_aaaa.1',tExist)
      if(.not.tExist) goto 123
      call f_Inquire('TwoRDM_abab.1',tExist)
      if(.not.tExist) goto 123
      call f_Inquire('TwoRDM_abba.1',tExist)
      if(.not.tExist) goto 123
      if(switch) then
        call f_Inquire('TwoRDM_bbbb.1',tExist)
        if(.not.tExist) goto 123
        call f_Inquire('TwoRDM_baba.1',tExist)
        if(.not.tExist) goto 123
        call f_Inquire('TwoRDM_baab.1',tExist)
        if(.not.tExist) goto 123
      end if

      D_alpha(:) = 0.0d0
      D_beta(:) = 0.0d0
      PSMAT(:) = 0.0d0
      PAMAT(:) = 0.0d0

      fac = merge(0.5d0, 1.0d0, switch)
      fcnacte = 1.0d0 / dble(nactel - 1)

*******************************************************************************************
*************************** Processing TwoRDM-AAAA ****************************************
*******************************************************************************************
      iUnit = IsFreeUnit(11)
      call Molcas_Open(iUnit, 'TwoRDM_aaaa.1')
      Rewind(iUnit)
      IF(IPRLEV >= DEBUG) THEN
        write(6,*) '    p     q     r     s    pq    rs   pqrs        ',
     &  'RDMval                PSMAT                   PAMAT'
        write(6,*) ' ********************** AAAA ****************** '
      end if
      do
******************* processing as PQRS ***********************
**************************************************************
        read(iUnit, "(4I6,G25.17)", iostat=iread) s, q, r, p, RDMval
        if(iread /= 0) exit
        pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
        PSMAT(pqrs) = PSMAT(pqrs) + fac * RDMval
        if (r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs) + fac * RDMval
        if (p > q.and.r > s) PAMAT(pqrs) = PAMAT(pqrs) + fac * RDMval
        if (p > q.and.r < s) PAMAT(pqrs) = PAMAT(pqrs) - fac * RDMval
        IF(IPRLEV >= DEBUG) THEN
           write(6,'(7I6,3G25.17)')
     &     p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
        END IF
******* Contribution to D_alpha (not final):
        if (p == q) D_alpha(rs) = D_alpha(rs) + RDMval
        if (r == s) D_alpha(pq) = D_alpha(pq) + RDMval
******************* processing as PSRQ ***********************
**************************************************************
        psrq = two_el_idx_flatten(p, s, r, q, ps, rq)
******* Contribution to PSMAT and PAMAT:
        if (r <= q) then
          PSMAT(psrq) = PSMAT(psrq) - fac * RDMval
          if (r /= q) PAMAT(psrq) = PAMAT(psrq) + fac * RDMval
        end if
        if (r > q) then
          PSMAT(psrq) = PSMAT(psrq) - fac * RDMval
          PAMAT(psrq) = PAMAT(psrq) - fac * RDMval
        end if
        IF(IPRLEV >= DEBUG) THEN
          write(6,'(7I6,3G25.17)')
     &      p,s,r,q,ps,rq,psrq, RDMval,PSMAT(psrq),PAMAT(psrq)
        END IF
******* Contribution to D_alpha (not final):
* The minus sign comes from the fact that in NECI these elements have opposite sign
* compared to the element in normal order, that is d_pqrs = -d_psrq.
        if (p == s) D_alpha(rq) = D_alpha(rq) - RDMval
        if (r == q) D_alpha(ps) = D_alpha(ps) - RDMval
      end do
      close(iunit)

*******************************************************************************************
*************************** Processing TwoRDM-BBBB ****************************************
*******************************************************************************************
      if (switch) then
        iUnit=IsFreeUnit(11)
        Call Molcas_Open(iUnit,'TwoRDM_bbbb.1')
        Rewind(iUnit)
        if (IPRLEV >= DEBUG) then
         write(6,*) '    p     q     r     s    pq    rs   pqrs       ',
     &   'RDMval                PSMAT                   PAMAT'
         write(6,*) ' ********************** BBBB ****************** '
        end if
        do
******************* processing as PQRS ***********************
**************************************************************
          read(iUnit,"(4I6,G25.17)",iostat=iread) s, q, r, p, RDMval
          if(iread /= 0) exit
          pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
          PSMAT(pqrs) = PSMAT(pqrs) + fac * RDMval
          if(r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
          if(p > q.and.r > s) PAMAT(pqrs) = PAMAT(pqrs) + fac * RDMval
          if(p > q.and.r < s) PAMAT(pqrs) = PAMAT(pqrs) - fac * RDMval
          IF(IPRLEV >= DEBUG) THEN
            write(6,'(7I6,3G25.17)')
     &         p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
          END IF
******* Contribution to D_beta (not final):
          if (p == q) D_beta(rs) = D_beta(rs) + RDMval
          if (r == s) D_beta(pq) = D_beta(pq) + RDMval
******************* processing as PSRQ ***********************
**************************************************************
          psrq = two_el_idx_flatten(p, s, r, q, ps, rq)
******* Contribution to PSMAT and PAMAT:
          if(r <= q) then
            PSMAT(psrq) = PSMAT(psrq) - fac*RDMval
            if(r /= q) PAMAT(psrq) = PAMAT(psrq) + fac*RDMval
          end if
          if(r > q) then
            PSMAT(psrq) = PSMAT(psrq) - fac*RDMval
            PAMAT(psrq) = PAMAT(psrq) - fac*RDMval
          end if
          IF(IPRLEV >= DEBUG) THEN
            write(6,'(7I6,3G25.17)')
     &        p,s,r,q,ps,rq,psrq, RDMval,PSMAT(psrq),PAMAT(psrq)
          END IF
******* Contribution to D_beta (not final):
* The minus sign comes from the fact that in NECI these elements have opposite sign
* compared to the element in normal order, that is d_pqrs = -d_psrq.
          if(p == s) D_beta(rq)=D_beta(rq)-RDMval
          if(r == q) D_beta(ps)=D_beta(ps)-RDMval
        end do
        close(iunit)
      end if ! End statement for spin-resolved RDMs.
*******************************************************************************************
*************************** Processing TwoRDM-ABAB ****************************************
*******************************************************************************************
      iUnit=IsFreeUnit(11)
      Call Molcas_Open(iUnit,'TwoRDM_abab.1')
      Rewind(iUnit)
      IF(IPRLEV >= DEBUG) THEN
       write(6,*) ' ********************** ABAB ****************** '
      END IF
      do
        read(iUnit,"(4I6,G25.17)",iostat=iread) s,q,r,p,RDMval
        if(iread /= 0) exit
        pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
        PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
        if(r > s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
        if(r < s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
        if(r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
        if (IPRLEV >= DEBUG) then
          write(6,'(7I6,3G25.17)')
     &        p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
        end if
******* Contribution to D_alpha and D_beta (not final):
        if (p == q) D_alpha(rs) = D_alpha(rs) + RDMval
        if (r == s .and. p /= r) D_beta(pq) = D_beta(pq) + RDMval
      end do

******* Copy D_beta to D_alpha and clean D_beta again for further use:
      if (.not. switch) then
        D_alpha(:) = D_beta(:) + D_alpha(:)
        D_beta(:) = 0.0d0
      end if
      close(iunit)
*******************************************************************************************
*************************** Processing TwoRDM-BABA ****************************************
*******************************************************************************************
      if (switch) then
        iUnit = IsFreeUnit(11)
        Call Molcas_Open(iUnit,'TwoRDM_baba.1')
        rewind(iUnit)
        if (IPRLEV >= DEBUG) then
           write(6,*) ' ********************** BABA ****************** '
        end if
        do
          read(iUnit,"(4I6,G25.17)",iostat=iread) s,q,r,p,RDMval
          if(iread /= 0) exit
          pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
          PSMAT(pqrs) = PSMAT(pqrs) + fac * RDMval
          if(r > s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
          if(r < s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
          if(r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
          IF(IPRLEV >= DEBUG) THEN
            write(6,'(7I6,3G25.17)')
     &          p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
          END IF
******* Contribution to D_alpha (not final):
          if(p == q) D_beta(rs) = D_beta(rs) + RDMval
          if(r == s.and.p /= r) D_alpha(pq) = D_alpha(pq)+RDMval
        end do
        close(iunit)
      end if ! End statement for spin-resolved RDMs.
*******************************************************************************************
*************************** Processing TwoRDM-ABBA ****************************************
*******************************************************************************************
      iUnit=IsFreeUnit(11)
      Call Molcas_Open(iUnit,'TwoRDM_abba.1')
      Rewind(iUnit)
      IF(IPRLEV >= DEBUG) THEN
        write(6,*) ' ********************** ABBA ****************** '
      END IF
      do
        read(iUnit,"(4I6,G25.17)",iostat=iread) q,s,r,p,RDMval
        if(iread /= 0) exit
        pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
        PSMAT(pqrs) = PSMAT(pqrs) - fac*RDMval
        if(r < s) then
          PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
        end if
        if(r > s) then
          PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
        end if
        IF(IPRLEV >= DEBUG) THEN
          write(6,'(7I6,3G25.17)')
     &        p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
        END IF
******* Contribution to D_alpha (not final):
          if(r == s) D_alpha(pq)=D_alpha(pq)-RDMval
          if(p == q) D_beta(rs)=D_beta(rs)-RDMval
      end do
      if (.not. switch) then
        D_alpha(:) = D_beta(:) + D_alpha(:)
        D_beta(:) = 0.0d0
      end if
      close(iunit)
*******************************************************************************************
*************************** Processing TwoRDM-BAAB ****************************************
*******************************************************************************************
      if(switch) then
        iUnit=IsFreeUnit(11)
        Call Molcas_Open(iUnit,'TwoRDM_baab.1')
        Rewind(iUnit)
        IF(IPRLEV >= DEBUG) THEN
          write(6,*) ' ********************** BAAB ****************** '
        END IF
        do
          read(iUnit,"(4I6,G25.17)",iostat=iread) q,s,r,p,RDMval
          if(iread /= 0) exit
          pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
          PSMAT(pqrs) = PSMAT(pqrs) - fac*RDMval
          if(r < s) then
            PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
          end if
          if(r > s) then
            PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
          end if
          IF(IPRLEV >= DEBUG) THEN
            write(6,'(7I6,3G25.17)')
     &          p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
          END IF
******* Contribution to D_alpha (not final):
          if(p == q) D_alpha(rs)=D_alpha(rs)-RDMval
          if(r == s) D_beta(pq)=D_beta(pq)-RDMval
        end do
        close(iunit)
      end if ! End statement for spin-resolved RDMs.
*******************************************************************************************
***************************   Final Updates to RDMs  **************************************
*******************************************************************************************
      if (.not.switch) D_beta(:) = D_alpha(:)
      D_alpha(:) = fcnacte * D_alpha(:)
      D_beta(:) = fcnacte * D_beta(:)
      DSPN(:) = D_Beta(:) - D_alpha(:)
      DMAT(:) = D_Beta(:) + D_alpha(:)

******* Clean evil non-positive semi-definite matrices. DMAT is input and output.
      call cleanMat(DMAT)

      IF(IPRLEV >= DEBUG) THEN
       norb  = (int(sqrt(dble(1 + 8 * size(DMAT)))) - 1) / 2
       call triprt('D_alpha in neci2molcas',' ',D_alpha,norb)
       call triprt('D_beta  in neci2molcas',' ',D_beta ,norb)
       call triprt('DMAT in neci2molcas',' ',DMAT,norb)
       call triprt('DSPN in neci2molcas',' ',DSPN,norb)
      END IF
      Return

123   continue
          write(6,*) 'RDM files not found!'
          write(6,*) 'Probably file not generated by NECI?'
          Call Abend()
      end subroutine read_neci_RDM

      subroutine bcast_2RDM(InFile)
        use filesystem, only : symlink_, strerror_, get_errno_
        character(len=*), intent(in) :: InFile
        character(len=1024) :: master
        integer :: lmaster1, err

        call prgmtranslate_master(InFile, master, lmaster1)
        call symlink_(trim(master), trim(InFile), err)
        if (err == 0) write(6, *) strerror_(get_errno_())
      end subroutine bcast_2RDM


      subroutine cleanup()
        ! Add your deallocations here.
        ! This routine will be called when exiting rasscf.
        continue
      end subroutine
      end module fciqmc_read_RDM
