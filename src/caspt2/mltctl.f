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
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1998  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MLTCTL(HEFF,EIGVEC,U0)
      use output_caspt2, only:iPrGlb,terse,usual,verbose
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      INTEGER LAXITY
      INTEGER  Cho_X_GetTol
      EXTERNAL Cho_X_GetTol
      CHARACTER(LEN=8) INLAB
      character(len=3) variant
      DIMENSION HEFF(NSTATE,NSTATE),EIGVEC(NSTATE,NSTATE)
      real(8) U0(Nstate,Nstate)
      real(8),allocatable :: Utmp(:,:)


      IF(IPRGLB.GE.TERSE) THEN
        CALL CollapseOutput(1,'Multi-State CASPT2 section:')
        WRITE(6,'(20A4)')('****',I=1,20)
        WRITE(6,*)' MULTI-STATE CASPT2 SECTION'
        IF(IPRGLB.GE.USUAL) THEN
          WRITE(6,'(20A4)')('----',I=1,20)
          WRITE(6,*)
        END IF
      END IF

C Write out the effective Hamiltonian, for use in e.g. RASSI:
      INLAB='HEFF'
      CALL put_darray(INLAB,HEFF,NSTATE**2)

C Analyze the effective Hamiltonian:
      DSHIFT=0.0D0
      IF(HEFF(1,1).LE.-100.0D0) THEN
        DSHIFT=-DBLE(INT(-HEFF(1,1)))
      END IF
      DO I=1,NSTATE
        HEFF(I,I)=HEFF(I,I)-DSHIFT
      END DO

      IF(IPRGLB.GE.TERSE .and. DSHIFT.NE.0.0D0) THEN
        WRITE(6,*)
     &  ' Output diagonal energies have been shifted. Add ',DSHIFT
       END IF

      IF((IPRGLB.GE.VERBOSE).OR.JMS) THEN
        WRITE(6,*)' Effective Hamiltonian matrix (Asymmetric):'
        DO ISTA=1,NSTATE,5
          IEND=MIN(ISTA+4,NSTATE)
          WRITE(6,*)
          WRITE(6,'(1x,5I16)')(MSTATE(I),I=ISTA,IEND)
          DO J=1,NSTATE
            WRITE(6,'(1x,I3,3X,5F16.8)')
     &            MSTATE(J),(HEFF(J,I),I=ISTA,IEND)
          END DO
        END DO
      END IF

C Diagonalize:
C Use a symmetrized matrix, in triangular storage:
      NUMAT=NSTATE**2
      NHTRI=(NUMAT+NSTATE)/2
      CALL GETMEM('UMAT','ALLO','REAL',LUMAT,NUMAT)
      CALL GETMEM('HTRI','ALLO','REAL',LHTRI,NHTRI)
      IJ=0
      DO I=1,NSTATE
        DO J=1,I
          IJ=IJ+1
          WORK(LHTRI-1+IJ)=0.5D0*(HEFF(I,J)+HEFF(J,I))
        END DO
      END DO
      IF(IPRGLB.GE.USUAL) THEN
        WRITE(6,*)
        WRITE(6,*)' Effective Hamiltonian matrix (Symmetric):'
        DO ISTA=1,NSTATE,5
          IEND=MIN(ISTA+4,NSTATE)
          WRITE(6,*)
          WRITE(6,'(1x,5I16)')(MSTATE(I),I=ISTA,IEND)
          DO I=ISTA,NSTATE
            II0=(I*(I-1))/2
            WRITE(6,'(1x,I3,3X,5F16.8)')
     &            MSTATE(I),(WORK(LHTRI-1+II0+J),J=ISTA,MIN(I,IEND))
          END DO
        END DO
      END IF
      CALL DCOPY_(NSTATE**2,[0.0D0],0,WORK(LUMAT),1)
      CALL DCOPY_(NSTATE,[1.0D0],0,WORK(LUMAT),NSTATE+1)
      CALL NIDiag(WORK(LHTRI),WORK(LUMAT),NSTATE,NSTATE)
      CALL JACORD(WORK(LHTRI),WORK(LUMAT),NSTATE,NSTATE)
      DO I=1,NSTATE
        ENERGY(I)=DSHIFT+WORK(LHTRI-1+(I*(I+1))/2)
        DO J=1,NSTATE
          EIGVEC(J,I)=WORK(LUMAT-1+J+NSTATE*(I-1))
        END DO
      END DO
      CALL GETMEM('UMAT','FREE','REAL',LUMAT,NUMAT)
      CALL GETMEM('HTRI','FREE','REAL',LHTRI,NHTRI)

      IF(IPRGLB.GE.TERSE) THEN
        If (IFRMS) Then
          variant = 'RMS'
        Else if (IFXMS.and.IFDW) then
          variant = 'XDW'
        Else if (IFXMS) then
          variant = 'XMS'
        Else if (IFDW) then
          variant = 'DW '
        Else
          variant = 'MS '
        End If
          WRITE(6,*)
          WRITE(6,'(6X,A,A)')' Total ',trim(variant)//
     &      '-CASPT2 energies:'
          DO I=1,NSTATE
            Call PrintResult(6,'(6x,A,I3,5X,A,F16.8)',trim(variant)//
     &      '-CASPT2 Root',I,'Total energy:',ENERGY(I),1)
          END DO
      END IF

      IF(IPRGLB.GE.USUAL) THEN
        WRITE(6,*)
        WRITE(6,'(6X,A)')' Eigenvectors:'
        DO ISTA=1,NSTATE,5
          IEND=MIN(ISTA+4,NSTATE)
          DO J=1,NSTATE
            WRITE(6,'(6x,5F16.8)')(EIGVEC(J,I),I=ISTA,IEND)
          END DO
          WRITE(6,*)
        END DO
        if (IFXMS.or.IFRMS) then
* Transform eigenvectors into the original input basis
          call mma_allocate(Utmp,Nstate,Nstate,Label='Utmp')
          call dgemm_('N','N',Nstate,Nstate,Nstate,
     &                1.0d0,U0,Nstate,eigvec,Nstate,
     &                0.0d0,Utmp,Nstate)
          WRITE(6,'(6X,A)')' In terms of the input states:'
          DO ISTA=1,NSTATE,5
            IEND=MIN(ISTA+4,NSTATE)
            DO J=1,NSTATE
              WRITE(6,'(6x,5F16.8)')(Utmp(J,I),I=ISTA,IEND)
            END DO
            WRITE(6,*)
          END DO
          call mma_deallocate(Utmp)
        end if
        CALL CollapseOutput(0,'Multi-State CASPT2 section:')
        WRITE(6,*)
      END IF

* Restore original effective Hamiltonian
      DO I=1,NSTATE
        HEFF(I,I)=HEFF(I,I)+DSHIFT
      END DO

* In automatic verification calculations, the precision is lower
* in case of Cholesky calculation.
      LAXITY=8
      IF(IfChol) LAXITY=Cho_X_GetTol(LAXITY)
      Call Add_Info('E_MSPT2',ENERGY,nState,LAXITY)

      RETURN
      END
