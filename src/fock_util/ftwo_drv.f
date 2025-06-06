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
      SUBROUTINE FTwo_Drv(nSym,nBas,nAsh,nSkipX,DI,D1A,FA,nTot1,
     &                    ExFac,nBMX,CMO)

      use Data_Structures, only: DSBA_Type, Allocate_DSBA,
     &                           Deallocate_DSBA
      Implicit Real*8 (A-H,O-Z)
#include "rasdim.fh"
#include "stdalloc.fh"
#include "real.fh"

      Integer nBas(8), nAsh(8), nSkipX(8)
      Real*8 CMO(*) , D1A(*) , DI(*), FA(*)
      Logical DoCholesky

      Type (DSBA_Type) :: WFSQ

      Real*8, Allocatable:: Temp(:)
#include "choras.fh"
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      SUBROUTINE CHORAS_DRV(nSym,nBas,nOcc,W_DSQ,W_DLT,W_FLT,ExFac,FSQ,
     &                      W_CMO)
      use Data_Structures, only: DSBA_Type
      Integer nSym, nBas(8)
      Integer, Target :: nOcc(nSym)
      Real*8 W_FLT(*), W_DSQ(*),W_DLT(*)
      Real*8 ExFac
      Type (DSBA_Type) FSQ
      Real*8 W_CMO(*)
      END SUBROUTINE CHORAS_DRV
      END Interface
*                                                                      *
************************************************************************
*

      Call DecideOncholesky(DoCholesky)

      IF (DoCholesky.and.ALGO.eq.2)THEN
*
* Building of the Fock matrix directly from Cholesky
* vectors
*
         Call Allocate_DSBA(WFSQ,nBas,nBas,nSym)
         WFSQ%A0(:)=Zero

         Call mma_allocate(Temp,nTot1,Label='nTot1')
         Temp(:)=Zero
*
         CALL CHORAS_DRV(nSym,nBas,nAsh,D1A,DI,Temp,ExFac,WFSQ,CMO)

         FA(1:nTot1) = FA(1:nTot1) + Temp(1:nTot1)
*
         Call mma_deallocate(Temp)
         Call Deallocate_DSBA(WFSQ)

      ELSE

*
* Standard building of the Fock matrix from Two-el integrals
*
         Call FockTwo_Drv(nSym,nBas,nAsh,nSkipX,DI,D1A,FA,nTot1,
     &                    ExFac,nBMX)

      ENDIF


      RETURN
      END
