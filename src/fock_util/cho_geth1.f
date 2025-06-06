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
      SUBROUTINE CHO_GETH1(nBtri,H1,RFpert,ERFNuc)

      IMPLICIT REAL*8 (A-H,O-Z)

#include "real.fh"
#include "stdalloc.fh"
      Integer nBTri
      Real*8 H1(nBTri)
      Logical RFpert
      Real*8 ERFNuc

      Character*8 OneLbl
      Real*8, Allocatable:: Tmp(:)

      iRc=-1
      iOpt=6
      iCmp=1
      iSyLab=1
      OneLbl='OneHam  '

      Call RdOne(iRc,iOpt,OneLbl,iCmp,H1,iSyLab)

      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'    *** ERROR IN SUBROUTINE  CHO_GETH1 *** '
        WRITE(6,*)'   BARE NUCLEI HAMILTONIAN IS NOT AVAILABLE'
        WRITE(6,*)
        CALL Abend
      ENDIF

      ERFNuc=Zero
      If ( RFpert ) then
         Call mma_allocate(Tmp,nBTri,Label='Tmp')
         Call Get_dScalar('RF Self Energy',ERFNuc)
         Call Get_dArray('Reaction field',Tmp,nBtri)

         H1(:) = H1(:) + Tmp(:)

         Call mma_deallocate(Tmp)
      End If


      RETURN
      END
