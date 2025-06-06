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
      SubRoutine Cho_SubScr_Dia(ChoVec,nVec,iSym,iLoc,DSPNorm)
C
C     Purpose: compute diagonal from Cholesky vectors and find norm
C              for each shell pair. Needed for screening of vector
C              subtraction. Which norm is used is determined by
C              string DSPNorm:
C
C              DSPNorm = 'Max' : max. element
C              DSPNorm = 'Fro' : Frobenius norm
C
C              Any other norm is taken to be 'Max'.
C
      use ChoSwp, only: nnBstRSh, iiBstRSh
      use ChoSubScr, only: DSubScr, DSPNm
#include "implicit.fh"
      Dimension ChoVec(*)
      Character*(*) DSPNorm
#include "cholesky.fh"

      Character*14 SecNam
      Parameter (SecNam = 'Cho_SubScr_Dia')

      Character*3 myDSPNorm

#if defined (_DEBUGPRINT_)
      If (iLoc.lt.1 .or. iLoc.gt.3) Then
         Call Cho_Quit('iLoc error in '//SecNam,104)
      End If
      If (iSym.lt.1 .or. iSym.gt.nSym) Then
         Call Cho_Quit('iSym error in '//SecNam,104)
      End If
      If (.not.Cho_SScreen) Then
         Call Cho_Quit('Cho_SScreen is .False. in '//SecNam,104)
      End If
#endif

C     Initialize and check for early return.
C     --------------------------------------

      Call Cho_dZero(DSubScr,nnBstR(iSym,iLoc))
      Call Cho_dZero(DSPNm,nnShl)
      If (nVec.lt.1 .or. nnBstR(iSym,iLoc).lt.1) return

C     Compute diagonal.
C     -----------------

      Do iVec = 1,nVec
         kOff = nnBstR(iSym,iLoc)*(iVec-1)
         Do iAB = 1,nnBstR(iSym,iLoc)
            DSubScr(iAB) = DSubScr(iAB)
     &                             + ChoVec(kOff+iAB)*ChoVec(kOff+iAB)
         End Do
      End Do

C     Find diagonal norm in each shell pair.
C     --------------------------------------

      lstr = len(DSPNorm)
      If (lstr .lt. 3) Then
         myDSPNorm = 'MAX'
      Else
         myDSPNorm = DSPNorm(1:3)
         Call UpCase(myDSPNorm)
      End If

#if defined (_DEBUGPRINT_)
      If (lstr .lt. 1) Then
         Write(Lupri,*) SecNam,': input norm: (null string)'
      Else If (lstr .lt. 3) Then
         Write(Lupri,*) SecNam,': input norm: ',DSPNorm,' (incomplete)'
      Else
         Write(Lupri,*) SecNam,': input norm: ',DSPNorm
      End If
      Write(Lupri,*) SecNam,': norm used : ',myDSPNorm
#endif

      If (myDSPNorm .eq. 'MAX') Then
         Do iSP = 1,nnShl
            iAB1 = iiBstRSh(iSym,iSP,iLoc) + 1
            iAB2 = iAB1 + nnBstRSh(iSym,iSP,iLoc) - 1
            Do iAB = iAB1,iAB2
               DSPNm(iSP) = max(DSPNm(iSP),DSubScr(iAB))
            End Do
         End Do
      Else If (myDSPNorm .eq. 'FRO') Then
         Do iSP = 1,nnShl
            iAB1 = iiBstRSh(iSym,iSP,iLoc) + 1
            iAB2 = iAB1 + nnBstRSh(iSym,iSP,iLoc) - 1
            Do iAB = iAB1,iAB2
               DSPNm(iSP) = DSPNm(iSP) + DSubScr(iAB)*DSubScr(iAB)
            End Do
            DSPNm(iSP) = sqrt(DSPNm(iSP))
         End Do
      Else
         Write(Lupri,*) SecNam,': WARNING: unkown norm: ',DSPNorm
         Write(Lupri,*) SecNam,': WARNING: using max element...'
         Do iSP = 1,nnShl
            iAB1 = iiBstRSh(iSym,iSP,iLoc) + 1
            iAB2 = iAB1 + nnBstRSh(iSym,iSP,iLoc) - 1
            Do iAB = iAB1,iAB2
               DSPNm(iSP) = max(DSPNm(iSP),DSubScr(iAB))
            End Do
         End Do
      End If

      End
