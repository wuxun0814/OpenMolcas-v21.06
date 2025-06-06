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
* Copyright (C) 1992, Markus P. Fuelscher                              *
************************************************************************
      Subroutine CRecPrt(Title,FmtIn,A,nRow,nCol,Type)
************************************************************************
* CRecPrt
*
*> @brief
*>   Write out a matrix on standard output
*> @author M. P. F&uuml;lscher, Lund, 1992
*>
*> @details
*> The matrix \p A of dimension \p nRow &times; \p nCol is printed
*> in output preceded by the character line \p Title. Format of the numerical
*> output is given by \p FmtIn. If \p FmtIn = ``''`` the utility will decide on format
*> for optimal output.
*>
*> @param[in] Title   Title card
*> @param[in] FmtIn   Format statement
*> @param[in] A       A matrix
*> @param[in] nRow    number of rows of \p A
*> @param[in] nCol    number of columns of \p A
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "standard_iounits.fh"
      Character*(*) Title
      Character*(*) FmtIn
      Character*1 Type
      Complex*16 A(nRow,nCol)
      Integer StrnLn
      Parameter (lPaper=120,lMaxTitle=60)
      Character*(lMaxTitle) Line
      Character*20 FMT
*----------------------------------------------------------------------*
      If (nRow*nCol.eq.0) Return
#ifdef _DEBUGPRINT_
C     Call CTrcPrt(Title,FmtIn,A,nRow,nCol,'R')
C     Call CTrcPrt(Title,FmtIn,A,nRow,nCol,'I')
      Return
#endif
*----------------------------------------------------------------------*
*     print the title                                                  *
*----------------------------------------------------------------------*
      lTitle=StrnLn(Title)
      If ( lTitle.gt.0 ) then
         Do 10 i=1,lMaxTitle
             Line(i:i)=' '
10       Continue
         lLeft=1
         Do 20 i=lTitle,1,-1
            If ( Title(i:i).ne.' ' ) lLeft=i
20       Continue
         lLeft=lLeft-1
         Do 25 i=1,lMaxTitle
            If ( i+lLeft.le.lTitle ) Line(i:i)=Title(i+lLeft:i+lLeft)
25       Continue
         Write(LuWr,*)
         If (Type.eq.'R') Then
            Write(LuWr,'(2X,A,A)') Line,' Real Component'
         Else
            Write(LuWr,'(2X,A,A)') Line,' Imaginary Component'
         End If
c         Do 30 i=1,StrnLn(Line)
c            Line(i:i)='-'
c30       Continue
c         Write(LuWr,'(2X,A)') Line
         Write(LuWr,'(2X,A,I5,A,I5)') 'mat. size = ',nRow,'x',nCol
      End If
*----------------------------------------------------------------------*
*     determine the printing format                                    *
*----------------------------------------------------------------------*
      lFmt=Strnln(FmtIn)
      If ( lFmt.ne.0 ) then
         FMT=FmtIn
      Else
         If (Type.eq.'R') Then
         Amax=DBLE(A(1,1))
         Amin=DBLE(A(1,1))
         Do j=1,nCol
            Do i=1,nRow
               Amax=Max(Amax,DBLE(A(i,j)))
               Amin=Min(Amin,DBLE(A(i,j)))
            End Do
         End Do
         Else
         Amax=DIMAG(A(1,1))
         Amin=DIMAG(A(1,1))
         Do j=1,nCol
            Do i=1,nRow
               Amax=Max(Amax,DIMAG(A(i,j)))
               Amin=Min(Amin,DIMAG(A(i,j)))
            End Do
         End Do
         End If
         Pmax=0.0D0
         If ( Abs(Amax).gt.1.0D-72 ) Pmax=Log10(Abs(Amax))
         iPmax=1+INT(Pmax)
         iPmax=Max(1,iPmax)
         Pmin=0.0D0
         If ( Abs(Amin).gt.1.0D-72 ) Pmin=Log10(Abs(Amin))
         iPmin=1+INT(Pmin)
         iPmin=Max(1,iPmin)
         nDigit=15
         nDecim=Min(9,nDigit-Max(iPmin,iPmax))
         nDecim=Max(nDecim,1)
         If ( Amax.lt.0.0D0 ) iPmax=iPmax+1
         If ( Amin.lt.0.0D0 ) iPmin=iPmin+1
         lNumbr=Max(iPmin,iPmax)+nDecim+2
         nCols=9
         lLine=nCols*lNumbr
         If ( lLine.gt.lPaper ) then
            If ( lLine.le.lPaper+nCols .and. nDecim.gt.1 ) then
               nDecim=nDecim-1
               lNumbr=Max(iPmin,iPmax)+nDecim
               lItem=Max(lNumbr,lPaper/nCols)
            Else
               nCols=5
               lItem=Max(lNumbr,lPaper/nCols)
            End If
         Else
            lItem=lNumbr
         End If
         Write(FMT,'(A,   I4.4,  A, I4.4,  A, I4.4,   A)')
     &             '(2X,',nCols,'F',lItem,'.',nDecim,')'
      End if
*----------------------------------------------------------------------*
*     print the data                                                   *
*----------------------------------------------------------------------*
c      Write(LuWr,*)
      If (Type.eq.'R') Then
      Do i=1,nRow
         Write(LuWr,FMT)(DBLE(A(i,j)),j=1,nCol)
      End Do
      Else
      Do i=1,nRow
         Write(LuWr,FMT)(DIMAG(A(i,j)),j=1,nCol)
      End Do
      End If
*----------------------------------------------------------------------*
*     End procedure                                                    *
*----------------------------------------------------------------------*
      Return
      End
