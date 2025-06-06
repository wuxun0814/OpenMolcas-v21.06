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
      Subroutine Square_A(Lu,nB,MaxMem_,Force_out_of_Core)
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
      Logical Force_out_of_Core

      Real*8, Allocatable :: Buf(:,:)
*
      If (nB.eq.0) Return
      MaxMem=MaxMem_
      nMem=nB**2
      If (Force_Out_of_Core) MaxMem=nMem/3
*
      If (nMem.le.MaxMem) Then
*
*        In-core case
*
         Call mma_allocate(Buf,nMem,1,Label='Buf')
         iAddr=0
         Call dDaFile(Lu,2,Buf(:,1),nMem,iAddr)
         Call In_place_Square(Buf(:,1),nB)
         iAddr=0
         Call dDaFile(Lu,1,Buf(:,1),nMem,iAddr)

      Else
*
*        Out-of-core case
*
         nBuff=MaxMem/2
         Call mma_allocate(Buf,nBuff,2,Label='Buf')
*
         Inc = nBuff/nB
         iAddr1=0
         Do iB = 1, nB, Inc
            mB=Min(Inc,nB-iB+1)
            iAddrs=iAddr1
            Call dDaFile(Lu,2,Buf(:,1),nB*mB,iAddr1)
*
            iAddr2=iAddr1
            Do jB = iB, nB, Inc
               kB=Min(Inc,nB-jB+1)
*
               If (jB.eq.iB) Then
                  Call In_place_Diag(Buf(:,1),nB,iB,iB+mB-1)
               Else
                  Call dDaFile(Lu,2,Buf(:,2),nB*kB,iAddr2)
                  Call Off_Diagonal(Buf(:,1),nB,iB,iB+mB-1,
     &                              Buf(:,2),   jB,jB+kB-1)
               End If
            End Do
*
            iAddr1=iAddrs
            Call dDaFile(Lu,1,Buf(:,1),nB*mB,iAddr1)
*
         End Do
      End If
      Call mma_deallocate(Buf)

      Return
      End
      Subroutine In_place_Square(Buff,nBuff)
      Implicit Real*8 (a-h,o-z)
      Real*8 Buff(nBuff,nBuff)
*
C     Call RecPrt('Buff',' ',Buff,nBuff,nBuff)
      Do j = 1, nBuff
         Do i = 1, j-1
            Buff(j,i)=Buff(i,j)
         End Do
      End Do
C     Call RecPrt('Buff',' ',Buff,nBuff,nBuff)
C     Write (6,'(10F10.3)') (Buff(i,i),i=1,nBuff)
*
      Return
      End
      Subroutine In_place_Diag(Buff,nBuff,iBs,iBe)
      Implicit Real*8 (a-h,o-z)
      Real*8 Buff(nBuff,iBs:iBe)
*
C     Call RecPrt('Buff',' ',Buff,nBuff,iBe-iBs+1)
      Do j = iBs, iBe
         Do i = iBs, j-1
            Buff(j,i)=Buff(i,j)
         End Do
      End Do
C     Call RecPrt('Buff',' ',Buff,nBuff,iBe-iBs+1)
*
      Return
      End
      Subroutine Off_Diagonal(B1,nB,iB1s,iB1e,B2,iB2s,iB2e)
      Implicit Real*8 (a-h,o-z)
      Real*8 B1(nB,iB1s:iB1e), B2(nB,iB2s:iB2e)
*
      Do j = iB2s, iB2e
         Do i = iB1s, iB1e
            B1(j,i)=B2(i,j)
         End Do
      End Do
*
      Return
      End
