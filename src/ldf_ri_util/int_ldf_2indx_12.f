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
      SubRoutine Int_LDF_2Indx_12(
     &                            iCmp,iShell,MapOrg,
     &                            iBas,jBas,kBas,lBas,kOp,
     &                            Shijij,IJeqKL,iAO,iAOst,ijkl,
     &                            AOInt,SOInt,nSOint,
     &                            iSOSym,nSkal,nSOs,
     &                            TInt,nTInt,itOffs,nSym)
*     calls the proper routines IndSft/PLF
*     if IntOrd_jikl==.TRUE. integral order within symblk: jikl
*                      else  integral order within symblk: ijkl
      Implicit Real*8 (a-h,o-z)
*
#include "localdf_int.fh"
*
      Character*16 SecNam
      Parameter (SecNam='Int_LDF_2Indx_12')
*
      Real*8 AOInt(*), SOInt(*), TInt(nTInt)
      Integer iCmp(4), iShell(4), iAO(4),
     &        iAOst(4), kOp(4), iSOSym(2,nSOs),
     &        itOffs(0:nSym-1,0:nSym-1,0:nSym-1), MapOrg(4)
      Logical Shijij,IJeqKL
*
      External LDF_nShell, LDF_nAuxShell
*
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
*
* call sorting routine
*
      If (nSym==1) Then
         nS_Val=LDF_nShell()
         nS_Aux=LDF_nAuxShell()
         iS_Dum=nS_Val+nS_Aux+1
         If (SHA.eq.iS_Dum .and.
     &       SHB.gt.nS_Val .and. SHB.lt.iS_Dum .and.
     &       SHC.le.nS_Val .and.
     &       SHD.le.nS_Val) Then
            !check that shells have not been re-ordered
            If (MapOrg(1).ne.1 .or. MapOrg(2).ne.2 .or.
     &          MapOrg(3).ne.3 .or. MapOrg(4).ne.4) Then
               Call WarningMessage(2,
     &      SecNam//': Shell reordering not implemented for this case!')
               Write(6,'(A,4I9)')
     &         'MapOrg.................',(MapOrg(i),i=1,4)
               Write(6,'(A,4I9)')
     &         'SHA,SHB,SHC,SHD........',SHA,SHB,SHC,SHD
               Write(6,'(A,2(9X,I9))')
     &         'SHAB,SHCD..............',iTri(SHA,SHB),iTri(SHC,SHD)
               Write(6,'(A,3I9)')
     &         'nS_Val,nS_Aux,iS_Dum...',nS_Val,nS_Aux,iS_Dum
               Call LDF_Quit(1)
            End If
            ! type (J|kl) where kl are two-center aux functions
            Call PLF_LDF_2Indx_12(TInt,nTInt,
     &                            AOInt,ijkl,
     &                            iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                            iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
         Else
            Call WarningMessage(2,
     &                  'Shell combination not implemented in '//SecNam)
            Write(6,'(A,4I9)')
     &      'SHA,SHB,SHC,SHD........',SHA,SHB,SHC,SHD
            Write(6,'(A,3I9)')
     &      'nS_Val,nS_Aux,iS_Dum...',nS_Val,nS_Aux,iS_Dum
            Call LDF_Quit(1)
         End If
      Else
         Call WarningMessage(2,'Symmetry not implemented in '//SecNam)
         Call LDF_Quit(1)
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iShell)
         Call Unused_logical(Shijij)
         Call Unused_logical(IJeqKL)
         Call Unused_real_array(SOInt)
         Call Unused_integer(nSOint)
         Call Unused_integer_array(iSOSym)
         Call Unused_integer(nSkal)
         Call Unused_integer_array(itOffs)
      End If
      End
