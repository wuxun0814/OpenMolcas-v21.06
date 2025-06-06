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
      SubRoutine Integral_RI_3(iCmp,iShell,MapOrg,
     &                         iBas,jBas,kBas,lBas,kOp,
     &                         Shijij,IJeqKL,iAO,iAOst,ijkl,
     &                         AOInt,SOInt,nSOint,
     &                         iSOSym,nSkal,nSOs,
     &                         TInt,nTInt,itOffs,nSym)
*     calls the proper routines IndSft/PLF
*     if IntOrd_jikl==.TRUE. integral order within symblk: jikl
*                      else  integral order within symblk: ijkl
      Use RICD_Info, only: LDF
      use j12
      Implicit Real*8 (A-H,O-Z)
*
      Real*8 AOInt(*), SOInt(*), TInt(nTInt)
      Integer iCmp(4), iShell(4), iAO(4),
     &        iAOst(4), kOp(4), iSOSym(2,nSOs),
     &        itOffs(0:nSym-1,0:nSym-1,0:nSym-1), MapOrg(4)
      Logical Shijij,IJeqKL
*                                                                      *
************************************************************************
*                                                                      *
      If (LDF) Then
*                                                                      *
************************************************************************
*                                                                      *
         If (nSym==1) Then
           Call PLF_LDF_3(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                   iShell,iAO,iAOst,Shijij.and.IJeqKL,
     &                   iBas,jBas,kBas,lBas,kOp,
     &                   TInt,nTInt,iTOffs,
     &                   ShlSO,nBasSh,
     &                   SOShl,nSO,nSkal_Valence,nSym,
     &                   iSSOff(0,0,klS))
         Else
           Call WarningMessage(2,'Not implemented yet!')
           Call Abend()
C          Call IndSft_RI_3(iCmp,iShell,
C    &                      iBas,jBas,kBas,lBas,Shijij,
C    &                      iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs,
C    &                      TInt,nTInt,iTOffs,
C    &                      ShlSO,nBasSh,
C    &                      SOShl,nSO,nSkal_Valence,nSym,
C    &                      iSSOff(:,:,klS))
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
         If (nSym==1) Then
           Call PLF_RI_3(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                   iShell,iAO,iAOst,Shijij.and.IJeqKL,
     &                   iBas,jBas,kBas,lBas,kOp,
     &                   TInt,nTInt,iTOffs,
     &                   ShlSO,nBasSh,
     &                   SOShl,nSO,nSkal_Valence,nSym,
     &                   iSSOff(0,0,klS))
      Else
           Call IndSft_RI_3(iCmp,iShell,
     &                      iBas,jBas,kBas,lBas,Shijij,
     &                      iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs,
     &                      TInt,nTInt,iTOffs,
     &                      ShlSO,nBasSh,
     &                      SOShl,nSO,nSkal_Valence,nSym,
     &                      iSSOff(:,:,klS))
         End If
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(MapOrg)
         Call Unused_integer(nSkal)
      End If
      End
