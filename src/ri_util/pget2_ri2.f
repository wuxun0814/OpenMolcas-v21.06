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
* Copyright (C) 1992,2007, Roland Lindh                                *
************************************************************************
      SubRoutine PGet2_RI2(iCmp,iBas,jBas,kBas,lBas,
     &                  Shijij, iAO, iAOst, nijkl,PSO,nPSO,
     &                  ExFac,CoulFac,PMax,V_K,mV_K,Z_p_K,nSA,
     &                  nZ_p_k)
************************************************************************
*  Object: to assemble the 2nd order density matrix of a SCF wave      *
*          function from the 1st order density matrix.                 *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             January '92.                                             *
*                                                                      *
*             Modified to RI-DFT, March 2007                           *
************************************************************************
      use Basis_Info, only: nBas, nBas_Aux
      use SOAO_Info, only: iAOtSO
      use pso_stuff, only: nnp, lPSO, lsa, DMdiag
      use Symmetry_Info, only: nIrrep
      use ExTerm, only: CijK, iMP2prpt, A
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "exterm.fh"
      Real*8 PSO(nijkl,nPSO), V_K(mV_K,nSA),Z_p_K(nZ_p_k,*)
      Integer iCmp(4), iAO(4), iAOst(4)
      Logical Shijij, Found
*     Local Array
      Integer jSym(0:7), lSym(0:7)
      Integer CumnnP(0:7),CumnnP2(0:7)

      Real*8, Pointer :: CiKj(:)=>Null(), CiKl(:)=>Null(),
     &                   V2(:)=>Null()
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt('V_K',' ',V_K,1,mV_K)
#endif

      Call CWTime(Cpu1,Wall1)
*                                                                      *
************************************************************************
*                                                                      *
      PMax=Zero
      iSO=1
*
      Call FZero(PSO,nijkl*nPSO)
*
      If (lPSO) Then
        CumnnP(0)=0
        CumnnP2(0)=0
        Do i=1,nIrrep-1
          nB = nBas_Aux(i-1)
          If (i.eq.1) nB = nB-1
          CumnnP(i)=CumnnP(i-1)+nnP(i-1)
          CumnnP2(i)=CumnnP2(i-1)+nnP(i-1)*nB
        End Do
      End If
*
      Call Qpg_iScalar('SCF mode',Found)
      If (Found) Then
         Call Get_iScalar('SCF mode',iUHF) ! either 0 or 1
      Else
         iUHF=0
      EndIf

*                                                                      *
************************************************************************
*                                                                      *
      Fac = One/Four
      MemSO2 = 0
*                                                                      *
************************************************************************
*                                                                      *
*     Pure DFT
*
      If (ExFac.eq.Zero) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do i2 = 1, iCmp(2)
         njSym = 0
         Do j = 0, nIrrep-1
            If (iAOtSO(iAO(2)+i2,j)>0) Then
               jSym(njSym) = j
               njSym = njSym + 1
            End If
         End Do
*
         Do i4 = 1, iCmp(4)
            nlSym = 0
            Do j = 0, nIrrep-1
               If (iAOtSO(iAO(4)+i4,j)>0) Then
                  lSym(nlSym) = j
                  nlSym = nlSym + 1
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over irreps which are spanned by the basis function.
*
            Do js = 0, njSym-1
               j2 = jSym(js)
               jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
*
               Do ls = 0, nlSym-1
                  j4 = lSym(ls)
                  If (j2/=j4) Cycle
                  lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                  MemSO2 = MemSO2 + 1
                  If (j2/=0) Cycle
*
                  mijkl = 0
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl - nBas(j4)
                     Do jAOj = 0, jBas-1
                        jSOj = jSO + jAOj - nBas(j2)
                        mijkl = mijkl + 1
*
*-----------------------Coulomb contribution
                        If (j2.eq.0) Then
*---------------------------j4.eq.0 also
                           temp=V_K(jSOj,1)*V_K(lSOl,1)*Coulfac
*                          temp=Zero
                        Else
                           temp = Zero
                        End If
*
                        PMax=Max(PMax,Abs(Temp))
                        PSO(mijkl,MemSO2) =  Fac * temp
*
                     End Do
                  End Do
*
               End Do
            End Do
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Hybrid DFT and HF
*
      Else If (iMP2prpt .ne. 2 .and. .Not.lPSO .and. iUHF.eq.0 ) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do i2 = 1, iCmp(2)
         njSym = 0
         Do j = 0, nIrrep-1
            If (iAOtSO(iAO(2)+i2,j)>0) Then
               jSym(njSym) = j
               njSym = njSym + 1
            End If
         End Do
*
         Do i4 = 1, iCmp(4)
            nlSym = 0
            Do j = 0, nIrrep-1
               If (iAOtSO(iAO(4)+i4,j)>0) Then
                  lSym(nlSym) = j
                  nlSym = nlSym + 1
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over irreps which are spanned by the basis function.
*
            Do js = 0, njSym-1
               j2 = jSym(js)
               jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
*
               Do ls = 0, nlSym-1
                  j4 = lSym(ls)
                  If (j2/=j4) Cycle
                  lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                  MemSO2 = MemSO2 + 1
*
                  A(1:jBas*lBas)=Zero
*
                  Do iSym = 1, nIrrep
                     kSym = iEor(j2,iSym-1)+1
                     nik = nIJ1(iSym,kSym,iSO)
*
                     If (nik==0) Cycle
*
                     iS = 1
                     iE = nik*jBas
                     CiKj(1:nik*jBas) => CijK(iS:iE)


                     jSOj= jSO-nBas(j2)
                     iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
                     Call dDaFile(LuCVector(j2+1,iSO),2,CikJ,nik*jBas,
     &                            iAdrJ)
*
                     If (lSO.ne.jSO) Then
                        iS = iE + 1
                        iE = iE + nik*lBas
                        CiKl(1:nik*lBas) => CijK(iS:iE)

                        lSOl=lSO-nBas(j4)
                        iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
                        Call dDaFile(LuCVector(j4+1,iSO),2,CiKl,
     &                               nik*lBas,iAdrL)
                        V2(1:) => CiKl(1:)
                     Else
                        V2(1:) => CiKj(1:)
                     End If
*
                     Fact=One
                     If (iSym.ne.kSym) Fact=Half
                     Call DGEMM_('T','N',jBas,lBas,nik,
     &                           Fact,CikJ,nik,
     &                                V2,nik,
     &                          1.0D0,A,jBas)
*
                  End Do

                  mijkl = 0
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl - nBas(j4)
                     Do jAOj = 0, jBas-1
                        jSOj = jSO + jAOj - nBas(j2)
                        mijkl = mijkl + 1
*
*-----------------------Coulomb contribution
                        If (j2.eq.0) Then
*---------------------------j4.eq.0 also
                           temp=V_K(jSOj,1)*V_K(lSOl,1)*Coulfac
*                          temp=Zero
                        Else
                           temp = Zero
                        End If
*
*-----------------------Exchange contribution
                        temp = temp - ExFac*A(mijkl)
*
                        PMax=Max(PMax,Abs(Temp))
                        PSO(mijkl,MemSO2) =  Fac * temp
*
                     End Do
                  End Do
*
               End Do
            End Do
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Hybrid UDFT and UHF
*
      Else If (iMP2prpt .ne. 2 .and. .Not.lPSO .and. iUHF.eq.1 ) Then
*
         Write (6,*) 'Pget2_RI2: UDFT/UHF not implemented yet.'
         Call Abend()
*                                                                      *
************************************************************************
*                                                                      *
*     CASSCF
*
      Else If (iMP2prpt .ne. 2 .and. lPSO .and. .Not. LSA) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do i2 = 1, iCmp(2)
         njSym = 0
         Do j = 0, nIrrep-1
            If (iAOtSO(iAO(2)+i2,j)>0) Then
               jSym(njSym) = j
               njSym = njSym + 1
            End If
         End Do
*
         Do i4 = 1, iCmp(4)
            nlSym = 0
            Do j = 0, nIrrep-1
               If (iAOtSO(iAO(4)+i4,j)>0) Then
                  lSym(nlSym) = j
                  nlSym = nlSym + 1
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over irreps which are spanned by the basis function.
*
            Do js = 0, njSym-1
               j2 = jSym(js)
               jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
*
               Do ls = 0, nlSym-1
                  j4 = lSym(ls)
                  If (j2/=j4) Cycle
                  lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                  MemSO2 = MemSO2 + 1
*
                  A(1:jBas*lBas)=Zero
*
                  Do iSym = 1, nIrrep
                     kSym = iEor(j2,iSym-1)+1
                     nik = nIJ1(iSym,kSym,iSO)
*
                     If (nik==0) Cycle

                     iS = 1
                     iE = nik*jBas
                     CiKj(1:nik*jBas) => CijK(iS:iE)
*
                     jSOj= jSO-nBas(j2)
                     iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
                     Call dDaFile(LuCVector(j2+1,iSO),2,CiKj,nik*jBas,
     &                            iAdrJ)
*
                     If (lSO.ne.jSO) Then
                        iS = iE + 1
                        iE = iE + nik*lBas
                        CiKl(1:nik*lBas) => CijK(iS:iE)

                        lSOl=lSO-nBas(j4)
                        iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
                        Call dDaFile(LuCVector(j4+1,iSO),2,CiKl,
     &                               nik*lBas,iAdrL)
                        V2(1:) => CiKl(1:)
                     Else
                        V2(1:) => CiKj(1:)
                     End If
*
                     Fact=One
                     If (iSym.ne.kSym) Fact=Half
                     Call DGEMM_('T','N',jBas,lBas,nik,
     &                           Fact,CikJ,nik,
     &                                V2,nik,
     &                          1.0D0,A,jBas)
*
                  End Do
*
                  mijkl = 0
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl - nBas(j4)
                     Do jAOj = 0, jBas-1
                        jSOj = jSO + jAOj - nBas(j2)
                        mijkl = mijkl + 1
*
*-----------------------Coulomb contribution
                        If (j2.eq.0) Then
*---------------------------j4.eq.0 also
                           temp=V_K(jSOj,1)*V_K(lSOl,1)*Coulfac
*                          temp=Zero
                        Else
                           temp = Zero
                        End If
*
*-----------------------Exchange contribution
                        temp = temp - ExFac*A(mijkl)
*
                        temp2=0.0d0
                        jpSOj=CumnnP2(j2)+(jSOj-1)*nnP(j2)
                        jpSOl=CumnnP2(j2)+(lSOl-1)*nnP(j2)
                        Do jp=1,nnP(j2)
                          temp2=temp2+sign(1.0d0,
     &                          DMdiag(CumnnP(j2)+jp,1))*
     &                          Z_p_K(jpSOj+jp,1)*Z_p_K(jpSOl+jp,1)
                        End Do
                        temp=temp+temp2
*
                        PMax=Max(PMax,Abs(Temp))
                        PSO(mijkl,MemSO2) =  Fac * temp
*
                     End Do
                  End Do
*
               End Do
            End Do
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     SA-CASSCF
*
      Else If ( iMP2prpt .ne. 2 .and. lPSO .and. lSA ) Then
*                                                                      *
************************************************************************
*                                                                      *
      Write (6,*) 'Pget2_ri2: SA-CASSCF not implemented yet'
      Call Abend()
*
      Do i2 = 1, iCmp(2)
         njSym = 0
         Do j = 0, nIrrep-1
            If (iAOtSO(iAO(2)+i2,j)>0) Then
               jSym(njSym) = j
               njSym = njSym + 1
            End If
         End Do
*
         Do i4 = 1, iCmp(4)
            nlSym = 0
            Do j = 0, nIrrep-1
               If (iAOtSO(iAO(4)+i4,j)>0) Then
                  lSym(nlSym) = j
                  nlSym = nlSym + 1
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over irreps which are spanned by the basis function.
*
            Do js = 0, njSym-1
               j2 = jSym(js)
               jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
*
               Do ls = 0, nlSym-1
                  j4 = lSym(ls)
                  If (j2/=j4) Cycle
                  lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                  MemSO2 = MemSO2 + 1
*
                  A(1:jBas*lBas) = Zero
*
                  Do iSym = 1, nIrrep
                     kSym = iEor(j2,iSym-1)+1
                     nik = nIJ1(iSym,kSym,iSO)
*
                     If (nik==0) Cycle

                     iS = 1
                     iE = nik*jBas
                     CiKj(1:nik*jBas) => CijK(iS:iE)
*
                     jSOj= jSO-nBas(j2)
                     iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
                     Call dDaFile(LuCVector(j2+1,iSO),2,CiKj,nik*jBas,
     &                            iAdrJ)
*
                     If (lSO.ne.jSO) Then
                        iS = iE + 1
                        iE = iE + nik*lBas
                        CiKl(1:nik*lBas) => CijK(iS:iE)

                        lSOl=lSO-nBas(j4)
                        iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
                        Call dDaFile(LuCVector(j4+1,iSO),2,CiKl,
     &                               nik*lBas,iAdrL)
                        V2(1:) => CiKl(1:)
                     Else
                        V2(1:) => CiKj(1:)
                     End If
*
                     Fact=One
                     If (iSym.ne.kSym) Fact=Half
                     Call DGEMM_('T','N',jBas,lBas,nik,
     &                           Fact,CiKJ,nik,
     &                                V2,nik,
     &                          1.0D0,A,jBas)
*
                  End Do

                  mijkl = 0
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl - nBas(j4)
                     Do jAOj = 0, jBas-1
                        jSOj = jSO + jAOj - nBas(j2)
                        mijkl = mijkl + 1
*
*-----------------------Coulomb contribution
                        If (j2.eq.0) Then
*---------------------------j4.eq.0 also
                           temp=CoulFac*(V_K(lSOl,1)*V_K(jSOj,2)+
     &                                   V_K(lSOl,2)*V_K(jSOj,1)+
     &                                   V_K(lSOl,3)*V_K(jSOj,4)+
     &                                   V_K(lSOl,4)*V_K(jSOj,3))
*                          temp=Zero
                        Else
                           temp = Zero
                        End If
*
*-----------------------Exchange contribution
                        temp = temp - ExFac*A(mijkl)
*
                        PMax=Max(PMax,Abs(Temp))
                        PSO(mijkl,MemSO2) =  Fac * temp
*
                     End Do
                  End Do
*
               End Do
            End Do
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     MP2
*
      Else
*                                                                      *
************************************************************************
*                                                                      *
      Write (6,*) 'Pget2_ri2: MP2 not implemented yet'
      Call Abend()
*
      Do i2 = 1, iCmp(2)
         njSym = 0
         Do j = 0, nIrrep-1
            If (iAOtSO(iAO(2)+i2,j)>0) Then
               jSym(njSym) = j
               njSym = njSym + 1
            End If
         End Do
*
         Do i4 = 1, iCmp(4)
            nlSym = 0
            Do j = 0, nIrrep-1
               If (iAOtSO(iAO(4)+i4,j)>0) Then
                  lSym(nlSym) = j
                  nlSym = nlSym + 1
               End If
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over irreps which are spanned by the basis function.
*
            Do js = 0, njSym-1
               j2 = jSym(js)
               jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
*
               Do ls = 0, nlSym-1
                  j4 = lSym(ls)
                  If (j2/=j4) Cycle
                  lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                  MemSO2 = MemSO2 + 1
*
                  A(1:jBas*lBas) = Zero
*
                  Do iSym = 1, nIrrep
                     kSym = iEor(j2,iSym-1)+1
                     nik = nIJ1(iSym,kSym,iSO)
*
                     If (nik==0) Cycle

                     iS = 1
                     iE = nik*jBas
                     CiKj(1:nik*jBas) => CijK(iS:iE)
*
                     jSOj= jSO-nBas(j2)
                     iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
                     Call dDaFile(LuCVector(j2+1,iSO),2,CiKj,nik*jBas,
     &                            iAdrJ)
*
                     If (lSO.ne.jSO) Then
                        iS = iE + 1
                        iE = iE + nik*lBas
                        CiKl(1:nik*lBas) => CijK(iS:iE)

                        lSOl=lSO-nBas(j4)
                        iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
                        Call dDaFile(LuCVector(j4+1,iSO),2,CiKl,
     &                               nik*lBas,iAdrL)
                        V2(1:) => CiKl(1:)
                     Else
                        V2(1:) => CiKj(1:)
                     End If
*
                     Fact=One
                     If (iSym.ne.kSym) Fact=Half
                     Call DGEMM_('T','N',jBas,lBas,nik,
     &                           Fact,CiKj,nik,
     &                                V2,nik,
     &                          1.0D0,A,jBas)
*
                  End Do
*
                  mijkl = 0
                  Do lAOl = 0, lBas-1
                     lSOl = lSO + lAOl - nBas(j4)
                     Do jAOj = 0, jBas-1
                        jSOj = jSO + jAOj - nBas(j2)
                        mijkl = mijkl + 1
*
*-----------------------Coulomb contribution
                        If (j2.eq.0) Then
*---------------------------j4.eq.0 also
                           temp=V_K(jSOj,1)*V_K(lSOl,1)*Coulfac
*                          temp=Zero
                        Else
                           temp = Zero
                        End If
*
*-----------------------Exchange contribution
                        temp = temp - ExFac*A(mijkl)
*
                        PMax=Max(PMax,Abs(Temp))
                        PSO(mijkl,MemSO2) =  Fac * temp
*
                     End Do
                  End Do
*
               End Do
            End Do
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If
      CiKj => Null()
      CiKl => Null()
      V2   => Null()
*                                                                      *
************************************************************************
*                                                                      *
      If (nPSO.ne.MemSO2) Then
        Write (6,*) ' PGet2: nPSO.ne.MemSO2'
        Write (6,*) nPSO, MemSO2
        Call Abend()
      End If
*
#ifdef _DEBUGPRINT_
      Call RecPrt(' In PGet2_RI2:PSO ',' ',PSO,nijkl,nPSO)
#endif
*
      Call CWTime(Cpu2,Wall2)
      Cpu = Cpu2 - Cpu1
      Wall = Wall2 - Wall1
      tavec(1) = tavec(1) + Cpu
      tavec(2) = tavec(2) + Wall
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iBas)
         Call Unused_integer(kBas)
         Call Unused_logical(Shijij)
      End If
      End
