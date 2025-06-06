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
      Subroutine Get_Fock_Occ(FockOcc,nFockOcc)
      Implicit Real*8 (A-H,O-Z)
      Real*8 FockOcc(nFockOcc)
#include "WrkSpc.fh"
#include "SysDef.fh"

      Character(LEN=24) Label
#ifdef _DEBUGPRINT_
#include "run_common.fh"
#endif
      Logical      Found

*
*
*...  Read the generalized Fock matrix
*                                                                      *
************************************************************************
*                                                                      *
      Label='FockOcc'
      Call qpg_dArray(Label,Found,mFockOcc)
      If(.not.Found .or. mFockOcc.eq.0) Then
         Call SysAbendMsg('get_fock_occ','Did not find:',Label)
      End If
      If(mFockOcc/=nFockOcc) Then
         Write (6,*) 'nFockOcc=',nFockOcc
         Write (6,*) 'mFockOcc=',mFockOcc
         Call SysAbendMsg('get_fock_occ','mFockOcc/=nFockOcc:',Label)
      End If
      Call Get_dArray(Label,FockOcc,nFockOcc)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      if(is_nSym.eq.0) then
       Call get_iScalar('nSym',nSym)
       is_nSym=1
      endif
      if(is_nBas.eq.0) then
       Call Get_iArray('nBas',nBas,nSym)
       is_nBas=1
      endif
      Write(6,*) 'Fock occ'
      ii=1
      Do iIrrep = 0, nSym - 1
         If (nBas(iIrrep).gt.0) Then
            Write(6,*) 'symmetry block',iIrrep
            Call TriPrt(' ',' ',FockOcc(ii),nBas(iIrrep))
            ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End If
      End Do
#endif
*
      Return
      End
