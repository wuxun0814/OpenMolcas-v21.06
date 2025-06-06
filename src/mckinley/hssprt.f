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
      SubRoutine HssPrt(Hess,nHess)
      use Symmetry_Info, only: nIrrep, lIrrep
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "disp2.fh"
#include "real.fh"
      Integer  nDisp(0:7)
      Character Label*39
      Real*8     Hess(nHess)
      Real*8, Allocatable:: Temp(:)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      Ind(idisp,jdisp)=idisp*(idisp-1)/2+jdisp
*                                                                      *
************************************************************************
*                                                                      *
 100  Format (A,A)
      iDisp=0
      Do iIrrep=0,nIrrep-1
           nDisp(iIrrep)=iDisp
           iDisp=iDisp+lDisp(iIrrep)
      End Do
*
      If (nirrep.eq.1) Then
         Write(Label,100) 'Hessian in Irrep ',lIrrep(0)
         Call TriPrt(Label,' ',Hess,ldisp(0))
      Else
         Call mma_allocate(Temp,nHess,Label='Temp')
         Do iIrrep=0,nIrrep-1
            Write(Label,100) 'Hessian in Irrep ',lIrrep(iIrrep)
            Do i=1,lDisp(iirrep)
               Do j=1,i
                  ii=ind(i,j)
                  jj=ind(ndisp(iirrep)+i,ndisp(iirrep)+j)
                  Temp(ii)=Hess(jj)
               End Do
            End Do
            Call TriPrt(Label,' ',Temp,ldisp(iirrep))
         End Do
         Call mma_deallocate(Temp)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
