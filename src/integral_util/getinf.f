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
* Copyright (C) 1992,2020, Roland Lindh                                *
************************************************************************
      SubRoutine GetInf(DoRys,nDiff)
************************************************************************
*                                                                      *
* Object: to read all input information on the file INFO.              *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January 1992                                             *
************************************************************************
      Use Iso_C_Binding
      use Real_Spherical
      use Her_RW
      use External_Centers
      use Temporary_Parameters, only: Test
      use DKH_Info, only: DKroll
      use Sizes_of_Seward, only: S
      Implicit Real*8 (A-H,O-Z)
#include "stdalloc.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "nq_info.fh"
#include "status.fh"
      Logical DoRys
      Integer iix(2)
      Real*8 rix(2)
#include "SysDef.fh"
      nbyte_i = iiloc(iix(2)) - iiloc(iix(1))
      nbyte_r = idloc(rix(2)) - idloc(rix(1))
*
      Call GetInf_Internal(cRFStrt,iRFStrt,lRFStrt,rRFStrt,
     &                     cQStrt,iQStrt,rQStrt)
*
*     This is to allow type punning without an explicit interface
      Contains
      SubRoutine GetInf_Internal(cRFStrt,iRFStrt,lRFStrt,rRFStrt,
     &                           cQStrt,iQStrt,rQStrt)
      Integer, Target :: cRFStrt,iRFStrt,lRFStrt,cQStrt,iQStrt
      Real*8, Target :: rRFStrt,rQStrt
      Integer, Pointer :: p_cRF(:),p_iRF(:),p_lRF(:),p_cQ(:),p_iQ(:)
      Real*8, Pointer :: p_rRF(:),p_rQ(:)
*
*     Load the dynamic input area.
*
      Call Get_Info_Dynamic()
*
*     Load the static input area.
*
      Call Get_Info_Static()
*                                                                      *
      Call External_Centers_Get()
*                                                                      *
************************************************************************
*                                                                      *
*     Reaction field parameters
*
      Len = ilLoc(lRFEnd)-ilLoc(lRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(lRFStrt),p_lRF,[Len])
      Call Get_iArray('RFlInfo',p_lRF,Len)
*
      Len = idLoc(rRFEnd)-idLoc(rRFStrt)
      Len = (Len+nByte_r)/nByte_r
      Call C_F_Pointer(C_Loc(rRFStrt),p_rRF,[Len])
      Call Get_dArray('RFrInfo',p_rRF,Len)
*
      Len = iiLoc(iRFEnd)-iiLoc(iRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(iRFStrt),p_iRF,[Len])
      Call Get_iArray('RFiInfo',p_iRF,Len)
*
      Len = iiLoc(cRFEnd)-iiLoc(cRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(cRFStrt),p_cRF,[Len])
      Call Get_iArray('RFcInfo',p_cRF,Len)
*
      Nullify(p_lRF,p_rRF,p_iRF,p_cRF)
*                                                                      *
************************************************************************
*                                                                      *
*     Numerical integration information and parameters
*
      Len = idLoc(rQEnd)-idLoc(rQStrt)
      Len = (Len+nByte_r)/nByte_r
      Call C_F_Pointer(C_Loc(rQStrt),p_rQ,[Len])
      Call Get_dArray('Quad_r',p_rQ,Len)
*
      Len = iiLoc(iQEnd)-iiLoc(iQStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(iQStrt),p_iQ,[Len])
      Call Get_iArray('Quad_i',p_iQ,Len)
*
      Len = iiLoc(cQEnd)-iiLoc(cQStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(cQStrt),p_cQ,[Len])
      Call Get_iArray('Quad_c',p_cQ,Len)
*
      Nullify(p_rQ,p_iQ,p_cQ)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate the transformation matrices
*
      If (S%iAngMx-1.ge.lMax) Then
         Call Sphere(S%iAngMx)
         lmax_internal=S%iAngMx
      Else
         Call Sphere(lMax)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Set the highest number of differentiations which will be
*     applied to the basis functions. In this case 2 + 1 ( the
*     kinetic energy operator and a differentiaion with respect
*     to a nuclear coordinate.
*
      nPrp=Max(lMax,3)
*
*     Setup of tables for coefficients of the Rys roots and weights.
*
      If (S%iAngMx.eq.0) nDiff=2
      If (DKroll.and.nOrdEF.gt.0) nDiff=nDiff+nOrdEF
      If (.Not.Test) Call Setup_RW(DoRys,nDiff)
*                                                                      *
************************************************************************
*                                                                      *
*     Set up for contracted calculation
*
      Call Flip_Flop(.False.)
*                                                                      *
************************************************************************
*                                                                      *
      Call Get_EFP()
*                                                                      *
************************************************************************
*
      Return
      End SubRoutine GetInf_Internal
*
      End
