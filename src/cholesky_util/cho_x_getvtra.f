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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
*  Cho_X_getVtra
*
*> @brief
*>    This routine performs a half-MO-transformation of Cholesky vectors stored in reduced storage
*> @author F. Aquilante
*>
*> @details
*> This routine performs a half-MO-transformation of Cholesky vectors stored in reduced
*> storage. For \p DoRead = ``.true.`` the vectors are read from
*> disk using array \p RedVec as scratch space, whereas for
*> \p DoRead = ``.false.`` the reduced vectors must be supplied in
*> array \p RedVec.
*>
*> Given (\p ChoT),the target arrays of type SBA_Type,
*> the routine performs a half-MO-transformation of \p NUMV Cholesky
*> vectors of compound symmetry \p ISYM starting with
*> vector \p IVEC1 and returns them in the target arrays.
*>
*> - \p iSwap = ``0``: \f$ L(k,b,J) \f$ is returned
*> - \p iSwap = ``1``: \f$ L(a,k,J) \f$ is returned
*> - \p iSwap = ``2``: \f$ L(k,J,b) \f$ is returned
*> - \p iSwap = ``3``: \f$ L(a,J,k) \f$ is returned
*>
*>
*> - \p IREDC: reduced set in core at the moment of the call to the routine.
*>             Can be set to ``-1`` (= unknown or undefined) by the calling routine.
*>
*> @param[out] irc     return code
*> @param[in]  RedVec  Vectors stored in reduced set(s) [\p DoRead option off] or scratch space for reading reduced vectors [\p DoRead option on]
*> @param[in]  lRedVec size of the \p RedVec
*> @param[in]  IVEC1   first vector to read
*> @param[in]  NUMV    number of vectors to transform starting from \p IVEC1
*> @param[in]  ISYM    compound symmetry of the Cholesky vectors
*> @param[in]  iSwap   type of the full storage for the half transformed Cholesky vectors
*> @param[in]  IREDC   reduced set in core
*> @param[in]  nDen    total number of densities to which MOs refer
*> @param[in]  kDen    first density for which the MO transformation has to be performed
*> @param[in]  MOs     the MOs coefficients stored in the data type DSBA_Type, i.e. symmetry blocked.
*> @param[in]  ChoT    the half transformed vectors, symmetry blocked as type SBA_Type
*> @param[in]  DoRead  flag for reading the reduced vectors
************************************************************************
      Subroutine Cho_X_getVtra(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,
     &                         iSwap,IREDC,nDen,kDen,MOs,ChoT,
     &                         DoRead)
      use Data_Structures, only: DSBA_Type, SBA_Type
      Implicit Real*8 (a-h,o-z)

      Type (DSBA_Type) MOs(nDen)
      Type (SBA_Type) Chot(nDen)

      Dimension RedVec(lRedVec)
      Integer   nDen,kDen
      Logical   DoRead
      Character(LEN=13), Parameter:: SECNAM = 'Cho_X_GetVtra'

#include "real.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      MXUSD = 0
      MUSED = 0

C zeroing the target arrays
C--------------------------
      Do jDen=kDen,nDen
         ChoT(jDen)%A0(:)=Zero
      End Do
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
      IF (DoRead) THEN
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
       JVEC1 = IVEC1             ! Absolute starting index
       IVEC2 = JVEC1 + NUMV - 1  ! Absolute ending index

       Do While (jVec1.le.iVec2)

        Call CHO_VECRD(RedVec,lRedVec,JVEC1,IVEC2,ISYM,JNUM,IREDC,MUSED)

        MXUSD = MAX(MXUSD,MUSED)

        If (JNUM.le.0 .or. JNUM.gt.(IVEC2-JVEC1+1)) then
           irc=77
           RETURN
        End If

        JVREF = JVEC1 - IVEC1 + 1 ! Relative index

        Call cho_vTra(irc,RedVec,lRedVec,JVREF,JVEC1,JNUM,NUMV,ISYM,
     &                IREDC,iSwap,nDen,kDen,MOs,ChoT)

        if (irc.ne.0) then
           return
        endif

        JVEC1 = jVec1 + JNUM

       End Do  ! end the while loop
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
      ELSE ! only MO transformation
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
       JNUM = NUMV
       JVREF= 1
       Call cho_vTra(irc,RedVec,lRedVec,JVREF,IVEC1,JNUM,NUMV,ISYM,
     &               IREDC,iSwap,nDen,kDen,MOs,ChoT)

       if (irc.ne.0) then
          return
       endif
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
      END IF
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
      irc=0

      RETURN
      END
