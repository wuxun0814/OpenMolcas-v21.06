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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine indxab_cvb(indxa,indxb,nstra,nstrb,nsa,nsb)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension indxa(nsa),indxb(nsb),nstra(mxirrep),nstrb(mxirrep)

      i1 = mstacki_cvb(norb+1)

      call indxab2_cvb(indxa,indxb,nstra,nstrb,
     >  iw(i1),nsa,nsb)
      call mfreei_cvb(i1)
      return
      end
