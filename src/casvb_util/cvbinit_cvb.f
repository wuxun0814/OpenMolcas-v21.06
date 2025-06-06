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
      subroutine cvbinit_cvb()
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "maxbfn.fh"
      parameter(iset=1)
      save is_set
      data is_set/0/

      if(is_set.eq.iset)return
      entry cvbfinit_cvb()

      mxaobf=maxbfn
      iprec=8
      iwidth=110
      call formats_cvb()
      call setidbl_cvb()
      call meminit_cvb()
      if(is_set.ne.iset)then
c  Initializations below are only carried out once :
        call io_init_cvb()
        call main_bdata_cvb()
      endif
      is_set=iset
      return
      end
