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
       subroutine vf (name,lun)
c
c      this routine open file vanisf file with a given name
c      name - name of the vanished file (I)
c      lun  - lun number with which file will be opened (I)
c
       character*8 name
       integer lun
c
       call molcas_open(lun,name)
c       open (unit=lun,file=name)
       write (lun,*) ' File scratched'
       close (lun)
c
       return
       end
