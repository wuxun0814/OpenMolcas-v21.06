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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine CmbnVe(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,Final,nComp,
     &                  Vxyz)
************************************************************************
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '91                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp),
     &       Zeta(nZeta), rKappa(nZeta),
     &       Rnxyz(nZeta,3,0:la,0:lb+1,0:lr),
     &       Vxyz(nZeta,3,0:la,0:lb)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
*     iRout = 161
*     iPrint = nPrint(iRout)
*     Call GetMem(' Enter CmbnVe','LIST','REAL',iDum,iDum)
*
      Do 10 ixa = 0, la
         iyaMax=la-ixa
      Do 11 ixb = 0, lb
         iybMax=lb-ixb
         Do 20 iya = 0, iyaMax
            iza = la-ixa-iya
            ipa= Ind(la,ixa,iza)
         Do 21 iyb = 0, iybMax
            izb = lb-ixb-iyb
            ipb= Ind(lb,ixb,izb)
*           If (iPrint.ge.99) Then
*              Write (*,*) ixa,iya,iza,ixb,iyb,izb
*              Write (*,*) ipa,ipb
*           End If
*
*           Combine integrals
*
            Do 30 iZeta = 1, nZeta
               Fact = rKappa(iZeta) * Zeta(iZeta)**(-Three/Two)
               Final(iZeta,ipa,ipb,1) = Fact *
     &               Vxyz(iZeta,1,ixa,ixb)*
     &               Rnxyz(iZeta,2,iya,iyb,0)*
     &               Rnxyz(iZeta,3,iza,izb,0)
               Final(iZeta,ipa,ipb,2) = Fact *
     &               Rnxyz(iZeta,1,ixa,ixb,0)*
     &               Vxyz(iZeta,2,iya,iyb)*
     &               Rnxyz(iZeta,3,iza,izb,0)
               Final(iZeta,ipa,ipb,3) = Fact *
     &               Rnxyz(iZeta,1,ixa,ixb,0)*
     &               Rnxyz(iZeta,2,iya,iyb,0)*
     &               Vxyz(iZeta,3,iza,izb)
 30         Continue
*
 21      Continue
 20      Continue
 11   Continue
 10   Continue
*
*     Call GetMem(' Exit CmbnVe','LIST','REAL',iDum,iDum)
      Return
      End
