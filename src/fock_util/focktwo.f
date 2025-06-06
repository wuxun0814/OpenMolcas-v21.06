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
* Copyright (C) Markus P. Fuelscher                                    *
*               1992,1994, Per Ake Malmqvist                           *
*               2002, Roland Lindh                                     *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE FOCKTWO(NSYM,NBAS,NFRO,KEEP,
     &                   W_DLT,W_DSQ,W_FLT,nFlt,W_FSQ,LBUF,X1,X2,ExFac)
      use Data_Structures, Only: DSBA_Type, Allocate_DSBA,
     &                           Deallocate_DSBA
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 W_FSQ(*),W_FLT(nFlt),W_DSQ(*),W_DLT(*),X1(*),X2(*)
      Integer KEEP(8),NBAS(8),NFRO(8)

      Type (DSBA_Type) DLT, FLT, DSQ, FSQ
c
c This routine has been nicked from the MOTRA package. It was
c originally written by Marcus Fuelscher, and has been slightly
c modified by P-A Malmqvist 1992-12-04.
c Further modifications by RL 2002-8-30.
c Purpose: Return FLT, which is the effective one-electron
c Hamiltonian when frozen orbitals are removed. It is the
c same as a closed-shell Fock matrix computed using the
c two-electron contribution only from frozen orbitals.
c FLT is returned as symmetry-blocked lower triangles. FSQ
c contains the same matrix, as symmetry-blocked square matrices.
c DSQ and DLT are computed in the calling routine.
c It is assumed that the SEWARD integral file was opened before
c call to this routine.
c
************************************************************************
*                                                                      *
      MUL(I,J)=IEOR(I-1,J-1)+1
*                                                                      *
************************************************************************
      Call Allocate_DSBA(DLT,nBas,nBas,nSym,Case='TRI',Ref=W_DLT)
      Call Allocate_DSBA(FLT,nBas,nBas,nSym,Case='TRI',Ref=W_FLT)
      Call Allocate_DSBA(DSQ,nBas,nBas,nSym,Ref=W_DSQ)
      Call Allocate_DSBA(FSQ,nBas,nBas,nSym,Ref=W_FSQ)

      DO 110 IS=1,NSYM
        IB=NBAS(IS)
        IK=KEEP(IS)
        NFI=NFRO(IS)
        DO 120 JS=1,IS
          JB=NBAS(JS)
          JK=KEEP(JS)
          NFJ=NFRO(JS)
          IJS=MUL(IS,JS)
          IJB=IB*JB
          IF( IS.EQ.JS ) IJB=(IB*(IB+1))/2
          DO 130 KS=1,IS
            KB=NBAS(KS)
            KK=KEEP(KS)
            NFK=NFRO(KS)
            LSMAX=KS
            IF ( KS.EQ.IS ) LSMAX=JS
            LS=MUL(IJS,KS)
            IF(LS.GT.LSMAX) GOTO 130
            LB=NBAS(LS)
            LK=KEEP(LS)
            NFL=NFRO(LS)
            KLB=KB*LB
            IF( KS.EQ.LS ) KLB=(KB*(KB+1))/2
C INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
            IF((IK+JK+KK+LK).NE.0) GOTO 130
C NO FROZEN ORBITALS?
            IF((NFI+NFJ+NFK+NFL).EQ.0) GOTO 130
C NO BASIS FUNCTIONS?
            IF((IJB*KLB).EQ.0 ) GOTO 130

            IF ( IS.EQ.JS .AND. IS.EQ.KS ) THEN
c CASE 1: Integrals are of symmetry type (II/II)
c Coulomb and exchange terms need to be accumulated
c Option code 1: Begin reading at first integral.
c NPQ: Nr of submatrices in buffer X1.
                IOPT=1
                LPQ=0
                IPQ=0
                NPQ=0
                DO 200 IP=1,IB
                  DO 201 JQ=1,IP
                    IPQ=IPQ+1
                    LPQ=LPQ+1
                    IF ( IPQ.GT.NPQ ) THEN
                      CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
                      IF(IRC.GT.1) GOTO 999
c Option code 2: Continue reading at next integral.
                      IOPT=2
                      IPQ=1
                    ENDIF
                    ISX=(IPQ-1)*KLB+1
                    TEMP=DDOT_(KLB,X1(ISX),1,DLT%SB(IS)%A1,1)
                    FLT%SB(IS)%A1(LPQ)=FLT%SB(IS)%A1(LPQ)+TEMP
                    CALL SQUARE (X1(ISX),X2(1),1,KB,LB)
                    CALL DGEMV_('N',KB,LB,
     &                        (-0.5D0*ExFac),X2(1),KB,
     &                                       DSQ%SB(IS)%A2(1:,IP),1,
     &                                 1.0D0,FSQ%SB(IS)%A2(1:,JQ),1)
                    IF ( IP.NE.JQ ) THEN
                      CALL DGEMV_('N',KB,LB,
     &                          (-0.5D0*ExFac),X2(1),KB,
     &                                         DSQ%SB(IS)%A2(1:,JQ),1,
     &                                   1.0D0,FSQ%SB(IS)%A2(1:,IP),1)
                    ENDIF
201               CONTINUE
200             CONTINUE
              ELSE IF ( IS.EQ.JS .AND. IS.NE.KS ) THEN
c CASE 2: Integrals are of symmetry type (II/JJ)
c Coulomb terms need to be accumulated only
                IOPT=1
                LPQ=0
                IPQ=0
                NPQ=0
                DO 210 IP=1,IB
                  DO 211 JQ=1,IP
                    IPQ=IPQ+1
                    LPQ=LPQ+1
                    IF ( IPQ.GT.NPQ ) THEN
                      CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
                      IF(IRC.GT.1) GOTO 999
                      IOPT=2
                      IPQ=1
                    ENDIF
                    ISX=(IPQ-1)*KLB+1
                    IF ( NFI.NE.0 ) THEN
                      TEMP=DLT%SB(IS)%A1(LPQ)
                      CALL DAXPY_(KLB,TEMP,X1(ISX),1,FLT%SB(KS)%A1,1)
                    ENDIF
                    IF ( NFK.NE.0 ) THEN
                      TEMP=DDOT_(KLB,X1(ISX),1,DLT%SB(KS)%A1,1)
                      FLT%SB(IS)%A1(LPQ)=FLT%SB(IS)%A1(LPQ)+TEMP
                    ENDIF
211               CONTINUE
210             CONTINUE
              ELSE IF ( IS.EQ.KS .AND. JS.EQ.LS ) THEN
c CASE 3: Integrals are of symmetry type (IJ/IJ)
c Exchange terms need to be accumulated only
                IOPT=1
                LPQ=0
                IPQ=0
                NPQ=0
                DO 220 IP=1,IB
                  DO 221 JQ=1,JB
                    IPQ=IPQ+1
                    LPQ=LPQ+1
                    IF ( IPQ.GT.NPQ ) THEN
                      CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
                      IF(IRC.GT.1) GOTO 999
                      IOPT=2
                      IPQ=1
                    ENDIF
                    ISX=(IPQ-1)*KLB+1
                    IF ( NFI.NE.0 ) THEN
                      CALL DGEMV_('N',LB,KB,
     &                          (-0.5D0*ExFac),X1(ISX),LB,
     &                                         DSQ%SB(IS)%A2(1:,IP),1,
     &                                   1.0D0,FSQ%SB(JS)%A2(1:,JQ),1)
                    ENDIF
                    IF ( NFJ.NE.0 ) THEN
                      CALL DGEMV_('T',LB,KB,
     &                          (-0.5D0*ExFac),X1(ISX),LB,
     &                                         DSQ%SB(JS)%A2(1:,JQ),1,
     &                                   1.0D0,FSQ%SB(IS)%A2(1:,IP),1)
                    ENDIF
221               CONTINUE
220             CONTINUE
            ENDIF
130       CONTINUE
120     CONTINUE
110   CONTINUE

c Accumulate the contributions
      DO 300 ISYM=1,NSYM
        NB=NBAS(ISYM)
        DO 310 IB=1,NB
          DO 315 JB=1,IB
            IJ = IB*(IB-1)/2 + JB
            FLT%SB(ISYM)%A1(IJ)=FLT%SB(ISYM)%A1(IJ)
     &                         +FSQ%SB(ISYM)%A2(IB,JB)
315       CONTINUE
310     CONTINUE
300   CONTINUE
*
      Call GADSum(W_FLT,nFlt)
*
#ifdef _DEBUGPRINT_
      WRITE(6,'(6X,A)')'FROZEN FOCK MATRIX IN AO BASIS:'
      DO ISYM=1,NSYM
        NB=NBAS(ISYM)
        IF ( NB.GT.0 ) THEN
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          CALL TRIPRT(' ',' ',FLT%SB(ISYM)%A1,NB)
        END IF
      END DO
#endif
      Call deallocate_DSBA(FSQ)
      Call deallocate_DSBA(DSQ)
      Call deallocate_DSBA(FLT)
      Call deallocate_DSBA(DLT)

      RETURN
 999  CONTINUE
      WRITE(6,*)' Error return code IRC=',IRC
      WRITE(6,*)' from RDORD call, in FTWOI.'
      CALL Abend
      END
