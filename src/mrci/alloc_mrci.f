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
      SUBROUTINE ALLOC_MRCI
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "warnings.fh"
#include "mrci.fh"
*------
* POW: This array is not used!
*     DIMENSION IPOF(9)
*------
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      NVSQ=0
      NVMAX=0
      DO 10 I=1,NSYM
        NVMAX=MAX(NVMAX,NVIR(I))
        NVSQ=NVSQ+NVIR(I)**2
10    CONTINUE
      NVT=(NVIRT*(NVIRT+1))/2
      if (NVIRT.eq.0) then
        Call SysAbendMsg('alloc_mrci.f:',
     &   'no virtual orbitals in the basis',' ')
      endif
*PAM04 MOVE ALLOCATION OF FOCK MATRIX TO SDCI.
*C FOCK MATRIX IN SORT,IIJJ,IJIJ,FIJ,AI.
*      LFOCK=LPERMA
*      LPERMX=LFOCK+NBTRI
*-----------------------------------------------
*PAM04 Change systematically:
*PAM04  (1) LPERMA, LPERMB, LPERMX set to 1
*PAM04  (2) MAXMEM changed to MEMWRK
      LPERMA=1
      LPERMB=1
      LPERMX=1
*-----------------------------------------------
CPAM97 C BUFFER FOR MOTRA INTEGRALS, TIBUF,  IN SORT, SORTA, SORTB.
CPAM97       LTIBUF=LPERMX
CPAM97       LPERMX=LTIBUF+NTIBUF
*PAM04      MEMX=MAXMEM-LPERMX+1
*PAM04      MEMX=MAXMEM
* PAM06: Evidently some miscounting margin needed:
      MEMX=INT(0.90D0*DBLE(MEMWRK))
C ALLOCATION FOR SORTA.
      NCHN1=LN*NVIRT+1
      IF(IFIRST.NE.0)NCHN1=1
      NBSIZ1=MEMX/NCHN1 - 1
      NBUFBI=KBUFF1
      NBSIZ1=MIN(NBSIZ1,MEMX-2*ISMAX-NBUFBI-1)
      NBSIZ1=MAX(NBSIZ1,256)
CPAM96 NBSIZ1=Size counted in real*8 words.
CPAM96 Must contain NBITM1 real*8 + NBITM1 integ + 2 integ:
CPAM96      NBITM1=2*((NBSIZ1-1)/3)
      NBITM1=(RTOI*NBSIZ1-2)/(RTOI+1)
      NBITM1=MIN(NBITM1,NVSQ)
      NBITM1=((NBITM1+2)/RTOI)*RTOI-2
CPAM96      NBSIZ1=(3*NBITM1+2)/2
      NBSIZ1=((RTOI+1)*NBITM1+2+(RTOI-1))/RTOI
C SORTING AREA, BUFOUT AND INDOUT
*PAM04      LBIN1=LPERMX
*PAM04      NBIN1=NBSIZ1*NCHN1
C BIAC
C NOTE: ONE SINGLE BIN IS IN USE TOGETHER WITH BIAC,BICA,BUFBI.
*PAM04      LBIAC1=LBIN1+NBSIZ1
*PAM04      LBIAC1=LPERMX
*PAM04C BICA
*PAM04      LBICA1=LBIAC1+ISMAX
*PAM04C BUFBI
*PAM04      LBUFBI=LBICA1+ISMAX
*PAM04*PAM04      LTOP1=MAX(LBUFBI+NBUFBI,LBIN1+NBIN1)-1
*PAM04      LTOP1=LBUFBI+NBUFBI-1
*PAM04      LTOP=LTOP1
      LTOP=LPERMX
C DYNAMIC ALLOCATION FOR SORTING ABCD
      NBITM2=1
*PAM04      LTOP2=0
      IF(IFIRST.EQ.0) THEN
        IPASS=0
110       IPASS=IPASS+1
          NCHN2=(NVT-1)/IPASS+1
          NBSIZ2=(MEMX-2*ISMAX-KBUFF1)/NCHN2
CPAM96          IF(2*NBSIZ2.GT.(3*NVSQ+2)) GOTO 120
          IF(RTOI*NBSIZ2.GT.((RTOI+1)*NVSQ+2)) GOTO 120
          IF(IPASS.EQ.5) GOTO 120
          IF(NBSIZ2.LT.1024) GOTO 110
120     CONTINUE
CPAM96        NBITM2=2*((NBSIZ2-1)/3)
        NBITM2=(RTOI*NBSIZ2-2)/(RTOI+1)
        NBITM2=MIN(NBITM2,NVSQ)
        NBITM2=((NBITM2+2)/RTOI)*RTOI-2
CPAM96        NBSIZ2=(3*NBITM2+2)/2
        NBSIZ2=((RTOI+1)*NBITM2+2+(RTOI-1))/RTOI
C SORTING BINS, BUFOUT AND INDOUT
*PAM04        LBIN2=LPERMX
*PAM04        NBIN2=NBSIZ2*NCHN2
C BFACBD, ACBDS, AND ACBDT:
*PAM04        LBACBD=LBIN2+NBIN2
*PAM04        LACBDS=LBACBD+KBUFF1
*PAM04        LACBDT=LACBDS+ISMAX
*PAM04        LTOP2=LACBDT+ISMAX-1
*PAM04        LTOP=MAX(LTOP,LTOP2)
      END IF
C DYNAMIC ALLOCATION FOR SORTING AIBJ, AND CREATING HDIAG:
      NOT2=IROW(LN+1)
      NCHN3=3*NOT2
*PAM04      LHDIAG=LPERMX
*PAM04      NHDIAG=MAX(NVT,IRC(1))
*PAM04      LIIJJ=LHDIAG+NHDIAG
*PAM04      LIJIJ=LIIJJ+NBTRI
*PAM04      NBSIZ3=(MAXMEM-LPERMX)/NCHN3
      NBSIZ3=(MEMWRK-LPERMX)/NCHN3
      NBSIZ3=MAX(NBSIZ3,256)
CPAM96      NBITM3=2*((NBSIZ3-1)/3)
      NBITM3=(RTOI*NBSIZ3-2)/(RTOI+1)
      NBITM3=MIN(NBITM3,NVSQ)
      NBITM3=((NBITM3+2)/RTOI)*RTOI-2
CPAM96      NBSIZ3=(3*NBITM3+2)/2
      NBSIZ3=((RTOI+1)*NBITM3+2+(RTOI-1))/RTOI
*PAM04      NBIN3=NBSIZ3*NCHN3
*PAM04      LTOP3=LBIN3+NBIN3
*PAM04      LTOP=MAX(LTOP,LTOP3)
C VECTORS PERMANENTLY IN CORE DURING CI ITERATIONS.
C LFOCK ALREADY ALLOCATED AT LPERMA. ALSO ALLOCATE DMO AND TDMO
C THERE, OVERLAYED, FOR FUTURE USE.
*PAM04 ALLOCATING DMO AND TDMO HAS BEEN MOVED TO SDCI.
*      LDMO=LPERMA
*      LTDMO=LPERMA
C AREF, EREF: EIGENVECTORS AND ENERGIES OF REFERENCE CI.
*PAM04 NOW AREF CAN START DIRECTLY AT LPERMA
*      LAREF=LPERMA+NBTRI
*      IF(ITRANS.EQ.1) LAREF=LPERMA+NBAST**2
*PAM04      LAREF=LPERMA
*PAM04      LEREF=LAREF+NREF**2
*PAM04      LPERMB=LEREF+NREF
C CALCULATE HOW MUCH SCRATCH WILL BE NEEDED FOR PERS PART:
C FIRST, SET ASIDE WHATS NEEDED FOR SIGMA GENERATION:
      NIJ=(LN*(LN+1))/2
      NIJKL=NIJ*(NIJ+1)/2
      NBMN=IAD10(1)
      IF(IFIRST.NE.0)NBMN=0
      NPER=5*NVSQ+NBSIZ3+2*NVMAX**2
CPAM97      NPER=MAX(NPER,NIJKL+1+KBUFF1/2)
      NPER=MAX(NPER,NIJKL)
      NPER=MAX(NPER,2*NBMN+2*ISMAX+KBUFF1)
      NPER=MAX(NPER,NBSIZ3+2*NVMAX**2+2*NVSQ)
C OVERLAY CI,((HREF,PLEN) & (SGM,PERS PART))
      NHREF=(NREF*(NREF+1))/2
      NPLEN=NREF
      NOVLY1=JSC(ILIM)+MAX(JSC(ILIM)+NPER,NHREF+NPLEN)
C THIS IS TO BE OVERLAYED WITH (CBUF,...,LSCR) IN MQCT. TWO ALT:
      NARR=11*NRROOT**2
*PAM04      MEMB=MAXMEM-LPERMB+1
      MEMB=MEMWRK-LPERMB+1
      MBUF1=MEMB-NOVLY1-NARR
      MBUF2=(MEMB-NARR-(3*NRROOT+2*MXVEC)*NSECT)/(3*MXVEC+2)
      MBUF=MIN(MBUF1,MBUF2,20249)
      MBUF=MAX(MBUF,1259)
C ICI, ONE BUFFER OF PACKED CI COEFFICIENTS:
*PAM04      LICI=LPERMB
C CI ARRAY:
*PAM04      LCI=LICI+MBUF
C SGM ARRAY:
*PAM04      LSGM=LCI+JSC(ILIM)
*PAM04      LPER=LSGM+JSC(ILIM)
*      LPER=LPERMB+JSC(ILIM)
*PAM04      LARR=LPER+NPER
C OVERLAY:
*PAM04      LHREF=LSGM
*PAM04      LPLEN=LHREF+NHREF
*PAM04      LARR=MAX(LARR,LPLEN+NPLEN)
C OVERLAY:
*PAM04      LMQ=LCI
*PAM04      LMQ=LPERMB
*PAM04      NMQ=(3*NRROOT+2*MXVEC)*NSECT+(3*MXVEC+1)*MBUF
*PAM04      LARR=MAX(LARR,LMQ+NMQ)
*PAM04      ltop4=larr+narr-1
*PAM04      ltop=max(ltop,ltop4)
C MORE DETAILED, PERS PART, FROM HERE -------------------------------
C DYNAMIC ALLOCATION FOR FAIBJ:
C MATRIX ABIJ
*PAM04      LABIJ=LPER
C MATRIX AIBJ
*PAM04      LAIBJ=LABIJ+NVSQ
C MATRIX AJBI
*PAM04      LAJBI=LAIBJ+NVSQ
C BUFIN, IBUFIN
*PAM04      LBFIN1=LAJBI+NVSQ
C A, SCRATCH AREA
*PAM04      LASCR1=LBFIN1+NBSIZ3
C B, SCRATCH AREA
*PAM04      LBSCR1=LASCR1+NVMAX**2
C F, SCRATCH AREA
*PAM04      LFSCR1=LBSCR1+NVMAX**2
C FSEC, SCRATCH AREA
*PAM04      LFSEC=LFSCR1+NVSQ
*PAM04      LTOP5=LFSEC+NVSQ-1
*PAM04      LTOP=MAX(LTOP,LTOP5)
C DYNAMIC ALLOCATION FOR IJKL
C FIJKL
CPAM04      LFIJKL=LPER
CPAM97C BUFIN, IBUFIN
CPAM97      LBFIN2=LFIJKL+NIJKL
CPAM97      NBFIN2=1+KBUFF1/2
CPAM97      LTOP6=LBFIN2+NBFIN2-1
CPAM97      LTOP=MAX(LTOP,LTOP6)
C DYNAMIC ALLOCATION FOR ABCI
C BMN
*PAM04      LBMN=LPER
C IBMN
*PAM04      LIBMN=LBMN+NBMN
C BIAC
*PAM04      LBIAC2=LIBMN+NBMN
C BICA
*PAM04      LBICA2=LBIAC2+ISMAX
C BUFIN
*PAM04      LBFIN3=LBICA2+ISMAX
*PAM04      NBFIN3=KBUFF1
*PAM04      LTOP7=LBFIN3+NBFIN3-1
*PAM04      LTOP=MAX(LTOP,LTOP7)
C DYNAMIC ALLOCATION FOR ABCD
C ACBDS
*PAM04      LAC1=LPER
C ACBDT
*PAM04      LAC2=LAC1+ISMAX
C BUFIN
*PAM04      LBFIN4=LAC2+ISMAX
*PAM04      NBFIN4=KBUFF1
*PAM04      LTOP8=LBFIN4+NBFIN4-1
*PAM04      LTOP=MAX(LTOP,LTOP8)
C DYNAMIC ALLOCATION FOR FIJ, AI AND AB
C BUFIN, IBUFIN
*PAM04      LBFIN5=LPER
C A, SCRATCH AREA
*PAM04      LASCR2=LBFIN5+NBSIZ3
C B, SCRATCH AREA
*PAM04      LBSCR2=LASCR2+NVMAX**2
C FK IN AI AND AB
*PAM04      LFSCR2=LBSCR2+NVMAX**2
C DBK
*PAM04      LDBK=LFSCR2+NVSQ
*PAM04      LTOP9=LDBK+NVSQ-1
*PAM04      LTOP=MAX(LTOP,LTOP9)
C ALLOCATION OF PERS PART ENDS HERE ------------------------------
C DYNAMIC ALLOCATION FOR NATURAL ORBITALS ETC.
*PAM04      NPRP=2*NCMO+NBAST+NBAST**2+3*NBTRI+MAX(NBTRI,NBMAX**2)
*PAM04      LPRP=LPERMB
*PAM04      LTOP10=LPRP+NPRP-1
*PAM04      LTOP=MAX(LTOP,LTOP10)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'      REQUIRED WORKSPACE SIZE:',LTOP
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'                    AVAILABLE:',MAXMEM
*PAM04      CALL XFLUSH(6)
*PAM04      IF((LTOP.LE.MAXMEM).AND.(IPRINT.LT.5)) RETURN
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)' DYNAMIC ALLOCATION INFORMATION:'
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  CSPCK:',LCSPCK
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' INTSYM:',LINTSY
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  INDX:',LINDX
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   ISAB:',LISAB
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  JREFX:',LJREFX
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  CISEL:',LCISEL
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
CPAM97      WRITE(6,'(A,I9)')'  TIBUF:',LTIBUF
CPAM97      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   FOCK:',LFOCK
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   BIN1:',LBIN1
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' NBSIZ1:',NBSIZ1
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  NCHN1:',NCHN1
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  BIAC1:',LBIAC1
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  BICA1:',LBICA1
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  BUFBI:',LBUFBI
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  LTOP1:',LTOP1
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   BIN2:',LBIN2
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' NBSIZ2:',NBSIZ2
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  NCHN2:',NCHN2
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' BFACBD:',LBACBD
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  ACBDS:',LACBDS
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  ACBDT:',LACBDT
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  LTOP2:',LTOP2
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  HDIAG:',LHDIAG
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  FIIJJ:',LIIJJ
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  FIJIJ:',LIJIJ
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   BIN3:',LBIN3
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' NBSIZ3:',NBSIZ3
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  NCHN3:',NCHN3
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  LTOP3:',LTOP3
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   AREF:',LAREF
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   EREF:',LEREF
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   HREF:',LHREF
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   PLEN:',LPLEN
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'    ICI:',LICI
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'     CI:',LCI
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'    SGM:',LSGM
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'    LMQ:',LMQ
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'    ARR:',LARR
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  LTOP4:',LTOP4
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   ABIJ:',LABIJ
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   AIBJ:',LAIBJ
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   AJBI:',LAJBI
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' BUFIN1:',LBFIN1
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  ASCR1:',LASCR1
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  BSCR1:',LBSCR1
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  FSCR1:',LFSCR1
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   FSEC:',LFSEC
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  LTOP5:',LTOP5
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
CPAM04      WRITE(6,'(A,I9)')'  FIJKL:',LFIJKL
CPAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'    BMN:',LBMN
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   IBMN:',LIBMN
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  BIAC2:',LBIAC2
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  BICA2:',LBICA2
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' BUFIN3:',LBFIN3
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  LTOP7:',LTOP7
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' ACBDS2:',LAC1
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' ACBDT2:',LAC2
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' BUFIN4:',LBFIN4
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  LTOP8:',LTOP8
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' BUFIN5:',LBFIN5
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  ASCR2:',LASCR2
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  BSCR2:',LBSCR2
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  FSCR2:',LFSCR2
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'    DBK:',LDBK
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'  LTOP9:',LTOP9
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'    DMO:',LDMO
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'   TDMO:',LTDMO
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'    PRP:',LPRP
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')' ltop10:',LTOP10
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,'(A,I9)')'    LTOP',LTOP
*PAM04      CALL XFLUSH(6)
*PAM04      WRITE(6,*)
*PAM04      CALL XFLUSH(6)
*PAM04      IF(LTOP.LE.MAXMEM) RETURN
      IF(LTOP.LE.MEMWRK) RETURN
      WRITE(6,*)'ALLOC Error: Too much workspace is needed.'
      WRITE(6,'(1X,A,2I10)')'      Needed LTOP=',LTOP
      WRITE(6,'(1X,A,2I10)')' Available MEMWRK=',MEMWRK
      CALL QUIT(_RC_GENERAL_ERROR_)
      END
