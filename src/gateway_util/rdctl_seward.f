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
      Subroutine RdCtl_Seward(LuRd,lOPTO,Do_OneEl)
      use SW_File
      use AMFI_Info
      use Basis_Info
      use Center_Info
      use Her_RW
      use Period
      use MpmC
      use EFP_Module
      use Real_Spherical, only: Sphere
      use fortran_strings, only : str
      use External_Centers
      use Symmetry_Info, only: Symmetry_Info_Setup, iSkip, nIrrep,
     &                         VarR, VarT
      use Temporary_Parameters
      use Integral_Parameters
      use Sizes_of_Seward, Only: S
      use Real_Info, only: ThrInt, Rtrnc, CutInt, PkAcc, Thrs, E1, E2,
     &                     RPQMin, SadStep, Shake, kVector, CoM
      use DKH_Info
      use RICD_Info, only: iRI_Type, LDF, Do_RI, Cholesky,
     &                     Do_acCD_Basis, Skip_High_AC, DiagCheck,
     &                     LocalDF, Do_nacCD_Basis, Thrshld_CD
      use Logical_Info
      use Gateway_Interfaces, only: GetBS
      use Gateway_global, only: Run_Mode, G_Mode, S_Mode
#ifdef _FDE_
      use Embedding_Global, only: embPot, embPotInBasis, embPotPath,
     &outGridPathGiven, embWriteDens, embWriteEsp, embWriteGrad,
     &embWriteHess
#endif
#ifndef _HAVE_EXTRA_
      use XYZ
#endif
      Use Para_Info, Only: MyRank
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: Is_Real_Par
#endif
      Implicit Real*8 (a-h,o-z)
      External NucExp
#include "Molcas.fh"
*
#include "angtp.fh"
#include "constants.fh"
#include "constants2.fh"
#include "SysDef.fh"
#include "notab.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "rctfld.fh"
#include "rmat.fh"
#include "real.fh"
#include "print.fh"
#ifdef _HAVE_EXTRA_
#include "hyper.fh"
#endif
*
      Real*8 Lambda
      Character Key*180, KWord*180, Oper(3)*3, BSLbl*80, Fname*256,
     &          DefNm*13, Ref(2)*180, ChSkip*80, AngTyp(0:iTabMx)*1,
     &          dbas*(LENIN),filename*180, KeepBasis*256, KeepGroup*180,
     &          Previous_Command*12, CtrLDK(10)*(LENIN),
     &          Directory*256, BasLib*256,ExtBasDir*256
      Character(LEN=72):: Header(2)=['','']
      Character(LEN=80):: Title(10)=['','','','','','','','','','']
      Character(LEN=14):: Vrsn='Gateway/Seward'
      Character(LEN=512):: Align_Weights='MASS'
#include "cgetl.fh"
      Character*180 Get_Ln
      External Get_Ln
      Logical lTtl, lSkip, lMltpl, DoRys, RF_read, Convert, IfTest,
     &        Exist,CutInt_UsrDef, ThrInt_UsrDef, MolWgh_UsrDef,
     &        CholeskyWasSet, GWInput, NoAMFI, lOPTO, Do_OneEl
      Logical do1CCD
      Logical:: CSPF=.False.
      Logical APThr_UsrDef, Write_BasLib
      Integer Cho_MolWgh, BasisTypes(4), BasisTypes_Save(4),
     &        iGeoInfo(2), iOpt_XYZ, RC
      Parameter (Cho_CutInt = 1.0D-40, Cho_ThrInt = 1.0D-40,
     &           Cho_MolWgh = 2)
*
      Real*8 NucExp, WellCff(3), WellExp(3), WellRad(3),
     &       OAMt(3), OMQt(3)
      Real*8, Allocatable :: RTmp(:,:), EFt(:,:),
     &                       DMSt(:,:), OrigTrans(:,:), OrigRot(:,:,:),
     &                       mIsot(:)
      Integer, Allocatable :: ITmp(:), nIsot(:,:), iScratch(:)
!     Temporary buffer
      Integer, Parameter:: nBuff=10000
      Real*8, Allocatable:: Buffer(:), Isotopes(:)
!
      Character*180, Allocatable :: STDINP(:)
      Character Basis_lib*256, CHAR4*4
      Character*256 Project, GeoDir, temp1, temp2
*
      Integer StrnLn
      External StrnLn

      Logical SymmSet
      Logical CoordSet,RPSet
      Logical BasisSet
      Logical GroupSet
      Logical DoneCoord
      Logical NoZMAT
      Logical ForceZMAT
      Logical NoDKroll
      Logical DoTinker
      Logical DoGromacs
      Logical OrigInput
      Logical OriginSet
      Logical FragSet
      Logical HyperParSet
      Logical WriteZMat
#ifdef _HAVE_EXTRA_
      Logical geoInput, oldZmat, zConstraints
#endif
      Logical EFgiven
      Logical Invert
      Real*8 HypParam(3), RandVect(3)
      Logical Vlct_, nmwarn, FOUND
*
      Logical Basis_test, lECP, lPP
      Logical :: lDMS=.FALSE., lOAM=.FALSE., lOMQ=.False.,
     &           lXF=.False., lFAIEMP=.False.
#include "embpcharg.fh"
*
#ifdef _GROMACS_
      Integer, Dimension(:), Allocatable :: CastMM
      Integer, Dimension(:,:), Allocatable :: DefLA
      Real*8, Dimension(:), Allocatable :: FactLA
#endif
      Character*256 Message
      parameter (MAX_XBAS=20)
      Character*12 xb_label(MAX_XBAS)
      Character*128 xb_bas(MAX_XBAS)
*
      Data WellCff/.35D0,0.25D0,5.2D0/
      Data WellExp/4.0D0,3.0D0,2.0D0/
      Data WellRad/-1.22D0,-3.20D0,-6.20D0/
*
#include "angstr.fh"
      Data DefNm/'basis_library'/
      Data IfTest/.False./
*                                                                      *
************************************************************************
*                                                                      *
      iRout=3
      iPrint = nPrint(iRout)
#ifdef _DEBUGPRINT_
      IfTest=.True.
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(Buffer,nBuff,Label='Buffer')
*                                                                      *
************************************************************************
*                                                                      *
      Do i=0,iTabMx
         AngTyp(i)=Angtp(i)
         Call UpCase(AngTyp(i))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call WhichMolcas(Basis_lib)
      If (Basis_lib(1:1).ne.' ') Then
         ib=index(Basis_lib,' ')-1
         if(ib.lt.1)
     &   Call SysAbendMsg('rdCtl','Too long PATH to MOLCAS',' ')
         BasLib=Basis_lib(1:ib)//'/basis_library'
      Else
         BasLib=DefNm
      End If
      Write_BasLib=.False.
*                                                                      *
************************************************************************
*                                                                      *
      Call Qpg_cArray('Align_Weights',Found,lAW)
      If (Found) Call Get_cArray('Align_Weights',Align_Weights,512)
*                                                                      *
************************************************************************
*                                                                      *
      CutInt_UsrDef=.False.
      ThrInt_UsrDef=.False.
      MolWgh_UsrDef=.False.
      APThr_UsrDef=.False.
      NoAMFI=.False.
*
      iChk_RI=0
      iChk_CH=0
      iChk_DC=0
      iOpt_XYZ=-1
*
      isXfield=0
      CholeskyThr=-9.99d9
*
      Basis_Test=.False.
*                                                                      *
************************************************************************
*                                                                      *
      LuWr=6
      LuFS=-1
      LuRdSave=-1
*
      nTemp=0
      lMltpl=.False.
*
      CholeskyWasSet=.False.
      do1CCD=.false.
      spanCD=-9.9d9
      lTtl = .False.
      RF_read=.False.
      lSkip=.False.
      NoDKroll=.false.
      SymmSet=.false.
      BasisSet=.false.
      GroupSet=.false.
      RPSet=.false.
      CoordSet=.false.
      DoneCoord=.false.
      KeepGroup='FULL'
      NoZMAT=.false.
      ForceZMAT=.false.
      DoTinker = .False.
      DoGromacs = .False.
      OrigInput = .False.
#ifdef _HAVE_EXTRA_
      origin_input = .False.
      geoInput = .False.
      OldZmat = .False.
      isHold=-1
      nCoord=0
      ZConstraints = .False.
#endif
      WriteZMat = .False.
      nFragment = 0
      iFrag = 0
      FragSet = .False.
      OriginSet = .False.
      HyperParSet = .False.
      stepFac1 = 60.0d0
      iOptimType = 1
      gradLim = 0.0d0
      Do_OneEl=.True.
      Vlct_=.False.
#ifdef _FDE_
      ! Embedding
      embPot=.false.
      embPotInBasis=.false.
      embPotPath='EMBPOT'
      outGridPathGiven=.false.
      embWriteDens=.false.
      embWriteEsp=.false.
      embWriteGrad=.false.
      embWriteHess=.false.
      Call Put_iScalar('embpot',0)
#endif
      DoEmPC=.False.
      EFgiven=.False.
      Invert=.False.
      Call Put_iScalar('agrad',0)
*
      ScaleFactor=1.0d0
      lSTDINP=0
      iCoord=0
      iBSSE=-1
      SymThr=0.01D0
      nTtl=0
*
      imix=0
      ifnr=-1
      ign=0
      itype=0
      ExtBasDir=' '
      isxbas=0
      nmwarn=.True.
*
*     Selective initialization
*
      If (Run_Mode.eq.S_Mode) Then
         iShll = S%Mx_Shll
         mdc = S%Mx_mdc
      Else
         iShll = 0
         mdc = 0
*
         ChSkip=' '
         Do i = 1,3
            Oper(i)=' '
         End Do
         nOper=0
      End If
*
      iDNG=0
      BasisTypes(:)=0
      KeepBasis=' '
cperiod
      lthCell = 0
      Cell_l = .FALSE.
      Call izero(ispread,3)
      Call fzero(VCell,9)
*     Set local DF variables (dummy)
      Call LDF_SetInc()
      Rewind(LuRd)
*     Count the number of calls to coord to set the number of fragments
*     for a geo/hyper-calculation
 400  Key= Get_Ln(LuRd)
      KWord = Key
      Call UpCase(KWord)
      If(KWord(1:4) .eq. 'COOR') nFragment=nFragment+1
      If((KWord(1:4) .eq. 'HYPE') .or.
     &   (KWord(1:4) .eq. 'GEO ')) Then
         call getenvf("Project",Project)
         call getenvf("GeoDir",GeoDir)
         temp1 = Project(1:index(Project,' ')-1)//
     &           '.Gateway.Input'
         temp2 = GeoDir(1:index(GeoDir,' ')-1)//'/'//
     &           Project(1:index(Project,' ')-1)//'.gwcopy.input'
         Call fCopy(temp1,temp2,ierr)
         if(ierr.ne.0) then
           write(6,*) '*** Detect Hyper input, but no GEO loop'
           call Quit_OnUserError()
         EndIf
      End If
      If (KWord(1:4).eq.'END ') Go To 401
      GoTo 400
 401  Call Put_iScalar('nCoordFiles',nFragment)
      Rewind(LuRd)
*                                                                      *
************************************************************************
*                                                                      *
      If (Run_Mode.eq.G_Mode) Then
         Call RdNLst(LuRd,'GATEWAY')
         Previous_Command='&Gateway'
      Else
         Call RdNLst(LuRd,'SEWARD')
         Previous_Command='&Seward'
      End If
*
*     Default setting of GWInput
*
      GWInput=Run_Mode.eq.G_Mode
*
*     GWInput is a logical flag which tells if a keyword is Gateway or
*     Seward specific. This is done down below for each keyword in the
*     case that the keyword is not Seward or Gateway specific enter the
*     line,
*
*     GWInput = Run_Mode.eq.G_Mode
*
*     in the section of the particular keyword
*                                                                      *
************************************************************************
*                                                                      *
*
*     KeyWord directed input
*
      Call mma_allocate(STDINP,MxAtom*2,label='STDINP')
*
      nDone=0
 998  lTtl = .False.
      If (Basis_Test.and.nDone.eq.1) Then
         nDone=0
         Basis_Test=.False.
      End If
 9988 Continue
      Key = Get_Ln(LuRd)
*
 9989 If (Run_Mode.eq.G_Mode.and..Not.GWInput) Then
         Call WarningMessage(2,'Gateway input error!')
         Write (LuWr,*) 'The keyword : "',Previous_Command,
     &                  '" is not allowed in the Gateway input!'
         Write (LuWr,*) 'This keyword most likely belongs in the '
     &                //'Seward input section!.'
         Write (LuWr,*)
         Call Quit_OnUserError()
      Else If (Run_Mode.eq.S_Mode.and.GWInput) Then
         Call WarningMessage(2,'Seward input error!')
         Write (LuWr,*) 'The keyword : "',Previous_Command,
     &                  '" is not allowed'
     &               //' in the Seward input when the Gateway is used!'
         Write (LuWr,*) ' Try putting the keyword in the Gateway '
     &                //'input section!'
         Write (LuWr,*)
         Call Quit_OnUserError()
      End If
      GWInput=.False.
*
      If (IfTest) Write (LuWr,*) ' RdCtl: Processing:',Key
      KWord = Key
      Call UpCase(KWord)
      Previous_Command=KWord(1:4)
      If (KWord(1:1).eq.'*') Go To 998
      If (KWord.eq.'')    Go To 998
      If (Basis_Test) nDone=1
*
*     KEYWORDs in ALPHABETIC ORDER!
*
      If (KWord(1:4).eq.'1C-C') Go To 9022
      If (KWord(1:4).eq.'1CCD') Go To 9022
      If (KWord(1:4).eq.'ACCD') Go To 8005
      If (KWord(1:4).eq.'ACD ') Go To 8004
      If (KWord(1:4).eq.'ALIG') Go To 8013
      If (KWord(1:4).eq.'AMFI') Go To 9761
      If (KWord(1:4).eq.'AMF1') Go To 8761
      If (KWord(1:4).eq.'AMF2') Go To 8762
      If (KWord(1:4).eq.'AMF3') Go To 8763
      If (KWord(1:4).eq.'AMPR') Go To 9951
      If (KWord(1:4).eq.'ANGM') Go To 995
      If (KWord(1:4).eq.'APTH') Go To 38
      If (KWord(1:4).eq.'AUXS') Go To 912
      If (KWord(1:4).eq.'BASD') Go To 7700
      If (KWord(1:4).eq.'BASI') Go To 920
      If (KWord(1:4).eq.'BASL') Go To 8029
      If (KWord(1:4).eq.'BSSE') Go To 6020
      If (KWord(1:4).eq.'BSSH') Go To 9121
      If (KWord(1:4).eq.'BSSM') Go To 9760
      If (KWord(1:4).eq.'CDTH') Go To 8001
      If (KWord(1:4).eq.'CELL') Go To 887
      If (KWord(1:4).eq.'CENT') Go To 973
      If (KWord(1:4).eq.'CHEC') Go To 39
      If (KWord(1:4).eq.'CLDF') Go To 42
      If (KWord(1:4).eq.'CHOL') Go To 9091
      If (KWord(1:4).eq.'CHOI') Go To 9092
      If (KWord(1:4).eq.'CLIG') Go To 9000
      If (KWord(1:4).eq.'CONS') Go To 8010
      If (KWord(1:4).eq.'COOR') Go To 6000
      If (KWord(1:4).eq.'CSPF') Go To 9110
      If (KWord(1:4).eq.'CUTO') Go To 942
      If (KWord(1:4).eq.'DCRN') Go To 958
      If (KWord(1:4).eq.'DIAG') Go To 9087
      If (KWord(1:4).eq.'DIRE') Go To 9770
      If (KWord(1:4).eq.'DK1H') Go To 9001
      If (KWord(1:4).eq.'DK2H') Go To 9002
      If (KWord(1:4).eq.'DK3F') Go To 9004
      If (KWord(1:4).eq.'DK3H') Go To 9003
      If (KWord(1:4).eq.'DOAN') Go To 8019
      If (KWord(1:4).eq.'DOFM') Go To 8006
      If (KWord(1:4).eq.'DOUG') Go To 976
      If (KWord(1:4).eq.'DOWN') Go To 1002
      If (KWord(1:4).eq.'DSHD') Go To 996
      If (KWord(1:4).eq.'ECPS') Go To 912
      If (KWord(1:4).eq.'EFLD') Go To 993
      If (KWord(1:4).eq.'EFP ') Go To 9088
#ifdef _FDE_
      If (KWord(1:4).eq.'EMBE') Go To 666
#endif
      If (KWord(1:4).eq.'EMFR') Go To 8035
      If (KWord(1:4).eq.'EMPC') Go To 974
      If (KWord(1:4).eq.'EPOT') Go To 9932
      If (KWord(1:4).eq.'EXPE') Go To 9771
      If (KWord(1:4).eq.'EXTR') Go To 9960
      If (KWord(1:4).eq.'FAKE') Go To 9759
      If (KWord(1:4).eq.'FAT-') Go To 8004
      If (KWord(1:4).eq.'FCD ') Go To 9091
      If (KWord(1:4).eq.'FILE') Go To 904
      If (KWord(1:4).eq.'FINI') Go To 9762
      If (KWord(1:4).eq.'FLDG') Go To 994
      If (KWord(1:4).eq.'FNMC') Go To 9086
      If (KWord(1:4).eq.'FOOC') Go To 8000
      If (KWord(1:4).eq.'FPCO') Go To 9764
      If (KWord(1:4).eq.'FPPR') Go To 9765
      If (KWord(1:4).eq.'FRGM') Go To 8025
      If (KWord(1:4).eq.'GEO ') Go To 8024
      If (KWord(1:4).eq.'GEOE') Go To 8020
      If (KWord(1:4).eq.'GIAO') Go To 9020
      If (KWord(1:4).eq.'GRID') Go To 9773
      If (KWord(1:4).eq.'GROM') Go To 8034
      If (KWord(1:4).eq.'GROU') Go To 6010
      If (KWord(1:4).eq.'HIGH') Go To 9096
      If (KWord(1:4).eq.'HYPE') Go To 8016
      If (KWord(1:4).eq.'ISOT') Go To 7654
      If (KWord(1:4).eq.'JMAX') Go To 971
      If (KWord(1:4).eq.'KHAC') Go To 8003
      If (KWord(1:4).eq.'LDF ') Go To 35
      If (KWord(1:4).eq.'LDF1') Go To 35
      If (KWord(1:4).eq.'LDF2') Go To 36
      If (KWord(1:4).eq.'RLOC') Go To 658
      If (KWord(1:4).eq.'LINK') Go To 8036
      If (KWord(1:4).eq.'LOCA') Go To 35
      If (KWord(1:4).eq.'LOW ') Go To 9094
      If (KWord(1:4).eq.'MEDI') Go To 9095
      If (KWord(1:4).eq.'MGAU') Go To 8009
      If (KWord(1:4).eq.'MOLC') Go To 957
      If (KWord(1:4).eq.'MOLE') Go To 960
      If (KWord(1:4).eq.'MOLP') Go To 959
      If (KWord(1:4).eq.'MOVE') Go To 6040
      If (KWord(1:4).eq.'MULT') Go To 972
      If (KWord(1:4).eq.'NACC') Go To 8030
      If (KWord(1:4).eq.'NEMO') Go To 800
      If (KWord(1:4).eq.'NGEX') Go To 501
      If (KWord(1:4).eq.'NOAL') Go To 7013
      If (KWord(1:4).eq.'NOAM') Go To 8007
      If (KWord(1:4).eq.'NOCD') Go To 9084
      If (KWord(1:4).eq.'NOCH') Go To 9089
      If (KWord(1:4).eq.'NODE') Go To 7070
      If (KWord(1:4).eq.'NODK') Go To 989
      If (KWord(1:4).eq.'NOGU') Go To 9100
      If (KWord(1:4).eq.'NOMO') Go To 6050
      If (KWord(1:4).eq.'NOON') Go To 8023
      If (KWord(1:4).eq.'NOPA') Go To 9910
      If (KWord(1:4).eq.'NOTA') Go To 980
      If (KWord(1:4).eq.'NOUN') Go To 46
      If (KWord(1:4).eq.'NUME') Go To 8031
      If (KWord(1:4).eq.'OLDD') Go To 8121
      If (KWord(1:4).eq.'OLDZ') Go To 8021
      If (KWord(1:4).eq.'OMQI') Go To 999
      If (KWord(1:4).eq.'ONEO') Go To 990
      If (KWord(1:4).eq.'OPTH') Go To 8022
      If (KWord(1:4).eq.'OPTO') Go To 940
      If (KWord(1:4).eq.'ORBA') Go To 913
      If (KWord(1:4).eq.'ORBC') Go To 906
      If (KWord(1:4).eq.'ORIG') Go To 8015
      If (KWord(1:4).eq.'OVER') Go To 41
      If (KWord(1:4).eq.'PAMF') Go To 8060
      If (KWord(1:4).eq.'PART') Go To 9763
      If (KWord(1:4).eq.'PKTH') Go To 9940
      If (KWord(1:4).eq.'MXTC') Go To 9023
      If (KWord(1:4).eq.'PRIN') Go To 930
c     If (KWord(1:1).eq.'R' .and.
c    &    (KWord(2:2).ge.'0' .and.
c    &     KWord(2:2).le.'9') .and.
c    &     (KWord(3:3).ge.'0' .and.
c    &      KWord(3:3).le.'9') .and.
c    &      (KWord(4:4).eq.'O' .or.
c    &       KWord(4:4).eq.'E' .or.
c    &       KWord(4:4).eq.'S' .or.
c    &       KWord(4:4).eq.'M' .or.
c    &       KWord(4:4).eq.'C') ) Go To 657
      If (KWord(1:4).eq.'RA0F') Go To 9012
      If (KWord(1:4).eq.'RA0H') Go To 9011
      If (KWord(1:4).eq.'RADI') Go To 909
      If (KWord(1:4).eq.'RAIH') Go To 9013
      If (KWord(1:4).eq.'RBSS') Go To 9015
      If (KWord(1:4).eq.'RELA') Go To 657
      If (KWord(1:4).eq.'RELI') Go To 962
      If (KWord(1:4).eq.'RESC') Go To 978
      If (KWord(1:4).eq.'RF-I') Go To 9970
      If (KWord(1:4).eq.'RLDF') Go To 8012
      If (KWord(1:4).eq.'RIC ') Go To 9097
      If (KWord(1:4).eq.'RICD') Go To 9080
      If (KWord(1:4).eq.'RIJ ') Go To 9098
      If (KWord(1:4).eq.'RIJK') Go To 9099
      If (KWord(1:4).eq.'RMAT') Go To 880
      If (KWord(1:4).eq.'RMBP') Go To 886
      If (KWord(1:4).eq.'RMDI') Go To 884
      If (KWord(1:4).eq.'RMEA') Go To 881
      If (KWord(1:4).eq.'RMER') Go To 882
      If (KWord(1:4).eq.'RMEQ') Go To 885
      If (KWord(1:4).eq.'RMQC') Go To 883
      If (KWord(1:4).eq.'ROT ') Go To 8027
      If (KWord(1:4).eq.'RP-C') Go To 9093
      If (KWord(1:4).eq.'RPQM') Go To 8008
      If (KWord(1:4).eq.'RTRN') Go To 950
      If (KWord(1:4).eq.'RX2C') Go To 9014
      If (KWord(1:4).eq.'SADD') Go To 9081
      If (KWord(1:4).eq.'SCAL') Go To 8018
      If (KWord(1:4).eq.'SHAK') Go To 8050
      If (KWord(1:4).eq.'SDEL') Go To 7071
      If (KWord(1:4).eq.'SDIP') Go To 992
      If (KWord(1:4).eq.'SHAC') Go To 8002
      If (KWord(1:4).eq.'SKIP') Go To 9950
      If (KWord(1:4).eq.'SLIM') Go To 8005
      If (KWord(1:4).eq.'SPAN') Go To 890
      If (KWord(1:4).eq.'SPRE') Go To 889
      If (KWord(1:4).eq.'STDO') Go To 9930
      If (KWord(1:4).eq.'SYMM') Go To 900
      If (KWord(1:4).eq.'SYMT') Go To 6060
      If (KWord(1:4).eq.'TARG') Go To 37
      If (KWord(1:4).eq.'TDEL') Go To 7072
      If (KWord(1:4).eq.'TEST') Go To 991
      If (KWord(1:4).eq.'THRC') Go To 9021
      If (KWord(1:4).eq.'THRE') Go To 941
      If (KWord(1:4).eq.'THRL') Go To 37
      If (KWord(1:4).eq.'THRS') Go To 907
      If (KWord(1:4).eq.'TINK') Go To 8014
      If (KWord(1:4).eq.'TITL') Go To 910
      If (KWord(1:4).eq.'TRAN') Go To 8026
      If (KWord(1:4).eq.'UNCO') Go To 43
      If (KWord(1:4).eq.'UNIQ') Go To 45
      If (KWord(1:4).eq.'UNNO') Go To 908
      If (KWord(1:4).eq.'UPON') Go To 1001
      If (KWord(1:4).eq.'VECT') Go To 905
      If (KWord(1:4).eq.'VART') Go To 8032
      If (KWord(1:4).eq.'VARR') Go To 8033
      If (KWord(1:4).eq.'VERB') Go To 9122
      If (KWord(1:4).eq.'VERI') Go To 40
      If (KWord(1:4).eq.'WEIG') Go To 7014
      If (KWord(1:4).eq.'WELL') Go To 986
      If (KWord(1:4).eq.'WRUC') Go To 44
      If (KWord(1:4).eq.'XBAS') Go To 1924
      If (KWord(1:4).eq.'XFIE') Go To 975
      If (KWord(1:4).eq.'XRIC') Go To 9085
      If (KWord(1:4).eq.'XYZ ') Go To 1917
      If (KWord(1:4).eq.'ZCON') Go To 8017
      If (KWord(1:4).eq.'ZMAT') Go To 1920
      If (KWord(1:4).eq.'ZONL') Go To 8028
*
      If (KWord(1:4).eq.'END ') Go To 997
      If (lTtl) Go To 911
*
      If (Basis_test) Then
*
*        So the Basis keyword was in the native format.
*        We have to back step until we find the command line!
*
         Backspace(LuRd)
         Backspace(LuRd)
         Read(LuRd,'(A)') Key
         Call UpCase(Key)
         Do While(Index(Key(1:4),'BASI').eq.0)
              Backspace(LuRd)
              Backspace(LuRd)
              Read(LuRd,'(A)') Key
              Call UpCase(Key)
         End Do
         Basis_test=.False.
         nDone=0
         Go To 9201
      End If
      iChrct=Len(KWord)
      Last=iCLast(KWord,iChrct)
      Write (LuWr,*)
      Call WarningMessage(2,KWord(1:Last)//' is not a keyword!,'//
     &               ' Error in keyword.')
      Call Quit_OnUserError()
*
c977  Call WarningMessage(2,' Premature end of input file.')
c     Call Quit_OnUserError()
*
#ifdef _FDE_
*                                                                      *
****** EMBE ************************************************************
*                                                                      *
*     Read in the name of a file for an embedding potential on a grid
*
 666  embPot=.true.
      Call Put_iScalar('embpot',1)
 667  KWord = Get_Ln(LuRd)
      if(KWord(1:4).eq."ENDE") then
       ! Sanity check
       if (embPotInBasis.and.(.not.outGridPathGiven).and.(embWriteDens
     &      .or.embWriteESP.or.embWriteGrad.or.embWriteHess)) then
        Call WarningMessage(2,' No grid to write output in embedding.')
        Call Quit_OnUserError()
       end if
       ! Write the embpot runfile
       Call EmbPotWrRun
       Go To 998
      end if
      ! Switch whether the potential is given in basis set
      ! representation instead of grid representation
      if (KWord(1:4).eq."BASI") then
       embPotInBasis = .true.
       write(LuWr,*) "Set embPotInBasis to ", embPotInBasis
       Go To 667
      end if
      ! Get the EMBInfile path containing an embedding pot. on a grid
      if(KWord(1:4).eq."EMBI") then
       KWord = Get_Ln(LuRd)
       Call Get_S(1,embPotPath,1)
       Go To 667
      end if
      ! If the output grid is different from the input grid a path to
      ! a grid file is specified here
      if(KWord(1:4).eq."OUTG") then
       KWord = Get_Ln(LuRd)
       Call Get_S(1,outGridPath,1)
       outGridPathGiven=.true.
       Go To 667
      end if
      ! If given, the final density is written out on a grid, the path
      ! to the file where to write it to must be given as well.
      if(KWord(1:4).eq."WRDE") then
       embWriteDens=.true.
       KWord = Get_Ln(LuRd)
       Call Get_S(1,embOutDensPath,1)
       Go To 667
      end if
      ! If given, the final electrostatic potential is written out on
      ! a grid, the path to the file where to write it to must be given
      ! as well.
      if(KWord(1:4).eq."WREP") then
       KWord = Get_Ln(LuRd)
       Call Get_S(1,embOutEspPath,1)
       embWriteESP=.true.
       Go To 667
      end if
      ! If given, the density gradient is written out on a grid, the
      ! path to the file where to write it to must be given as well.
      if(KWord(1:4).eq."WRGR") then
       KWord = Get_Ln(LuRd)
       Call Get_S(1,embOutGradPath,1)
       embWriteGrad=.true.
       Go To 667
      end if
      ! If given, the Hessian of the density is written out on a grid,
      ! the path to the file where to write it to must be given as
      ! well.
      if(KWord(1:4).eq."WRHE") then
       embWriteHess=.true.
       KWord = Get_Ln(LuRd)
       Call Get_S(1,embOutHessPath,1)
       Go To 667
      end if
      ! Should not be reached.
       Call WarningMessage(2,
     &             'Error in input of EMBEdding block!')
       Call Quit_OnUserError()
#endif
*                                                                      *
****** SYMM ************************************************************
*                                                                      *
*     Read distinct symmetry operators apart for the unity operator
*
 900  SymmSet=.true.
      KWord = Get_Ln(LuRd)
      GWInput=.True.
      If (.not.DoneCoord.and.iCoord.ne.0) Then
       Call WarningMessage(2,
     &             'SYMMETRY keyword is not compatible with COORD')
       Call Quit_OnUserError()
      End If
      If (Run_Mode.eq.S_Mode) Then
         Call WarningMessage(2,'Seward input error!')
         Write (LuWr,'(A,A,A)') 'The command : "',Previous_Command,
     &                          '" is not allowed'
     &                        //' in Seward input when Gateway is used!'
         Write (LuWr,*)
         Call Quit_OnUserError()
      End If
      Call UpCase(KWord)
      iChrct=Len(KWord)
 901  Last=iCLast(KWord,iChrct)
      iFrst=iCFrst(KWord,iChrct)
      If (iFrst.le.Last) Then
         nOper = nOper + 1
         Oper(nOper)(1:1) = KWord(iFrst:iFrst)
         KWord(iFrst:iFrst) = ' '
         iFrst = iFrst + 1
         If (KWord(iFrst:iFrst).eq.' ') Go To 901
         Oper(nOper)(2:2) = KWord(iFrst:iFrst)
         KWord(iFrst:iFrst) = ' '
         iFrst = iFrst + 1
         If (KWord(iFrst:iFrst).eq.' ') Go To 901
         Oper(nOper)(3:3) = KWord(iFrst:iFrst)
         KWord(iFrst:iFrst) = ' '
         Go to 901
      Else
         Go To 998
      End If
*                                                                      *
****** FILE ************************************************************
*                                                                      *
*     Specify filename for input orbitals
*
 904  Line = Get_Ln(LuRd)
      Call FileOrb(Line,SW_FileOrb)
      Go To 998
*                                                                      *
****** VECT ************************************************************
*                                                                      *
*     Change mode of Seward to property calculations.
*
 905  Prprt = .True.
      Go To 998
*                                                                      *
****** ORBC ************************************************************
*                                                                      *
*     Request property output with explicit listing of orbital
*     contributions.
*
 906  Short = .False.
      Go To 998
*                                                                      *
****** THRS ************************************************************
*                                                                      *
*     Change default for non-zero occupation numbers
*
 907  KWord = Get_Ln(LuRd)
      Call Get_F1(1,Thrs)
      Thrs = Abs(Thrs)
      Go To 998
*                                                                      *
****** UNNO ************************************************************
*                                                                      *
*     Change default to unnormalized integrals
*
 908  UnNorm = .True.
      Go To 998
*                                                                      *
****** RADI ************************************************************
*                                                                      *
*     Integral cutoff to be based on radial overlap integrals
*
 909  lSchw = .False.
      Go To 998
*                                                                      *
****** TITL ************************************************************
*                                                                      *
*     Read the Title card
*
 910  Key = Get_Ln(LuRd)
 911  Continue
      GWInput = Run_Mode.eq.G_Mode
      lTtl = .True.
      nTtl = nTtl + 1
      If (nTtl.gt.10) Then
         Call WarningMessage(2,' Too many title cards')
         Call Quit_OnUserError()
      End If
      i1=iCFrst(Key,80)
      i2=iCLast(Key,80)
      nc = 80-(i2-i1+1)
      nc2=nc/2
      Title(nTtl)=''
      Title(nTtl)(nc2+1:nc2+i2-i1+1)=Key(i1:i2)
      Go To 9988
*                                                                      *
****** ECPS **** or ****** AUXS ****************************************
*                                                                      *
*     Allow printing of ECP data
*
 912  nPrint(2)=Max(10,nPrint(2))
      GWInput=.True.
      Go To 998
*                                                                      *
****** BSSH ************************************************************
*                                                                      *
*     Allow printing of basis set data
*
 9121 nPrint(2)=Max(6,nPrint(2))
      GWInput=.True.
      Go To 998
*                                                                      *
****** VERB ************************************************************
*                                                                      *
*     Verbose printing
*
 9122 nPrint(2)=Max(10,nPrint(2))
      nPrint(117)=6
      nPrint(80)=6
      nPrint( 1)=6
      GWInput = Run_Mode.eq.G_Mode
      Go To 998
*                                                                      *
****** ORBA ************************************************************
*                                                                      *
*     Request property output with explicit listing of properties of
*     all orbitals, including all unoccupied (ignoring THRS), and not
*     weighted by occupation numbers. (S.S.Dong, 2018)
*
 913  ifallorb = .True.
      Go To 998
*                                                                      *
****** ZMAT ************************************************************
*                                                                      *
*     Read Basis Sets & Coordinates in Z-Matrix format
*
1920  Continue
      if(isxbas.eq.0) Call Quit_OnUserError()
      Call ZMatrixConverter(LuRd,LuWr,mxAtom,STDINP,lSTDINP,
     &   iglobal,nxbas,xb_label,xb_bas,iErr)
      If (iErr.ne.0) Call Quit_OnUserError()
      GWInput=.True.
      Call StdSewInput(LuRd,ifnr,mdc,iShll,BasisTypes,
     &                 STDINP,lSTDINP,iErr)
      If (iErr.ne.0) Call Quit_OnUserError()
      Go To 998
*                                                                      *
****** XBAS ************************************************************
*                                                                      *
1924  Continue
      call read_xbas(LuRd,iglobal,nxbas,xb_label,xb_bas,ierr)
      GWInput=.True.
      isxbas=1
      if(ierr.eq.1) Call Quit_OnUserError()
      goto 998
*                                                                      *
****** XYZ  ************************************************************
*                                                                      *
1917  Continue
      if(isxbas.eq.0) Call Quit_OnUserError()
      Call XMatrixConverter(LuRd,LuWr,mxAtom,STDINP,lSTDINP,
     &   iglobal,nxbas,xb_label,xb_bas,iErr)
      If (iErr.ne.0) Call Quit_OnUserError()
      GWInput=.True.
      Call StdSewInput(LuRd,ifnr,mdc,iShll,BasisTypes,
     &                 STDINP,lSTDINP,iErr)
      If (iErr.ne.0) Call Quit_OnUserError()
*      If (SymmSet) Then
*         Call WarningMessage(2,
*     &                 'SYMMETRY keyword is not compatible with XYZ')
*         Call Quit_OnUserError()
*      End If
      If (GroupSet) Then
         Call WarningMessage(2,
     &                 'GROUP keyword is not compatible with XYZ')
         Call Quit_OnUserError()
      End If
      Go To 998

*                                                                      *
****** COOR ************************************************************
*                                                                      *
*     Read Basis Sets & Coordinates in xyz format
*
6000  Continue
      If (SymmSet) Then
         Call WarningMessage(2,
     &                 'SYMMETRY keyword is not compatible with COORD')
         Call Quit_OnUserError()
      End If
      If (iOpt_XYZ.eq.0) Then
         Call WarningMessage(2,
     &            'COORD keyword is not compatible with non-XYZ format')
         Call Quit_OnUserError()
      End If
c      If (CoordSet) Then
c         Call WarningMessage(2,'There is more than one COORD keyword!')
c*        Should this be an error and stop? Only if no HYPER or GEO?
c      End If
      if(isxbas.eq.1) Call Quit_OnUserError()
      iCoord=iCoord+1

      CoordSet=.true.
      GWInput=.True.
#ifdef _HAVE_EXTRA_
      Call XYZread(LuRd,ForceZMAT,nCoord, iErr)
      If (iErr.ne.0) Call Quit_OnUserError()
      Call XYZcollect(iCoord,nCoord,OrigTrans,OrigRot,nFragment)
#else
      Call Read_XYZ(LuRd,OrigRot,OrigTrans)
#endif
      Go To 998
*                                                                      *
****** GROUP ***********************************************************
*                                                                      *
*     Read information for a group
*
6010  Continue
      If (SymmSet) Then
         Call WarningMessage(2,
     &                 'SYMMETRY keyword is not compatible with GROUP')
         Call Quit_OnUserError()
      End If
      If (.not.CoordSet) Then
         Call WarningMessage(2,'COORD keyword is not found')
         Call Quit_OnUserError()
      End If
      KeepGroup=Get_Ln(LuRd)
c Simplistic validity check for value
      temp1=KeepGroup
      Call UpCase(temp1)
      If (temp1(1:4).eq.'FULL') Goto 6015
      If (temp1(1:1).eq.'E') Goto 6015
      If (temp1(1:2).eq.'C1') Goto 6015
      If (temp1(1:5).eq.'NOSYM') Goto 6015
      Do i=1,StrnLn(temp1)
        If (temp1(i:i).ne.'X'.and.
     &      temp1(i:i).ne.'Y'.and.
     &      temp1(i:i).ne.'Z'.and.
     &      temp1(i:i).ne.' ') Then
          Call WarningMessage(2,
     &    'Illegal symmetry group or operator: '//temp1(:StrnLn(temp1)))
          Call Quit_OnUserError()
        End If
      End Do
6015  Continue
      GroupSet=.true.
      GWInput=.True.
      DoneCoord=.True.
      goto 998
*                                                                      *
****** BSSE ************************************************************
*                                                                      *
*     Read information for BSSE
*
6020  Continue
      if(.not.CoordSet) Then
         Call WarningMessage(2,'COORD keyword is not found')
         Call Quit_OnUserError()
      End If
      KWord=Get_Ln(LuRd)
      read(Kword,*,end=6666,err=6666) iBSSE
      GWInput=.True.
      goto 998
*                                                                      *
****** MOVE ************************************************************
*                                                                      *
*     allow to MOVE coordinates
*
6040  Continue
      if(.not.CoordSet) Then
         Call WarningMessage(2,'COORD keyword is not found')
         Call Quit_OnUserError()
      End If
#ifdef _HAVE_EXTRA_
      isHold=0
#endif
      GWInput=.True.
      goto 998
*                                                                      *
****** NOMOVE **********************************************************
*                                                                      *
*     Do NOT allow to MOVE coordinates
*
6050  Continue
      if(.not.CoordSet) Then
         Call WarningMessage(2,'COORD keyword is not found')
         Call Quit_OnUserError()
      End If
#ifdef _HAVE_EXTRA_
      isHold=1
#endif
      GWInput=.True.
      goto 998
*                                                                      *
****** SYMT ************************************************************
*                                                                      *
*     Threshold for findsym
*
6060  Continue
      if(.not.CoordSet) Then
         Call WarningMessage(2,'COORD keyword is not found')
         Call Quit_OnUserError()
      End If
      KWord=Get_Ln(LuRd)
      read(Kword,*,end=6666, err=6666) SymThr
      GWInput=.True.
      goto 998
*                                                                      *
****** NODE ************************************************************
*                                                                      *
*     Set global delete parameters to disable orbital deleting
*     Keywords: NODElete at 7070
*               SDELete  at 7071
*               TDELete  at 7072
*
7070  continue
      Call Put_dScalar('S delete thr',0.0d0)
      Call Put_dScalar('T delete thr',1.0d15)
      goto 998
7071  continue
      KWord = Get_Ln(LuRd)
      Call Get_F1(1,sDel)
      Call Put_dScalar('S delete thr',sDel)
      goto 998
7072  continue
      KWord = Get_Ln(LuRd)
      Call Get_F1(1,tDel)
      Call Put_dScalar('T delete thr',tDel)
      goto 998
7700  continue
      KWord = Get_Ln(LuRd)
      ExtBasDir=KWord
      GWInput=.True.
      goto 998
*                                                                      *
****** BASI ************************************************************
*                                                                      *
*     Read information for a basis set
*
 920  continue
*
*     Check if the format is old or new style. Damn the person who used
*     the same keword for two different styles of input and making the
*     input require a specific order of the keyword. Comrade 55?
*
      Basis_Test=.True.
*
      GWInput=.True.
      Key = Get_Ln(LuRd)
      BSLbl = Key(1:80)
      If (BasisSet) Then
         KeepBasis=KeepBasis(1:index(KeepBasis,' '))//','//BSLbl
      Else
         KeepBasis=BSLbl
         BasisSet=.True.
      Endif
      temp1=KeepBasis
      Call UpCase(temp1)
*     If (INDEX(temp1,'INLINE').ne.0) then
*        Write(LuWr,*)
*    &        'XYZ input and Inline basis set are not compatible'
*        Write(LuWr,*)
*    &        'Consult the manual how to change inline basis set'
*        Write(LuWr,*) ' into basis set library'
*        Call Quit_OnUserError()
*     End If
      iOpt_XYZ=1
      Goto 998
*
 9201 Continue
      iOpt_XYZ=0
      GWInput=.True.
      nCnttp = nCnttp + 1
      If (Run_Mode.eq.S_Mode) Then
         Call WarningMessage(2,'Seward input error!')
         Write (LuWr,*) 'The command : "',Previous_Command,
     &               '" is not allowed'
     &              //' in the Seward input when the Gateway is used!'
         Write (LuWr,*)
         Call Quit_OnUserError()
      End If
      If (nCnttp.gt.Mxdbsc) Then
         Call WarningMessage(2,' Increase Mxdbsc')
         Call Quit_OnUserError()
      End If
      NoZMAT=.True.
*
*     Read the basis set label
*
      Key = Get_Ln(LuRd)
      BSLbl = Key(1:80)
*
*     If dummy atom point at the ANO-RCC set where the specification is.
*
      Call UpCase(BSLbl)
      iDummy_basis=0
      Call ICopy(4,BasisTypes,1,BasisTypes_save,1)
      If (BSLbl(1:2).eq.'X.'.and.Index(BSLbl,'INLINE').eq.0.and.
     &    Index(BSLbl,'RYDBERG').eq.0) Then
         BSLbl='X.ANO-RCC.'
         Do i=11,80
           BSLbl(i:i)='.'
         End do
         iDummy_basis=1
      End If
*     Call UpCase(BSLbl)
      LenBSL=Len(BSLbl)
      Last=iCLast(BSLbl,LenBSL)
      GWInput=.True.
      Indx=Index(BSLbl,'/')
      If (Indx.eq.0) Then
         Fname=BasLib
         Indx = Last+1
         dbsc(nCnttp)%Bsl=BSLbl
      Else
         Fname= BSLbl(Indx+2:Last)
         If (Fname.eq.' ') Then
            Call WarningMessage(2,
     &                     ' No basis set library specified for'
     &                   //'BSLbl='//BSLbl(1:Indx-1)//' Fname='//Fname)
            Call Quit_OnUserError()
         End If
 1919    If (Fname(1:1).eq.' ') Then
            Fname(1:79)=Fname(2:80)
            Fname(80:80) = ' '
            Go To 1919
         End If
         dbsc(nCnttp)%Bsl=BSLbl(1:Indx-1)
      End If
*
      n=INDEX(dbsc(nCnttp)%Bsl,' ')
      If (n.eq.0) n=81
      Do i=n,80
        dbsc(nCnttp)%Bsl(i:i)='.'
      End Do
*
      If ((Show.and.nPrint(2).ge.6) .or.
     &    Write_BasLib) Then
         Write (LuWr,*)
         Write (LuWr,*)
         Write(LuWr,'(1X,A,I5,A,A)')
     &           'Basis Set ',nCnttp,' Label: ', BSLbl(1:Indx-1)
         Write(LuWr,'(1X,A,A)') 'Basis set is read from library:',
     *         Fname(1:index(Fname,' '))
      End if
*
      jShll = iShll
      dbsc(nCnttp)%Bsl_old=dbsc(nCnttp)%Bsl
      dbsc(nCnttp)%mdci=mdc
      Call GetBS(Fname,dbsc(nCnttp)%Bsl,iShll,Ref,UnNorm,LuRd,
     &           BasisTypes,STDINP,lSTDINP,.False.,Expert,ExtBasDir)
*
      Do_FckInt = Do_FckInt .and. dbsc(nCnttp)%FOp .and.
     &            dbsc(nCnttp)%AtmNr.le.96
#ifdef _DEMO_
      Do_GuessOrb = .False.
#else
      Do_GuessOrb = Do_GuessOrb .and. dbsc(nCnttp)%AtmNr.le.96
#endif
*
      If (iDummy_Basis.eq.1) Call ICopy(4,BasisTypes_Save,1,
     &                                    BasisTypes,1)
      If (itype.eq.0) Then
         If (BasisTypes(3).eq.1 .or. BasisTypes(3).eq.2 .or.
     &       BasisTypes(3).eq.14)
     &       iType=BasisTypes(3)
      Else
         If (BasisTypes(3).eq.1 .or. BasisTypes(3).eq.2 .or.
     &       BasisTypes(3).eq.14) Then
            If (BasisTypes(3).ne.iType) Then
               imix=1
               BasisTypes(3)=-1
            End If
            iType=BasisTypes(3)
         End If
      End If
      If (itype.eq.1) ifnr=1
      If (itype.eq.2 .or. itype.eq.14) ifnr=0
*
      If (ign.eq.0) Then
         ign=BasisTypes(4)
      Else If (Abs(BasisTypes(4)).ne.Abs(ign)) Then
         If (nmwarn) Then
            Call WarningMessage(1,
     &        'SEWARD found basis sets of mixed nuclear charge model. '
     &      //'The most advanced one will be used.')
         End If
         nmwarn=.False.
         ign=Max(ign,BasisTypes(4))
         BasisTypes(4)=ign
      End If
*
      If (Show.and.nPrint(2).ge.6 .and.
     &   Ref(1).ne.'' .and. Ref(2).ne.'') Then
         Write (LuWr,'(1x,a)')  'Basis Set Reference(s):'
         If (Ref(1).ne.'') Write (LuWr,'(5x,a)') Trim(Ref(1))
         If (Ref(2).ne.'') Write (LuWr,'(5x,a)') Trim(Ref(2))
         Write (LuWr,*)
         Write (LuWr,*)
      End If
      dbsc(nCnttp)%ECP=(dbsc(nCnttp)%nPP
     &                 +dbsc(nCnttp)%nPrj
     &                 +dbsc(nCnttp)%nSRO
     &                 +dbsc(nCnttp)%nSOC
     &                 +dbsc(nCnttp)%nM1
     &                 +dbsc(nCnttp)%nM2) .NE. 0
*
      lAng=Max(dbsc(nCnttp)%nVal,
     &         dbsc(nCnttp)%nSRO,
     &         dbsc(nCnttp)%nPrj)-1
      S%iAngMx=Max(S%iAngMx,lAng)
*     No transformation needed for s and p shells
      Shells(jShll+1)%Transf=.False.
      Shells(jShll+1)%Prjct =.False.
      Shells(jShll+2)%Transf=.False.
      Shells(jShll+2)%Prjct =.False.
      dbsc(nCnttp)%nShells = dbsc(nCnttp)%nVal
     &                     + dbsc(nCnttp)%nPrj
     &                     + dbsc(nCnttp)%nSRO
     &                     + dbsc(nCnttp)%nSOC
     &                     + dbsc(nCnttp)%nPP
      nCnt = 0
      If (dbsc(nCnttp)%Aux) Then
         Do iSh = jShll+1, iShll
            Shells(iSh)%Aux=.True.
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Set Cartesian functions if specified by the basis type
*
      If (BasisTypes(1).eq.9) Then
         Do iSh = jShll+3, iShll
            Shells(iSh)%Transf=.False.
            Shells(iSh)%Prjct =.False.
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Automatic onset of muonic charge if the basis type is muonic.
*     This will also automatically activate finite nuclear mass
*     correction.
*
      KWord=''
      KWord(1:Indx-1)=BSLbl(1:Indx-1)
      Call UpCase(KWord)
      If (INDEX(KWord,'MUONIC').ne.0) Then
         dbsc(nCnttp)%fMass=
     &    CONST_MUON_MASS_IN_SI_ / CONST_ELECTRON_MASS_IN_SI_
         FNMC=.True.
         tDel=1.0D50
         Call Put_dScalar('T delete thr',tDel)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Update BasisTypes
*
      Do i=1,4
         If (BasisTypes_save(i).eq.0) Cycle
         If (BasisTypes(i).ne.BasisTypes_save(i)) BasisTypes(i)=-1
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*
 777  KWord = Get_Ln(LuRd)
      Call UpCase(KWord)
      Call LeftAd(KWord)
      If (KWord(1:4).eq.'PSEU') Then
         dbsc(nCnttp)%pChrg=.True.
         dbsc(nCnttp)%Fixed=.True.
         Go To 777
      End If
      If (KWord(1:4).eq.'ACDT') Then
         KWord = Get_Ln(LuRd)
         Call Get_F1(1,dbsc(nCnttp)%aCD_Thr)
         Go To 777
      End If
      If (KWord(1:4).eq.'MUON') Then
         dbsc(nCnttp)%fMass=
     &    CONST_MUON_MASS_IN_SI_ / CONST_ELECTRON_MASS_IN_SI_
         Go To 777
      End If
      If (KWord(1:4).eq.'NUCL') Then
         KWord = Get_Ln(LuRd)
         Call Get_F1(1,dbsc(nCnttp)%ExpNuc)
         Go To 777
      End If
      If (KWord(1:4).eq.'FIXE') Then
         dbsc(nCnttp)%Fixed=.True.
         Go To 777
      End If
      If (KWord(1:4).eq.'SPHE') Then
         If (Index(KWord,'ALL').ne.0) Then
            Do iSh = jShll+3, iShll
               Shells(iSh)%Transf=.True.
               Shells(iSh)%Prjct =.True.
            End Do
            Go To 777
         End If
         ist = index(KWord,' ')
         iAng = 2
         Do iSh = jShll+3, iShll
            If (Index(KWord(ist:80),AngTyp(iAng)).ne.0) Then
               Shells(iSh)%Transf = .True.
               Shells(iSh)%Prjct  = .True.
            End If
            iAng = iAng + 1
         End Do
         Go To 777
      End If
      If (KWord(1:4).eq.'CART') Then
         If (Index(KWord,'ALL').ne.0) Then
            Do iSh = jShll+1, iShll
               Shells(iSh)%Transf=.False.
               Shells(iSh)%Prjct =.False.
            End Do
            Go To 777
         End If
         ist = index(KWord,' ')
         iAng = 0
         Do iSh = jShll+1, iShll
            If (Index(KWord(ist:80),AngTyp(iAng)).ne.0) Then
               Shells(iSh)%Transf = .False.
               Shells(iSh)%Prjct  = .False.
            End If
            iAng = iAng + 1
         End Do
         Go To 777
      End If
      If (KWord(1:4).eq.'CONT') Then
         If (Index(KWord,'ALL').ne.0) Then
            Do iSh = jShll+1, iShll
               Shells(iSh)%Prjct  = .False.
            End Do
            Go To 777
         End If
         ist = index(KWord,' ')
         iAng = 0
         Do iSh = jShll+1, iShll
            If (Index(KWord(ist:80),AngTyp(iAng)).ne.0)
     &          Shells(iSh)%Prjct  = .False.
            iAng = iAng + 1
         End Do
         Go To 777
      End If
      If (KWord(1:4).eq.'CHAR') Then
         KWord = Get_Ln(LuRd)
         Call UpCase(KWord)
         Call Get_F1(1,dbsc(nCnttp)%Charge)
         ist = index(KWord,' ')
         If (dbsc(nCnttp)%IsMM.ne.0) Then
            Call WarningMessage(1,
     &         ' Found a charge associated with a MM atom. Ignore it')
            dbsc(nCnttp)%Charge = Zero
         End If
         Go To 777
      End If
      If (KWord(1:4).eq.'FRAG') Then
         dbsc(nCnttp)%pChrg=.True.
         dbsc(nCnttp)%Fixed=.True.
         lFAIEMP=.True.
         Go To 777
      End If
      If (KWord(1:4).eq.'END ') Then
         If (nCnt.eq.0) Then
            Call WarningMessage(2,' Input error, no center specified!')
            Call Quit_OnUserError()
         End If
         dbsc(nCnttp)%nCntr = nCnt
         mdc = mdc + nCnt
!        Now allocate the array for the coordinates and copy them over.
!        Call Allocate(dbsc(nCnttp)%Coor(1:3,1:nCnt)
         Call mma_Allocate(dbsc(nCnttp)%Coor_Hidden,3,nCnt,
     &                     Label='dbsc:C')
         dbsc(nCnttp)%Coor => dbsc(nCnttp)%Coor_Hidden(:,:)
         Call DCopy_(3*nCnt,Buffer,1,dbsc(nCnttp)%Coor,1)
!
         Go To 998
      End If
*
*     Read Coordinates
*
      nCnt = nCnt + 1
      n_dc=max(mdc+nCnt,n_dc)
      If (mdc+nCnt.gt.MxAtom) Then
         Call WarningMessage(2,' RdCtl: Increase MxAtom')
         Write (LuWr,*) '        MxAtom=',MxAtom
         Call Quit_OnUserError()
      End If
      jend=Index(KWord,' ')
      If (jEnd.gt.LENIN+1) Then
         Write (6,*) 'Warning: the label ', KWord(1:jEnd),
     &               ' will be truncated to ',LENIN,' characters!'
      End If
      dc(mdc+nCnt)%LblCnt = KWord(1:Min(LENIN,jend-1))
      dbas=dc(mdc+nCnt)%LblCnt(1:LENIN)
      Call Upcase(dbas)
      If (dbas.eq.'DBAS') Then
         RMat_On=.True.
      End If
      If (mdc+nCnt.gt.1) then
        Call Chk_LblCnt(dc(mdc+nCnt)%LblCnt,mdc+nCnt-1)
      endif
      iOff=1+(nCnt-1)*3
      Call Get_F(2,Buffer(iOff),3)
      If (Index(KWord,'ANGSTROM').ne.0) Then
         Do i = 0, 2
            Buffer(iOff+i) = Buffer(iOff+i)/angstr
         End Do
      End If
*
      If(Cell_l) Then
      nCnt0 = nCnt
      iOff0=iOff
      lthCell = lthCell + 1
      AdCell(lthCell) = mdc+nCnt0  ! the sequence atom No
      ii = 0
      Do n1=-ispread(1),ispread(1)
         Do n2=-ispread(2),ispread(2)
            Do n3=-ispread(3),ispread(3)
               If((n1.Eq.0).And.(n2.Eq.0).And.(n3.Eq.0)) GoTo 110

                  ii = ii + 1

                  If(ii.GE.10000) Then
                     Call WarningMessage(2,' Too many atoms in Seward')
                     Call Quit_OnUserError()
                  Else
                     If(ii.LT.1000) Then
                        CHAR4 = '_'//str(ii)
                     Else
                        CHAR4 = str(ii)
                     End If
                  End If

                  nCnt = nCnt + 1

                  n_dc=max(mdc+nCnt,n_dc)
                  If (mdc+nCnt.gt.MxAtom) Then
                     Call WarningMessage(2,' RdCtl: Increase MxAtom')
                     Write (LuWr,*) '        MxAtom=',MxAtom
                     Call Quit_OnUserError()
                  End If

                  jend=Index(KWord,' ')
                  If (jEnd.gt.5) Then
                     Write (6,*) 'Warning: the label ', KWord(1:jEnd),
     &               ' will be truncated to ',LENIN,' characters!'
                  End If
                  dc(mdc+nCnt)%LblCnt = KWord(1:Min(LENIN,jend-1))//
     &              CHAR4

                  Call Chk_LblCnt(dc(mdc+nCnt)%LblCnt,mdc+nCnt-1)

                  iOff=1+(nCnt-1)*3

*                 Copy old coordinate  first
                  CALL DCOPY_(3,Buffer(iOff0),1,Buffer(iOff),1)
                  CALL DAXPY_(3,DBLE(n1),VCell(1,1),1,Buffer(iOff),1)
                  CALL DAXPY_(3,DBLE(n2),VCell(1,2),1,Buffer(iOff),1)
                  CALL DAXPY_(3,DBLE(n3),VCell(1,3),1,Buffer(iOff),1)
*
  110          Continue
*
            End Do
         End Do
      End Do
      End If
      Go To 777
*                                                                      *
****** PRIN ************************************************************
*                                                                      *
*     Print level
*
 930  KWord = Get_Ln(LuRd)
      GWInput = Run_Mode.eq.G_Mode
      Call Get_I1(1,n)
      Do i = 1, n
         KWord = Get_Ln(LuRd)
         Call Get_I1(1,jRout)
         Call Get_I1(2,iPrint)
         nPrint(jRout)=iPrint
      End Do
      Go To 998
*                                                                      *
****** OPTO ************************************************************
*                                                                      *
*     Reduce the output for optimizations
*
 940  lOPTO = .True.
      Go To 998
*                                                                      *
****** THRE ************************************************************
*                                                                      *
*     Threshold for writing integrals to disk
*
 941  KWord = Get_Ln(LuRd)
      Call Get_F1(1,ThrInt)
      ThrInt = Abs(ThrInt)
      ThrInt_UsrDef = .True.
      Go To 998
*                                                                      *
****** CUTO ************************************************************
*                                                                      *
*     Cutoff for computing primitive integrals [a0|c0]
*
 942  KWord = Get_Ln(LuRd)
      Call Get_F1(1,CutInt)
      CutInt = Abs(CutInt)
      CutInt_UsrDef = .True.
      Go To 998
*                                                                      *
****** RTRN ************************************************************
*                                                                      *
*     Defining max bond distance for bonds, angles and dihedrals
*     Define max number of atoms to list for
*
 950  KWord = Get_Ln(LuRd)
      Call Upcase(KWord)
      Call Get_I1(1,S%Max_Center)
      Call Get_F1(2,rtrnc)
      If (Index(KWord,'ANGSTROM').ne.0)
     &    Rtrnc = Rtrnc/angstr
      GWInput=.True.
      Go To 998
*                                                                      *
****** DIRE ************************************************************
*                                                                      *
*     Force direct calculations & disable two-electron integrals
*
 9770 DirInt = .True.
      Onenly = .True.
      iChk_DC=1
      If ( (iChk_RI+iChk_CH) .gt. 0) Then
         Call WarningMessage(2,
     &           'Direct is incompatible with RI and Cholesky keywords')
         Call Quit_OnUserError()
      End If
      Do_RI=.False.
      Go To 998
*                                                                      *
****** CSPF ************************************************************
*                                                                      *
*     Turn on the use of Condon-Shortley phase factors
*
 9110 CSPF=.True.
      GWInput = Run_Mode.eq.G_Mode
      Go To 998
*                                                                      *
****** EXPE ************************************************************
*                                                                      *
*     Expert mode
*
 9771 Expert = .True.
      GWInput = Run_Mode.eq.G_Mode
      Call WarningMessage(1,
     &   ' EXPERT option is ON!')
      Go To 998
*                                                                      *
****** MOLC or DCRN ****************************************************
*                                                                      *
*     Weight for DCR summation (this is the default for conventional
*     calculations)
*
 957  Continue
 958  MolWgh=0
      MolWgh_UsrDef=.true.
      Go To 998
*                                                                      *
****** MOLP ************************************************************
*                                                                      *
*     Weight for DCR summation modified to MOLPRO format
*
 959  MolWgh=2
      MolWgh_UsrDef=.true.
      Go To 998
*                                                                      *
****** MOLE ************************************************************
*                                                                      *
*     Weight for DCR summation modified to MOLECULE format
*
 960  MolWgh=1
      MolWgh_UsrDef=.true.
      Go To 998
*                                                                      *
****** RELI ************************************************************
*                                                                      *
*     Compute integrals for first order relativistic corrections
*     of the energy, i.e. the mass-velocity integrals and the
*     one-electron Darwin contract term integrals.
*
 962  lRel=.True.
      Go To 998
*                                                                      *
****** JMAX ************************************************************
*                                                                      *
*     Change max j quantum number for the rigid rotor analysis
*
 971  KWord = Get_Ln(LuRd)
      Call Get_I1(1,S%jMax)
      Go To 998
*                                                                      *
****** MULT ************************************************************
*                                                                      *
*     Read order of highest multipole to be computed
*
 972  KWord = Get_Ln(LuRd)
      Call Get_I1(1,S%nMltpl)
      Go To 998
*                                                                      *
****** CENT ************************************************************
*                                                                      *
*     User specified centers of multipole moment operators.
*
 973  KWord = Get_Ln(LuRd)
      Call Get_I1(1,nTemp)
      If (lMltpl) Then
         Call WarningMessage(2,
     &               ' Abend: User specified centers already defined;'
     &             //' Abend: Correct input error!')
         Call Quit_OnUserError()
      Else
         lMltpl = .TRUE.
      End If
*     Allocate temporary memory for user defined centers
      Call mma_allocate(RTmp,3,nTemp,label='RTmp')
      Call mma_allocate(ITmp,nTemp,label='ITmp')
      Do 1502 i = 1, nTemp
         KWord = Get_Ln(LuRd)
         Call Upcase(KWord)
         Call Get_I1(1,iMltpl)
         Call Get_F(2,RTmp(1,i),3)
         If (Index(KWord,'ANGSTROM').ne.0)
     &       Call DScal_(3,One/angstr,RTmp(1,i),1)
         ITmp(i) = iMltpl
 1502 Continue
      Go To 998
*                                                                      *
****** EMPC ************************************************************
*                                                                      *
*     Compute Orbital-Free Embedding integrals from Point Charges
*             specified in XFIEld
*
 974  DoEmPC=.True.
#ifdef _HAVE_EXTRA_
      isHold=1 ! avoid coordinate moving
#endif
      GWInput=.True.
      Go To 998
*                                                                      *
****** XFIE ************************************************************
*                                                                      *
*     User specified external field
*
 975  lXF=.True.
      GWInput=.True.
      KWord = Get_Ln(LuRd)
*     Open external file if the line does not start with an integer
*     Note that the "Err" signal cannot be completely trusted, since
*     a slash (i.e., an absolute path) marks end of input and gives
*     no error
      LuRd_saved=LuRd
      ibla = -1
      Read(KWord,*,Err=9751) ibla
      If (ibla.lt.0) Goto 9751
      Goto 9752
9751  LuRd=1
      Call Get_S(1,filename,1)
9753  Call OpnFl(filename(1:(Index(filename,' ')-1)),LuRd,Exist)
      If (.not.Exist) Then
         Call WarningMessage(2,
     &                'Error! File not found: '//
     &        filename(1:(Index(filename,' ')-1)))
         Call Quit_OnUserError()
      End If
      write(LuWr,*)'Reading external field from file: ',
     &        filename(1:(Index(filename,' ')-1))
      KWord = Get_Ln(LuRd)
9752  Call Get_I1(1,nXF)
      Convert=.False.
      Call Upcase(KWord)
      If (Index(KWord,'ANGSTROM').ne.0) Then
          Convert=.True.
          ix=Index(KWord,'ANGSTROM')
          KWord(ix:ix+7)='        '
      End If
*
      KWord(170:180)='-2 -2 -2 -2'
      Call Put_Ln(KWord)
      Call Get_I1(2,nOrd_XF)
      Call Get_I1(3,iXPolType)
      Call Get_I1(4,nXMolnr)
      Call Get_I1(5,nReadEle)
*
*     Set defaults: ch+dip, no polarisabilities,
*                   exclude only its own multipole,
*                   no element read
*
      if (nOrd_XF  .eq.-2) nOrd_XF  =1
      if (iXPolType.eq.-2) iXPolType=0
      if (nXMolnr  .eq.-2) nXMolnr  =0
      if (nReadEle .eq.-2) nReadEle =0

      If ((nOrd_XF.gt.2).or.(nOrd_XF.lt.-1)) Then
         Call WarningMessage(2,'Error! Illegal value of nOrd_XF')
         Write(LuWr,*)'nOrd_XF= ', nOrd_XF
         Call Quit_OnUserError()
      EndIf
      If ((iXPolType.gt.2).or.(iXPolType.lt.0)) Then
         Call WarningMessage(2,'Error! Illegal value of iXPolType')
         Write(LuWr,*)'iXPolType= ', iXPolType
         Call Quit_OnUserError()
      EndIf
      If ((nXMolnr.gt.100).or.(nXMolnr.lt.0)) Then
         Call WarningMessage(2,'Error! Illegal value of nXMolnr')
         Write(LuWr,*)'nXMolnr= ', nXMolnr
         Call Quit_OnUserError()
      EndIf
      If ((nReadEle.gt.1).or.(nReadEle.lt.0)) Then
         Call WarningMessage(2,'Error! Illegal value of nReadEle')
         Write(LuWr,*)'nReadEle= ', nReadEle
         Call Quit_OnUserError()
      EndIf
*
      nData_XF=3
      Do iOrd_XF = 0, nOrd_XF
         nData_XF = nData_XF +  (iOrd_XF+1)*(iOrd_XF+2)/2
      End Do

      if(iXPolType.gt.0) then
         nData_XF=nData_XF+6
         lRF=.true.  ! Polarisabilities treated as Langevin
      endif
*
      If(iXPolType.eq.1) Then
         nDataRead=nData_XF-5   !Read only one pol value
      Else
         nDataRead=nData_XF
      Endif
*
      Call mma_allocate(XF,nData_XF,nXF,Label='XF')
      Call mma_allocate(XMolnr,nXMolnr,nXF,Label='XMolnr')
      Call mma_allocate(XEle,nXF,Label='XEle')
*
      Call Upcase(KWord)
*
      Do iXF = 1, nXF
         XEle(iXF)=0   ! default: no element spec.
*
*        If reading from external file, use free format to allow
*        long lines of input. On the other hand, comments are
*        not allowed in external files.
*
         If (LuRd.ne.LuRd_saved) then
            Call mma_Allocate(iScratch,nXMolnr+nReadEle,
     &                        Label='iScratch')
            Read(LuRd,*)(iScratch(k),k=1,nXMolnr),
     &                  (iScratch(nXMolnr+k),k=1,nReadEle),
     &           (XF(k,iXF),k=1,nDataRead)
            Do i = 1, nXMolnr
               XMolnr(i,iXF)=iScratch(i)
            End Do
            Do i = 1, nReadEle
               XEle(iXF+(i-1))=iScratch(nXMolnr+i)
            End Do
            Call mma_deallocate(iScratch)
         Else
            KWord = Get_Ln(LuRd)
            KWord(170:180)=' 0.0 0.0 0.0'
            Call Put_Ln(KWord)

            Do i = 1, nXMolnr
               Call Get_I1(i,iTemp)
               XMolnr(i,iXF)=iTemp
            End Do
            Do i = 1, nReadEle
               Call Get_I1(nXMolnr+i,iTemp)
               XEle(iXF+(i-1))=iTemp
            End Do
            Call Get_F(nXMolnr+nReadEle+1,XF(1,iXF),nDataRead)
         EndIf
*
         XF(1:3,iXF) = XF(1:3,iXF) * ScaleFactor
         If (Convert) XF(1:3,iXF) = XF(1:3,iXF) / angstr
*
      End Do
*
*---- Close file and reset LuRd if external file was used
      If(LuRd.ne.LuRd_saved) then
         Close(LuRd)
         LuRd = LuRd_saved
      EndIf
      If (isXfield.eq.1) Then
         goto 9755
      End If
      Go To 998
*                                                                      *
****** DOUG ************************************************************
*                                                                      *
*     Full Douglas-Kroll operator
*
 976  continue
      IRELAE=0
      DKroll=.True.
      GWInput=.True.
      go to 998
*
*     Full Douglas-Kroll (DK1) operator
*
 9001 continue
      IRELAE=1
      DKroll=.True.
      GWInput=.True.
      go to 998
*
*     Full Douglas-Kroll (DK2) operator
*
 9002 continue
      IRELAE=2
      DKroll=.True.
      GWInput=.True.
      go to 998
*
*     Douglas-Kroll (DK3) operator
*
 9003 continue
      IRELAE=3
      DKroll=.True.
      GWInput=.True.
      go to 998
*
*     Full Douglas-Kroll (DK3) operator
*
 9004 continue
      IRELAE=4
      DKroll=.True.
      GWInput=.True.
      go to 998
*
*     Full RESC operator
*
 978  continue
      IRELAE=11
      DKroll=.True.
      GWInput=.True.
      go to 998
*
*     Full ZORA
*
 9011 continue
      IRELAE=21
      DKroll=.True.
      GWInput=.True.
      go to 998
*
*     Full ZORA-FP
*
 9012 continue
      IRELAE=22
      DKroll=.True.
      GWInput=.True.
      go to 998
*
*     Full IORA
*
 9013 continue
      IRELAE=23
      DKroll=.True.
      GWInput=.True.
      go to 998
*
*     Exact decoupling X2C method
*
 9014 continue
      IRELAE=101
      DKroll=.True.
      GWInput=.True.
      go to 998
*
*     Exact decoupling BSS method
*
 9015 continue
      IRELAE=102
      DKroll=.True.
      GWInput=.True.
      go to 998
*
*     Switch on the old DKH routine
*
 8121 continue
      IRFLAG1=1
      GWInput=.True.
      goto 998
*                                                                      *
****** BSSM ************************************************************
*                                                                      *
*     BSS method
 9760 continue
      BSS=.True.
      GWInput=.True.
      go to 976
*                                                                      *
****** AMFI ************************************************************
*                                                                      *
*     AMFI integrals
*
 9761 Continue
      lAMFI=.True.
      GWInput=.True.
      Go To 998
*                                                                      *
****** AMF1 ************************************************************
*                                                                      *
*     AMFI integrals
*
 8761 Continue
      lAMFI=.True.
      GWInput=.True.
      Go To 998
*                                                                      *
****** AMF2 ************************************************************
*                                                                      *
*     AMFI integrals (including 2nd-order)
*
 8762 Continue
      lAMFI=.True.
      GWInput=.True.
      Go To 998
*                                                                      *
****** AMF3 ************************************************************
*                                                                      *
*     AMFI integrals (including 3rd-order)
*
 8763 Continue
      lAMFI=.True.
      GWInput=.True.
      Go To 998
*                                                                      *
****** FAKE ************************************************************
*                                                                      *
*     Fake run : conventional, RI or CD ERIs not computed
*                but some info is set to runfile (e.g. CD thrs)
*                Not the same as ONEOnly !!
*
 9759 Continue
      Fake_ERIs = .true.
      Go To 998
*                                                                      *
****** FINI ************************************************************
*                                                                      *
*     Finite nuclei - Gaussian type
*
 9762 Continue
      Nuclear_Model=Gaussian_Type
      GWInput=.True.
      Go To 998
*                                                                      *
****** MGAU ************************************************************
*                                                                      *
*     Finite nuclei - modified Gaussian type
*
 8009 Continue
      Nuclear_Model=mGaussian_Type
      GWInput=.True.
      Go To 998
*                                                                      *
****** PART ************************************************************
*                                                                      *
*     Show partitioning statistics
*
 9763 Continue
      nPrint(10)=6
      Go To 998
*                                                                      *
****** FPCO ************************************************************
*                                                                      *
*     Force partitioning for contracted functions
*
 9764 Continue
      force_part_c=.True.
      Go To 998
*                                                                      *
****** FPPR ************************************************************
*                                                                      *
*     Force partitioning for primitive functions
*
 9765 Continue
      force_part_p=.True.
      Go To 998
*                                                                      *
****** NOTA ************************************************************
*                                                                      *
*     Do not use tables for the roots and weights of the
*     Rys polynomials.
*
 980  NoTab = .TRUE.
      Go To 998
*                                                                      *
****** WELL ************************************************************
*                                                                      *
*-----Read radius and exponents for spherical well integrals
*     and coefficient
*
 986  KWord = Get_Ln(LuRd)
      GWInput=.True.
      Call Get_I1(1,nWel)
      If (nWel.le.0) Then
*--------Use automatic set up for well integrals
         nWel=3
         Call mma_allocate(Wel_Info,3,nWel,Label='Wel_Info')
         Do iWel = 1, nWel
            Wel_Info(3,iWel)=WellCff(iWel)
            Wel_Info(2,iWel)=WellExp(iWel)
            Wel_Info(1,iWel)=WellRad(iWel)
         End Do
      Else
         Call mma_allocate(Wel_Info,3,nWel,Label='Wel_Info')
         Do iWel = 1, nWel
*---------- Read the Coefficient, Exponent, and Radius
            KWord = Get_Ln(LuRd)
            call Get_F1(1,Wel_Info(3,iWel))
            call Get_F1(2,Wel_Info(2,iWel))
            call Get_F1(3,Wel_Info(1,iWel))
            Call Upcase(KWord)
            If (Index(KWord,'ANGSTROM').ne.0) Then
               Wel_Info(1,iWel)=Wel_Info(1,iWel)/angstr
               Wel_Info(2,iWel)=Wel_Info(2,iWel)*angstr
            End If
         End Do
      End If
      Go To 998
*                                                                      *
****** NODK ************************************************************
*                                                                      *
*     Do not compute Douglas-Kroll integrals.
*
 989  NoDKroll = .TRUE.
      Go To 998
*                                                                      *
****** ONEO ************************************************************
*                                                                      *
*     Do not compute two electron integrals.
*
 990  Onenly = .TRUE.
      Go To 998
*                                                                      *
****** TEST ************************************************************
*                                                                      *
*     Process only the input.
*
 991  Test = .TRUE.
      GWInput=Run_Mode.eq.G_Mode
      Go To 998
*                                                                      *
****** SDIP ************************************************************
*                                                                      *
*     Compute integrals for transition dipole moment
*
 992  Vlct_ = .TRUE.
      GWInput=.True.
      Go To 998
*                                                                      *
****** EPOT ************************************************************
*                                                                      *
*     Compute the electric potential for a number of points.  If nEF is
*     set to 0 this will cause the points to coincide with the unique
*     centers.
*
 9932 nOrdEF=Max(nOrdEF,0)
      GWInput=.True.
      IF (EFgiven) Then
        Call WarningMessage(2,'Only one of EPOT,EFLD,FLDG may be given')
        Call Quit_OnUserError()
      End If
      EFgiven=.True.
      Go To 9931
*                                                                      *
****** EFLD ************************************************************
*                                                                      *
*     Compute the electric potential and electric field for a number of
*     points. If nEF is set to 0 this will cause the points to coincide
*     with the unique centers.
*
 993  nOrdEF=Max(nOrdEF,1)
      GWInput=.True.
      IF (EFgiven) Then
        Call WarningMessage(2,'Only one of EPOT,EFLD,FLDG may be given')
        Call Quit_OnUserError()
      End If
      EFgiven=.True.
      Go To 9931
*                                                                      *
****** FLDG ************************************************************
*                                                                      *
*     Compute the electric potential, electric field, and electric field
*     gradient for a number of points. If nEF is set to 0 this will
*     cause the points to coincide with the Unique centers.
*
 994  nOrdEF=Max(nOrdEF,2)
      GWInput=.True.
      IF (EFgiven) Then
        Call WarningMessage(2,'Only one of EPOT,EFLD,FLDG may be given')
        Call Quit_OnUserError()
      End If
      EFgiven=.True.
      Go To 9931
*
 9931 KWord = Get_Ln(LuRd)
      Call Get_I1(1,nEF)
      If (nEF.lt.0) nEF = 0
      If (nEF.eq.0) Go To 998
      Call mma_allocate(EFt,3,nEF,label='nEF')
      Do iEF = 1, nEF
         KWord = Get_Ln(LuRd)
* Check whether a label is specified instead of a coordinate
         Call Get_S(1,Key,1)
         Call LeftAd(Key)
         Call Upcase(Key)
         jTmp = iChar(Key(1:1))
         If (jTmp .ge. 65 .AND. jTmp .le. 90) Then
            jEnd=Index(Key,' ')-1
            iOff = 0
            iFound_Label = 0
            Do iCnttp = 1, nCnttp
               Do iCnt = iOff+1, iOff+dbsc(iCnttp)%nCntr
                  If (Key(1:jEnd) .Eq. dc(iCnt)%LblCnt(1:jEnd)) Then
                     iFound_Label = 1
                     EFt(1:3,iEF)=dbsc(iCnttp)%Coor(1:3,iCnt-iOff)
                  End If
               End Do
               iOff = iOff + dbsc(iCnttp)%nCntr
            End Do
            If (iFound_Label .Eq. 0) Then
               Call WarningMessage(2,';'
     &                     //' Error in processing the keyword FLDG.;'
     &                     //' The label '''//Key(1:jEnd)
     &                  //''' could not be found among the centers.;'
     &                     //' Remember to specify the atom center'
     &                     //' before specifying the FLDG keyword.')
               Call Quit_OnUserError()
            End If
         Else
            Call Get_F(1,EFt(1,iEF),3)
            Call Upcase(KWord)
            If (Index(KWord,'ANGSTROM').ne.0)
     &         Call DScal_(3,One/angstr,EFt(1,iEF),1)
         End If
      End Do
      Go To 998
*                                                                      *
****** ANGM ************************************************************
*                                                                      *
*     Orbital angular momentum
*
 995  lOAM = .True.
      GWInput=.True.
      KWord = Get_Ln(LuRd)
      Call Upcase(KWord)
      Call Get_F(1,OAMt,3)
      If (Index(KWord,'ANGSTROM').ne.0)
     &    Call DScal_(3,One/angstr,OAMt,1)
      Go To 998
*                                                                      *
****** ANGM derivative restriction *************************************
*                                                                      *
*     Orbital angular momentum restriction
*
 1001 lUPONLY = .True.
      Go To 998
*
 1002 lDOWNONLY = .True.
      Go To 998
*                                                                      *
****** OMQI ************************************************************
*                                                                      *
*     Orbital magnetic quadrupole
*
 999  lOMQ = .True.
      GWInput=.True.
      KWord = Get_Ln(LuRd)
      Call Upcase(KWord)
      Call Get_F(1,OMQt,3)
      If (Index(KWord,'ANGSTROM').ne.0)
     &    Call DScal_(3,One/angstr,OMQt,1)
      Go To 998
*                                                                      *
****** AMPR ************************************************************
*                                                                      *
*     Angular momentum products
*
 9951 GWInput=.True.
      If (Run_Mode.eq.S_Mode.and.GWInput) Go To 9989
      Call mma_allocate(AMP_Center,3,Label='AMP_Center')
      KWord = Get_Ln(LuRd)
      Call Upcase(KWord)
      Call Get_F(1,AMP_Center,3)
      If (Index(KWord,'ANGSTROM').ne.0)
     &   AMP_Center(:)=(One/angstr)*AMP_Center(:)
      Go To 998
*                                                                      *
****** DSHD ************************************************************
*                                                                      *
*     Compute the diamagnetic shielding for a number of points. If nDMS
*     is set to 0 this will cause the points to coincide with the
*     unique centers.
*
 996  lDMS = .True.
      GWInput=.True.
      KWord = Get_Ln(LuRd)
      Call Get_F(1,Dxyz,3)
      KWord = Get_Ln(LuRd)
      Call Get_I1(1,nDMS)
      If (nDMS.lt.0) nDMS = 0
      If (nDMS.eq.0) Go To 998
      Call mma_allocate(DMSt,3,nDMS,label='DMSt')
      Do iDMS = 1, nDMS
         KWord = Get_Ln(LuRd)
         Call Upcase(KWord)
         Call Get_F(1,DMSt(1,iDMS),3)
         If (Index(KWord,'ANGSTROM').ne.0)
     &        Call DScal_(3,One/angstr,DMSt(1,iDMS),1)
      End Do
      Go To 998
*                                                                      *
****** NOPA ************************************************************
*                                                                      *
*     Set integral packing flag
*     Note      : this flag is only active if iWRopt=0
*     iPack=0   : pack 2el integrals (= Default)
*     iPack=1   : do not pack 2el integrals
*
 9910 iPack=1
      Go To 998
*                                                                      *
****** STDO ************************************************************
*                                                                      *
*     Set integral write option for 2 el integrals
*     iWRopt=0  : 2 el integrals are written in the MOLCAS2 format,
*                 i.e., in canonical order, no labels and packed format
*                 (= Default)
*     iWRopt=1  : 2 el integrals are written in a format identical
*                 to MOLECULE, i.e., values and labels
*
 9930 iWRopt=1
      Go To 998
*                                                                      *
****** PKTH ************************************************************
*                                                                      *
*     Read desired packing accuracy ( Default = 1.0D-14 )
*     Note      : this flag is only active if iWRopt=0
*
 9940 KWord = Get_Ln(LuRd)
      Call Get_F1(1,PkAcc)
      PkAcc = Abs(PkAcc)
      Go To 998
*                                                                      *
****** SKIP ************************************************************
*                                                                      *
*     Read skip parameters,i.e.,
*     if 2el integral symmetry blocks containing a given symmetry
*     will not be needed in subsequent calculations their computation
*     ans storage can be omitted.
*     ( Default = 0,0,0,0,0,0,0,0 )
*     Note      : this flag is only activ if iWRopt=0
*
 9950 KWord = Get_Ln(LuRd)
      lSkip = .True.
      ChSkip = KWord(1:80)
      Go To 998
*                                                                      *
****** EXTR ************************************************************
*                                                                      *
*     Put the program name and the time stamp onto the extract file
*
c23456789012345678901234567890123456789012345678901234567890123456789012
 9960 Write (LuWr,*)'RdCtl: keyword EXTRACT is obsolete and is ignored!'
      Go To 998
*                                                                      *
****** REAC ************************************************************
*                                                                      *
*     Read reaction field input.
*
 9970 Continue
      GWInput=.True.
      If (.Not.RF_read) Then
         Call InpRct(LuRd)
         If (lLangevin.or.PCM) Go To 9971
*
*------- Add a center corresponding to the center of the RF cavity.
*
         RF_read=.True.
 9971    Continue
*
      Else
         Call WarningMessage(2,
     &                  'RdCtl: A second RF-input block discovered!')
         Call Quit_OnUserError()
      End If
      Go To 998
*                                                                      *
***** GRID *************************************************************
*                                                                      *
 9773 Continue
      Call Funi_input(LuRd)
      Go To 998
*                                                                      *
****** CLIG ************************************************************
*                                                                      *
*     Speed of light (in au)
*
 9000 KWord = Get_Ln(LuRd)
      Call Get_F1(1,CLightAU)
      CLightAU = Abs(CLightAU)
      write(LuWr,*)'The speed of light in this calculation =', CLightAU
      Go To 998
*                                                                      *
***** NEMO *************************************************************
*                                                                      *
 800  Continue
      NEMO=.True.
      Go To 998
*                                                                      *
***** RMAT *************************************************************
*                                                                      *
*     RmatR    : radius of the R-matrix sphere (bohr)
*
 880  KWord = Get_Ln(LuRd)
      Call Get_F1(1,RMatR)
      Go To 998
*                                                                      *
***** RMEA *************************************************************
*                                                                      *
*     Epsabs   : absolute precision of numerical radial integration
*
 881  KWord = Get_Ln(LuRd)
      Call Get_F1(1,Epsabs)
      Go To 998
*                                                                      *
***** RMER *************************************************************
*                                                                      *
*     Epsrel   : relative precision of numerical radial integration
*
 882  KWord = Get_Ln(LuRd)
      Call Get_F1(1,Epsrel)
      Go To 998
*                                                                      *
***** RMQC *************************************************************
*                                                                      *
*     qCoul    : effective charge of the target molecule
*
 883  KWord = Get_Ln(LuRd)
      Call Get_F1(1,qCoul)
      Go To 998
*                                                                      *
***** RMDI *************************************************************
*                                                                      *
*     dipol(3) : effective dipole moment of the target molecule
*
 884  KWord = Get_Ln(LuRd)
      Call Get_F(1,dipol,3)
      dipol1=Abs(dipol(1))+Abs(dipol(2))+Abs(dipol(3))
      Go To 998
*                                                                      *
***** RMEQ *************************************************************
*                                                                      *
*     epsq     : minimal value of qCoul and/or dipol1 to be considered
*
 885  KWord = Get_Ln(LuRd)
      Call Get_F1(1,epsq)
      Go To 998
*                                                                      *
***** RMBP *************************************************************
*                                                                      *
*     bParm    : Bloch term parameter
*
 886  KWord = Get_Ln(LuRd)
      Call Get_F1(1,bParm)
      Go To 998
*                                                                      *
***** GIAO *************************************************************
*                                                                      *
*     Enable GIAO integrals.
*
 9020 GIAO=.True.
      Go To 998
*                                                                      *
***** NOCH *************************************************************
*                                                                      *
*     Deactivate Cholesky decomposition.
*
 9089 Continue
      Cholesky=.False.
      CholeskyWasSet=.True.
      Do_RI=.False.
      Go To 998
*                                                                      *
***** CHOL *************************************************************
*                                                                      *
*     Activate Cholesky decomposition with default settings.
*     This section can only be executed once.
*
 9091 Continue
      Do_RI=.False.
      If (.not. CholeskyWasSet) Then
         CholeskyWasSet=.True.
         Cholesky=.True.
         Do_RI=.False.
         DirInt = .True.
         Call Cho_Inp(.True.,-1,6)
         iChk_CH=1
      End If
      If ( (iChk_RI+iChk_DC) .gt. 0) Then
         Call WarningMessage(2,
     &           'Cholesky is incompatible with RI and Direct keywords')
         Call Quit_OnUserError()
      EndIf
      GWInput = .False.
      Go To 998
*                                                                      *
***** THRC *************************************************************
*                                                                      *
*     Set Cholesky decomposition threshold to specified value.
*
 9021 Continue
      KWord=Get_Ln(LuRd)
      Call Get_F1(1,CholeskyThr)
      Go To 998
*                                                                      *
***** 1CCD *************************************************************
*                                                                      *
*     Use one-center Cholesky.
*
 9022 Continue
      do1CCD=.true.
      Do_RI=.False.
      iChk_Ch=1
      If ( (iChk_RI+iChk_DC) .gt. 0) Then
         Call WarningMessage(2,
     &           'Cholesky is incompatible with RI and Direct keywords')
         Call Quit_OnUserError()
      EndIf
      Go To 998
*                                                                      *
***** CHOI *************************************************************
*                                                                      *
*     Activate Cholesky decomposition with user-defined settings.
*     This section can be executed any number of times.
*     User-defined settings will be preferred to defaults also if the
*     keywords appear in "wrong" order,
*
*     ChoInput
*     ...
*     End ChoInput
*     Cholesky
*
 9092 Continue
      Do_RI=.False.
      CholeskyWasSet=.True.
      Cholesky=.True.
      DirInt = .True.
      Call Cho_Inp(.False.,LuRd,6)
      iChk_CH=1
      If ( (iChk_RI+iChk_DC) .gt. 0) Then
         Call WarningMessage(2,
     &           'Cholesky is incompatible with RI and Direct keywords')
         Call Quit_OnUserError()
      EndIf
      Go To 998
*                                                                      *
***** RP-C *************************************************************
*                                                                      *
 9093 lRP=.True.
      KWord = Get_Ln(LuRd)
      jTmp=0
      nRP_prev=-1
      Read(KWord,*,err=9082) nRP
      isnumber=1
      ifile=index(KWord,' ')
      do i=1,ifile
        ii=index(' 0123456789',KWord(i:i))
        if(ii.eq.0) then
           isnumber=0
        endif
      enddo
*
      if(isnumber.eq.0) goto 9082
      nRP=3*nRP
      Call mma_allocate(RP_Centers,3,nRP/3,2,Label='RP_Centers')
*
**    Inline input
*
      Call UpCase(KWord)
      If (Index(KWord,'ANGSTROM').ne.0) Then
         Fact=One/Angstr
      Else
         Fact=One
      End If
*
*
      KWord = Get_Ln(LuRd)
      Call Get_F1(1,E1)
      Call Read_v(LuRd,RP_Centers(1,1,1),1,nRP,1,iErr)
      KWord = Get_Ln(LuRd)
      Call Get_F1(1,E2)
      Call Read_v(LuRd,RP_Centers(1,1,2),1,nRP,1,iErr)
      RP_Centers(:,:,:)=Fact*RP_Centers(:,:,:)
      GWInput = .True.
      Go To 998
*
**    Files
*
 9082 Continue
      RPSet=.true.
      jTmp=jTmp+1
      ifile=index(KWord,' ')
      If(KWord(1:1).eq.'/') then
        call f_inquire(KWord(1:ifile-1),Exist)
        Key=KWord
      Else
        call getenvf('MOLCAS_SUBMIT_DIR',Directory)
        if(Directory(1:1).ne.' ') then
          i=index(Directory,' ')
          Key=Directory(1:i-1)//'/'//KWord(1:ifile-1)
          ifile=i+ifile
          call f_inquire(Key(1:iFile-1),Exist)
        Else
          Exist=.false.
        End If
        If (.not.Exist) Then
          Key=Key(i+1:iFile-1)
          ifile=ifile-i
          call f_inquire(Key(1:iFile-1),Exist)
        End If
      End If
      If (.not.Exist) Then
        Call WarningMessage(2,'File '//
     &                            Key(1:ifile)//' is not found')
        Call Quit_OnUserError()
      End If
      LuIn=8
      LuIn=isFreeUnit(LuIn)
      call molcas_open(LuIn,Key(1:iFile-1))
*
      KWord = Get_Ln(LuIn)
      Read(KWord,*,err=9083) nRP
      If (nRP_prev.ge.0 .and. nRP.ne.nRP_prev) Then
        Call WarningMessage(2,'The numbers of atoms in the two RP'//
     &                        ' structures do not match.')
        Call Quit_OnUserError()
      End IF
      nRP_prev=nRP
      nRP=3*nRP
      If (.NOT.Allocated(RP_Centers))
     &   Call mma_allocate(RP_Centers,3,nRP/3,2,Label='RP_Centers')
      KWord = Get_Ln(LuIn)
      Call UpCase(KWord)
      If (Index(KWord,'BOHR').ne.0) Then
         Fact=One
      Else
         Fact=One/Angstr
      End If
      If (jTmp.eq.1) Then
         LuRP=10
         LuRP=isFreeUnit(LuRP)
         call molcas_open(LuRP,'findsym.RP1')
         Read(KWord,*,err=9083) E1
*
**  write a separate file for findsym
*
         Write(LuRP,*) nRP/3
#ifdef _HAVE_EXTRA_
         Write(LuRP,'(a)')
#else
         Write(LuRP,'(a)') 'bohr'
#endif
         Do i=1,nRP/3
            KWord = Get_Ln(LuIn)
            Read(KWord,*,err=9083) Key,(RP_Centers(j,i,1),j=1,3)
            Write(LuRP,'(A,3F20.12)') Key(1:LENIN),
     &                   (RP_Centers(j,i,1)*Fact,j=1,3)
         End Do
         KWord = Get_Ln(LuRd)
         close(LuIn)
         close(LuRP)
         Go To 9082
      Else
         LuRP=10
         LuRP=isFreeUnit(LuRP)
         call molcas_open(LuRP,'findsym.RP2')
         Write(LuRP,*) nRP/3
#ifdef _HAVE_EXTRA_
         Write(LuRP,'(a)')
#else
         Write(LuRP,'(a)') 'bohr'
#endif
         Read(KWord,*,err=9083) E2
         Do i=1,nRP/3
            KWord = Get_Ln(LuIn)
            Read(KWord,*,err=9083) Key,(RP_Centers(j,i,2),j=1,3)
            Write(LuRP,'(A,3F20.12)') Key(1:LENIN),
     &            (RP_Centers(j,i,2)*Fact,j=1,3)
         End Do
         close(LuRP)
      End If
      RP_Centers(:,:,:)= Fact* RP_Centers(:,:,:)
*
      close(LuIn)
      GWInput = Run_Mode.eq.G_Mode
      Go To 998
*
**    Error
*
 9083 Continue
      Write(6,'(a,a)') 'Error reading from file ',Key(1:iFile-1)
      Write(6,'(a,a)') 'unable to process line: ',KWord
      Call Quit_OnUserError()
*                                                                      *
***** SADD *************************************************************
*                                                                      *
*     Saddle options
 9081 Key = Get_Ln(LuRd)
      Call Get_F1(1,SadStep)
      GWInput = .True.
      Go To 998
*                                                                      *
***** CELL *************************************************************
*                                                                      *
*     VCell(3,3)    : the vectors of the cell
*
 887  Key = Get_Ln(LuRd)
      Call Upcase(Key)
      If (Index(Key,'ANGSTROM').ne.0) KWord = Get_Ln(LuRd)
      Call Get_F(1,VCell(1,1),3)
      KWord = Get_Ln(LuRd)
      Call Get_F(1,VCell(1,2),3)
      KWord = Get_Ln(LuRd)
      Call Get_F(1,VCell(1,3),3)
      If (Index(Key,'ANGSTROM').ne.0) Call DScal_(9,One/angstr,VCell,1)
      Cell_l = .TRUE.
      Call mma_allocate(AdCell,MxAtom)
      Go To 998
*                                                                      *
***** SPAN *************************************************************
*                                                                      *
*     Set span factor in Cholesky decomposition (0 < span < 1).
*     The span decides the smallest diagonal element that can be
*     treated as span*max(Diag). Span=1 thus implies full pivoting.
*
 890  Continue
      KWord = Get_Ln(LuRd)
      Call Get_F1(1,spanCD)
      spanCD=abs(spanCD)
      Go To 998
*                                                                      *
***** SPREAD ***********************************************************
*                                                                      *
*     ispread(3)    : the number of cells to spread in different directions
*
 889  KWord = Get_Ln(LuRd)
      Call Get_I(1,ispread,3)
      Go To 998
*                                                                      *
***** LOW  *************************************************************
*                                                                      *
*     Activate low-accuracy Cholesky decomposition.
*
 9094 Continue
      Do_RI=.False.
      If (.not. CholeskyWasSet) Then
         CholeskyWasSet=.True.
         Cholesky=.True.
         DirInt = .True.
         Call Cho_Inp(.True.,-1,6)
         Call Cho_InpMod('LOW ')
         Thrshld_CD=1.0D-4
      End If
      Go To 998
*                                                                      *
***** MEDI *************************************************************
*                                                                      *
*     Activate medium-accuracy Cholesky decomposition.
*
 9095 Continue
      Do_RI=.False.
      If (.not. CholeskyWasSet) Then
         CholeskyWasSet=.True.
         Cholesky=.True.
         DirInt = .True.
         Call Cho_Inp(.True.,-1,6)
         Call Cho_InpMod('MEDI')
         Thrshld_CD=1.0D-6
      End If
      Go To 998
*                                                                      *
***** HIGH *************************************************************
*                                                                      *
*     Activate high-accuracy Cholesky decomposition.
*
 9096 Continue
      Do_RI=.False.
      If (.not. CholeskyWasSet) Then
         CholeskyWasSet=.True.
         Cholesky=.True.
         DirInt = .True.
         Call Cho_Inp(.True.,-1,6)
         Call Cho_InpMod('HIGH')
         Thrshld_CD=1.0D-8
      End If
      Go To 998
*                                                                      *
***** DIAG *************************************************************
*                                                                      *
 9087 Continue
      DiagCheck=.true.
      Go To 998
*                                                                      *
***** RI   *************************************************************
*                                                                      *
*     Activate RI approach
*
 9097 Continue
      Do_RI=.True.
      GWInput=.True.
      iRI_Type=3
      iChk_RI=1
      If ( (iChk_DC+iChk_CH) .gt. 0) Then
         Call WarningMessage(2,
     &          'RI is incompatible with Direct and Cholesky keywords')
         Call Quit_OnUserError()
      End If
      Go To 998
 9098 Continue
      Do_RI=.True.
      GWInput=.True.
      iRI_Type=1
      iChk_RI=1
      If ( (iChk_DC+iChk_CH) .gt. 0) Then
         Call WarningMessage(2,
     &          'RI is incompatible with Direct and Cholesky keywords')
         Call Quit_OnUserError()
      End If
      Go To 998
 9099 Continue
      Do_RI=.True.
      GWInput=.True.
      iRI_Type=2
      iChk_RI=1
      If ( (iChk_DC+iChk_CH) .gt. 0) Then
         Call WarningMessage(2,
     &          'RI is incompatible with Direct and Cholesky keywords')
         Call Quit_OnUserError()
      End If
      Go To 998
 9080 Continue
      Do_RI=.True.
      GWInput=.True.
      iRI_Type=4
      iChk_RI=1
      If ( (iChk_DC+iChk_CH) .gt. 0) Then
         Call WarningMessage(2,
     &          'RI is incompatible with Direct and Cholesky keywords')
         Call Quit_OnUserError()
      End If
      Go To 998
 9085 Continue
      Do_RI=.True.
      GWInput=.True.
      iRI_Type=5
      iChk_RI=1
      If ( (iChk_DC+iChk_CH) .gt. 0) Then
         Call WarningMessage(2,
     &          'RI is incompatible with Direct and Cholesky keywords')
         Call Quit_OnUserError()
      End If
      Go To 998
*                                                                      *
***** NOGU *************************************************************
*                                                                      *
*     Disable atomatic execution of GuessOrb
*
 9100 Do_GuessOrb=.FALSE.
      Go To 998
*                                                                      *
***** RELA *************************************************************
*                                                                      *
*     DKH option: order and parameterization.
*     xx: order of Hamiltonian
*      y: parameterization
*     zz: order of properties
*
 657  Continue
      kWord = Get_Ln(LuRd)
      If (KWord(1:1).eq.'R' .and.
     &    (KWord(2:2).ge.'0' .and.
     &     KWord(2:2).le.'9') .and.
     &     (KWord(3:3).ge.'0' .and.
     &      KWord(3:3).le.'9') .and.
     &      (KWord(4:4).eq.'O' .or.
     &       KWord(4:4).eq.'E' .or.
     &       KWord(4:4).eq.'S' .or.
     &       KWord(4:4).eq.'M' .or.
     &       KWord(4:4).eq.'C') ) Then
      Else
         Call WarningMessage(2,'Error in RELA keyword')
         Call Quit_OnUserError()
      End If
      DKroll=.True.
*
*     DKH order in the Hamiltonian
*
      read( KWord(2:3), * ) idk_ord
      IRELAE = 1000 + idk_ord * 10
*
*     Method of parametrization
*
      If (kWord(4:4).eq.'O') IRELAE=IRELAE+1
      If (kWord(4:4).eq.'E') IRELAE=IRELAE+2
      If (kWord(4:4).eq.'S') IRELAE=IRELAE+3
      If (kWord(4:4).eq.'M') IRELAE=IRELAE+4
      If (kWord(4:4).eq.'C') IRELAE=IRELAE+5
*
*     DKH order in the property integrals
*
      If ( KWord(5:5).ge.'0' .and.
     &     KWord(5:5).le.'9' .and.
     &     KWord(6:6).ge.'0' .and.
     &     KWord(6:6).le.'9' ) then
        read( KWord(5:6), * ) iprop_ord
          IRELAE = IRELAE + iprop_ord * 10000
      Else
          IRELAE = IRELAE + idk_ord * 10000
      End If
*
      Go To 998
*                                                                      *
***** LDKH *************************************************************
*                                                                      *
*     Local Douglas-Kroll-Hess/X2C/BSS
*
 658  LDKroll=.True.
*     GWInput=.True.
      nCtrLD=0
      radiLD=5.5d0
*
      KWord = Get_Ln(LuRd)
      Call Upcase(KWord)
      If (KWord(1:3).eq.'DLU') Go To 998
      If (KWord(1:3).eq.'DLH') Then
        radiLD=0.0d0
        Go To 998
      End If
      read(KWord,*,end=6582, err=6582) nCtrLD, radiLD
      If (nCtrLD.gt.10) Then
         Call WarningMessage(2,
     &           'The number of centers for LDKH is limited to 10')
         call Quit_OnUserError()
      End If
      If (Index(KWord,'ANGSTROM').ne.0) Then
          radiLD = radiLD/angstr
      End If
*
      KWord = Get_Ln(LuRd)
      Call Upcase(KWord)
      read(KWord,*,end=6666, err=6581) (iCtrLD(i),i=1,nCtrLD)
      Go To 998
 6581 read(Kword,*,end=6666, err=6666) (CtrLDK(i),i=1,nCtrLD)
      Call Get_nAtoms_all(nAtom)
      k=0
      Do i=1,nAtom
        Do j=1,nCtrLD
          if (CtrLDK(j).eq.dc(i)%LblCnt(1:LENIN)) Then
             iCtrLD(j)=i
             k=k+1
          End If
        End Do
      End Do
      If (k.ne.nCtrLD) Then
         Call WarningMessage(2,'Error in LDKH Centers definitions')
         Call Quit_OnUserError()
      End If
      Go To 998
*
**    Automatic choice: all heavy elements (from K 19)
*
 6582 Continue
CDP      If (nCtrLD.eq.0) radiLD=0.0d0
      Key=KWord
      Go To 9989
*                                                                      *
***** FOOC *************************************************************
*                                                                      *
*     Force the use of the out-of-core RI algorithm.
*
 8000 Force_Out_of_Core=.True.
      Go To 998
*                                                                      *
***** CDTH *************************************************************
*                                                                      *
*     Threshold for CD to generate RICD auxiliary basis sets
*
 8001 Key = Get_Ln(LuRd)
      Call Get_F1(1,Thrshld_CD)
      GWInput=.True.
      Go To 998
*                                                                      *
***** SHAC *************************************************************
*                                                                      *
*     Skip high angular combinations when constructing RICD aux basis.
*
 8002 Skip_High_AC=.True.
      GWInput=.True.
      Go To 998
*                                                                      *
***** KHAC *************************************************************
*                                                                      *
*     Keep high angular combinations when constructing RICD aux basis.
*
 8003 Skip_High_AC=.False.
      GWInput=.True.
      Go To 998
*                                                                      *
***** ACD  *************************************************************
*                                                                      *
*     Generate a aCD basis.
*
 8004 Do_acCD_Basis=.False.
      GWInput=.True.
      Go To 998
*                                                                      *
***** ACCD *************************************************************
*                                                                      *
*     Generate a acCD basis.
*
 8005 Do_acCD_Basis=.True.
      GWInput=.True.
      Go To 998
*                                                                      *
***** NACC *************************************************************
*                                                                      *
*     Generate a nacCD basis.
*
 8030 Do_acCD_Basis=.False.
      Do_nacCD_Basis=.True.
      GWInput=.True.
      Go To 998
*                                                                      *
***** DOFM *************************************************************
*                                                                      *
*     DoFMM: activate FMM option
*
 8006 DoFMM = .True.
      Go To 998
*                                                                      *
***** NOAM *************************************************************
*                                                                      *
*     No computation of AMFI integrals
*
 8007 NoAMFI=.True.
      GWInput=.True.
      Go To 998
*                                                                      *
***** RPQM *************************************************************
*                                                                      *
*     Set RPQMin for FMM option
*
 8008 Key = Get_Ln(LuRd)
      Call Get_F1(1,RPQMin)
      Go To 998
*                                                                      *
***** CONS *************************************************************
*                                                                      *
*     Have the Gateway read the constraints for Slapaf
*
 8010 Continue
      GWInput=.True.
      Lu_UDC=97
      Lu_UDC = IsFreeUnit(Lu_UDC)
      Call Molcas_Open(Lu_UDC,'UDC.Gateway')
 8011 Continue
         Key=Get_Ln(LuRd)
         Call UpCase(Key)
         Write (Lu_UDC,'(A)') Trim(Key)
         If (Key(1:4).ne.'END ') Go To 8011
*     This rather obscure feature seems to be needed to to make Intel
*     compilers behave like the others when detecting EOF
      End File(Lu_UDC)
      Close(Lu_UDC)
      Go To 998
*                                                                      *
***** NGEX *************************************************************
*                                                                      *
*     Have the Gateway read the constraints for Numerical_gradient
*
  501 Continue
      GWInput=.True.
      Lu_UDC=97
      Lu_UDC = IsFreeUnit(Lu_UDC)
      Call Molcas_Open(Lu_UDC,'UDC.NG')
  502 Continue
         Key=Get_Ln(LuRd)
         Call UpCase(Key)
         If (AdjustL(Key).eq.'INVERT') Then
            Invert=.True.
            Go To 502
         End If
         Write (Lu_UDC,'(A)') Trim(Key)
         If (Key(1:4).ne.'END ') Go To 502
*     This rather obscure feature seems to be needed to to make Intel
*     compilers behave like the others when detecting EOF
      End File(Lu_UDC)
      Close(Lu_UDC)
      Go To 998
*                                                                      *
***** LOCA or LDF1 or LDF **********************************************
*                                                                      *
*     Activate Local Density Fitting.
*
   35 Continue
         LocalDF=.True.
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** LDF2 *************************************************************
*                                                                      *
*     Activate Local Density Fitting with 2-center functions included
*     when needed to achieve target accuracy.
*
   36 Continue
         LocalDF=.True.
         Call LDF_SetLDF2(.True.)
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** TARG or THRL *****************************************************
*                                                                      *
*     Set target accuracy for Local Density Fitting.
*     This implies inclusion of 2-center functions (the only way we can
*     affect accuracy).
*
   37 Continue
         Key=Get_Ln(LuRd)
         Call Get_F1(1,Target_Accuracy)
         Call LDF_SetThrs(Target_Accuracy)
         LocalDF=.True.
         Call LDF_SetLDF2(.True.)
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** APTH *************************************************************
*                                                                      *
*     Set screening threshold for LDF - i.e. threshold for defining
*     significant atom pairs.
*
   38 Continue
         Key=Get_Ln(LuRd)
         Call Get_F1(1,APThr)
         Call LDF_SetPrescreen(APThr)
         LocalDF=.True.
         APThr_UsrDef=.True.
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** CHEC *************************************************************
*                                                                      *
*     LDF debug option: check pair integrals.
*
   39 Continue
         Call LDF_SetOptionFlag('CHEC',.True.)
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** VERI *************************************************************
*                                                                      *
*     LDF debug option: verify fit for each atom pair.
*
   40 Continue
         Call LDF_SetOptionFlag('VERI',.True.)
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** OVER *************************************************************
*                                                                      *
*     LDF debug option: check overlap integrals (i.e. charge)
*
   41 Continue
         Call LDF_SetOptionFlag('OVER',.True.)
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** CLDF *************************************************************
*                                                                      *
*     Constrained LDF - read constraint order
*     order=-1 --- unconstrained
*     order=0  --- charge (i.e. overlap)
*
   42 Continue
         Key=Get_Ln(LuRd)
         Call Get_I1(1,iCLDF)
         Call LDF_AddConstraint(iCLDF)
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** UNCO *************************************************************
*                                                                      *
*     Unconstrained LDF (same as CLDF=-1)
*
   43 Continue
         Call LDF_AddConstraint(-1)
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** WRUC *************************************************************
*                                                                      *
*     Write unconstrained coefficients to disk.
*     Only meaningful along with constrained fitting.
*     For debugging purposes: enables constrained fit verification in
*     modules other than Seward.
*
   44 Continue
         Call LDF_SetOptionFlag('WRUC',.True.)
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** UNIQ *************************************************************
*                                                                      *
*     LDF: use unique atom pairs.
*
   45 Continue
         Call LDF_SetOptionFlag('UNIQ',.True.)
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** NOUN *************************************************************
*                                                                      *
*     LDF: do not use unique atom pairs.
*
   46 Continue
         Call LDF_SetOptionFlag('UNIQ',.False.)
         GWInput=.False. ! Only in Seward
      Go To 998
*                                                                      *
***** RLDF *************************************************************
*                                                                      *
*     Activate local DF/RI, Roland's original LDF test implementation
*
 8012 LDF=.True.
      GWInput=.True.
      Go To 998
*                                                                      *
***** NOAL *************************************************************
*                                                                      *
*     Do not align reactants and products
*
 7013 Do_Align=.False.
      If (Align_Only) Then
         Call WarningMessage(2,
     &       'Keywords ALIG and NOAL are not compatible')
         Call Quit_OnUserError()
      End If
      GWInput=.True.
      Go To 998
*                                                                      *
***** WEIG *************************************************************
*                                                                      *
*     Weights for alignment of reactants and products
*
 7014 Align_Weights=Get_Ln(LuRd)
      Call UpCase(Align_Weights)
      GWInput=.True.
      Go To 998
*                                                                      *
***** ALIG *************************************************************
*                                                                      *
*     Align reactants and products
*
 8013 Align_Only=.True.
      If (.not.Do_Align) Then
         Call WarningMessage(2,
     &       'Keywords ALIG and NOAL are not compatible')
         Call Quit_OnUserError()
      End If
      GWInput=.True.
      Go To 998
*                                                                      *
****** TINK ************************************************************
*                                                                      *
*     Read Coordinates in Tinker's xyz format
*
8014  If (SymmSet) Then
         Call WarningMessage(2,
     &                 'SYMMETRY keyword is not compatible with TINKER')
         Call Quit_OnUserError()
      End If
      DoTinker = .True.
      ITkQMMM = 1
      If (MyRank.eq.0) Then
        ITkQMMM = IsFreeUnit(ITkQMMM)
        Call Molcas_Open (ITkQMMM,'QMMM')
        Write(ITkQMMM,'(A)') 'Molcas -1 0'
        Close (ITkQMMM)
c
        Call Getenvf('TINKER ',Key)
        mLine = Len(Key)
        iLast = iCLast(Key,mLine)
        if (iLast.eq.0) Then
          Call Getenvf('MOLCAS',Key)
          mLine = Len(Key)
          iLast = iCLast(Key,mLine)
          Key = Key(1:iLast)//'/tinker/bin'
        End If
        iLast = iCLast(Key,mLine)
        Call Getenvf('Project',Project)
        mLine = Len(Project)
        jLast = iCLast(Project,mLine)
        Key = Key(1:iLast)//'/tkr2qm_s '//Project(1:jLast)//'.xyz'//
     &                 '>'//Project(1:jLast)//'.Tinker.log'
        mLine = Len(Key)
        iLast = iCLast(Key,mLine)
        Write(6,*) 'TINKER keyword found, run ',Key(1:iLast)
        Call StatusLine(' Gateway:',' Read input from Tinker')
        RC=0
        Call Systemf(Key(1:iLast),RC)
        If (RC.ne.0) Then
          Key='RdCtl_Seward: Tinker call terminated abnormally'
          Call WarningMessage(2,Key)
          Call Abend()
        End If
      End If
#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
         Call GA_Sync()
         Call PFGet_ASCII('QMMM')
         Call GA_Sync()
      End If
#endif
      iCoord=iCoord+1
      CoordSet=.True.
      ITkQMMM = IsFreeUnit(ITkQMMM)
      Call Molcas_Open (ITkQMMM,'QMMM')
#ifdef _HAVE_EXTRA_
      If (Expert) Then
         If (iCoord.gt.1) Then
            Call WarningMessage(1,
     &         'TINKER and COORD keywords cannot be combined '//
     &         'with molcas_extra')
         End If
      End If
      Call XYZread(ITkQMMM,ForceZMAT,nCoord,iErr)
      If (iErr.ne.0) Then
        Key='RdCtl_Seward: Tinker+XYZread failed:'//
     &          ' check Tinker input files'
        Call WarningMessage(2,Key)
        Call Abend()
      End If
      Call XYZcollect(iCoord,nCoord,OrigTrans,OrigRot,nFragment)
#else
      If (Expert) Then
         If (iCoord.gt.1) Then
            Call WarningMessage(1,
     &          'TINKER coordinates replacing COORD')
         End If
         Call Read_XYZ(ITkQMMM,OrigRot,OrigTrans,Replace=(iCoord.gt.1))
      Else
         Call Read_XYZ(ITkQMMM,OrigRot,OrigTrans)
      End If
#endif
      Close(ITkQMMM)
      GWInput = .True.
      Go To 998
*                                                                      *
***** ORIG *************************************************************
*                                                                      *
*     Defines translation and rotation for each xyz-file
*
 8015 OrigInput = .True.
#ifdef _HAVE_EXTRA_
      Origin_input = .True.
#endif
      If(FragSet) Then
         Write(6,*) 'Keywords FRGM and ORIG are mutually exclusive!'
         Call Quit_OnUserError()
      End If
      If(.not. OriginSet) Then
         Call mma_allocate(OrigTrans,3,nFragment,label='OrigTrans')
         Call mma_allocate(OrigRot,3,3,nFragment,label='OrigRot')
         OriginSet = .True.
      End If
      Do iFrag = 1, nFragment
         KWord = Get_Ln(LuRd)
         Call Get_F(1,OrigTrans(1,iFrag),3)
         KWord = Get_Ln(LuRd)
         Call Get_F(1,OrigRot(1,1,iFrag),9)
      End Do
      GWinput = .True.
      Go To 998
*                                                                      *
***** HYPE *************************************************************
*                                                                      *
 8016 KWord = Get_Ln(LuRd)
#ifdef _HAVE_EXTRA_
      geoInput = .true.
#endif
      writeZmat = .true.
      Call Get_F(1,HypParam,3)
      GWinput = .True.
      HyperParSet = .True.
      Go To 998
*                                                                      *
***** ZCON *************************************************************
*                                                                      *
 8017 writeZMat = .true.
#ifdef _HAVE_EXTRA_
      ZConstraints = .true.
#endif
      GWinput = .True.
      Go To 998
*                                                                      *
***** SCAL *************************************************************
*                                                                      *
 8018 KWord = Get_Ln(LuRd)
      Call Get_F1(1,ScaleFactor)
      GWinput = .True.
      if(.not.CoordSet) then
         Call WarningMessage(2,'Scale can be used only with xyz input')
         Call Quit_OnUserError()
      endif
      Go To 998
*                                                                      *
***** DOAN *************************************************************
*                                                                      *
 8019 Call Put_iScalar('agrad',1)
      Go To 998
*                                                                      *
***** GEOE *************************************************************
*                                                                      *
 8020 Kword = Get_Ln(LuRd)
      Call Get_I1(1,iGeoInfo(2))
      GWinput = .True.
      iGeoInfo(1) = 1
      Call Put_iArray('GeoInfo',iGeoInfo,2)
      Go To 998
*                                                                      *
***** OLDZ *************************************************************
*                                                                      *
 8021 GWinput = .True.
#ifdef _HAVE_EXTRA_
      oldZmat = .True.
#endif
      Go To 998
*                                                                      *
***** OPTH *************************************************************
*                                                                      *
 8022 GWinput = .True.
      Kword = Get_Ln(LuRd)
      Call Get_I1(1,iOptimType)
      Kword = Get_Ln(LuRd)
      Call Get_F1(1,StepFac1)
      if(iOptimType .eq. 2) Then
         KWord = Get_Ln(LuRd)
         Call Get_F1(1,gradLim)
      end if
      Go To 998
*                                                                      *
***** NOON *************************************************************
*                                                                      *
 8023 Do_OneEl = .False.
      Go To 998
*                                                                      *
***** GEO  *************************************************************
*                                                                      *
 8024 Continue
#ifdef _HAVE_EXTRA_
      geoInput = .true.
#endif
      writeZMat = .true.
*     Parameters for the gridsize is set to default-values if geo is
*     used instead of hyper
      If(.not. HyperParSet) Then
            HypParam(1) = 0.15d0
            HypParam(2) = 2.5d0
            HypParam(3) = 2.5d0
      End If
      GWinput = .True.
      Go To 998
*                                                                      *
***** GEN1INT **********************************************************
*                                                                      *
*        GEN1INT integrals
 9023 IF(IRELAE.EQ.101) Then
        lMXTC=.true.
      ELSE
       Write(6,*) 'Keyword MXTC must be preceded by keyword RX2C!'
       Call Quit_OnUserError()
      ENDIF
      Go To 998
*                                                                      *
***** FRGM *************************************************************
*                                                                      *
 8025 OrigInput = .True.
#ifdef _HAVE_EXTRA_
      Origin_input = .True.
#endif
      GWinput = .True.
      If(OriginSet) Then
         Write(6,*) 'Keywords FRGM and ORIG are mutually exclusive!'
         Call Quit_OnUserError()
      End If
      If(.not.FragSet) then
         Call mma_allocate(OrigTrans,3,nFragment,label='OrigTrans')
         Call mma_allocate(OrigRot,3,3,nFragment,label='OrogRot')
*     Set up no translation and no rotation as default
         Call FZero(OrigTrans,3*nFragment)
         Call FZero(OrigRot,9*nFragment)
         Do i = 1, nFragment
            OrigRot(1,1,i)   = 1.0d0
            OrigRot(2,2,i)   = 1.0d0
            OrigRot(3,3,i)   = 1.0d0
         End Do
         FragSet = .True.
      End If
      Kword = Get_Ln(LuRd)
      Call Get_I1(1,iFrag)
      Go To 998
*                                                                      *
***** TRAN *************************************************************
*                                                                      *
 8026 If(.not. FragSet) Then
         Write(6,*) 'Keyword TRANS must be preceded by keyword FRAG!'
         Call Quit_OnUserError()
      End If
      GWinput = .True.
      Kword = Get_Ln(LuRd)
      Call Get_F(1,OrigTrans(1,iFrag),3)
      Go To 998
*                                                                      *
****** ROT  ************************************************************
*                                                                      *
 8027 If(.not. FragSet) Then
         Write(6,*) 'Keyword ROT must be preceded by keyword FRAG!'

      End If
      GWinput = .True.
      Kword = Get_Ln(LuRd)
      Call Get_F(1,OrigRot(1,1,iFrag),9)
      Go To 998
*                                                                      *
******* ZONL ***********************************************************
*                                                                      *
 8028 GWinput = .True.
      WriteZMat = .True.
      Go To 998
*                                                                      *
******* BASL ***********************************************************
*                                                                      *
 8029 GWinput = .True.
      BasLib=Get_Ln(LuRd)
      Write_BasLib=.True.
      Go To 998
*                                                                      *
******* NUME ***********************************************************
*                                                                      *
 8031 GWinput = .True.
      iDNG=1
      Go To 998
*                                                                      *
******* VART ***********************************************************
*                                                                      *
 8032 GWinput = .True.
      VarT=.True.
      Go To 998
*                                                                      *
******* VARR ***********************************************************
*                                                                      *
 8033 GWinput = .True.
      VarR=.True.
      Go To 998
*                                                                      *
******* SHAK ***********************************************************
*                                                                      *
 8050 Continue
      GWinput = .True.
      KWord = Get_Ln(LuRd)
      Call Upcase(KWord)
      Call Get_F1(1,Shake)
      If (Index(KWord,'ANGSTROM').ne.0) Shake = Shake/angstr
      Go To 998
*                                                                      *
****** PAMF ************************************************************
*                                                                      *
*     Disable AMFI for an atom type
*
 8060 KWord = Get_Ln(LuRd)
      Call Get_I1(1,iAtom_Number)
      No_AMFI(iAtom_Number)=.True.
      Go To 998
*                                                                      *
******* GROM ***********************************************************
*                                                                      *
*     Import definition of QMMM system from Gromacs
*
 8034 Continue
#ifdef _GROMACS_
      If (SymmSet) Then
         Message = 'SYMMETRY keyword is not compatible with GROMACS'
         Call WarningMessage(2,Message)
         Call Quit_OnUserError()
      End If
      DoGromacs = .True.
      GWInput = .True.
      VarT = .True.
      VarR = .True.
*
* Check for options
      KWord = Get_Ln(LuRd)
      Call UpCase(KWord)
      If (KWord(1:4).Eq.'SIMP') Then
         nCastMM = 0
         Call mma_allocate(CastMM,nCastMM)
      Else If (KWord(1:4).Eq.'CAST') Then
         KWord = Get_Ln(LuRd)
         Call Get_I(1,nCastMM,1)
         If (nCastMM.LE.0) Then
            Message = 'nCastMM is zero or negative'
            Call WarningMessage(2,Message)
            Call Quit_OnUserError()
         End If
         Call mma_allocate(CastMM,nCastMM)
         KWord = Get_Ln(LuRd)
         Call Get_I(1,CastMM,nCastMM)
         Do iCastMM = 1,nCastMM
            If (CastMM(iCastMM).LE.0) Then
               Message = 'Impossible, MM index < 1'
               Call WarningMessage(2,Message)
               Call Quit_OnUserError()
            End If
         End Do
      Else
         Message='GROMACS keyword found, but no valid option'
         Call WarningMessage(2,Message)
         Call Quit_OnUserError()
      End If
*
* After the call to Fetch_QMMM, the inner subsystem is in a temporary
* xyz file and the outer subsystem is on the runfile
      Call Fetch_QMMM(CastMM,nCastMM)
*
      Call mma_deallocate(CastMM)
*
* Let Molcas read the xyz file
      iCoord = iCoord+1
      CoordSet = .True.
      LuXYZ = 1
      LuXYZ = isFreeUnit(LuXYZ)
      Call molcas_open(LuXYZ,'GMX.XYZ')
#ifdef _HAVE_EXTRA_
      If (Expert) Then
         If (iCoord.gt.1) Then
            Call WarningMessage(1,
     &         'GROMACS and COORD keywords cannot be combined '//
     &         'with molcas_extra')
         End If
      End If
      Call XYZread(LuXYZ,ForceZMAT,nCoord,iErr)
      If (iErr.NE.0) Then
         Message='RdCtl_Seward: XYZread returned non-zero error code'
         Call WarningMessage(2,Message)
         Call Abend()
      End If
      Call XYZcollect(iCoord,nCoord,OrigTrans,OrigRot,nFragment)
#else
      If (Expert) Then
         If (iCoord.gt.1) Then
            Call WarningMessage(1,
     &          'TINKER coordinates replacing COORD')
         End If
         Call Read_XYZ(LuXYZ,OrigRot,OrigTrans,Replace=(iCoord.gt.1))
      Else
         Call Read_XYZ(LuXYZ,OrigRot,OrigTrans)
      End If
#endif
      Close(LuXYZ)
#else
      Message = 'Interface to Gromacs not installed'
      Call WarningMessage(2,Message)
      Call Quit_OnUserError()
#endif
      Go To 998
*                                                                      *
******* LINK ***********************************************************
*                                                                      *
*     Define link atoms for a Molcas/Gromacs run
*
 8036 Continue
#ifdef _GROMACS_
      GWInput = .True.
      KWord = Get_Ln(LuRd)
      Call Get_I(1,nLA,1)
      If (nLA.LE.0) Then
         Message = 'LA definition: nLA is zero or negative'
         Call WarningMessage(2,Message)
         Call Quit_OnUserError()
      End If
#ifdef _DEBUGPRINT_
      Write(LuWr,'(/,a)') ' Link atoms (Gromacs numbering):'
      Write(LuWr,'(/,a)') '      LA     QM     MM     Scaling factor'
#endif
      Call mma_allocate(DefLA,3,nLA)
      Call mma_allocate(FactLA,nLA)
      Do iLA = 1,nLA
         KWord = Get_Ln(LuRd)
         Call Get_I(1,DefLA(1,iLA),3)
         Call Get_F(4,FactLA(iLA),1)
#ifdef _DEBUGPRINT_
         Write(LuWr,'(i8,2i7,F19.8)') (DefLA(i,iLA),i=1,3),FactLA(iLA)
#endif
         If (DefLA(1,iLA).LE.0) Then
            Call WarningMessage(2,'LA definition: index of LA atom < 1')
            Call Quit_OnUserError()
         Else If (DefLA(2,iLA).LE.0) Then
            Call WarningMessage(2,'LA definition: index of QM atom < 1')
            Call Quit_OnUserError()
         Else If (DefLA(3,iLA).LE.0) Then
            Call WarningMessage(2,'LA definition: index of MM atom < 1')
            Call Quit_OnUserError()
         Else If (FactLA(iLA).LE.Zero.Or.FactLA(iLA).GE.One) Then
            Call WarningMessage(2,'LA definition: bad scaling factor')
            Call Quit_OnUserError()
         End If
      End Do
      Call Put_iArray('LA Def',DefLA,3*nLA)
      Call Put_dArray('LA Fact',FactLA,nLA)
      Call mma_deallocate(DefLA)
      Call mma_deallocate(FactLA)
#else
      Message = 'Interface to Gromacs not installed'
      Call WarningMessage(2,Message)
      Call Quit_OnUserError()
#endif
      Go To 998
*                                                                      *
****** EMFR ************************************************************
*                                                                      *
 8035 GWinput = .True.
      Kword = Get_Ln(LuRd)
      Call Upcase(KWord)
      EMFR=.True.
      Call Get_F(1,KVector,3)
      Temp=Sqrt(KVector(1)**2+KVector(2)**2+KVector(3)**2)
      KVector(1)=KVector(1)/Temp
      KVector(2)=KVector(2)/Temp
      KVector(3)=KVector(3)/Temp
*     Get the wavelength in atomic units.
      Call Get_F1(4,Lambda)
      If (Index(KWord,'ANGSTROM').ne.0) Lambda  = Lambda/angstr
      If (Index(KWord,'NANOMETER').ne.0) Then
         Lambda  = Ten*Lambda/angstr
      ENd If
      KVector(1)=((Two*Pi)/Lambda)*KVector(1)
      KVector(2)=((Two*Pi)/Lambda)*KVector(2)
      KVector(3)=((Two*Pi)/Lambda)*KVector(3)
      Go To 998
*                                                                      *
****** NOCD ************************************************************
*                                                                      *
 9084 GWinput = .True.
      If (.NOT.CholeskyWasSet) Then
         Do_RI=.False.
         iRI_Type=0
         Cholesky=.False.
         CholeskyWasSet=.True.
      End If
      Go To 998
*                                                                      *
****** FNMC ************************************************************
*                                                                      *
 9086 GWinput = .True.
      FNMC=.True.
      Go To 998
*                                                                      *
****** ISOT ************************************************************
*                                                                      *
 7654 GWinput = .True.
      KWord = Get_Ln(LuRd)
      Call Upcase(KWord)
      Call Get_I1(1,nIsotopes)
      Call mma_allocate(nIsot,nIsotopes,2)
      Call mma_allocate(mIsot,nIsotopes)
      Do i=1,nIsotopes
         KWord = Get_Ln(LuRd)
         Call Upcase(KWord)
         Call Get_I1(1,iAt)
         nIsot(i,1)=iAt
         If (Index(KWord,'DALTON').ne.0) Then
            Call Get_F1(2,dMass)
            nIsot(i,2) = -1
            mIsot(i) = dMass*UToAU
         Else
            Call Get_I1(2,iIso)
            nIsot(i,2) = iIso
            mIsot(i) = -One
         End If
      End Do
      Go To 998
*                                                                      *
****** EFP  ************************************************************
*                                                                      *
 9088 GWinput = .True.
      Kword = Get_Ln(LuRd)
      Call Get_I1(1,nEFP_fragments)
      Allocate(FRAG_TYPE(nEFP_fragments))
      Allocate(ABC(3,nEFP_fragments))
      Kword = Get_Ln(LuRd)
      Call Upcase(kWord)
      If (KWord.eq.'XYZABC') Then
         Coor_Type=XYZABC_type
         nEFP_Coor=6
         Allocate(EFP_COORS(nEFP_Coor,nEFP_fragments))
         Write (LuWr,*) 'XYZABC option to be implemented'
         Call Abend()
      Else If (KWord.eq.'POINTS') Then
         Coor_Type=POINTS_type
         nEFP_Coor=9
         Allocate(EFP_COORS(nEFP_Coor,nEFP_fragments))
         Do iFrag = 1, nEFP_fragments
            KWord = Get_Ln(LuRd)
            FRAG_Type(iFrag)=KWord
            Do i = 1, 3
               KWord = Get_Ln(LuRd)
               jend=Index(KWord,' ')
               If (jEnd.gt.LENIN+1) Then
                  Write (LuWr,*) 'Warning: the label ', KWord(1:jEnd),
     &                        ' will be truncated to ',LENIN,
     &                        ' characters!'
               End If
               ABC(i,iFrag) = KWord(1:Min(LENIN,jend-1))
               Call Get_F(2,EFP_COORS((i-1)*3+1,iFrag),3)
            End Do
         End Do
      Else If (KWord.eq.'ROTMAT') Then
         Coor_Type=ROTMAT_type
         nEFP_Coor=12
         Allocate(EFP_COORS(nEFP_Coor,nEFP_fragments))
         Write (LuWr,*) 'ROTMAT option to be implemented'
         Call Abend()
      Else
         Write (LuWr,*) 'Illegal EFP format :',KWord
         Write (LuWr,*)
         Write (LuWr,*) 'Allowed format: XYZABC,'
         Write (LuWr,*) '                POINTS, and'
         Write (LuWr,*) '                ROTMAT'
      End If
      lEFP=.True.
      Go To 998
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     P O S T   P R O C E S S I N G                                    *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
 997  Continue
c     Postprocessing for COORD
c      ik=index(KeepBasis,'....')
c      if(ik.ne.0) then
c        KeepBasis=KeepBasis(1:ik-1)
c      endif
      If (CoordSet) Then
         Do ik=len(KeepBasis),1,-1
            If (KeepBasis(ik:ik).ne.' '.and.KeepBasis(ik:ik).ne.'.')
     &         Go To 1997
         End Do
1997     Continue
         KeepBasis=KeepBasis(1:ik)
#ifdef _HAVE_EXTRA_
         Call ProcessXYZ(BasisSet, KeepBasis, KeepGroup,iBSSE,
     &                   SymThr,isHold,ScaleFactor,HyperParSet,
     &                   isXfield)
#else
         Call Parse_Basis(KeepBasis)
         Call Parse_Group(KeepGroup, SymThr)
         Call Write_SewInp('COORD',[iBSSE])
#endif
         if(writeZmat) then
#ifdef _HAVE_EXTRA_
            stepFactor = stepFac1/(hypParam(1)*hypParam(1))
            Call Geo_Setup_Drv(ishold,oldZMat,zConstraints,
     &                         geoInput,hypParam,nFragment,iOptimType,
     &                         stepFactor,gradLim)
#else
            Call WarningMessage(2,'molcas-extra not installed')
            Call Quit_OnUserError()
#endif
         end if
         DoneCoord=.true.
         if(isXfield.eq.1) then
            LuRd_saved=LuRd
            filename='findsym.xfield'
            lXF=.True.
            goto 9753
         endif
      endif
*
9755  continue
      If (CoordSet) Then
         CoordSet=.false.
         LuRdSave=LuRd
         LuFS=IsFreeUnit(1)
#ifdef _HAVE_EXTRA_
         Call Molcas_Open(LuFS,'FS.std')
#else
         Call Molcas_Open(LuFS,'COORD')
#endif
         LuRd=LuFS
         GWInput=.True.
         Go To 998
      Else
         If (DoneCoord) Then
           Close(LuFS)
           LuRd=LuRdSave
#ifndef _HAVE_EXTRA_
           Call Clear_XYZ()
#endif
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Put_lScalar('CSPF',CSPF)
      Call Put_cArray('Align_Weights',Align_Weights,512)
*                                                                      *
************************************************************************
*                                                                      *
*     Isotopic specifications
*
      If (.not.Allocated(nIsot)) Call mma_allocate(nIsot,0,2)

      If (Run_Mode.ne.S_Mode) Then
*        Loop over unique centers
         iUnique = 0
         Call mma_allocate(Isotopes,nCnttp,Label='Isotopes')
         Do iCnttp = 1, nCnttp
            nCnt = dbsc(iCnttp)%nCntr
            Do iCnt = 1, nCnt
               iUnique = iUnique+1
*              Get the mass for this center
               dm = rMass(dbsc(iCnttp)%AtmNr)
               Do j = 1, Size(nIsot, 1)
                  If (nIsot(j,1).eq.iUnique) Then
                     If (nIsot(j,2).ge.0) Then
                        dm = rMassx(dbsc(iCnttp)%AtmNr,nIsot(j,2))
                     Else
                        dm = mIsot(j)
                     End If
                     Exit
                  End If
               End Do
               If (iCnt.eq.1) Then
                  dbsc(iCnttp)%CntMass = dm
               Else
                  If (dm.ne.dbsc(iCnttp)%CntMass) Then
                     Call WarningMessage(2,
     &                 'Error: All centers of the same type must '//
     &                 'have the same mass')
                     Call Quit_OnUserError()
                  End If
               End If
            End Do
            Isotopes(iCnttp)=dbsc(iCnttp)%CntMass
         End Do ! iCnttp
         Call Put_dArray('Isotopes',Isotopes,nCnttp)
         Call mma_deallocate(Isotopes)

*        Find errors
         Do j = 1, Size(nIsot, 1)
            If (nIsot(j,1).gt.iUnique) Then
               Call WarningMessage(2,
     &           'Error: Isotope specification index larger than the '//
     &           'number of unique centers')
               Call Quit_OnUserError()
            End If
         End Do
      End If

*     Deallocate
      If (Allocated(nIsot)) Call mma_deallocate(nIsot)
      If (Allocated(mIsot)) Call mma_deallocate(mIsot)
*                                                                      *
************************************************************************
*                                                                      *
**    post-processing for RP-Coord
*
      If (lRP.and.RPset) Call processRP(KeepGroup,SymThr)
*
**
*
      lAMFI=lAMFI .and. .Not. NoAMFI
*
*     Disable the RI flag if only one-electron integrals are
*     requested
*
      Do_RI = .Not.Onenly .and. Do_RI
      If (Do_RI) Then
         If (LDF .and. LocalDF) Then
            Call WarningMessage(2,'LDF and LocalDF are incompatible')
            Call Quit_OnUserError()
         End If
      End If
*
      iPrint = nPrint(iRout)
*
      S%Mx_Shll = iShll + 1
      Max_Shells=S%Mx_Shll
*
      If (nCnttp.eq.0) then
         Call WarningMessage(2,'Input does not contain any basis sets')
         Call Quit_OnUserError()
      End If
      If (mdc.eq.0) Then
         Call WarningMessage(2,'Input does not contain coordinates')
         Call Quit_OnUserError()
      End If
      If (S%iAngMx.lt.0) Then
         Call WarningMessage(2,
     &     ' There is an error somewhere in the input!;S%iAngMx.lt.0')
         Call Quit_OnUserError()
      End If
      If (S%iAngMx.gt.iTabMx) Then
         Call WarningMessage(2,' Too High angular momentum !!!')
         Call Quit_OnUserError()
      End If
      If (DoTinker.and.DoGromacs) Then
         Call WarningMessage(2,
     &      'TINKER and GROMACS keywords cannot be used together')
         Call Quit_OnUserError()
      End If
      If (DoTinker.and.iCoord.gt.1) Then
         If (.Not.Expert) Then
            Call WarningMessage(2,
     &         'TINKER and COORD keywords cannot be used together')
            Call Quit_OnUserError()
         End If
      End If
      If (DoGromacs.and.iCoord.gt.1) Then
         If (.Not.Expert) Then
            Call WarningMessage(2,
     &         'GROMACS and COORD keywords cannot be used together')
            Call Quit_OnUserError()
         End If
      End If
*
      If (Test) Then
         Do_GuessOrb=.False.
         Do_FckInt=.False.
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Automatic DK and AMFI option for relativistic basis sets
*                                                                      *
************************************************************************
*                                                                      *
      IF (BSS.AND..Not.DKroll) Then
         Call WarningMessage(2,
     &           ';BSSM GOES ALWAYS WITH DOUGLAS.'//
     &           'THE OPPOSITE IS NOT TRUE')
         Call Abend()
      End If
*
      lECP = .False.
      lPP  = .False.
      Do i = 1, nCnttp
         lECP = lECP .or. dbsc(i)%ECP
         lPP  = lPP  .or. dbsc(i)%nPP.ne.0
      End Do
      If ((lECP.or.lPP).and.DKroll.and..Not.Expert) Then
         Call WarningMessage(2,
     &               ' ECP option not compatible with Douglas-Kroll!')
         Call Quit_OnUserError()
      End If
*
      If (imix.eq.1) Then
         if (Expert) then
           Call WarningMessage(1,
     &        ' input is inconsistent!;'
     &      //'SEWARD found basis sets of mixed relativistic'
     &      //' (or non-relativistic) types!;'
     &      //'No relativistic option will be automatically enabled')
         else
           Call WarningMessage(2,
     &        ' input is inconsistent!;'
     &      //'SEWARD found basis sets of mixed relativistic'
     &      //' (or non-relativistic) types!')
           Call Quit_OnUserError()
         endif
      End If
      If (ifnr.eq.1.and..not.Expert) Then
         If (DKroll) Then
         Call WarningMessage(1,
     *    ';you requested the DK-option for;'
     *   //'a non-relativistic basis.;'
     *   //'This request will be ignored')
         End If
         DKroll=.False.
      Else If (ifnr.eq.0) Then
         lAMFI=.True. .and. .NOT. NoAMFI
         If (.Not.DKroll) Then
            DKroll=.True.
C           If (iRELAE.eq.-1) IRELAE=201022
            If (iRELAE.eq.-1) Then
               If (itype.eq.2) Then
                  IRELAE=  1022
               Else If (itype.eq.14) Then
                  IRELAE=   101
               End If
            End If
         End If
         If (MolWgh.ne.0 .and. MolWgh.ne.2) MolWgh=2
      End If
*
      If (NoDKroll) DKroll=.false.
#ifdef _HAVE_EXTRA_
      If (DoEMPC) isHold=1
#endif
*
      If ((lECP.or.lPP).and.lAMFI.and..Not.Expert) Then
         Call WarningMessage(2,
     &               ' ECP option not compatible with AMFI!')
         Call Quit_OnUserError()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Activate Finite Nucleus parameters
*
      If (Nuclear_Model.eq.Point_Charge) Then
         If (ign.eq.2) Nuclear_Model=Gaussian_Type
         If (ign.eq.3) Nuclear_Model=mGaussian_Type
      End If
*
      Do iCnttp = 1, nCnttp
         If (Nuclear_Model.eq.Gaussian_Type) Then
*
*           If ExpNuc not explicitly defined use default value.
*
            nMass = nInt(dbsc(iCnttp)%CntMass/UToAU)
            If (dbsc(iCnttp)%ExpNuc.lt.Zero)
     &          dbsc(iCnttp)%ExpNuc=NucExp(nMass)
         Else If (Nuclear_Model.eq.mGaussian_Type) Then
*
*           Get parameters for the Modified Gaussian Nuclear
*           charge distribution.
*
            jAtmNr=dbsc(iCnttp)%AtmNr
            nMass = nInt(dbsc(iCnttp)%CntMass/UToAU)
            Call ModGauss(DBLE(jAtmNr),nMass,dbsc(iCnttp)%ExpNuc,
     &                    dbsc(iCnttp)%w_mGauss)
*
         Else
*
*           Nothing to do for point charges!
*
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Cholesky-specific postprocessing:
*     0) if 1C-CD is requested, do it.
*     1) reset integral prescreening thresholds (if not user-defined).
*     2) use default Cholesky normalization (if not user-defined);
*        for AMFI or Douglas-Kroll, use Molcas normalization (again,
*        if not user-defined).
*     3) Turn off Onenly flag, which might be set through the Direct
*        keyword. Thus, specifying Cholesky will force 2-el. int.
*        processing even with Direct specifed as well!
*     4) Turn off Dist flag (makes no sense with Cholesky).
*     5) Integral format on disk is irrelevant, except that Aces II is
*        not allowed. So, reset iWrOpt or quit (for Aces II).
*     6) if Cholesky threshold is specified, use it.
*     7) if span factor is specified, use it.
*
      If (do1CCD) Then
         If (.not.Cholesky) Then
            DirInt = .True.
            Call Cho_Inp(.True.,-1,6)
         End If
         Cholesky=.True.
         Call Cho_InpMod('1CCD')
      End If
      If (Cholesky) Then
         If (Onenly) Then
            Cholesky=.false. ! we gotta be lazy in such cases
         Else
            If (.not. CutInt_UsrDef) CutInt = Cho_CutInt
            If (.not. ThrInt_UsrDef) ThrInt = Cho_ThrInt
            If (.not. MolWgh_UsrDef) Then
               If (lAMFI .or. DKroll) Then
                  MolWgh = 2
               Else
                  MolWgh = Cho_MolWgh
               End If
            End If
            If (iWrOpt .eq. 2) Then
               Write(LuWr,*)
     &         'Acess II format not allowed with Cholesky!!'
               Call Quit_OnUserError()
            Else If (iWrOpt.ne.0 .and. iWrOpt.ne.3) Then
               iWrOpt = 0
            End If
            If (CholeskyThr.ge.0.0d0) Then
               Thrshld_CD=CholeskyThr
               Call Cho_SetDecompositionThreshold(Thrshld_CD)
               Call Put_Thr_Cho(Thrshld_CD)
            End If
            If (spanCD.ge.0.0d0) Then
               v=min(spanCD,1.0d0)
               Call Cho_SetSpan(v)
            End If
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (Prprt) Then
         Onenly = .True.
         Vlct   = .False.
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Post processing for FAIEMP fragment data
*
      If (lFAIEMP.and.Run_Mode.ne.S_Mode) Call FragExpand(LuRd)
*                                                                      *
************************************************************************
*                                                                      *
*     Post processing for RI and RI/CD option
*
      If ( Do_RI .and. Run_Mode.ne.S_Mode ) Then
         If (iRI_Type.eq.4) Then
*
*           Generate on-the-fly aCD or aTrue.cCD auxiliary basis sets.
*
            Call Mk_RICD_Shells()
*
         Else
*
*           Pick up an externally defined auxiliary basis set.
*
            Call Mk_RI_Shells(LuRd)
*
         End If
      End If
      If (Do_RI .and. LocalDF .and. Run_Mode.eq.S_Mode) Then
         Call SetTargetAccuracy_LDF()
         If (CutInt_UsrDef .and. .not.APThr_UsrDef) Then
            Call LDF_SetPrescreen(CutInt)
            Call LDF_CheckThrs()
         End If
         Call LDF_CheckConfig()
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Post processing for Well integrals
*
      Do iWel = 1, nWel
         If (Wel_Info(1,iWel).lt.Zero) Then
            If (.Not.lRF) Then
               Call WarningMessage(2,
     &                        '; Input inconsistency!; ;'
     &                      //'Relative positions of well integrals'
     &                      //' can only be used if the cavity radius'
     &                      //' has been specified!')
               Call Quit_OnUserError()
            End If
            Wel_Info(1,iWel)=rds+Abs(Wel_Info(1,iWel))
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Put up list for point at which the electric field will be
*     evaluated. If nEF=0 the default points will be the unique
*     centers.
*
      If (nOrdEF.ge.0.and. .NOT.(Run_Mode.eq.S_Mode)) Then
         If (nEF.ne.0) Then
            Call mma_allocate(EF_Centers,3,nEF,Label='EF_Centers')
            EF_Centers(:,:) = EFt(:,:)
            Call mma_deallocate(EFt)
         Else
            nEF = 0
            Do iCnttp = 1, nCnttp
               If (.NOT.dbsc(iCnttp)%Aux .and. .NOT.dbsc(iCnttp)%Frag)
     &         nEF = nEF + dbsc(iCnttp)%nCntr
            End Do
            Call mma_allocate(EF_Centers,3,nEF,Label='EF_Centers')
*
            iEF = 1
            Do iCnttp = 1, nCnttp
               If (.NOT.dbsc(iCnttp)%Aux .and.
     &             .NOT.dbsc(iCnttp)%Frag) Then
                  call dcopy_(3*dbsc(iCnttp)%nCntr,
     &                                        dbsc(iCnttp)%Coor,1,
     &                                        EF_Centers(1,iEF),1)
                  iEF = iEF + dbsc(iCnttp)%nCntr
               End If
            End Do
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Put up list for point at which the diamagnetic shielding will
*     be evaluated. If nDMS=0 the default points will be the unique
*     centers.
*
      If (lDMS.and. .NOT.(Run_Mode.eq.S_Mode)) Then
         If (nDMS.ne.0) Then
            Call mma_allocate(DMS_Centers,3,nDMS,Label='DMS_Centers')
            DMS_Centers(:,:)=DMSt(:,:)
            call mma_deallocate(DMSt)
         Else
            nDMS = 0
            Do iCnttp = 1, nCnttp
               nDMS = nDMS + dbsc(iCnttp)%nCntr
            End Do
            Call mma_allocate(DMS_Centers,3,nDMS,Label='DMS_Centers')
            iDMS = 1
            Do iCnttp = 1, nCnttp
               call dcopy_(3*dbsc(iCnttp)%nCntr,
     &                                     dbsc(iCnttp)%Coor,1,
     &                                     DMS_Centers(1,iDMS),1)
               iDMS = iDMS + dbsc(iCnttp)%nCntr
            End Do
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
*     If no multipole moment integrals are requested turn also of the
*     computation of the velocity integrals.
*
      If (S%nMltpl.eq.0) Vlct=.False.
*
*     But turn it on again if explicitly requested
*
      If (Vlct_) Vlct=.True.
*
*     This is the highest order of any property operator.
*     The default value of 4 is due to the mass-velocity operator
*     which is computed by default.
*
      nPrp = Max(4,S%nMltpl)
*
*     Setup of tables for coefficients of the Rys roots and weights.
*
      nDiff=0
      If (S%iAngMx.eq.0) nDiff=2
      DoRys=.True.
      If (DKroll.and.nOrdEF.gt.0) nDiff=nDiff+nOrdEF
      If (.Not.Test.and.Run_Mode.ne.S_Mode) Call SetUp_RW(DoRys,nDiff)
*                                                                      *
************************************************************************
*                                                                      *
*     Store information for the Douglas-Kroll code.
*
      If (DKroll.or.NEMO) Call Fill_rInfo1()
*                                                                      *
************************************************************************
*                                                                      *
*     Compute kOffAO and lOffAO
*
      Call Setup_OffAO()
*                                                                      *
************************************************************************
*                                                                      *
*---- Generate labels for Cartesian and spherical basis sets.
*     Generate the transformation matrix for cartesian to sphericals
*     and contaminants. This has to be done after adding auxiliary or
*     fragment basis sets.
*
      Call Sphere(S%iAngMx)
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
*     Set up Symmetry_Info
*
      Call Symmetry_Info_Setup(nOper,Oper,Max(S%iAngMx,3))

      If (lSkip) then
         Call Put_Ln(ChSkip)
         Call Get_I(1,iSkip,nIrrep)
         Do_GuessOrb=.FALSE.
      End If
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
*     Generate list of Stabilizers , Stabilizer Index
*     and distinct cosets
*
      S%mCentr=0
      S%mCentr_Aux=0
      S%mCentr_Frag=0
      Do iCnttp = 1, nCnttp
         nCnt = dbsc(iCnttp)%nCntr
         Do iCnt = 1, nCnt
            mdc = iCnt + dbsc(iCnttp)%mdci
            S%Mx_mdc = Max(S%Mx_mdc,mdc)
            n_dc=max(mdc,n_dc)
            If (mdc.gt.MxAtom) Then
               Call WarningMessage(2,' mdc.gt.MxAtom!;'
     &                      //' Increase MxAtom in info.fh.')
               Write (LuWr,*) ' MxAtom=',MxAtom
               Call Abend()
            End If
*
*           The symmetry operators of the fragment's atoms should
*           always be identical to that of the fragment's
*           pseudocenter/placeholder
*
            If (dbsc(iCnttp)%Frag) Then
*              Check the FragExpand routine!
               If (Abs(dbsc(iCnttp)%nFragCoor)>mdc) Then
                  Write (6,*) 'rdctl_seward: incorrect mdc index'
                  Call Abend()
               End If
               iChxyz = dc(Abs(dbsc(iCnttp)%nFragCoor))%iChCnt
            Else
*
*------------- To assign the character of a center we need to find
*              the cartesian components that are permutable. We
*              will only need to loop over the generators of the
*              group. We will use the three first bits to indicate if
*              the cartesian component is affected by any symmetry
*              operation.
*
               iChxyz=iChAtm(dbsc(iCnttp)%Coor(:,iCnt))
            End If
            dc(mdc)%iChCnt = iChxyz
            Call Stblz(iChxyz,dc(mdc)%nStab,dc(mdc)%iStab,
     &                 nIrrep,dc(mdc)%iCoSet)
*
*           Perturb the initial geometry if the SHAKE keyword was given,
*           but maintain the symmetry
*
            If ((Shake.gt.Zero).And.
     &          .Not.(dbsc(iCnttp)%pChrg.Or.
     &                dbsc(iCnttp)%Frag.Or.
     &                dbsc(iCnttp)%Aux)) Then
               jTmp=0
               Do j=1,dc(mdc)%nStab-1
                  jTmp=iOr(jTmp,dc(mdc)%iStab(j))
               End Do
               S%nDim=0
               Do j=0,2
                  If (iAnd(jTmp,2**j).eq.0) S%nDim=S%nDim+1
               End Do
               If (S%nDim.gt.0) Then
                  Call Random_Vector(S%nDim,RandVect(1:S%nDim),.False.)
                  jDim=0
                  Do j=0,2
                     If (iAnd(jTmp,2**j).eq.0) Then
                        jDim=jDim+1
                        dbsc(iCnttp)%Coor(j+1,iCnt)=
     &                      dbsc(iCnttp)%Coor(j+1,iCnt)
     &                     +Shake*RandVect(jDim)
                  End If
                  End Do
               End If
            End If
            If (dbsc(iCnttp)%Frag) Then
               S%mCentr_Frag = S%mCentr_Frag + nIrrep/dc(mdc)%nStab
            Else If (dbsc(iCnttp)%Aux) Then
               S%mCentr_Aux = S%mCentr_Aux + nIrrep/dc(mdc)%nStab
            Else
               S%mCentr = S%mCentr + nIrrep/dc(mdc)%nStab
            End If
         End Do
      End Do
      If (S%mCentr.gt.MxAtom) Then
         Call WarningMessage(2,'RdCtl: S%mCentr.gt.MxAtom')
         Write (6,*) 'S%mCentr=',S%mCentr
         Write (6,*) 'Edit src/Include/Molcas.fh'
         Write (6,*) 'Set MxAtom to the value of S%mCentr.'
         Write (6,*) 'Recompile MOLCAS and try again!'
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If ((SymmSet.or.nIrrep.gt.1).and.LDKroll) Then
         Call WarningMessage(2,
     &      'Local DKH approach is not yet implemented with SYMMETRY')
         Call Quit_OnUserError()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (Do_GuessOrb.and.Run_Mode.ne.S_Mode) Call Fix_FockOp(LuRd)
*                                                                      *
************************************************************************
*                                                                      *
      If ((iXPolType.ne.0).and.(nIrrep.ne.1)) Then
         Call WarningMessage(2,
     &                'Polarizabilities are not compatible'
     &              //' with symmetry.')
         Call Quit_OnUserError()
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
      If (nTtl.ne.0.and.Run_Mode.eq.G_Mode) Then
         If (iPrint.ge.6) Then
            Write (LuWr,*)
            Write (LuWr,'(15X,88A)') ('*',i=1,88)
            Write (LuWr,'(15X,88A)') '*', (' ',i=1,86), '*'
            Do iTtl = 1, nTtl
               Write (LuWr,'(15X,A,A,A)') '*   ',Title(iTtl),'   *'
            End Do
            Write (LuWr,'(15X,88A)') '*', (' ',i=1,86), '*'
            Write (LuWr,'(15X,88A)') ('*',i=1,88)
         Else
            Write (LuWr,*)
            Write (LuWr,'(A)') ' Title:'
            Do iTtl = 1, nTtl
               Write (LuWr,'(8X,A)') Title(iTtl)
            End Do
            Write (LuWr,*)
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Process the weights used for alignment and distance measurement
*
      Call Process_Weights(iPrint)
*
************************************************************************
*                                                                      *
*     Set structures for TS optimization according to the Saddle
*     method.
*
      If (Run_Mode.ne.G_Mode) Then
         Call Saddle()
*                                                                      *
************************************************************************
*                                                                      *
*---- Read coordinates from run file (if any), ditto for external
*     field. Do not do this in the Gateway!
*
         Call GeoNew(Show)
         If (lXF) Call GeoNew_PC()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Gen_GeoList()
*                                                                      *
************************************************************************
*                                                                      *
*     Generate list of centers for multipole operators
*
      Call SetMltplCenters()
*
*     Put in user specified centers if any
*
      If (lMltpl) Then
         Do i = 1, nTemp
            iMltpl = ITmp(i)
            If (iMltpl.le.S%nMltpl) call dcopy_(3,RTmp(1,i),1,
     &                                         Coor_MPM(1,iMltpl+1),1)
         End Do
         Call mma_deallocate(RTmp)
         Call mma_deallocate(ITmp)
      End If
#ifdef _DEBUGPRINT_
       Call RecPrt(' Multipole centers',' ',Coor_MPM,3,S%nMltpl+1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Put up list for point at which the orbital angular momentum
*     will be computed.
*
      If (lOAM .and. .NOT.(Run_Mode.eq.S_Mode)) Then
         Call mma_allocate(OAM_Center,3,Label='OAM_Center')
         call dcopy_(3,OAMt,1,OAM_Center,1)
      Else If (.NOT.allocated(OAM_Center)) Then
         Call mma_allocate(OAM_Center,3,Label='OAM_Center')
         call dcopy_(3,CoM,1,OAM_Center,1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Put up list for point at which the orbital magnetic quadrupole
*     will be computed.
*
      If (lOMQ .and. .NOT.(Run_Mode.eq.S_Mode)) Then
         Call mma_allocate(OMQ_Center,3,Label='OMQ_Center')
         Call DCopy_(3,OMQt,1,OMQ_Center,1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocate fields from keyword ORIGIN
*
      If(OrigInput) Then
         Call mma_deallocate(OrigRot)
         Call mma_deallocate(OrigTrans)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (NoZMAT .and. .NOT.ForceZMAT) Call Put_iScalar('N ZMAT',0)
      If (Run_Mode.ne.S_Mode) Call Put_iArray('BasType',BasisTypes,4)
*                                                                      *
************************************************************************
*                                                                      *
      If (Run_Mode.eq.G_Mode)
     &   Call Put_lScalar('Invert constraints',Invert)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Buffer)
*                                                                      *
************************************************************************
*                                                                      *
      Call Datimx(KWord)
      Header(1)=Title(1)(5:76)
      Write (Header(2),'(4A)')
     &          ' Integrals generated by ',
     &            Vrsn,', ',KWord(1:24)
      Call Put_cArray('Seward Title',Header(1),144)
      If (nTtl>0) Call Put_cArray('SewardXTitle',Title(1),nTtl*80)
*                                                                      *
************************************************************************
*                                                                      *
      If (Run_Mode.eq.G_Mode) Call Put_iScalar('DNG',iDNG)
*                                                                      *
************************************************************************
*                                                                      *
*
      Call mma_deallocate(STDINP)
      Return
6666  Call WarningMessage(2,'Unable to read data from '//KWord)
      call Quit_OnUserError()
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine SetTargetAccuracy_LDF()
      Implicit None
#include "localdf.fh"
      If (Thr_Accuracy.lt.0.0d0) Call LDF_SetDefaultThrs()
      Return
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CheckThrs()
      Implicit None
#include "localdf.fh"
      If (Thr_Accuracy.lt.0.0d0) Then
         Call WarningMessage(2,'LDF: Thr_Accuracy<0')
         Call Quit_OnUserError()
      End If
      If (Thr_Prescreen.lt.0.0d0) Then
         Call WarningMessage(2,'LDF: Thr_Prescreen<0')
         Call Quit_OnUserError()
      End If
      Thr_Prescreen=min(Thr_Prescreen,Thr_Accuracy)
      Return
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_SetOptionFlag(Option,Value)
      Implicit None
      Character*4 Option
      Logical Value
#include "localdf.fh"
      If (Option.eq.'LDF2') Then
         LDF2=Value
      Else If (Option.eq.'CHEC') Then
         CheckPairIntegrals=Value
      Else If (Option.eq.'VERI') Then
         VerifyFit=Value
      Else If (Option.eq.'OVER') Then
         CheckOverlapIntegrals=Value
      Else If (Option.eq.'WRUC') Then
         WriteUnconstrainedC=Value
      Else If (Option.eq.'UNIQ') Then
         UseUniqueAtomPairs=Value
      Else
         Call WarningMessage(2,'LDF_SetOptionFlag: unknown Option')
         Write(6,'(A,A)') 'Option=',Option
         Write(6,'(A,L1)') 'Value=',Value
         Call LDF_Quit(1)
      End If
      Return
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CheckConfig()
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit None
#include "localdf.fh"
      ! Debug write of unconstrained coefficients:
      ! 1) makes no sense for unconstrained LDF => reset
      ! 2) not implemented in parallel => error
      WriteUnconstrainedC=WriteUnconstrainedC .and. LDF_Constraint.ne.-1
#ifdef _MOLCAS_MPP_
      If (WriteUnconstrainedC) Then
         If (nProcs.gt.1 .and. Is_Real_Par()) Then
            Call WarningMessage(2,
     &   'Write unconstrained coefficients not implemented in parallel')
            Call Quit_OnUserError()
         End If
      End If
#endif
      ! Using unique atom pairs is buggy, warn!
      If (UseUniqueAtomPairs) Then
         Call WarningMessage(1,
     &   'WARNING: using unique atom pairs may cause erroneous results')
         Call xFlush(6)
      End If
      Return
      End
