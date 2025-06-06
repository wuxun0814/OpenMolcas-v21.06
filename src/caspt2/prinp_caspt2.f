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
* Copyright (C) 1994, Markus P. Fuelscher                              *
*               1994, Per Ake Malmqvist                                *
************************************************************************
      Subroutine PrInp_CASPT2
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     - echo the input parameters                                      *
*                                                                      *
*     calling parameters: none                                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and P.-AA. Malmqvist                              *
*     University of Lund, Sweden, 1994                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      use output_caspt2, only:iPrGlb,terse,usual,verbose
      Implicit Real*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
      Character(Len=8)   Fmt1,Fmt2
      Character(Len=120)  Line
      Character(Len=3) lIrrep(8)
      Character(Len=20) calctype,FockOpType
*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      Line=' '
      lLine=Len(Line)
      lPaper=132
      left=(lPaper-lLine)/2
      WRITE(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
      WRITE(Fmt2,'(A,I3.3,A)') '(',left,'X,'
*----------------------------------------------------------------------*
*     Print the ONEINT file identifier                                 *
*----------------------------------------------------------------------*
      IF(IPRGLB.GE.VERBOSE) THEN
      WRITE(6,*)
      WRITE(6,Fmt1) 'Header of the ONEINT file:'
      WRITE(6,Fmt1) '--------------------------'
      WRITE(Line,'(36A2)') (Header(i),i=1,36)
      Call LeftAd(Line)
      WRITE(6,Fmt1)  Line(:mylen(Line))
      WRITE(Line,'(36A2)') (Header(i),i=37,72)
      Call LeftAd(Line)
      WRITE(6,Fmt1)  Line(:mylen(Line))
      WRITE(6,*)
      END IF
*----------------------------------------------------------------------*
*     Print cartesian coordinates of the system                        *
*----------------------------------------------------------------------*
      IF(IPRGLB.GE.VERBOSE) THEN
      Call PrCoor
      END IF
*----------------------------------------------------------------------*
*     Print orbital and wavefunction specifications                    *
*----------------------------------------------------------------------*
      IF(IPRGLB.GE.USUAL) THEN
      WRITE(6,*)
      Line=' '
      WRITE(Line(left-2:),'(A)') 'Wave function specifications:'
      CALL CollapseOutput(1,Line)
      WRITE(6,Fmt1)'-----------------------------'
      WRITE(6,*)
      WRITE(6,Fmt2//'A,T45,I6)')'Number of closed shell electrons',
     &                           2*NISHT
      WRITE(6,Fmt2//'A,T45,I6)')'Number of electrons in active shells',
     &                           NACTEL ! 活性电子数
      WRITE(6,Fmt2//'A,T45,I6)')'Max number of holes in RAS1 space',
     &                           NHOLE1 ! RASSCF RAS1空间中最大空穴数，CASSCF中都是0 
      WRITE(6,Fmt2//'A,T45,I6)')'Max number of electrons in RAS3 space',
     &                           NELE3 ! RASSCF RAS3空间中最大电子数，CASSCF中都是0 
      WRITE(6,Fmt2//'A,T45,I6)')'Number of inactive orbitals',
     &                           NISHT ! 没冻结的非活性轨道的个数
      WRITE(6,Fmt2//'A,T45,I6)')'Number of active orbitals',
     &                           NASHT ! 活性轨道的个数
      WRITE(6,Fmt2//'A,T45,I6)')'Number of secondary orbitals',
     &                           NSSHT ! 虚轨道的个数
      WRITE(6,Fmt2//'A,T45,F6.1)')'Spin quantum number',
     &                           0.5D0*DBLE(ISPIN-1) ! ISPIN是多重度
      WRITE(6,Fmt2//'A,T45,I6)')'State symmetry',
     &                           STSYM
      WRITE(6,Fmt2//'A,T40,I11)')'Number of CSFs',
     &                           NCONF ! CSF的个数
      WRITE(6,Fmt2//'A,T45,I6)')'Number of root(s) available',
     &                           NROOTS ! 根的个数
      WRITE(6,Fmt2//'A,T45,I6)')'Root passed to geometry opt.',
     &                           iRlxRoot
      IF(IFMIX) THEN
        WRITE(6,Fmt2//'A,T45,10I3)')'A file JOBMIX will be created'
      END IF
      IF(NSTATE.GT.1) THEN
        WRITE(6,Fmt1) 'This is a MULTI-STATE CASSCF reference'
        WRITE(6,Fmt2//'A,T45,I6)')'Number of CI roots used',
     &                           NSTATE
        WRITE(6,Fmt2//'A,(T47,10I4))')'These are:',
     &                           (MSTATE(I),I=1,NSTATE)
        IF(IFMSCOUP) THEN
           WRITE(6,Fmt1) 'Off-diagonal elements of Heff are computed'
        ELSE
           WRITE(6,Fmt1) 'Heff is assumed diagonal'
        ENDIF
      ELSE
        If ( ISCF.eq.0 ) then
           WRITE(6,Fmt1) 'This is a CASSCF or RASSCF reference function'
#ifdef _ENABLE_BLOCK_DMRG_
           If (DoCumulant) then
              write(6,Fmt1) 'Using 4-RDM cumulant approximation,' //
     &                      ' activated by 3RDM keyword in RASSCF'
           End If
#elif _ENABLE_CHEMPS2_DMRG_
           If (DoCumulant) then
            write(6,Fmt1) 'This is a DMRG reference with exact 4-RDM,'//
     &                    ' activated by 3RDM keyword in RASSCF'
           End If
#endif
        Else If ( ISCF.eq.1 ) then
           WRITE(6,Fmt1) 'This is a closed shell RHF reference function'
        Else
           WRITE(6,Fmt1)
     &     'This is a high spin open shell RHF reference function'
        End If
      END IF
      CALL CollapseOutput(0,'Wave function specifications:')
      END IF
*
      Call Get_cArray('Irreps',lIrrep,24)
      Do iSym = 1, nSym
         Call RightAd(lIrrep(iSym))
      End Do
*
      IF(IPRGLB.GE.USUAL  ) THEN
      WRITE(6,*)
      Line=' '
      WRITE(Line(left-2:),'(A)') 'Orbital specifications:'
      CALL CollapseOutput(1,Line)
      WRITe(6,Fmt1)'-----------------------'
      WRITE(6,*)
      WRITE(6,Fmt2//'A,T47,8I4)') 'Symmetry species',
     &                            (iSym,iSym=1,nSym)
      WRITE(6,Fmt2//'A,T47,8(1X,A))') '                ',
     &                            (lIrrep(iSym),iSym=1,nSym) ! lIrrep(iSym)是第iSym个对称性的不可约表示
      WRITE(6,Fmt2//'A,T47,8I4)') 'Frozen orbitals',
     &                            (nFro(iSym),iSym=1,nSym) ! nFro(iSym)是第iSym个对称性的冻结轨道的个数
      WRITE(6,Fmt2//'A,T47,8I4)') 'Inactive orbitals',
     &                            (nIsh(iSym),iSym=1,nSym) ! nIsh(iSym)是第iSym个对称性的没冻结的非活性轨道的个数
      WRITE(6,Fmt2//'A,T47,8I4)') 'Active orbitals',
     &                            (nAsh(iSym),iSym=1,nSym) ! nAsh(iSym)是第iSym个对称性的活性轨道个数
      WRITE(6,Fmt2//'A,T47,8I4)') 'Secondary orbitals',
     &                            (nSsh(iSym),iSym=1,nSym) ! nSsh(iSym)是第iSym个对称性的虚轨道个数
      WRITE(6,Fmt2//'A,T47,8I4)') 'Deleted orbitals',
     &                            (nDel(iSym),iSym=1,nSym)
      WRITE(6,Fmt2//'A,T47,8I4)') 'Number of basis functions',
     &                            (nBas(iSym),iSym=1,nSym) ! nBas(iSym)是第iSym个对称性的基函数个数
      CALL CollapseOutput(0,'Orbital specifications:')
      END IF
*----------------------------------------------------------------------*
*     Print routing information                                        *
*----------------------------------------------------------------------*
      IF (IPRGLB.GE.TERSE) THEN
        If ( RFpert ) then
          WRITE(6,*)
          WRITE(6,Fmt1)'Reaction field specifications:'
          WRITE(6,Fmt1)'------------------------------'
          WRITE(6,*)
          WRITE(6,'(6X,A)')'An external reaction field was determined'//
     &        ' previously and added to the one-electron Hamiltonian'
          WRITE(6,'(6X,A)')'It will not be reevaluated even though'//
     &               ' dynamic correlation may change the density.'
          WRITE(6,*)
        End If

        write(6,*)
        Line=' '
        write(Line(left-2:),'(A)') 'CASPT2 specifications:'
        call CollapseOutput(1,Line)
        write(6,Fmt1)'----------------------'
        write(6,*)

        if (IFMSCOUP) then
          if (IFDW) then
            FockOpType='dynamically weighted'
            if (IFXMS) then
              calctype='XDW-CASPT2'
            else
              calctype='DW-CASPT2'
            end if
          else if (IFRMS) then
              FockOpType='state-specific'
              calctype='RMS-CASPT2'
          else
            if (IFXMS) then
              FockOpType='state-average'
              calctype='XMS-CASPT2'
            else
              FockOpType='state-specific'
              calctype='MS-CASPT2'
            end if
          end if
        else
          FockOpType='state-specific'
          calctype='SS-CASPT2'
        end if

        write(6,Fmt2//'A,T50,A)')'Type of calculation',trim(calctype)

        write(6,Fmt2//'A,T50,A)')'Fock operator',trim(FockOpType)
        if (IFDW) then
          write(6,Fmt2//'A,T45,I6)')'DW Type',DWType
          if (zeta.ge.0) then
            write(6,Fmt2//'A,T50,E10.4)')'DW exponent',zeta
          else
            write(6,Fmt2//'A,T50,A)')'DW exponent','infinity'
          end if
        end if

        if (Hzero.ne.'STANDARD') then
          write(6,Fmt2//'A,T50,A)')'0th-order Hamiltonian',trim(Hzero)
        end if

        write(6,Fmt2//'A,T45,F9.2)')'IPEA shift',BSHIFT ! IPEA shift
        write(6,Fmt2//'A,T45,F9.2)')'Real shift',SHIFT ! Real shift
        write(6,Fmt2//'A,T45,F9.2)')'Imaginary shift',SHIFTI ! Imaginary shift

        if (ORBIN.eq.'TRANSFOR') then
          write(6,Fmt1)'The input orbitals will be transformed'//
     &                 ' to quasi-canonical'
        else
          write(6,Fmt1)'The input orbitals will not be transformed'//
     &                 ' to quasi-canonical'
        end if

        if (IFXMS .or. IFRMS) then
          write(6,Fmt1)'The input states will be rotated to '//
     &     'diagonalize the Fock operator'
        end if

        call CollapseOutput(0,'CASPT2 specifications:')
        write(6,*)
      end if

C Compute necessary quantities for subsequent gradient calc:
      IF (IPRGLB.GE.USUAL) THEN
        IF(IFDENS) THEN
        WRITE(6,*)
          WRITE(6,*)' The wave functions (P_CAS) H W |Psi_0> and '//
     &         '(P_CAS) W H |Psi_0>'
          WRITE(6,*)' will be computed and written out for subsequent'
          WRITE(6,*)' use in a gradient calculation.'
        END IF
      END IF

      Return
      End
