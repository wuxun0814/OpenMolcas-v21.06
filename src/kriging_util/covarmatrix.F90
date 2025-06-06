!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************

subroutine covarMatrix()

use kriging_mod, only: eps, eps2, full_R, Index_PGEK, l, m_t, nD, nInter, nInter_Eff, nPoints, x
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i, j, i0, i1, j0, j1, k, i_eff, j_eff
real(kind=wp), allocatable :: diffx_j(:,:), diffx_i(:,:), matFder(:,:), matSder(:,:), r(:,:,:), d(:,:)

!**********************************************************************
!
! Allocate temporary memory

call mma_Allocate(diffx_j,nPoints,nPoints,label='diffx_j')
call mma_Allocate(diffx_i,nPoints,nPoints,label='diffx_i')
call mma_Allocate(matFder,nPoints,nPoints,label='matFder')
call mma_Allocate(matSder,nPoints,nPoints,label='matSder')
call mma_Allocate(r,nPoints,nPoints,nInter,label='r')
call mma_Allocate(d,nPoints,nPoints,label='d')

!**********************************************************************
!
! Compute the distances between the sample points using
! the characteristic length.

full_R(:,:) = Zero
d(:,:) = Zero
diffx_j(:,:) = Zero
diffx_i(:,:) = 0

do i=1,nInter

  do k=1,nPoints
    do j=1,nPoints
      r(k,j,i) = (x(i,k)-x(i,j))/l(i)
    end do
  end do

  ! Accumulate contributions to the square of the individual distances.

  d(:,:) = d(:,:)+r(:,:,i)**2

# ifdef _DEBUGPRINT_
  call RecPrt('r',' ',r(:,:,i),nPoints,nPoints)
# endif

end do

!**********************************************************************

#ifdef _DEBUGPRINT_
call RecPrt('l',' ',l,1,nInter)
call RecPrt('x',' ',x,nInter,nPoints)
call RecPrt('d',' ',d,nPoints,nPoints)
#endif

!**********************************************************************
!**********************************************************************
!
! Now evaluate the covariance function over all the distances. For GEK
! we will need gradients and 2nd order derivatives of the covariance
! function too.
!
! Currently we use the 5/2 Matern function as the covariance function
!
!**********************************************************************
!
! Note that we will evaluate the derivative of the covariance function
! w.r.t d. For the full derivative this has to be complemented by
! the derivative of d w.r.t the individual components of the coordinates.
!
!**********************************************************************
! 1) Evaluate the covariance function for all the distances.

call matern(d,full_R(1:nPoints,1:nPoints),nPoints,nPoints)

! Writing the covariant matrix in GEK (eq 2 of doi:10.1007/s00366-015-0397)
!
!**********************************************************************
!
! 2) Evaluate first derivatives of the covariance function with respect to d at all distances.

call matderiv(1,d,MatFder,nPoints,nPoints)

! Covariant matrix in Gradient Enhanced Kriging (eq 2 of doi:10.1007/s00366-015-0397):
!
! First line and first column derivative in Psi matrix

do i_eff=1,nInter_Eff      ! Loop over component of the coordinate to differentiate
  i = Index_PGEK(i_eff)

  ! Compute the range of the block in the covariance matrix.

  i0 = nPoints+1+(i_eff-1)*(nPoints-nD)
  i1 = i0+(nPoints-nD)-1

  ! Do an on-the-fly evaluation of the dervative w.r.t x_i
  diffx_i(1:nPoints,1+nD:nPoints) = -Two*r(1:nPoints,1+nD:nPoints,i)/l(i)

  ! Writing the 1st row of 1st derivatives with respect the coordinates

  full_R(1:nPoints,i0:i1) = matFDer(1:nPoints,1+nD:nPoints)*diffx_i(1:nPoints,1+nD:nPoints)

end do
! Complete by filling in the transpose blocks

full_R(nPoints+1:m_t,1:nPoints) = transpose(Full_R(1:nPoints,nPoints+1:m_t))

!**********************************************************************
!
! 3) Evaluate the second derivatives.
!
! Matern second derivative with respect to d

call matderiv(2,d,matSder,nPoints,nPoints)

! Second derivatives
do i_Eff=1,nInter_Eff
  i = Index_PGEK(i_Eff)
  i0 = nPoints+1+(i_Eff-1)*(nPoints-nD)
  i1 = i0+(nPoints-nD)-1

  diffx_i(1+nD:nPoints,1+nD:nPoints) = -Two*r(1+nD:nPoints,1+nD:nPoints,i)/l(i)

  do j_Eff=i_Eff,nInter_Eff
    j = Index_PGEK(j_Eff)
    j0 = nPoints+1+(j_Eff-1)*(nPoints-nD)
    j1 = j0+(nPoints-nD)-1

    diffx_j(1+nD:nPoints,1+nD:nPoints) = Two*r(1+nD:nPoints,1+nD:nPoints,j)/l(j)

    ! if differentiating twice on the same dimension
    full_R(i0:i1,j0:j1) = matSder(1+nD:nPoints,1+nD:nPoints)*diffx_j(1+nD:nPoints,1+nD:nPoints)*diffx_i(1+nD:nPoints,1+nD:nPoints)

    if (i == j) full_R(i0:i1,j0:j1) = full_R(i0:i1,j0:j1)-matFder(1+nD:nPoints,1+nD:nPoints)*(Two/(l(i)*l(j)))

    ! Writing the second derivatives in eq(2)
    if (i /= j) full_R(j0:j1,i0:i1) = transpose(Full_r(i0:i1,j0:j1))
  end do

end do

! Add constants to reflect the error in the energy and the
! gradient, respectively.

do j=1,m_t
  if (j <= nPoints) then
    Full_R(j,j) = Full_R(j,j)+eps
  else
    Full_R(j,j) = Full_R(j,j)+eps2
  end if
end do

#ifdef _DEBUGPRINT_
call RecPrt('The covariance matrix:','(12(2x,E9.3))',full_R,m_t,m_t)
#endif

call mma_deallocate(diffx_j)
call mma_deallocate(diffx_i)
call mma_deallocate(matFder)
call mma_deallocate(matSder)
call mma_deallocate(r)
call mma_deallocate(d)

end subroutine covarMatrix
