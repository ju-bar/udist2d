!**********************************************************************!
!**********************************************************************!
!
! file: "interpol.f90"
!
! author: Juri Barthel
!         juribarthel@gmx.de
!         October 2008
!
! purpose: Supply of interpolation core routines
!          interpolation of real*4, real*8, complex*8 data
!                        in 1d 2d 3d
!                        nearest neighbor, linear, cubic FFT
!
! LINK: SFFTs.f, interpolprm.fi
!
!**********************************************************************!
!**********************************************************************!




!**********************************************************************
!                                                                      
!   Date: 2018-11-15                                                   
!                                                                      
!   Author: Juri Barthel                                               
!           Ernst Ruska-Centre                                         
!           Forschungszentrum Jülich GmbH, 52425 Jülich, Germany       
!                                                                      
!----------------------------------------------------------------------
!                                                                      
! This program is free software: you can redistribute it and/or modify 
! it under the terms of the GNU General Public License as published by 
! the Free Software Foundation, either version 3 of the License, or    
! (at your option) any later version.                                  
!                                                                      
! This program is distributed in the hope that it will be useful,      
! but WITHOUT ANY WARRANTY; without even the implied warranty of       
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        
! GNU General Public License for more details.                         
!                                                                      
! You should have received a copy of the GNU General Public License    
! along with this program. If not, see <http://www.gnu.org/licenses/>. 
!                                                                      
!----------------------------------------------------------------------


















!**********************************************************************!
!**********************************************************************!
!
!
!
! *** CUBIC INTERPOLATION KERNELS
!
!
!
!**********************************************************************!
!**********************************************************************!


!**********************************************************************!
!
! real*4 function cubickernelr4
!
! defines cubic interpolation kernel
!
! INPUT:
!   real*4 :: t             = offset of supporting point from
!                             interpolation position
!   real*4 :: a             = interpolation smoothness parameter
!
real*4 function cubickernelr4(t,a)
  implicit none
  real*4, intent(in) :: t, a
  real*4 :: tabs
  tabs = abs(t)
  cubickernelr4 = 0.0
  if (tabs<1.0) then
    cubickernelr4 = 1.0+tabs*tabs*((a+2.0)*tabs-(a+3.0))
  else if (tabs<2.0 .or. tabs>=1.0) then
    cubickernelr4 = a*(tabs*(8.0-tabs*(5.0-tabs))-4.0)
  end if
  return
end function cubickernelr4

!**********************************************************************!
!
! real*8 function cubickernelr8
!
! defines cubic interpolation kernel
!
! INPUT:
!   real*8 :: t             = offset of supporting point from
!                             interpolation position
!   real*8 :: a             = interpolation smoothness parameter
!
real*8 function cubickernelr8(t,a)
  implicit none
  real*8, intent(in) :: t, a
  real*8 :: tabs
  tabs = dabs(t)
  cubickernelr8 = 0.0D+0
  if (tabs<1.0D+0) then
    cubickernelr8 = 1.0D+0+tabs*tabs*((a+2.0D+0)*tabs-(a+3.0D+0))
  else if (tabs<2.0D+0 .or. tabs>=1.0D+0) then
    cubickernelr8 = a*(tabs*(8.0D+0-tabs*(5.0D+0-tabs))-4.0D+0)
  end if
  return
end function cubickernelr8










































!**********************************************************************!
!**********************************************************************!
!
!
!
! *** 1D INTERPOLATION FUNCTIONS
!
!     return an interpolated value
!
!     real*4, external :: cubicr4wrap
!     real*8, external :: cubicr8wrap
!     complex*8, external :: cubicc8wrap
!     real*4, external :: cubicr4
!     real*8, external :: cubicr8
!     complex*8, external :: cubicc8
!     real*4, external :: linr4wrap
!     real*8, external :: linr8wrap
!     complex*8, external :: linc8wrap
!     real*4, external :: linr4
!     real*8, external :: linr8
!     complex*8, external :: linc8
!
!     real*4, external :: r4interpol1d
!     real*8, external :: r8interpol1d
!     complex*8, external :: c8interpol1d
!     external :: r4interpol1dplanes
!
!**********************************************************************!
!**********************************************************************!


!**********************************************************************!
!
! function cubicr4wrap
!
! returns cubic interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n)          = input data array
!   integer*4 :: n          = number of data samples
!   real*4 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*4 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation uses wrap around
!
real*4 function cubicr4wrap(a,n,bip,x)
  implicit none
  real*4, intent(in) :: a(n), bip, x
  integer*4, intent(in) :: n
  integer*4 :: i                    ! iterators
  integer*4 :: i1                   ! 1st support pixel
  integer*4 :: i2                   ! current pixel index
  real*4 :: f1                      ! positional fractions
  real*4 :: rn, bcx                 ! temp values
  real*4 :: rval                    ! temp values
  real*4 :: mx                      ! wrapped coordinate
  real*4, external :: cubickernelr4
  rval = 0.0
  rn = real(n)
  mx = modulo(x-1.0,rn)+1.0
  i1 = int(mx)
  do i=i1-1, i1+2
    f1 = mx-real(i)
    bcx = cubickernelr4(f1,bip)
    i2 = 1+modulo(i-1,n)
    rval = rval + a(i2)*bcx
  end do
  cubicr4wrap = rval
end function cubicr4wrap


!**********************************************************************!
!
! function cubicr8wrap
!
! returns cubic interpolation value (real*8)
!
! INPUT:
!   real*8 :: a(n)          = input data array
!   integer*4 :: n          = number of data samples
!   real*8 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*8 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation uses wrap around
!
real*8 function cubicr8wrap(a,n,bip,x)
  implicit none
  real*8, intent(in) :: a(n), bip, x
  integer*4, intent(in) :: n
  integer*4 :: i                    ! iterators
  integer*4 :: i1                   ! 1st support pixel
  integer*4 :: i2                   ! current pixel index
  real*8 :: f1                      ! positional fractions
  real*8 :: rn, bcx                 ! temp values
  real*8 :: rval                    ! temp values
  real*8 :: mx                      ! wrapped coordinate
  real*8, external :: cubickernelr8
  rval = 0.0D+0
  rn = dble(n)
  mx = modulo(x-1.0D+0,rn)+1.0D+0
  i1 = int(mx)
  do i=i1-1, i1+2
    f1 = mx-dble(i)
    bcx = cubickernelr8(f1,bip)
    i2 = 1+modulo(i-1,n)
    rval = rval + a(i2)*bcx
  end do
  cubicr8wrap = rval
end function cubicr8wrap


!**********************************************************************!
!
! function cubicc8wrap
!
! returns cubic interpolation value (complex*8)
!
! INPUT:
!   complex*8 :: a(n)       = input data array
!   integer*4 :: n          = number of data samples
!   real*4 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*4 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation uses wrap around
!
complex*8 function cubicc8wrap(a,n,bip,x)
  implicit none
  complex*8, intent(in) :: a(n)
  real*4, intent(in) :: bip, x
  integer*4, intent(in) :: n
  integer*4 :: i                    ! iterators
  integer*4 :: i1                   ! 1st support pixel
  integer*4 :: i2                   ! current pixel index
  real*4 :: f1                      ! positional fractions
  real*4 :: rn, bcx                 ! temp values
  complex*8 :: rval                 ! temp values
  real*4 :: mx                      ! wrapped coordinate
  real*4, external :: cubickernelr4
  rval = cmplx(0.0,0.0)
  rn = real(n)
  mx = modulo(x-1.0,rn)+1.0
  i1 = int(mx)
  do i=i1-1, i1+2
    f1 = mx-real(i)
    bcx = cubickernelr4(f1,bip)
    i2 = 1+modulo(i-1,n)
    rval = rval + a(i2)*bcx
  end do
  cubicc8wrap = rval
end function cubicc8wrap


!**********************************************************************!
!
! function cubicr4
!
! returns cubic interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n)          = input data array
!   integer*4 :: n          = number of data samples
!   real*4 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*4 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation without wrap around
!
real*4 function cubicr4(a,n,bip,x)
  implicit none
  real*4, intent(in) :: a(n), bip, x
  integer*4, intent(in) :: n
  integer*4 :: i                    ! iterators
  integer*4 :: i1                   ! 1st support pixel
  integer*4 :: i2                   ! current pixel index
  real*4 :: f1                      ! positional fractions
  real*4 :: rn, bcx                 ! temp values
  real*4 :: rval                    ! temp values
  real*4 :: mx                      ! wrapped coordinate
  real*4, external :: cubickernelr4
  rval = 0.0
  rn = real(n)
  if ((x>=1.0).and.(x<rn+1.0)) then ! avoid wrap
  mx = modulo(x-1.0,rn)+1.0
  i1 = int(mx)
  do i=i1-1, i1+2
    f1 = mx-real(i)
    bcx = cubickernelr4(f1,bip)
    i2 = 1+modulo(i-1,n)
    rval = rval + a(i2)*bcx
  end do
  end if ! avoid wrap
  cubicr4 = rval
end function cubicr4


!**********************************************************************!
!
! function cubicr8
!
! returns cubic interpolation value (real*8)
!
! INPUT:
!   real*8 :: a(n)          = input data array
!   integer*4 :: n          = number of data samples
!   real*8 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*8 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation without wrap around
!
real*8 function cubicr8(a,n,bip,x)
  implicit none
  real*8, intent(in) :: a(n), bip, x
  integer*4, intent(in) :: n
  integer*4 :: i                    ! iterators
  integer*4 :: i1                   ! 1st support pixel
  integer*4 :: i2                   ! current pixel index
  real*8 :: f1                      ! positional fractions
  real*8 :: rn, bcx                 ! temp values
  real*8 :: rval                    ! temp values
  real*8 :: mx                      ! wrapped coordinate
  real*8, external :: cubickernelr8
  rval = 0.0D+0
  rn = dble(n)
  if ((x>=1.0D+0).and.(x<rn+1.0D+0)) then ! avoid wrap
  mx = modulo(x-1.0D+0,rn)+1.0D+0
  i1 = int(mx)
  do i=i1-1, i1+2
    f1 = mx-dble(i)
    bcx = cubickernelr8(f1,bip)
    i2 = 1+modulo(i-1,n)
    rval = rval + a(i2)*bcx
  end do
  end if ! avoid wrap
  cubicr8 = rval
end function cubicr8


!**********************************************************************!
!
! function cubicc8
!
! returns cubic interpolation value (complex*8)
!
! INPUT:
!   complex*8 :: a(n)       = input data array
!   integer*4 :: n          = number of data samples
!   real*4 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*4 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation without wrap around
!
complex*8 function cubicc8(a,n,bip,x)
  implicit none
  complex*8, intent(in) :: a(n)
  real*4, intent(in) :: bip, x
  integer*4, intent(in) :: n
  integer*4 :: i                    ! iterators
  integer*4 :: i1                   ! 1st support pixel
  integer*4 :: i2                   ! current pixel index
  real*4 :: f1                      ! positional fractions
  real*4 :: rn, bcx                 ! temp values
  complex*8 :: rval                 ! temp values
  real*4 :: mx                      ! wrapped coordinate
  real*4, external :: cubickernelr4
  rval = cmplx(0.0,0.0)
  rn = real(n)
  if ((x>=1.0).and.(x<rn+1.0)) then ! avoid wrap
  mx = modulo(x-1.0,rn)+1.0
  i1 = int(mx)
  do i=i1-1, i1+2
    f1 = mx-real(i)
    bcx = cubickernelr4(f1,bip)
    i2 = 1+modulo(i-1,n)
    rval = rval + a(i2)*bcx
  end do
  end if ! avoid wrap
  cubicc8 = rval
end function cubicc8


!**********************************************************************!
!
! function linr4wrap
!
! returns linear interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n)          = input data array
!   integer*4 :: n          = number of data samples
!   real*4 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation uses wrap around
!
real*4 function linr4wrap(a,n,x)
  implicit none
  real*4, intent(in) :: a(n)
  real*4 :: x
  integer*4, intent(in) :: n
  integer*4 :: i1, i2               ! surrounding pixels
  real*4 :: f1                      ! positional fractions
  real*4 :: rn                      ! temp value
  real*4 :: rval                    ! temp value
  rval = 0.0
  rn = real(n)
  i1 = int(modulo(x-1.0,rn))+1
  i2 = int(modulo(x,rn))+1
  f1 = modulo(x-1.0,rn)+1.0-i1
  rval = rval + (1.0-f1)*a(i1)
  rval = rval + f1      *a(i2)
  linr4wrap = rval
end function linr4wrap


!**********************************************************************!
!
! function linr8wrap
!
! returns linear interpolation value (real*8)
!
! INPUT:
!   real*8 :: a(n)          = input data array
!   integer*4 :: n          = number of data samples
!   real*8 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation uses wrap around
!
real*8 function linr8wrap(a,n,x)
  implicit none
  real*8, intent(in) :: a(n)
  real*8 :: x
  integer*4, intent(in) :: n
  integer*4 :: i1, i2               ! surrounding pixels
  real*8 :: f1                      ! positional fractions
  real*8 :: rn                      ! temp value
  real*8 :: rval                    ! temp value
  rval = 0.0D+0
  rn = dble(n)
  i1 = int(modulo(x-1.0D+0,rn))+1
  i2 = int(modulo(x,rn))+1
  f1 = modulo(x-1.0D+0,rn)+1.0D+0-i1
  rval = rval + (1.0D+0-f1)*a(i1)
  rval = rval + f1         *a(i2)
  linr8wrap = rval
end function linr8wrap


!**********************************************************************!
!
! function linc8wrap
!
! returns linear interpolation value (complex*8)
!
! INPUT:
!   complex*8 :: a(n)       = input data array
!   integer*4 :: n          = number of data samples
!   real*4 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation uses wrap around
!
complex*8 function linc8wrap(a,n,x)
  implicit none
  complex*8, intent(in) :: a(n)
  real*4 :: x
  integer*4, intent(in) :: n
  integer*4 :: i1, i2               ! surrounding pixels
  real*4 :: f1                      ! positional fractions
  real*4 :: rn                      ! temp value
  complex*8 :: rval                 ! temp value
  rval = cmplx(0.0,0.0)
  rn = real(n)
  i1 = int(modulo(x-1.0,rn))+1
  i2 = int(modulo(x,rn))+1
  f1 = modulo(x-1.0,rn)+1.0-i1
  rval = rval + (1.0-f1)*a(i1)
  rval = rval + f1      *a(i2)
  linc8wrap = rval
end function linc8wrap


!**********************************************************************!
!
! function linr4
!
! returns linear interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n)          = input data array
!   integer*4 :: n          = number of data samples
!   real*4 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation without wrap around
!
real*4 function linr4(a,n,x)
  implicit none
  real*4, intent(in) :: a(n)
  real*4 :: x
  integer*4, intent(in) :: n
  integer*4 :: i1, i2               ! surrounding pixels
  real*4 :: f1                      ! positional fractions
  real*4 :: rn                      ! temp value
  real*4 :: rval                    ! temp value
  rval = 0.0
  rn = real(n)
  if ((x>=1.0).and.(x<rn+1.0)) then ! avoid wrap
  i1 = int(modulo(x-1.0,rn))+1
  i2 = int(modulo(x,rn))+1
  f1 = modulo(x-1.0,rn)+1.0-i1
  rval = rval + (1.0-f1)*a(i1)
  rval = rval + f1      *a(i2)
  end if ! avoid wrap
  linr4 = rval
end function linr4


!**********************************************************************!
!
! function linr8
!
! returns linear interpolation value (real*8)
!
! INPUT:
!   real*8 :: a(n)          = input data array
!   integer*4 :: n          = number of data samples
!   real*8 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation without wrap around
!
real*8 function linr8(a,n,x)
  implicit none
  real*8, intent(in) :: a(n)
  real*8 :: x
  integer*4, intent(in) :: n
  integer*4 :: i1, i2               ! surrounding pixels
  real*8 :: f1                      ! positional fractions
  real*8 :: rn                      ! temp value
  real*8 :: rval                    ! temp value
  rval = 0.0D+0
  rn = dble(n)
  if ((x>=1.0D+0).and.(x<rn+1.0D+0)) then ! avoid wrap
  i1 = int(modulo(x-1.0D+0,rn))+1
  i2 = int(modulo(x,rn))+1
  f1 = modulo(x-1.0D+0,rn)+1.0D+0-i1
  rval = rval + (1.0D+0-f1)*a(i1)
  rval = rval + f1         *a(i2)
  end if ! avoid wrap
  linr8 = rval
end function linr8


!**********************************************************************!
!
! function linc8
!
! returns linear interpolation value (complex*8)
!
! INPUT:
!   complex*8 :: a(n)       = input data array
!   integer*4 :: n          = number of data samples
!   real*4 :: x             = interpolation position
!
! REMARKS:
!   position indices run from 1 to n
!   interpolation without wrap around
!
complex*8 function linc8(a,n,x)
  implicit none
  complex*8, intent(in) :: a(n)
  real*4 :: x
  integer*4, intent(in) :: n
  integer*4 :: i1, i2               ! surrounding pixels
  real*4 :: f1                      ! positional fractions
  real*4 :: rn                      ! temp value
  complex*8 :: rval                 ! temp value
  rval = cmplx(0.0,0.0)
  rn = real(n)
  if ((x>=1.0).and.(x<rn+1.0)) then ! avoid wrap
  i1 = int(modulo(x-1.0,rn))+1
  i2 = int(modulo(x,rn))+1
  f1 = modulo(x-1.0,rn)+1.0-i1
  rval = rval + (1.0-f1)*a(i1)
  rval = rval + f1      *a(i2)
  end if ! avoid wrap
  linc8 = rval
end function linc8


!**********************************************************************!
!
! function r4interpol1d
!
! returns interpolation value (real*4) on 1d arrays
! allows selection of interpolation method
! allows selection of wrap around option
!
! INPUT:
!   real*4 :: a(n)          = input data array
!   integer*4 :: n          = number of data samples
!   real*4 :: x             = interpolation position
!   integer*4 :: im         = interpolation method
!                             0 = nearest neighbour
!                             1 = linear
!                             2 = cubic
!   integer*4 :: wa         = wrap around flag
!                             0 = without wrap around, zeroes returned offsides
!                             1 = with wrap around
!
! REMARKS:
!   position indices run from 1 to n
!   cubic interpolation kernel parameter is defined in include 'interpolprm.fi'
!
real*4 function r4interpol1d(a, n, p, im, wa)

  implicit none
  
  include 'interpolprm.fi'
  
  real*4, dimension(n), intent(in) :: a
  integer*4, intent(in) :: n, im, wa
  real*4, intent(in) :: p
  
  integer*4 :: i
  real*4 :: bip, xm
  real*4 :: rn
  real*4 :: rval
  real*4, external :: linr4, linr4wrap, cubicr4, cubicr4wrap
  
  rval = 0.0
  bip = prm_interpol_cubic
  rn = real(n)
  
  ! distinguish interpolation method
  select case (im)
  
  case (0) ! nearest neighbor
    xm = modulo(p-1.0,rn)+1.0
    if ( ((p>=0.5).and.(p<rn+0.5)).or.(wa==1) ) then ! either (inside array) OR (wraparound asked)
      i = nint(xm)
      rval = a(i)
    end if
  
  case (1) ! linear interpolation
    if (wa==0) then
      rval = linr4(a,n,p)
    else
      rval = linr4wrap(a,n,p)
    end if
  
  case (2) ! cubic interpolation
    if (wa==0) then
      rval = cubicr4(a,n,bip,p)
    else
      rval = cubicr4wrap(a,n,bip,p)
    end if
  
  end select ! case (im)
  
  r4interpol1d = rval
  
  return
  
end function r4interpol1d


!**********************************************************************!
!
! function r8interpol1d
!
! returns interpolation value (real*8) on 1d arrays
! allows selection of interpolation method
! allows selection of wrap around option
!
! INPUT:
!   real*8 :: a(n)          = input data array
!   integer*4 :: n          = number of data samples
!   real*8 :: x             = interpolation position
!   integer*4 :: im         = interpolation method
!                             0 = nearest neighbour
!                             1 = linear
!                             2 = cubic
!   integer*4 :: wa         = wrap around flag
!                             0 = without wrap around, zeroes returned offsides
!                             1 = with wrap around
!
! REMARKS:
!   position indices run from 1 to n
!   cubic interpolation kernel parameter is defined in include 'interpolprm.fi'
!
real*8 function r8interpol1d(a, n, p, im, wa)

  implicit none
  
  include 'interpolprm.fi'
  
  real*8, dimension(n), intent(in) :: a
  integer*4, intent(in) :: n, im, wa
  real*8, intent(in) :: p
  
  integer*4 :: i
  real*8 :: bip, xm
  real*8 :: rn
  real*8 :: rval
  real*8, external :: linr8, linr8wrap, cubicr8, cubicr8wrap
  
  rval = 0.0D+0
  bip = dble(prm_interpol_cubic)
  rn = dble(n)
  
  ! distinguish interpolation method
  select case (im)
  
  case (0) ! nearest neighbor
    xm = modulo(p-1.0D+0,rn)+1.0D+0
    if ( ((p>=0.5D+0).and.(p<rn+0.5D+0)).or.(wa==1) ) then ! either (inside array) OR (wraparound asked)
      i = nint(xm)
      rval = a(i)
    end if
  
  case (1) ! linear interpolation
    if (wa==0) then
      rval = linr8(a,n,p)
    else
      rval = linr8wrap(a,n,p)
    end if
  
  case (2) ! cubic interpolation
    if (wa==0) then
      rval = cubicr8(a,n,bip,p)
    else
      rval = cubicr8wrap(a,n,bip,p)
    end if
  
  end select ! case (im)
  
  r8interpol1d = rval
  
  return
  
end function r8interpol1d


!**********************************************************************!
!
! function c8interpol1d
!
! returns interpolation value (real*4) on 1d arrays
! allows selection of interpolation method
! allows selection of wrap around option
!
! INPUT:
!   complex*8 :: a(n)       = input data array
!   integer*4 :: n          = number of data samples
!   real*4 :: x             = interpolation position
!   integer*4 :: im         = interpolation method
!                             0 = nearest neighbour
!                             1 = linear
!                             2 = cubic
!   integer*4 :: wa         = wrap around flag
!                             0 = without wrap around, zeroes returned offsides
!                             1 = with wrap around
!
! REMARKS:
!   position indices run from 1 to n
!   cubic interpolation kernel parameter is defined in include 'interpolprm.fi'
!
complex*8 function c8interpol1d(a, n, p, im, wa)

  implicit none
  
  include 'interpolprm.fi'
  
  complex*8, dimension(n), intent(in) :: a
  integer*4, intent(in) :: n, im, wa
  real*4, intent(in) :: p
  
  integer*4 :: i
  real*4 :: bip, xm
  real*4 :: rn
  complex*8 :: rval
  complex*8, external :: linc8, linc8wrap, cubicc8, cubicc8wrap
  
  rval = cmplx(0.0,0.0)
  bip = prm_interpol_cubic
  rn = real(n)
  
  ! distinguish interpolation method
  select case (im)
  
  case (0) ! nearest neighbor
    xm = modulo(p-1.0,rn)+1.0
    if ( ((p>=0.5).and.(p<rn+0.5)).or.(wa==1) ) then ! either (inside array) OR (wraparound asked)
      i = nint(xm)
      rval = a(i)
    end if
  
  case (1) ! linear interpolation
    if (wa==0) then
      rval = linc8(a,n,p)
    else
      rval = linc8wrap(a,n,p)
    end if
  
  case (2) ! cubic interpolation
    if (wa==0) then
      rval = cubicc8(a,n,bip,p)
    else
      rval = cubicc8wrap(a,n,bip,p)
    end if
  
  end select ! case (im)
  
  c8interpol1d = rval
  
  return
  
end function c8interpol1d



!**********************************************************************!
!
! function r4interpol1dplanes
!
! interpolation between x-y planes of 3d (real*4) arrays
! allows selection of interpolation method
! allows selection of wrap around option
!
! INPUT:
!   real*4 :: a(n,m,p)      = input data array
!   real*4 :: b(n,m)        = output data array
!   integer*4 :: n, m, p    = number of data samples
!   real*4 :: z             = interpolation position
!   integer*4 :: im         = interpolation method
!                             0 = nearest neighbour
!                             1 = linear
!                             2 = cubic
!   integer*4 :: wa         = wrap around flag
!                             0 = without wrap around, zeroes returned offsides
!                             1 = with wrap around
!
! REMARKS:
!   position indices run from 1,1,1 to n,m,p
!   cubic interpolation kernel parameter is defined in include 'interpolprm.fi'
!   KEEP IN MIND FOR LATER UPDATES:
!   The interpolation is implemented inline again, no calls of other routines
!
subroutine r4interpol1dplanes(a, b, n, m, p, z, im, wa)

  implicit none
  
  include 'interpolprm.fi'
  
  integer*4, intent(in) :: n, m, p, im, wa
  real*4, dimension(n,m,p), intent(in) :: a
  real*4, dimension(n,m), intent(out) :: b
  real*4, intent(in) :: z
  
  integer*4 :: i, i1, i2
  real*4 :: bip, mz, bcz
  real*4 :: rp, f1
  real*4 :: rval
  real*4, external :: cubickernelr4
  
  rval = 0.0
  ! b = 0.0 ! commented out in favour of computation speed
  bip = prm_interpol_cubic
  rp = real(p)
  mz = modulo(z-1.0,rp)+1.0
  
  ! handle wrap-around
  ! if interpolation position is outsides the input array index range of dimension 3
  ! 1) nothing is done if (wa==0), i.e. wrap-around is deactivated
  ! 2) coordinates are wrapped around if (wa==1)
  if (wa==1) then ! (wraparound asked)
  
  ! distinguish interpolation method
  select case (im)
  
    case (0) ! nearest neighbor
      i = nint(mz)
      b(1:n,1:m) = a(1:n,1:m,i)
  
    case (1) ! linear interpolation
      i1 = int(mz)
      i2 = modulo(i1,p)+1
      f1 = mz-i1
      b(1:n,1:m) = b(1:n,1:m) + (1.0-f1)*a(1:n,1:m,i1)
      b(1:n,1:m) = b(1:n,1:m) + f1      *a(1:n,1:m,i2)
  
    case (2) ! cubic interpolation
  
      i1 = int(mz)
      do i=i1-1, i1+2
        f1 = mz-real(i)
        bcz = cubickernelr4(f1,bip)
        i2 = 1+modulo(i-1,p)
        b(1:n,1:m) = b(1:n,1:m) + a(1:n,1:m,i2)*bcz
      end do
  
  end select ! case (im)
  
  else ! avoid wrap, also internally
  
  ! distinguish interpolation method
  select case (im)
  
    case (0) ! nearest neighbor
      i = nint(mz)
      if (i>=1.and.i<=p) b(1:n,1:m) = a(1:n,1:m,i)
  
    case (1) ! linear interpolation
      i1 = int(mz)
      i2 = i1+1
      f1 = mz-i1
      if (i1>=1.and.i1<=p) b(1:n,1:m) = b(1:n,1:m) + (1.0-f1)*a(1:n,1:m,i1)
      if (i2>=1.and.i2<=p) b(1:n,1:m) = b(1:n,1:m) + f1      *a(1:n,1:m,i2)
  
    case (2) ! cubic interpolation
      i1 = int(mz)
      do i=i1-1, i1+2
        f1 = mz-real(i)
        bcz = cubickernelr4(f1,bip)
        if (i>=1.and.i<=p) b(1:n,1:m) = b(1:n,1:m) + a(1:n,1:m,i)*bcz
      end do
  
  end select ! case (im)
  
  end if ! avoid wrap
  
  return
  
end subroutine r4interpol1dplanes









































!**********************************************************************!
!**********************************************************************!
!
!
!
! *** 2D INTERPOLATION FUNCTIONS
! *** WITH AND WITHOUT WRAP-AROUND
!
!     return an interpolated value
!
!     real*4, external :: bicubicr4wrap
!     real*8, external :: bicubicr8wrap
!     complex*8, external :: bicubicc8wrap
!     real*4, external :: bicubicr4
!     real*8, external :: bicubicr8
!     complex*8, external :: bicubicc8
!     real*4, external :: bilinr4wrap
!     real*8, external :: bilinr8wrap
!     complex*8, external :: bilinc8wrap
!     real*4, external :: bilinr4
!     real*8, external :: bilinr8
!     complex*8, external :: bilinc8
!
!     real*4, external :: r4interpol2d
!     real*4, external :: r8interpol2d
!     real*4, external :: c8interpol2d
!
!**********************************************************************!
!**********************************************************************!


!**********************************************************************!
!
! function bicubicr4wrap
!
! returns bicubic interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*4 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation uses wrap around
!
real*4 function bicubicr4wrap(a,n,m,bip,x,y)
  implicit none
  real*4, intent(in) :: a(n,m), bip, x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i, j                 ! iterators
  integer*4 :: i1, j1               ! 1st support pixel
  integer*4 :: i2, j2               ! current pixel index
  real*4 :: f1, f2                  ! positional fractions
  real*4 :: rn, rm, rval, bcx, bcy  ! temp values
  real*4 :: mx, my                  ! wrapped coordinates
  real*4, external :: cubickernelr4
  rval = 0.0
  rn = real(n)
  rm = real(m)
  mx = modulo(x-1.0,rn)+1.0
  my = modulo(y-1.0,rm)+1.0
  i1 = int(mx)
  j1 = int(my)
  do j=j1-1, j1+2
    f2 = my-real(j)
    bcy = cubickernelr4(f2,bip)
    j2 = 1+modulo(j-1,m)
    do i=i1-1, i1+2
      f1 = mx-real(i)
      bcx = cubickernelr4(f1,bip)
      i2 = 1+modulo(i-1,n)
      rval = rval + a(i2,j2)*bcx*bcy
    end do
  end do
  bicubicr4wrap = rval
end function bicubicr4wrap


!**********************************************************************!
!
! function bicubicr4
!
! returns bicubic interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*4 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation does not use wrap around
!
real*4 function bicubicr4(a,n,m,bip,x,y)
  implicit none
  real*4, intent(in) :: a(n,m), bip, x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i, j                 ! iterators
  integer*4 :: i1, j1               ! 1st support pixel
  integer*4 :: i2, j2               ! current pixel index
  real*4 :: f1, f2                  ! positional fractions
  real*4 :: rn, rm, rval, bcx, bcy  ! temp values
  real*4, external :: cubickernelr4
  rval = 0.0
  rn = real(n)
  rm = real(m)
  if ((x>=1.0).and.(x<rn+1.0).and.(y>=1.0).and.(y<rm+1.0)) then ! wrap avoid
  i1 = int(modulo(x-1.0,rn))+1
  j1 = int(modulo(y-1.0,rm))+1
  do j=j1-1, j1+2
    f2 = y-real(j)
    bcy = cubickernelr4(f2,bip)
    j2 = 1+modulo(j-1,m)
    do i=i1-1, i1+2
      f1 = x-real(i)
      bcx = cubickernelr4(f1,bip)
      i2 = 1+modulo(i-1,n)
      rval = rval + a(i2,j2)*bcx*bcy
    end do
  end do
  end if ! rwap avoid
  bicubicr4 = rval
end function bicubicr4


!**********************************************************************!
!
! function bicubicr8wrap
!
! returns bicubic interpolation value (real*8)
!
! INPUT:
!   real*8 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*8 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*8 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation uses wrap around
!
real*8 function bicubicr8wrap(a,n,m,bip,x,y)
  implicit none
  real*8, intent(in) :: a(n,m), bip, x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i, j                 ! iterators
  integer*4 :: i1, j1               ! 1st support pixel
  integer*4 :: i2, j2               ! current pixel index
  real*8 :: f1, f2                  ! positional fractions
  real*8 :: rn, rm, rval, bcx, bcy  ! temp values
  real*8 :: mx, my                  ! wrapped coordinates
  real*8, external :: cubickernelr8
  rval = 0.0
  rn = dble(n)
  rm = dble(m)
  mx = modulo(x-1.0D+0,rn)+1.0D+0
  my = modulo(y-1.0D+0,rm)+1.0D+0
  i1 = int(mx)
  j1 = int(my)
  do j=j1-1, j1+2
    f2 = my-dble(j)
    bcy = cubickernelr8(f2,bip)
    j2 = 1+modulo(j-1,m)
    do i=i1-1, i1+2
      f1 = mx-dble(i)
      bcx = cubickernelr8(f1,bip)
      i2 = 1+modulo(i-1,n)
      rval = rval + a(i2,j2)*bcx*bcy
    end do
  end do
  bicubicr8wrap = rval
end function bicubicr8wrap


!**********************************************************************!
!
! function bicubicc8wrap
!
! returns bicubic interpolation value (complex*8)
!
! INPUT:
!   complex*8 :: a(n,m)     = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*4 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation uses wrap around
!
complex*8 function bicubicc8wrap(a,n,m,bip,x,y)
  implicit none
  complex*8, intent(in) :: a(n,m)
  real*4 :: bip, x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i, j                 ! iterators
  integer*4 :: i1, j1               ! 1st support pixel
  integer*4 :: i2, j2               ! current pixel index
  real*4 :: f1, f2                  ! positional fractions
  real*4 :: rn, rm, bcx, bcy        ! temp values
  real*4 :: mx, my                  ! wrapped coordinates
  complex*8 :: cval                 ! temp values
  real*4, external :: cubickernelr4
  cval = cmplx(0.0,0.0)
  rn = real(n)
  rm = real(m)
  mx = modulo(x-1.0,rn)+1.0
  my = modulo(y-1.0,rm)+1.0
  i1 = int(mx)
  j1 = int(my)
  do j=j1-1, j1+2
    f2 = my-real(j)
    bcy = cubickernelr4(f2,bip)
    j2 = 1+modulo(j-1,m)
    do i=i1-1, i1+2
      f1 = mx-real(i)
      bcx = cubickernelr4(f1,bip)
      i2 = 1+modulo(i-1,n)
      cval = cval + a(i2,j2)*bcx*bcy
    end do
  end do
  bicubicc8wrap = cval
end function bicubicc8wrap





!**********************************************************************!
!
! function bicubicr8
!
! returns bicubic interpolation value (real*8)
!
! INPUT:
!   real*8 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*8 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*8 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation does not use wrap around
!
real*8 function bicubicr8(a,n,m,bip,x,y)
  implicit none
  real*8, intent(in) :: a(n,m), bip, x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i, j                 ! iterators
  integer*4 :: i1, j1               ! 1st support pixel
  integer*4 :: i2, j2               ! current pixel index
  real*8 :: f1, f2                  ! positional fractions
  real*8 :: rn, rm, rval, bcx, bcy  ! temp values
  real*8, external :: cubickernelr8
  rval = 0.0
  rn = dble(n)
  rm = dble(m)
  i1 = int(modulo(x-1.0D+0,rn))+1
  j1 = int(modulo(y-1.0D+0,rm))+1
  if ((x>=1.0D+0).and.(x<rn+1.0D+0).and.(y>=1.0D+0).and.(y<rm+1.0D+0)) then ! wrap avoid
  do j=j1-1, j1+2
    f2 = y-dble(j)
    bcy = cubickernelr8(f2,bip)
    j2 = 1+modulo(j-1,m)
    do i=i1-1, i1+2
      f1 = x-dble(i)
      bcx = cubickernelr8(f1,bip)
      i2 = 1+modulo(i-1,n)
      rval = rval + a(i2,j2)*bcx*bcy
    end do
  end do
  end if ! wrap avoid
  bicubicr8 = rval
end function bicubicr8


!**********************************************************************!
!
! function bicubicc8
!
! returns bicubic interpolation value (complex*8)
!
! INPUT:
!   complex*8 :: a(n,m)     = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: bip           = bicubic interpolation parameter (-0.5 ... -1.0)
!   real*4 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation does not use wrap around
!
complex*8 function bicubicc8(a,n,m,bip,x,y)
  implicit none
  complex*8, intent(in) :: a(n,m)
  real*4 :: bip, x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i, j                 ! iterators
  integer*4 :: i1, j1               ! 1st support pixel
  integer*4 :: i2, j2               ! current pixel index
  real*4 :: f1, f2                  ! positional fractions
  real*4 :: rn, rm, bcx, bcy        ! temp values
  complex*8 :: cval                 ! temp values
  real*4, external :: cubickernelr4
  cval = cmplx(0.0,0.0)
  rn = real(n)
  rm = real(m)
  i1 = int(modulo(x-1.0,rn))+1
  j1 = int(modulo(y-1.0,rm))+1
  if ((x>=1.0).and.(x<rn+1.0).and.(y>=1.0).and.(y<rm+1.0)) then ! wrap avoid
  do j=j1-1, j1+2
    f2 = y-real(j)
    bcy = cubickernelr4(f2,bip)
    j2 = 1+modulo(j-1,m)
    do i=i1-1, i1+2
      f1 = x-real(i)
      bcx = cubickernelr4(f1,bip)
      i2 = 1+modulo(i-1,n)
      cval = cval + a(i2,j2)*bcx*bcy
    end do
  end do
  end if ! wrap avoid
  bicubicc8 = cval
end function bicubicc8


!**********************************************************************!
!
! function bilinr4wrap
!
! returns bilinear interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation uses wrap around
!
real*4 function bilinr4wrap(a,n,m,x,y)
  implicit none
  real*4, intent(in) :: a(n,m), x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i1, j1, i2, j2       ! surrounding pixels
  real*4 :: f1, f2                  ! positional fractions
  real*4 :: rn, rm, rval            ! temp values
  rval = 0.0
  rn = real(n)
  rm = real(m)
  i1 = int(modulo(x-1.0,rn))+1
  j1 = int(modulo(y-1.0,rm))+1
  i2 = int(modulo(x,rn))+1
  j2 = int(modulo(y,rm))+1
  f1 = modulo(x-1.0,rn)+1.0-i1
  f2 = modulo(y-1.0,rm)+1.0-j1
  rval = rval + (1.0-f1)*(1.0-f2)*a(i1,j1)
  rval = rval + f1      *(1.0-f2)*a(i2,j1)
  rval = rval + (1.0-f1)*f2      *a(i1,j2)
  rval = rval + f1      *f2      *a(i2,j2)
  bilinr4wrap = rval
end function bilinr4wrap


!**********************************************************************!
!
! function bilinr8wrap
!
! returns bilinear interpolation value (real*8)
!
! INPUT:
!   real*8 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*8 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation uses wrap around
!
real*8 function bilinr8wrap(a,n,m,x,y)
  implicit none
  real*8, intent(in) :: a(n,m), x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i1, j1, i2, j2       ! surrounding pixels
  real*8 :: f1, f2                  ! positional fractions
  real*8 :: rn, rm, rval            ! temp values
  rval = 0.0
  rn = dble(n)
  rm = dble(m)
  i1 = int(modulo(x-1.0D+0,rn))+1
  j1 = int(modulo(y-1.0D+0,rm))+1
  i2 = int(modulo(x,rn))+1
  j2 = int(modulo(y,rm))+1
  f1 = modulo(x-1.0D+0,rn)+1.0D+0-i1
  f2 = modulo(y-1.0D+0,rm)+1.0D+0-j1
  rval = rval + (1.0D+0-f1)*(1.0D+0-f2)*a(i1,j1)
  rval = rval + f1         *(1.0D+0-f2)*a(i2,j1)
  rval = rval + (1.0D+0-f1)*f2      *a(i1,j2)
  rval = rval + f1         *f2      *a(i2,j2)
  bilinr8wrap = rval
end function bilinr8wrap


!**********************************************************************!
!
! function bilinc8wrap
!
! returns bilinear interpolation value (complex*8)
!
! INPUT:
!   complex*8 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation uses wrap around
!
complex*8 function bilinc8wrap(a,n,m,x,y)
  implicit none
  complex*8, intent(in) :: a(n,m)
  real*4, intent(in) :: x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i1, j1, i2, j2       ! surrounding pixels
  real*4 :: f1, f2                  ! positional fractions
  real*4 :: rn, rm                  ! temp values
  complex*8 :: cval                 ! temp value
  cval = cmplx(0.0,0.0)
  rn = real(n)
  rm = real(m)
  i1 = int(modulo(x-1.0,rn))+1
  j1 = int(modulo(y-1.0,rm))+1
  i2 = int(modulo(x,rn))+1
  j2 = int(modulo(y,rm))+1
  f1 = modulo(x-1.0,rn)+1.0-i1
  f2 = modulo(y-1.0,rm)+1.0-j1
  cval = cval + (1.0-f1)*(1.0-f2)*a(i1,j1)
  cval = cval + f1      *(1.0-f2)*a(i2,j1)
  cval = cval + (1.0-f1)*f2      *a(i1,j2)
  cval = cval + f1      *f2      *a(i2,j2)
  bilinc8wrap = cval
end function bilinc8wrap



!**********************************************************************!
!
! function bilinr4
!
! returns bilinear interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation does not use wrap around
!
real*4 function bilinr4(a,n,m,x,y)
  implicit none
  real*4, intent(in) :: a(n,m), x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i1, j1, i2, j2       ! surrounding pixels
  real*4 :: f1, f2                  ! positional fractions
  real*4 :: rn, rm, rval            ! temp values
  rval = 0.0
  rn = real(n)
  rm = real(m)
  if ((x>=1.0).and.(x<rn+1.0).and.(y>=1.0).and.(y<rm+1.0)) then ! wrap avoid
  i1 = int(modulo(x-1.0,rn))+1
  j1 = int(modulo(y-1.0,rm))+1
  i2 = int(modulo(x,rn))+1
  j2 = int(modulo(y,rm))+1
  f1 = modulo(x-1.0,rn)+1.0-i1
  f2 = modulo(y-1.0,rm)+1.0-j1
  rval = rval + (1.0-f1)*(1.0-f2)*a(i1,j1)
  rval = rval + f1      *(1.0-f2)*a(i2,j1)
  rval = rval + (1.0-f1)*f2      *a(i1,j2)
  rval = rval + f1      *f2      *a(i2,j2)
  end if ! wrap avoid
  bilinr4 = rval
end function bilinr4


!**********************************************************************!
!
! function bilinr8
!
! returns bilinear interpolation value (real*8)
!
! INPUT:
!   real*8 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*8 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation uses wrap around
!
real*8 function bilinr8(a,n,m,x,y)
  implicit none
  real*8, intent(in) :: a(n,m), x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i1, j1, i2, j2       ! surrounding pixels
  real*8 :: f1, f2                  ! positional fractions
  real*8 :: rn, rm, rval            ! temp values
  rval = 0.0
  rn = dble(n)
  rm = dble(m)
  if ((x>=1.0D+0).and.(x<rn+1.0D+0).and.(y>=1.0D+0).and.(y<rm+1.0D+0)) then ! wrap avoid
  i1 = int(modulo(x-1.0D+0,rn))+1
  j1 = int(modulo(y-1.0D+0,rm))+1
  i2 = int(modulo(x,rn))+1
  j2 = int(modulo(y,rm))+1
  f1 = modulo(x-1.0D+0,rn)+1.0D+0-i1
  f2 = modulo(y-1.0D+0,rm)+1.0D+0-j1
  rval = rval + (1.0D+0-f1)*(1.0D+0-f2)*a(i1,j1)
  rval = rval + f1         *(1.0D+0-f2)*a(i2,j1)
  rval = rval + (1.0D+0-f1)*f2         *a(i1,j2)
  rval = rval + f1         *f2         *a(i2,j2)
  end if
  bilinr8 = rval
end function bilinr8


!**********************************************************************!
!
! function bilinc8
!
! returns bilinear interpolation value (complex*8)
!
! INPUT:
!   complex*8 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: x, y          = interpolation position
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   interpolation does not use wrap around
!
complex*8 function bilinc8(a,n,m,x,y)
  implicit none
  complex*8, intent(in) :: a(n,m)
  real*4, intent(in) :: x, y
  integer*4, intent(in) :: n, m
  integer*4 :: i1, j1, i2, j2       ! surrounding pixels
  real*4 :: f1, f2                  ! positional fractions
  real*4 :: rn, rm                  ! temp values
  complex*8 :: cval                 ! temp value
  cval = cmplx(0.0,0.0)
  rn = real(n)
  rm = real(m)
  if ((x>=1.0).and.(x<rn+1.0).and.(y>=1.0).and.(y<rm+1.0)) then ! wrap avoid
  i1 = int(modulo(x-1.0,rn))+1
  j1 = int(modulo(y-1.0,rm))+1
  i2 = int(modulo(x,rn))+1
  j2 = int(modulo(y,rm))+1
  f1 = modulo(x-1.0,rn)+1.0-i1
  f2 = modulo(y-1.0,rm)+1.0-j1
  cval = cval + (1.0-f1)*(1.0-f2)*a(i1,j1)
  cval = cval + f1      *(1.0-f2)*a(i2,j1)
  cval = cval + (1.0-f1)*f2      *a(i1,j2)
  cval = cval + f1      *f2      *a(i2,j2)
  end if ! wrap avoid
  bilinc8 = cval
end function bilinc8




!**********************************************************************!
!
! function r4interpol2d
!
! returns interpolation value (real*4) on 2d arrays
! allows selection of interpolation method
! allows selection of wrap around option
!
! INPUT:
!   real*4 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: x, y          = interpolation position
!   integer*4 :: im         = interpolation method
!                             0 = nearest neighbour
!                             1 = linear
!                             2 = cubic
!   integer*4 :: wa         = wrap around flag
!                             0 = without wrap around, zeroes returned offsides
!                             1 = with wrap around
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   cubic interpolation kernel parameter is defined in include 'interpolprm.fi'
!
real*4 function r4interpol2d(a, n, m, x, y, im, wa)

  implicit none
  
  include 'interpolprm.fi'
  
  real*4, dimension(n,m), intent(in) :: a
  integer*4, intent(in) :: n, m, im, wa
  real*4, intent(in) :: x, y
  
  integer*4 :: i, j
  real*4 :: bip, xm, ym
  real*4 :: rn, rm
  real*4 :: rval
  real*4, external :: bilinr4, bilinr4wrap, bicubicr4, bicubicr4wrap
  
  rval = 0.0
  bip = prm_interpol_cubic
  rn = real(n)
  rm = real(m)
  
  ! distinguish interpolation method
  select case (im)
  
  case (0) ! nearest neighbor
    xm = modulo(x-1.0,rn)+1.0
    ym = modulo(y-1.0,rm)+1.0
    if ( ((x>=0.5).and.(x<rn+0.5).and.(y>=0.5).and.(y<rm+0.5)).or. &
     &   ( wa==1 ) ) then ! either (inside array) OR (wraparound asked)
      i = nint(xm)
      j = nint(ym)
      rval = a(i,j)
    end if
  
  case (1) ! linear interpolation
    if (wa==0) then
      rval = bilinr4(a,n,m,x,y)
    else
      rval = bilinr4wrap(a,n,m,x,y)
    end if
  
  case (2) ! cubic interpolation
    if (wa==0) then
      rval = bicubicr4(a,n,m,bip,x,y)
    else
      rval = bicubicr4wrap(a,n,m,bip,x,y)
    end if
  
  end select ! case (im)
  
  r4interpol2d = rval
  
  return
  
end function r4interpol2d



!**********************************************************************!
!
! function r8interpol2d
!
! returns interpolation value (real*8) on 2d arrays
! allows selection of interpolation method
! allows selection of wrap around option
!
! INPUT:
!   real*8 :: a(n,m)        = input data array
!   integer*4 :: n, m       = number of data samples
!   real*8 :: x, y          = interpolation position
!   integer*4 :: im         = interpolation method
!                             0 = nearest neighbour
!                             1 = linear
!                             2 = cubic
!   integer*4 :: wa         = wrap around flag
!                             0 = without wrap around, zeroes returned offsides
!                             1 = with wrap around
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   cubic interpolation kernel parameter is defined in include 'interpolprm.fi'
!
real*8 function r8interpol2d(a, n, m, x, y, im, wa)

  implicit none
  
  include 'interpolprm.fi'
  
  real*8, dimension(n,m), intent(in) :: a
  integer*4, intent(in) :: n, m, im, wa
  real*8, intent(in) :: x, y
  
  integer*4 :: i, j
  real*8 :: bip, xm, ym
  real*8 :: rn, rm
  real*8 :: rval
  real*8, external :: bilinr8, bilinr8wrap, bicubicr8, bicubicr8wrap
  
  rval = 0.0D+0
  bip = dble(prm_interpol_cubic)
  rn = dble(n)
  rm = dble(m)
  
  ! distinguish interpolation method
  select case (im)
  
  case (0) ! nearest neighbor
    xm = modulo(x-1.0D+0,rn)+1.0D+0
    ym = modulo(y-1.0D+0,rm)+1.0D+0
    if ( ((x>=0.5D+0).and.(x<rn+0.5D+0).and.(y>=0.5D+0).and.(y<rm+0.5D+0)).or. &
     &   ( wa==1 ) ) then ! either (inside array) OR (wraparound asked)
      i = nint(xm)
      j = nint(ym)
      rval = a(i,j)
    end if
  
  case (1) ! linear interpolation
    if (wa==0) then
      rval = bilinr8(a,n,m,x,y)
    else
      rval = bilinr8wrap(a,n,m,x,y)
    end if
  
  case (2) ! cubic interpolation
    if (wa==0) then
      rval = bicubicr8(a,n,m,bip,x,y)
    else
      rval = bicubicr8wrap(a,n,m,bip,x,y)
    end if
  
  end select ! case (im)
  
  r8interpol2d = rval
  
  return
  
end function r8interpol2d



!**********************************************************************!
!
! function c8interpol2d
!
! returns interpolation value (complex*8) on 2d arrays
! allows selection of interpolation method
! allows selection of wrap around option
!
! INPUT:
!   complex*8 :: a(n,m)     = input data array
!   integer*4 :: n, m       = number of data samples
!   real*4 :: x, y          = interpolation position
!   integer*4 :: im         = interpolation method
!                             0 = nearest neighbour
!                             1 = linear
!                             2 = cubic
!   integer*4 :: wa         = wrap around flag
!                             0 = without wrap around, zeroes returned offsides
!                             1 = with wrap around
!
! REMARKS:
!   position indices run from 1,1 to n,m
!   cubic interpolation kernel parameter is defined in include 'interpolprm.fi'
!
complex*8 function c8interpol2d(a, n, m, x, y, im, wa)

  implicit none
  
  include 'interpolprm.fi'
  
  complex*8, dimension(n,m), intent(in) :: a
  integer*4, intent(in) :: n, m, im, wa
  real*4, intent(in) :: x, y
  
  integer*4 :: i, j
  real*4 :: bip, xm, ym
  real*4 :: rn, rm
  complex*8 :: rval
  complex*8, external :: bilinc8, bilinc8wrap, bicubicc8, bicubicc8wrap
  
  rval = cmplx(0.0,0.0)
  bip = prm_interpol_cubic
  rn = real(n)
  rm = real(m)
  
  ! distinguish interpolation method
  select case (im)
  
  case (0) ! nearest neighbor
    xm = modulo(x-1.0,rn)+1.0
    ym = modulo(y-1.0,rm)+1.0
    if ( ((x>=0.5).and.(x<rn+0.5).and.(y>=0.5).and.(y<rm+0.5)).or. &
     &   ( wa==1 ) ) then ! either (inside array) OR (wraparound asked)
      i = nint(xm)
      j = nint(ym)
      rval = a(i,j)
    end if
  
  case (1) ! linear interpolation
    if (wa==0) then
      rval = bilinc8(a,n,m,x,y)
    else
      rval = bilinc8wrap(a,n,m,x,y)
    end if
  
  case (2) ! cubic interpolation
    if (wa==0) then
      rval = bicubicc8(a,n,m,bip,x,y)
    else
      rval = bicubicc8wrap(a,n,m,bip,x,y)
    end if
  
  end select ! case (im)
  
  c8interpol2d = rval
  
  return
  
end function c8interpol2d







































!**********************************************************************!
!**********************************************************************!
!
!
!
! *** 3D INTERPOLATION FUNCTIONS
!
!     return an interpolated value
!
!     real*4, external :: tricubicr4wrap
!     real*4, external :: tricubicr4
!     real*4, external :: trilinr4wrap
!     real*4, external :: trilinr4
!     real*4, external :: trinnr4wrap
!     real*4, external :: trinnr4
!
!
!
!**********************************************************************!
!**********************************************************************!

!**********************************************************************!
!
! function tricubicr4wrap
!
! returns tricubic interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m,p)      = input data array
!   integer*4 :: n, m, p    = number of data samples
!   real*4 :: bip           = bicubic interpolation parameter (-0.5 ... -2.0)
!   real*4 :: x, y, z       = interpolation position
!
! REMARKS:
!   position indices run from 1,1,1 to n,m,p
!   interpolation uses wrap around
!
real*4 function tricubicr4wrap(a,n,m,p,bip,x,y,z)
  implicit none
  integer*4, intent(in) :: n, m, p
  real*4, intent(in) :: a(n,m,p), bip, x, y, z
  integer*4 :: i, j, k              ! iterators
  integer*4 :: i1, j1, k1           ! 1st support pixel
  integer*4 :: i2, j2, k2           ! current pixel index
  real*4 :: f1, f2, f3              ! positional fractions
  real*4 :: rn, rm, rp, rval, bcx, bcy, bcz  ! temp values
  real*4 :: mx, my, mz              ! wrapped position
  real*4, external :: cubickernelr4
  rval = 0.0
  rn = real(n)
  rm = real(m)
  rp = real(p)
  mx = modulo(x-1.0,rn)+1.0
  my = modulo(y-1.0,rm)+1.0
  mz = modulo(z-1.0,rp)+1.0
  i1 = int(mx)
  j1 = int(my)
  k1 = int(mz)
  do k=k1-1, k1+2
    f3 = mz-real(k)
    bcz = cubickernelr4(f3,bip)
    k2 = 1+modulo(k-1,p)
    do j=j1-1, j1+2
      f2 = my-real(j)
      bcy = cubickernelr4(f2,bip)
      j2 = 1+modulo(j-1,m)
      do i=i1-1, i1+2
        f1 = mx-real(i)
        bcx = cubickernelr4(f1,bip)
        i2 = 1+modulo(i-1,n)
        rval = rval + a(i2,j2,k2)*bcx*bcy*bcz
      end do
    end do
  end do
  tricubicr4wrap = rval
end function tricubicr4wrap


!**********************************************************************!
!
! function tricubicr4
!
! returns tricubic interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m,p)      = input data array
!   integer*4 :: n, m, p    = number of data samples
!   real*4 :: bip           = bicubic interpolation parameter (-0.5 ... -2.0)
!   real*4 :: x, y, z       = interpolation position
!
! REMARKS:
!   position indices run from 1,1,1 to n,m,p
!   interpolation does not use wrap around
!
real*4 function tricubicr4(a,n,m,p,bip,x,y,z)
  implicit none
  integer*4, intent(in) :: n, m, p
  real*4, intent(in) :: a(n,m,p), bip, x, y, z
  integer*4 :: i, j, k              ! iterators
  integer*4 :: i1, j1, k1           ! 1st support pixel
  integer*4 :: i2, j2, k2           ! current pixel index
  real*4 :: f1, f2, f3              ! positional fractions
  real*4 :: rn, rm, rp, rval, bcx, bcy, bcz  ! temp values
  real*4 :: mx, my, mz              ! wrapped position
  real*4, external :: cubickernelr4
  rval = 0.0
  rn = real(n)
  rm = real(m)
  rp = real(p)
  if ( (x>=1.0).and.(x<rn+1.0).and. &
     & (y>=1.0).and.(y<rm+1.0).and. &
     & (z>=1.0).and.(z<rp+1.0) ) then ! avoid wrap
  mx = modulo(x-1.0,rn)+1.0
  my = modulo(y-1.0,rm)+1.0
  mz = modulo(z-1.0,rp)+1.0
  i1 = int(mx)
  j1 = int(my)
  k1 = int(mz)
  do k=k1-1, k1+2
    f3 = mz-real(k)
    bcz = cubickernelr4(f3,bip)
    k2 = 1+modulo(k-1,p)
    do j=j1-1, j1+2
      f2 = my-real(j)
      bcy = cubickernelr4(f2,bip)
      j2 = 1+modulo(j-1,m)
      do i=i1-1, i1+2
        f1 = mx-real(i)
        bcx = cubickernelr4(f1,bip)
        i2 = 1+modulo(i-1,n)
        rval = rval + a(i2,j2,k2)*bcx*bcy*bcz
      end do
    end do
  end do
  end if ! avoid wrap
  tricubicr4 = rval
end function tricubicr4



!**********************************************************************!
!
! function trilinr4
!
! returns trilinear interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m,p)      = input data array
!   integer*4 :: n, m, p    = number of data samples
!   real*4 :: x, y, z       = interpolation position
!
! REMARKS:
!   position indices run from 1,1,1 to n,m,p
!   interpolation does not use wrap around
!
real*4 function trilinr4(a,n,m,p,x,y,z)
  implicit none
  integer*4, intent(in) :: n, m, p
  real*4, intent(in) :: a(n,m,p), x, y, z
  integer*4 :: i1, j1, k1, i2, j2, k2   ! surrounding pixels
  real*4 :: f1, f2, f3                  ! positional fractions
  real*4 :: mx, my, mz                  ! wrapped positions
  real*4 :: rn, rm, rp, rval            ! temp values
  rval = 0.0
  rn = real(n)
  rm = real(m)
  rp = real(p)
  if ( (x>=1.0).and.(x<rn+1.0).and. &
     & (y>=1.0).and.(y<rm+1.0).and. &
     & (z>=1.0).and.(z<rp+1.0)) then ! wrap avoid
  mx = modulo(x-1.0,rn)+1.0
  my = modulo(y-1.0,rm)+1.0
  mz = modulo(z-1.0,rp)+1.0
  i1 = int(mx)
  j1 = int(my)
  k1 = int(mz)
  i2 = modulo(i1,n)+1
  j2 = modulo(j1,m)+1
  k2 = modulo(k1,p)+1
  f1 = modulo(mx-1.0,rn)+1.0-i1
  f2 = modulo(my-1.0,rm)+1.0-j1
  f3 = modulo(mz-1.0,rp)+1.0-k1
  rval = rval + (1.0-f1)*(1.0-f2)*(1.0-f3)*a(i1,j1,k1)
  rval = rval + f1      *(1.0-f2)*(1.0-f3)*a(i2,j1,k1)
  rval = rval + (1.0-f1)*f2      *(1.0-f3)*a(i1,j2,k1)
  rval = rval + f1      *f2      *(1.0-f3)*a(i2,j2,k1)
  rval = rval + (1.0-f1)*(1.0-f2)*f3      *a(i1,j1,k2)
  rval = rval + f1      *(1.0-f2)*f3      *a(i2,j1,k2)
  rval = rval + (1.0-f1)*f2      *f3      *a(i1,j2,k2)
  rval = rval + f1      *f2      *f3      *a(i2,j2,k2)
  end if ! wrap avoid
  trilinr4 = rval
end function trilinr4


!**********************************************************************!
!
! function trilinr4wrap
!
! returns trilinear interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m,p)      = input data array
!   integer*4 :: n, m, p    = number of data samples
!   real*4 :: x, y, z       = interpolation position
!
! REMARKS:
!   position indices run from 1,1,1 to n,m,p
!   interpolation uses wrap around
!
real*4 function trilinr4wrap(a,n,m,p,x,y,z)
  implicit none
  integer*4, intent(in) :: n, m, p
  real*4, intent(in) :: a(n,m,p), x, y, z
  integer*4 :: i1, j1, k1, i2, j2, k2   ! surrounding pixels
  real*4 :: f1, f2, f3                  ! positional fractions
  real*4 :: mx, my, mz                  ! wrapped positions
  real*4 :: rn, rm, rp, rval            ! temp values
  rval = 0.0
  rn = real(n)
  rm = real(m)
  rp = real(p)
  mx = modulo(x-1.0,rn)+1.0
  my = modulo(y-1.0,rm)+1.0
  mz = modulo(z-1.0,rp)+1.0
  i1 = int(mx)
  j1 = int(my)
  k1 = int(mz)
  i2 = modulo(i1,n)+1
  j2 = modulo(j1,m)+1
  k2 = modulo(k1,p)+1
  f1 = modulo(mx-1.0,rn)+1.0-i1
  f2 = modulo(my-1.0,rm)+1.0-j1
  f3 = modulo(my-1.0,rm)+1.0-k1
  rval = rval + (1.0-f1)*(1.0-f2)*(1.0-f3)*a(i1,j1,k1)
  rval = rval + f1      *(1.0-f2)*(1.0-f3)*a(i2,j1,k1)
  rval = rval + (1.0-f1)*f2      *(1.0-f3)*a(i1,j2,k1)
  rval = rval + f1      *f2      *(1.0-f3)*a(i2,j2,k1)
  rval = rval + (1.0-f1)*(1.0-f2)*f3      *a(i1,j1,k2)
  rval = rval + f1      *(1.0-f2)*f3      *a(i2,j1,k2)
  rval = rval + (1.0-f1)*f2      *f3      *a(i1,j2,k2)
  rval = rval + f1      *f2      *f3      *a(i2,j2,k2)
  trilinr4wrap = rval
end function trilinr4wrap



!**********************************************************************!
!
! function trinnr4
!
! returns 3d nearest neighbour interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m,p)      = input data array
!   integer*4 :: n, m, p    = number of data samples
!   real*4 :: x, y, z       = interpolation position
!
! REMARKS:
!   position indices run from 1,1,1 to n,m,p
!   interpolation does not use wrap around
!
real*4 function trinnr4(a,n,m,p,x,y,z)
  implicit none
  integer*4, intent(in) :: n, m, p
  real*4, intent(in) :: a(n,m,p), x, y, z
  integer*4 :: i1, j1, k1               ! surrounding pixels
  real*4 :: mx, my, mz                  ! wrapped positions
  real*4 :: rn, rm, rp, rval            ! temp values
  rval = 0.0
  rn = real(n)
  rm = real(m)
  rp = real(p)
  if ( (x>=1.0).and.(x<rn+1.0).and. &
     & (y>=1.0).and.(y<rm+1.0).and. &
     & (z>=1.0).and.(y<rp+1.0)) then ! wrap avoid
  mx = modulo(x-1.0,rn)+1.0
  my = modulo(y-1.0,rm)+1.0
  mz = modulo(z-1.0,rp)+1.0
  i1 = nint(mx)
  j1 = nint(my)
  k1 = nint(mz)
  rval = a(i1,j1,k1)
  end if ! wrap avoid
  trinnr4 = rval
end function trinnr4


!**********************************************************************!
!
! function trinnr4wrap
!
! returns 3d nearest neighbour interpolation value (real*4)
!
! INPUT:
!   real*4 :: a(n,m,p)      = input data array
!   integer*4 :: n, m, p    = number of data samples
!   real*4 :: x, y, z       = interpolation position
!
! REMARKS:
!   position indices run from 1,1,1 to n,m,p
!   interpolation uses wrap around
!
real*4 function trinnr4wrap(a,n,m,p,x,y,z)
  implicit none
  integer*4, intent(in) :: n, m, p
  real*4, intent(in) :: a(n,m,p), x, y, z
  integer*4 :: i1, j1, k1               ! surrounding pixels
  real*4 :: mx, my, mz                  ! wrapped positions
  real*4 :: rn, rm, rp, rval            ! temp values
  rval = 0.0
  rn = real(n)
  rm = real(m)
  rp = real(p)
  mx = modulo(x-1.0,rn)+1.0
  my = modulo(y-1.0,rm)+1.0
  mz = modulo(z-1.0,rp)+1.0
  i1 = nint(mx)
  j1 = nint(my)
  k1 = nint(mz)
  rval = a(i1,j1,k1)
  trinnr4wrap = rval
end function trinnr4wrap



!**********************************************************************!
!
! function r4interpol3d
!
! returns interpolation value (real*4) on 3d arrays
! allows selection of interpolation method
! allows selection of wrap around option
!
! INPUT:
!   real*4 :: a(n,m,p)      = input data array
!   integer*4 :: n, m, p    = number of data samples
!   real*4 :: x, y, z       = interpolation position
!   integer*4 :: im         = interpolation method
!                             0 = nearest neighbour
!                             1 = linear
!                             2 = cubic
!   integer*4 :: wa         = wrap around flag
!                             0 = without wrap around, zeroes returned offsides
!                             1 = with wrap around
!
! REMARKS:
!   position indices run from 1,1,1 to n,m,p
!   cubic interpolation kernel parameter is defined in include 'interpolprm.fi'
!
real*4 function r4interpol3d(a, n, m, p, x, y, z, im, wa)

  implicit none
  
  include 'interpolprm.fi'
  
  integer*4, intent(in) :: n, m, p, im, wa
  real*4, dimension(n,m,p), intent(in) :: a
  real*4, intent(in) :: x, y, z
  
  integer*4 :: i, j
  real*4 :: bip
  real*4 :: rval
  real*4, external :: trinnr4, trinnr4wrap
  real*4, external :: trilinr4, trilinr4wrap
  real*4, external :: tricubicr4, tricubicr4wrap
  
  rval = 0.0
  bip = prm_interpol_cubic
  
  ! distinguish interpolation method
  select case (im)
  
  case (0) ! nearest neighbor
    if (wa==0) then
      rval = trinnr4(a,n,m,p,x,y,z)
    else
      rval = trinnr4wrap(a,n,m,p,x,y,z)
    end if
  
  case (1) ! linear interpolation
    if (wa==0) then
      rval = trilinr4(a,n,m,p,x,y,z)
    else
      rval = trilinr4wrap(a,n,m,p,x,y,z)
    end if
  
  case (2) ! cubic interpolation
    if (wa==0) then
      rval = tricubicr4(a,n,m,p,bip,x,y,z)
    else
      rval = tricubicr4wrap(a,n,m,p,bip,x,y,z)
    end if
  
  end select ! case (im)
  
  r4interpol3d = rval
  
  return
  
end function r4interpol3d