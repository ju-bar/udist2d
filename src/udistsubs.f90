!**********************************************************************!
!**********************************************************************!
!
! FILE: "udistsubs.f90"
!
! AUTHOR: Dr. J. Barthel
!         Forschungszentrum Jülich
!         Jülich, Germany
!
! PURPOSE: Implementations for image undistortions (udist2d)
!
! VERSION: 0.11, J.B., 15.11.2018
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
!
! subroutine CriticalError
!
! posts an error message to console and halts the program.
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine CriticalError(smessage)

  implicit none
  
  include "global.fi"
  
  character*(*) :: smessage

  nume = nume + 1
  
  call PostMessage("")
  call PostMessage("Error: "//trim(smessage))
  call PostMessage("Critical error. Halting program.")
  call PostMessage("")
  call Outroduce()
  stop

  return

end subroutine CriticalError
!**********************************************************************!


!**********************************************************************!
!
! subroutine PostError
!
! posts an error message to console, continues program
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine PostError(smessage)

  implicit none
  
  include "global.fi"
  
  character*(*) :: smessage

  nume = nume + 1

  call PostMessage("Error: "//trim(smessage))

  return

end subroutine PostError
!**********************************************************************!



!**********************************************************************!
!
! subroutine PostWarning
!
! posts a warning message to console.
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine PostWarning(smessage)

  implicit none
  
  include "global.fi"

  character*(*) :: smessage
  
  numw = numw + 1

  call PostMessage("Warning: "//trim(smessage))

  return

end subroutine PostWarning


!**********************************************************************!
!
! subroutine PostMessage
!
! posts a normal message to console.
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine PostMessage(smessage)

  implicit none
  
  include "global.fi"

  character*(*) :: smessage

  if (nsil==0) then ! do not post if silence flag is set in global options
    write (*,*) " > "//trim(smessage)
  end if

  return

end subroutine PostMessage


!**********************************************************************!
!
! subroutine PostDBGMessage
!
! posts a normal message to console.
!
! INPUT:
!   character(len=*) :: smessage            = the error meassage as string
!
! IN/OUTPUT: none
!
subroutine PostDBGMessage(smessage)

  implicit none
  
  include "global.fi"

  character*(*) :: smessage

  if (ndbg/=0) then ! post if debug flag is set in global options
    write (*,*) " debug > "//trim(smessage)
  end if

  return

end subroutine PostDBGMessage



!**********************************************************************!
!
! subroutine Introduce
!
! writes introduction info to console
!
! INPUT: none
!
! IN/OUTPUT: none
!
subroutine Introduce
  
  implicit none
  
  call PostMessage("")
  call PostMessage(" +---------------------------------------------------+")
  call PostMessage(" | Program [udist2d]                                 |")
  call PostMessage(" | Version: 0.11b (20181115JB)                       |")
  call PostMessage(" | Author : Dr. J. Barthel, ju.barthel@fz-juelich.de |")
  call PostMessage(" |          Forschungszentrum Juelich GmbH           |")
  call PostMessage(" |          Germany 2010 - 2018                      |")
  call PostMessage(" | License: GNU GPL 3 <http://www.gnu.org/licenses/> |")
  call PostMessage(" +---------------------------------------------------+")
  call PostMessage("")
  call PostMessage("")
  
  return

end subroutine Introduce


!**********************************************************************!
!
! subroutine Outroduce
!
! writes introduction info to console
!
! INPUT: none
!
! IN/OUTPUT: none
!
subroutine Outroduce
  
  implicit none
  
  include "global.fi"
  
  character(len=400) :: smsg
  
  write(unit=smsg,fmt='(A,I3,A,I3,A)') "Number of warnings: ",numw,",  number of errors: ",nume,"."
  call PostMessage("")
  call PostMessage(trim(smsg))
  call PostMessage("")
  
  return

end subroutine Outroduce





!**********************************************************************!
!**********************************************************************!
FUNCTION factorial(n)
! function: calculates the factorial of n -> n!
! -------------------------------------------------------------------- !
! parameter: integer*4 :: n
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4 :: factorial
  integer*4, intent(in) :: n
  integer*4 :: i
! ------------

! ------------
! init
!  write(unit=*,fmt=*) " > factorial: INIT."
  factorial = 0 ! precheck default -> this means ERROR!
  if (n<0) return
  factorial = 1 ! postcheck default -> this means NO ERROR!
  i=2
! ------------

! ------------
  do while (n>=i)
    factorial = factorial * i
    i = i + 1
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > factorial: EXIT."
  return

END FUNCTION factorial
!**********************************************************************!



!**********************************************************************!
!**********************************************************************!
FUNCTION binomial(n,k)
! function: calculates the binomial coefficient of (n over k), which is
!           equal to (n!)/( (n-k)! * k! )
! -------------------------------------------------------------------- !
! parameter: integer*4 :: n,k
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  integer*4 :: binomial
  integer*4, intent(in) :: n, k
  integer*4, external :: factorial
! ------------

! ------------
! init
!  write(unit=*,fmt=*) " > binomial: INIT."
  binomial = 0 ! precheck default -> this means ERROR!
  if (n<0.or.k<0.or.n<k) return
  binomial = factorial(n)/( factorial(n-k) * factorial(k) )
! ------------

! ------------
!  write(unit=*,fmt=*) " > binomial: EXIT."
  return

END FUNCTION binomial
!**********************************************************************!





!**********************************************************************!
!**********************************************************************!
FUNCTION sigmoid(x,x0,dx)
! function: 0.5*(tanh((x-x0)/dx)+1)
! -------------------------------------------------------------------- !
! parameter: all real*4
!            
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  real*4 :: sigmoid
  real*4, intent(in) :: x, x0, dx
! ------------

! ------------
! init
!  write(unit=*,fmt=*) " > sigmoid: INIT."
  sigmoid = 0.5*(tanh((x-x0)/dx)+1.0)
! ------------



! ------------
!  write(unit=*,fmt=*) " > sigmoid: EXIT."
  return

END FUNCTION sigmoid
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE SetVarString(carray,n,string)
! function: copy data from string to integer*1 array
! -------------------------------------------------------------------- !
! parameter: carray - integer*1(n) array
!            n      - array size
!            string - character(len=*)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  character(len=*) :: string
  integer*1, dimension(1:n) :: carray
  integer*4 :: i, n, slen, alen, mlen
! ------------

! ------------
! init
!  write(unit=*,fmt=*) " > SetVarString: INIT."
  alen = size(carray,dim=1)
  slen = len(string)
  mlen = min(alen,slen)
  if (mlen<=0) return
! ------------

! ------------
! copy
  do i=1, mlen
    carray(i)=mod(ichar(string(i:i)),256)
  end do
  if (alen>mlen) then
    do i=1+mlen, alen
      carray(i) = 0
    end do
  end if
! ------------

! ------------
!  write(unit=*,fmt=*) " > SetVarString: EXIT."
  return

END SUBROUTINE SetVarString
!**********************************************************************!


!**********************************************************************!
!**********************************************************************!
SUBROUTINE GetVarString(string,carray,n)
! function: copy n data from character array to string
! -------------------------------------------------------------------- !
! parameter: carray - integer*1(n) array
!            n      - size of carray
!            string - character(len=*)
! -------------------------------------------------------------------- !

  implicit none

! ------------
! declaration
  character(len=n) :: string
  integer*1, dimension(1:n) :: carray
  integer*4 :: n, i, slen, alen, mlen
! ------------

! ------------
! init
!  write(unit=*,fmt=*) " > GetVarString: INIT."
  alen = size(carray,dim=1)
  slen = len(string)
  mlen = min(alen,slen)
  if (mlen<=0) return
! ------------

! ------------
! copy
  do i=1, mlen
    string(i:i) = achar(carray(i))
  end do
! ------------

! ------------
!  write(unit=*,fmt=*) " > GetVarString: EXIT."
  return

END SUBROUTINE GetVarString
!**********************************************************************!









!**********************************************************************!
!
! subroutine ExplainUsage
!
! posts usage info to console.
!
! INPUT: none
!
! IN/OUTPUT: none
!
subroutine ExplainUsage()

  implicit none

  call PostMessage("")
  call PostMessage("Usage of udist2d in command line:")
  call PostMessage("UDIST2D [-prm 'parameter file name', e.g. 'udist.prm']")
  call PostMessage("        [-in  'input image data (binary) file name', e.g. 'img_001.dat']")
  call PostMessage("        [-out 'output image data (binary) file name', e.g. 'img_001ud.dat']")
  call PostMessage("        [/sil deactivates console output]")
  call PostMessage("        [/dbg activates extra debug console output]")
  call PostMessage("")

  return

end subroutine ExplainUsage


!**********************************************************************!
!
! subroutine ParseCommandLine
!
! parses the command line options and sets global variables.
! also performs parameter checks
!
! INPUT: none
!
! IN/OUTPUT: none
!
subroutine ParseCommandLine()

  implicit none
  
  include "global.fi"

  character*512 :: buffer, cmd
  logical :: fex
  integer*4 :: lfu
  integer*4 :: i, cnt, status, len, nfound
  integer*4 :: nprm, nout, nin
  integer*4, external :: getfreelun

! ------------
! initialize
!  write(unit=*,fmt=*) " > ParseCommandLine: INIT."
  i = 0
  cnt = command_argument_count()
  if (cnt==0) then
    call ExplainUsage()
    call CriticalError("No arguments found.")
  end if
  nprm = 0
  nout = 0
  nin = 0

! ------------
! LOOP OVER ALL GIVEN ARGUMENTS
  do
    i = i + 1
    if (i>cnt) exit
    
    call get_command_argument (i, buffer, len, status)
    if (status/=0) then
      call ExplainUsage()
      call CriticalError("Command line parsing error.")
    end if
    
    ! CHECK COMMAND
    nfound = 0
    cmd = buffer(1:len)
    CHECK_COMMAND: select case (cmd(1:len))
    
    ! THE PARAMETER FILE
    case ("-prm")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status)
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      write(unit = sprmfile, fmt='(A)') buffer(1:len)
      inquire(file=trim(sprmfile),exist=fex)
      if (.not.fex) then
        call CriticalError("Invalid argument: Specified parameter file ["//trim(sprmfile)//"] does not exist.")
      end if
      nprm = 1
      
    ! THE INPUT FILE
    case ("-in")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status)
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      write(unit = simgfile, fmt='(A)') buffer(1:len)
      inquire(file=trim(simgfile),exist=fex)
      if (.not.fex) then
        call CriticalError("Invalid argument: Specified input image data file ["//trim(simgfile)//"] does not exist.")
      end if
      nin = 1
      
    ! THE INPUT FILE
    case ("-out")
      nfound = 1
      i = i + 1
      if (i>cnt) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      call get_command_argument (i, buffer, len, status)
      if (status/=0) then
        call ExplainUsage()
        call CriticalError("Command line parsing error.")
      end if
      write(unit = soutfile, fmt='(A)') buffer(1:len)
      if (len_trim(soutfile)==0) then
        call CriticalError("Invalid argument: Specified output file name is invalid.")
      end if
      nout = 1
      
    case ("/sil")
      nfound = 1
      nsil = 1
      
    case ("/dbg")
      nfound = 1
      ndbg = 1
    
    end select CHECK_COMMAND
    
    if (nfound == 0) then
      call ExplainUsage()
      call CriticalError("Command line parsing error. Unknown command ["//cmd(1:len)//"].")
    end if
  
  end do

! ------------
! final option existence checks
  if (nprm==0) then
    call ExplainUsage()
    call CriticalError("Command line error. Parameter file not specified")
  end if
  
  if (nin==0) then
    call ExplainUsage()
    call CriticalError("Command line error. Input file not specified")
  end if
  
  if (nout==0) then
    call ExplainUsage()
    call CriticalError("Command line error. Output file not specified")
  end if

! ------------
  return

END SUBROUTINE ParseCommandLine



!**********************************************************************!
!
! subroutine loadprm
!
! loads the program parameters from file
!
! Format of the parameter file: TEXT ASCII
!
! line 01: 4                                ! input data type option (0=float 32-bit, 1=float 64-bit, 2=int 8-bit, 3=int 16-bit, 5=int 32-bit)
! line 02: 0                                ! byte offset of image data in input file
! line 03: 2048, 2048                       ! dimension of input data array (cols and rows)
! line 04: -0.45574E-04, 0.91449E-02, -0.89793E-02, 0.22102E-03    ! sampling of input data (xi, xj, yi, yj) [nm/pix] (x=xi*i+xj*j, y=yi*i+yj*j)
! line 05: 1025, 1025                       ! reference pixel (fix position) in input data
! line 06: 2048, 2048                       ! dimension of output data array (cols and rows)
! line 07: 0.008                            ! sampling of output data [nm/pix]
! line 08: 1025, 1025                       ! reference pixel (fix position) in output data
! line 09: 0                                ! flag: do periodic wrap around
! line 10: 1                                ! flag: do over-write existing output files
!
subroutine loadprm

  implicit none
  
  include "global.fi"
  
  integer*4 :: lun, na, ai, i
  integer*4, external :: getfreelun
  character(len=600) :: smsg
  
  
  lun = getfreelun()
  if (lun<=0) goto 13
  
  write(unit=sdbgmsg,fmt=*) "Opening file [",trim(sprmfile),"]."
  call PostDBGMessage(trim(sdbgmsg))
  open( unit=lun, file=trim(sprmfile), iostat=nerr, &
     &  action='read', status='old', err=14)
  
  read(unit=lun, fmt=*, err=16) ndtype
  read(unit=lun, fmt=*, err=16) ndoff
  read(unit=lun, fmt=*, err=16) nix, niy
  read(unit=lun, fmt=*, err=16) sixi, sixj, siyi, siyj
  read(unit=lun, fmt=*, err=16) rpii, rpij
  read(unit=lun, fmt=*, err=16) nox, noy
  read(unit=lun, fmt=*, err=16) sout
  read(unit=lun, fmt=*, err=16) rpoi, rpoj
  read(unit=lun, fmt=*, err=16) nwrap
  read(unit=lun, fmt=*, err=16) noverw
  
  close( unit=lun, iostat=nerr, err=15)
  write(unit=sdbgmsg,fmt=*) "Closing file [",trim(sprmfile),"]."
  call PostDBGMessage(trim(sdbgmsg))
  
  if (ndtype<0) ndtype = 0
  if (ndtype>4) ndtype = 4
  if (ndoff<0) ndoff = 0
  if (nwrap<0) nwrap = 0
  if (nwrap>0) nwrap = 1
  if (noverw<0) noverw = 0
  if (noverw>0) noverw = 1
  
  call PostMessage("Finished loading of parameters.")
  
  write(unit=sdbgmsg,fmt=*) "prm load: data type: ", ndtype
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: data offset: ", ndoff," bytes."
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: input data size: ", nix, niy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: input data sampling: x =", sixi,"*i +",sixj,"*j"
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: input reference pixel: ", rpii, rpij
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: input data sampling: y =", siyi,"*i +",siyj,"*j"
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: output data size: ", nox, noy
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: output data sampling: ", sout
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: output reference pixel: ", rpoi, rpoj
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: do wrap around: ", nwrap
  call PostDBGMessage(trim(sdbgmsg))
  write(unit=sdbgmsg,fmt=*) "prm load: do over-write existing output files: ", noverw
  call PostDBGMessage(trim(sdbgmsg))
  
  return
  
  ! handle errors
13 nerr = 1
  call CriticalError("Failed to acquire free logical file unit.")
14 nerr = 2
  call CriticalError("Failed to open parameter file.")
15 nerr = 3
  call CriticalError("Failed to close parameter file.")
16 nerr = 3
  call CriticalError("Failed to read data from parameter file.")

end subroutine loadprm


!**********************************************************************!
!
! subroutine checkprm
!
! checks input parameters and transfers data to additional variables
!
subroutine checkprm

  implicit none
  
  include "global.fi"
  
  logical :: fex
  
  
  ! check file existence
  ! input
  inquire(file=trim(simgfile),exist=fex)
  if (.not.fex) call CriticalError("Input data file does not exist.")
  ! output
  if (noverw==0) then
    inquire(file=trim(soutfile),exist=fex)
    if (fex) call CriticalError("Output file exists and should not be overwritten.")
  end if
  
  ! check input image size
  if (nix<ndmin .or. nix>ndmax) then
    write(unit=sdbgmsg,fmt='(A,I5,A,I5,A)') "Input data x-dim is out of range (",ndmin," - ",ndmax,")."
    call CriticalError(trim(sdbgmsg))
  end if
  if (niy<ndmin .or. niy>ndmax) then
    write(unit=sdbgmsg,fmt='(A,I5,A,I5,A)') "Input data y-dim is out of range (",ndmin," - ",ndmax,")."
    call CriticalError(trim(sdbgmsg))
  end if
  
  ! check output image size
  if (nox<ndmin .or. nox>ndmax) then
    write(unit=sdbgmsg,fmt='(A,I5,A,I5,A)') "Output data x-dim is out of range (",ndmin," - ",ndmax,")."
    call CriticalError(trim(sdbgmsg))
  end if
  if (noy<ndmin .or. noy>ndmax) then
    write(unit=sdbgmsg,fmt='(A,I5,A,I5,A)') "Output data y-dim is out of range (",ndmin," - ",ndmax,")."
    call CriticalError(trim(sdbgmsg))
  end if
  
  ! check sampling determinant
  if (0.0==( sixj*siyi - sixi*siyj )) then
    call CriticalError("Determinant of input sampling matrix is zero.")
  end if
  
  ! check input reference pixel numbers 
  if (rpii<1 .or. rpii>nix) then
    write(unit=sdbgmsg,fmt='(A,I5,A)') "Input reference pixel x-coordinate is not inside the input image, range (1 .. ",nix,")."
    call CriticalError(trim(sdbgmsg))
  end if
  if (rpij<1 .or. rpij>niy) then
    write(unit=sdbgmsg,fmt='(A,I5,A)') "Input reference pixel y-coordinate is not inside the input image, range (1 .. ",niy,")."
    call CriticalError(trim(sdbgmsg))
  end if
  ! check output reference pixel numbers 
  
  return

end subroutine checkprm


!**********************************************************************!
!
! subroutine Undistort2d
!
! performs the undistort operations by inversion of input sampling
!
! input:
!   include 'global.fi'
!     real*4, parameter :: d2r = 0.0174533      ! degree to radian factor
!     real*4, parameter :: pi = 3.14159265      ! Pi
!     real*4, parameter :: cip = -0.8           ! cubic interpolation parameter   
!     integer*4 :: nix, niy                     ! input discretization
!     real*4 :: sixi, sixj, siyi, siyj          ! input sampling [nm/pix]
!     integer*4 :: nox, noy                     ! output discretization
!     real*4 :: sout                            ! output sampling [nm/pix]
!     integer*4 :: nwrap                        ! do periodic wrap around
!   real*4 :: rin(nix,niy)
!
! output:
!   real*4 :: rout(nox,noy)
!   
!
subroutine Undistort2d(rin,rout)

  implicit none
  
  include 'global.fi'
  
  real*4, intent(in) :: rin(nix,niy)
  real*4, intent(out) :: rout(nox,noy)
  
  integer*4 :: i, j                             ! iterators
  integer*4 :: oi, oj                           ! output pixel indices
  integer*4 :: o0i, o0j                         ! output offset pixel
  real*4 :: ii, ij                              ! float input pixel indices for interpolation
  real*4 :: i0i, i0j                            ! float input offset pixel indices for interpolation
  real*4 :: x, y                                ! real data positions
  real*4 :: isidet                              ! inverse of sampling determinant
  
  real*4, external :: bicubicr4 ! (a,n,m,bip,x,y)
  real*4, external :: bicubicr4wrap ! (a,n,m,bip,x,y)
  
  ! init
  nerr = 0
  o0i = rpoi
  o0j = rpoj
  i0i = real(rpii)
  i0j = real(rpij)
  
  isidet = 1.0 / ( sixj*siyi - sixi*siyj )
  
  ! loop through output
  if (nwrap==1) then ! switch wrap option
    ! use wrap around
    do j=1, noy
      oj = j-o0j
      y = real(oj)*sout
      do i=1, nox
        oi = i-o0i
        x = real(oi)*sout
      
        ! inverse transform to pixels
        ii = i0i + ( sixj*y - siyj*x ) * isidet
        ij = i0j + ( siyi*x - sixi*y ) * isidet
      
        rout(i,j) = bicubicr4wrap(rin,nix,niy,cip,ii,ij)
      
      end do
    end do
  else
    ! no wrap around
    do j=1, noy
      oj = j-o0j
      y = real(oj)*sout
      do i=1, nox
        oi = i-o0i
        x = real(oi)*sout
      
        ! inverse transform to pixels
        ii = i0i + ( sixj*y - siyj*x ) * isidet
        ij = i0j + ( siyi*x - sixi*y ) * isidet
      
        rout(i,j) = bicubicr4(rin,nix,niy,cip,ii,ij)
      
      end do
    end do
  end if
  
  
  return

end subroutine Undistort2d