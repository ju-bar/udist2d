!****************************************************************************
!
!  file "udist2d.f90"
!
!****************************************************************************
    


!****************************************************************************
!
!  PROGRAM: udist2d
!
!  PURPOSE: Entry point for the console application.
!           Image undistortion with bicubic interpolation
!
!  LINK:    global.fi
!           udistsubs.f90
!           binio2.f90
!           interpol.f90
!           csprog.fi
!           ConsoleProgressBar.f90
!
!  HISTORY:
!           100121: J.B.
!                   birth, version 0.1b
!           181115: J.B.
!                   aligment with new I/O basis code.
!
!****************************************************************************



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


program udist2d

  implicit none
  
  include "global.fi"
  
  character(len=600) :: stmp
  
  real*4, allocatable :: rin(:,:)               ! input data
  real*4, allocatable :: rout(:,:)              ! output data
  
  
  ! initialization
  call PostMessage("")
  call Introduce()
  nerr = 0
  nume = 0
  numw = 0
  nsil = 0
  ndbg = 0
  
  ! *** 
  ! input
  call PostDBGMessage("Parsing command line.")
  call ParseCommandLine()
  if (nerr/=0) call CriticalError("Failed to determine input parameters.")
  
  ! read parameter file
  call PostDBGMessage("Loading parameters from file ["//trim(sprmfile)//"].")
  call loadprm()
  
  ! check parameters
  call PostDBGMessage("Checking parameters.")
  call checkprm()
  call PostDBGMessage("Parameters checked successfully.")
  
  ! allocate
  write(unit=sdbgmsg,fmt=*) "Allocating memory: rin: 32-bit * ",nix," * ",niy
  call PostDBGMessage(trim(sdbgmsg))
  allocate(rin(nix,niy), stat=nerr)
  if (nerr/=0) then
    write(unit=stmp,fmt=*) "Failed to allocate input data. Code(",nerr,")."
    call CriticalError(trim(stmp))
  end if
  write(unit=sdbgmsg,fmt=*) "Allocating memory: rout: 32-bit * ",nox," * ",noy
  call PostDBGMessage(trim(sdbgmsg))
  allocate(rout(nox,noy), stat=nerr)
  if (nerr/=0) then
    write(unit=stmp,fmt=*) "Failed to allocate output data. Code(",nerr,")."
    call CriticalError(trim(stmp))
  end if
  
  
  ! load input data
  call PostDBGMessage("Loading input image data from file ["//trim(simgfile)//"].")
  call loaddata(trim(simgfile), ndtype, nix * niy, ndoff, rin, nerr)
  if (nerr/=0) then
    write(unit=stmp,fmt=*) "Failed to read input data from file. Code(",nerr,")."
    call PostError(trim(stmp))
    select case (nerr)
    case (-6)
      call CriticalError("Unknown data type.")
    case (-5)
      call CriticalError("File access error while reading.")
    case (-4)
      call CriticalError("File access error while seeking offset.")
    case (-3)
      call CriticalError("Failed to open file for reading.")
    case (-2)
      call CriticalError("File is not existing.")
    case (-1)
      call CriticalError("No free logical file unit available.")
    case (1)
      call CriticalError("File access error while seeking offset.")
    case (2)
      call CriticalError("Failed to allocate temporary array memory.")
    case default
      call CriticalError("Failed to deallocate temporary array memory.")
    end select
  end if
  
  ! zero output data
  rout = 0.0
  call PostMessage("Undistorting image.")
  call Undistort2d(rin,rout)
  
  ! save data
  call PostDBGMessage("Saving output data to file ["//trim(soutfile)//"].")
  call savedata(trim(soutfile), nox * noy, rout, nerr)
  
  ! deallocate
  call PostDBGMessage("Undistorting image.")
  deallocate(rin, stat=nerr)
  deallocate(rout, stat=nerr)
  
  ! *** 
  ! quit
  call PostMessage("Finished.")
  call Outroduce()
  
end program udist2d