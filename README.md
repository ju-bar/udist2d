# udist2d

udist2d is a program to undistort 2-dimensional maps of scalar data.

## Usage

udist2d [options]

## List of options

    -prm 'parameter file name' (required)
	     Sets the control parameters of the program.
	
	-in  'input data file name' (required)
	     Specifies the input data file by name.
	
	-out 'output data file name' (required)
	     Specifies the ouput data file by name.
	
	/sil
	    Suppresses console output of the program
	
	/dbg
	    Activates additional output to the console.

## Example

udist2d -prm "udist2d.prm" -in "distorted_map.dat" -out "undistorted_map.dat"

## Parameter file

The parameter file is a fix text list of program parameter values controlling the
file I/O and undistortion. Multiple parameters in one line can be seprated by
spaces, comma, or colon. Combination like ", ", "; ", or multiple spaces are 
interpreted as single separator.

Following is a documented example parameter file:

    3                                // line 01: input data type (0=float 32-bit, 1=float 64-bit, 2=int 8-bit, 3=int 16-bit, 4=int 32-bit)
	0                                // line 02: byte offset of image data in input file
	2048, 2048                       // line 03: dimension of input data array (cols, rows) (ni1,ni2)
	-0.4815E-04, 0.1289E-01, -0.1276E-01, 0.2809E-03    // line 04: sampling of input data (xi, xj, yi, yj) [nm/pix]
	1025, 1025                       // line 05: reference pixel (fix position) in input data (iref1,jref1)
	2048, 2048                       // line 06: dimension of output data array (cols, rows) (no1,no2)
	0.0125                           // line 07: isotropic target sampling of output data [nm/pix] (starg)
	1025, 1025                       // line 08: reference pixel (fix position) in output data (iref2,jref2)
	0                                // line 09: flag: do periodic wrap around
	1                                // line 10: flag: do over-write existing output files

... any content below line 10 will be ignored.

## Distortion data

The distortion of the input 2d data map is parameterized by four numbers: xi, xj, yi, yj.
These numbers represent a linear distortion matrix, translating from a pixel coordinate (i, j) to
a physical coordinate (x, y) via 

   (x, y) = ( xi * (i-iref1) + xj * (j-jref1), yi * (i-iref1) + yj * (j-jref1) ).

The program applies an "undistortion" by inverting the above assingment and reading intensities
from the input map using bi-cubic interpolation. This is done for each output pixel using the assignment

   (x, y) = starg * (i-iref2, j-ref2).

The resulting map may contain 0 at the bordes in case of no periodic wrap-around (line 09).
