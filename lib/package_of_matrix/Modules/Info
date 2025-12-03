!*******************************************************************************
!*                                                                             *
!*  precdef.mod.f90 - type definitions and constants                           *
!*                                                                             *
!*  Thomas C. Lang <lang@physik.uni-wuerzburg.de>, 2005-2007                   *
!*                                                                             *
!*******************************************************************************

!===============================================================================
   MODULE precdef
!-------------------------------------------------------------------------------
   IMPLICIT NONE

!  for intel and lahey fortran compilers on x86-architecture:
!
!  kind(byte)   = 1;  8 bit  -128 ... 127
!  kind(short)  = 2; 16 bit  -32768 ... 32767
!  kind(long)   = 4; 32 bit  -2147483648 ... 2147483647
!  kind(single) = 4; 32 bit, 6 digits, 10^(+/- 37) range
!  kind(double) = 8; 64 bit, 15 digits, 10^(+/- 307) range

   INTEGER, PARAMETER :: &
      byte   = selected_int_kind(2),           &
      short  = selected_int_kind(4),           &
      long   = selected_int_kind(9),           &
      single = selected_real_kind(p=6,r=37),   & ! kind(1.0),
      double = selected_real_kind(p=15,r=307)    ! selected_real_kind(2*precision(1.0_double))
   
   INTEGER, PARAMETER :: &
      sgl = single, &
      dbl = double

   REAL(KIND=double), PARAMETER :: &
      zero    =  0.0_dbl, &
      half    =  0.5_dbl, &
      onehalf =  half,    &
      one     =  1.0_dbl, &
      two     =  2.0_dbl, &
      three   =  3.0_dbl, &
      four    =  4.0_dbl, &
      pi      =  3.1415926535897932384626433832795028_dbl, &  ! = acos(-one)
      E       =  2.7182818284590452353602874713526624_dbl, &  ! = exp(1)
      kB      =  1.38065_dbl*10_dbl**(-23_dbl),            &
      eps     =  epsilon(1.0_dbl)

   COMPLEX(KIND=double), PARAMETER :: &
      ci     = (0.0_dbl,1.0_dbl), &
      cpi    = (3.1415926535897932384626433832795028_dbl,0.0_dbl), &
      ipi    = (0.0_dbl,3.1415926535897932384626433832795028_dbl), &
      iE     = (2.7182818284590452353602874713526624_dbl,0.0_dbl), &
      czero  = (0.0_dbl,0.0_dbl), &
      cone   = (1.0_dbl,0.0_dbl), &
      ctwo   = (2.0_dbl,0.0_dbl), &
      cthree = (3.0_dbl,0.0_dbl), &
      cfour  = (4.0_dbl,0.0_dbl)
      
   END MODULE precdef
