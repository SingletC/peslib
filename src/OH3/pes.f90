!********************************************************************** 
!   System:                     OH3 (lowest 3 state)
!   Functional form:            permutationally invariant polynomials
!   Common name:                OH3 three states PES
!   Number of derivatives:      3
!   Number of bodies:           4
!   Number of electronic surfaces: 3
!   Interface: Section-2
!
!   References:  Yinan Shu, Joanna Kryven, Antonio Gustavo Sampaio de Oliveira Filho,
!                Steven L. Mielke, Linyao Zhang, Guo-Liang Song, Shaohong L. Li, 
!                Rubén Meana-Pañeda, Bina Fu, Joel M. Bowman, and Donald G. Truhlar
!                to be published
!
!   Notes:    -PES of OH3
!             -Including diabatic states, diabatic couplings, and transformed adiabatic states
!
!**********************************************************************

!***************************************************************************
! PES program for OH3 
!
! Authors: Yinan Shu, Antonio Gustavo Sampaio de Oliveira Filho
! Dec. 2018
!
! Reference for this potential energy surface:
!   "Direct  Diabatization and Analytic Representation of Coupled 
!    Potential Energy Surfaces and Couplings for the Reactive 
!    Quenching of the Excited 2Σ+ State of OH by Molecular Hydrogen"
!   by Yinan Shu, Joanna Kryven, Antonio Gustavo Sampaio de Oliveira Filho,
!      Steven L. Mielke, Guo-Liang Song, Shaohong L. Li, Rubén Meana-Pañeda,
!      Bina Fu, Joel Bowman, and Donald G. Truhlar
!   to be published
!
!***************************************************************************

!===========================================================================
! *def oh3_pes
!  1. Main PES program for calculating:
!     energies and gradients for diabatic states and diabatic couplings
!     energies, gradients and non-adiabatic couplings for adiabatic states 
!  2. The input file: ****.xyz (**** is a file name, for example, OH3;) 
!     !!!THE UNIT SHOULD BE ANGSTROM !!!
!     the format of ****.xyz should be a standard xyz format that is readbale
!     by Molden with 3 hydrogen atoms first, the oxygen at last, for example: 
!     4
!    
!     H  X Y Z 
!     H  X Y Z
!     H  X Y Z
!     O  X Y Z
!  3. Executing the program with the input, the program will generate 2 output 
!     they correspond to: ****.out and ****.tmp
!     ****.out has a detailed print out of energies, gradients and nonadiabatic 
!     couplings
!     ****.tmp has a three line simple  print of diabatic states, diabatic couplings and
!     adiabatic states in the following format
!     u(1,1) u(2,2) u(3,3)
!     u(1,2) u(1,3) u(2,3)
!     a(1)   a(2)   a(3)
!  4. The zero of energy is: -76.80042447 hartree
!
!  ==========================definition of variables========================
!
!  4. xyz: the xyz coordinate
!          unit: in angstroms
!          xyz is a 4 by 3 matrix
!  5.  u : the diabatic matrix
!          unit: in hartrees
!          u is a 3 by 3 matrix
!          notice u(i,j)=u(j,i) for i not equal to j 
!          the diabatic states are diagonal elements of u 
!          the diabatic couplings are off-diagonal elements of u 
!  6.  gu: the gradients of diabatic matrix w.r.t. the xyz coordinates
!          unit: in hartrees/angstrom
!          gu is a 4 by 3 by 3 by 3 matrix
!          "gu(:,:,i,j)" is the gradient of u(i,j)
!  7.  a : the adiabatic states
!          unit: in hartrees
!          a is a 3 element vector 
!  8.  ga: the gradient of adiabatic states w.r.t. the xyz coordinates
!          unit: in hartrees/angstrom
!          ga is a 4 by 3 by 3 matrix
!          "ga(:,:,1)" is the gradient of the ground adiabatic state
!          "ga(:,:,2)" is the gradient of the first excited adiabatic state
!          "ga(:,:,3)" is the gradient of the second excited adiabatic state
!  9.  f : the non-adiabatic couplings
!          unit: in angstrom**(-1)
!          f is a 4 by 3 by 3 by 3 matrix
!          "f(:,:,i,j)" is the non-adiabatic coupling of state i and j
!          f(:,:,i,j)=-f(:,:,j,i), and f(:,:,i,i)=0 
!  10. rv: the diabatic to adiabatic rotation matrix 
!===========================================================================

subroutine oh3_pes(xyz,a,ga)
 
  integer, parameter:: dp=kind(0.d0)                   ! double precision
  character (len=1) :: symb
  integer :: i,j,k,l
  integer :: natm
  real(dp) :: u11, u22, u33
  real(dp) :: u12, u13, u23
  real(dp) :: s12, s23
  real(dp) :: xyz(4,3)
  character(len=32)::fname, bname
  real(dp) :: du11dx(4,3), du22dx(4,3), du33dx(4,3)
  real(dp) :: du12dx(4,3), du13dx(4,3), du23dx(4,3)
  real(dp) :: gu(4,3,3,3)
  real(dp) :: u(3,3),rv(3,3),a(3)
  real(dp) :: ga(4,3,3),f(4,3,3,3)
  real(dp) :: tmpgu(3,3),tmpga(3,3)



!compute diabatic states
  call u11_prepot_8
  call u11_pot_8(xyz,u11,du11dx)
  call u22_prepot_8
  call u22_pot_8(xyz,u22,du22dx)
  call u33_prepot_8
  call u33_pot_8(xyz,u33,du33dx)

!  write(13,*) "diabatic energies:"
!  write(13,*) u11, u22, u33

!  write(13,*) "analytical gradiens of u11"
!  do i=1,4
!  write(13,*) du11dx(i,1), du11dx(i,2), du11dx(i,3)
!  end do

!  write(13,*) "analytical gradiens of u22"
!  do i=1,4
!  write(13,*) du22dx(i,1), du22dx(i,2), du22dx(i,3)
!  end do

!  write(13,*) "analytical gradiens of u33"
!  do i=1,4
!  write(13,*) du33dx(i,1), du33dx(i,2), du33dx(i,3)
!  end do

!compute diabatic couplings
  call u12_order10
  call u13_order12
  call u23_order10
  call u12_dc_10(xyz,u12,s12)
  call u13_dc_12(xyz,u13)
  call u23_dc_10(xyz,u23,s23)
  
!  write(13,*) "diabatic couplings:"
!  write(13,*) u12, u13, u23

!compute diabatic couplings derivatives
  call u12_dvdx(xyz,u12,s12,du12dx)
  call u13_dvdx(xyz,u13,du13dx)
  call u23_dvdx(xyz,u23,s23,du23dx)

!  write(13,*) "analytical gradiens of u12"
!  do i=1,4
!  write(13,*) du12dx(i,1), du12dx(i,2), du12dx(i,3)
!  end do

!  write(13,*) "analytical gradiens of u13"
!  do i=1,4
!  write(13,*) du13dx(i,1), du13dx(i,2), du13dx(i,3)
!  end do

!  write(13,*) "analytical gradiens of u23"
!  do i=1,4
!  write(13,*) du23dx(i,1), du23dx(i,2), du23dx(i,3)
!  end do

!compute the adiabatic states and eigenvectors
  u(1,1)=u11
  u(1,2)=u12
  u(1,3)=u13
  u(2,1)=u12
  u(2,2)=u22
  u(2,3)=u23
  u(3,1)=u13
  u(3,2)=u23
  u(3,3)=u33
  call dsyevj3(u,rv,a)
  call order(a,rv)
!  write(13,*) "adiabatic energies:"
!  write(13,*) a(1), a(2), a(3)
!  write(13,*) "adiabatic eigenvectors:"
!  write(13,*) "v1:",rv(1,1), rv(2,1), rv(3,1)
!  write(13,*) "v2:",rv(1,2), rv(2,2), rv(3,2)
!  write(13,*) "v3:",rv(1,3), rv(2,3), rv(3,3)

!assign dudx to gu
  do i=1,4
  do j=1,3
  gu(i,j,1,1)=du11dx(i,j)
  gu(i,j,1,2)=du12dx(i,j)
  gu(i,j,1,3)=du13dx(i,j)
  gu(i,j,2,1)=du12dx(i,j)
  gu(i,j,2,2)=du22dx(i,j)
  gu(i,j,2,3)=du23dx(i,j)
  gu(i,j,3,1)=du13dx(i,j)
  gu(i,j,3,2)=du23dx(i,j)
  gu(i,j,3,3)=du33dx(i,j)
  enddo
  enddo

!initialize the ga and f matrices
  do i=1,4
  do j=1,3
  do k=1,3
  do l=1,3
  ga(i,j,k)=0.000
  f(i,j,k,l)=0.000
  enddo
  enddo
  enddo
  enddo

!compute the adiabatic states gradients,nonadiabatic coupling vectors
  do i=1,4
  do j=1,3
  tmpgu(1,1)=du11dx(i,j)
  tmpgu(1,2)=du12dx(i,j)
  tmpgu(1,3)=du13dx(i,j)
  tmpgu(2,1)=du12dx(i,j)
  tmpgu(2,2)=du22dx(i,j)
  tmpgu(2,3)=du23dx(i,j)
  tmpgu(3,1)=du13dx(i,j)
  tmpgu(3,2)=du23dx(i,j)
  tmpgu(3,3)=du33dx(i,j)
  !C^T*dU*C
  tmpga = matmul(transpose(rv), matmul(tmpgu, rv))
  !dV_i = (C^T*dU*C)_ii
  do k=1,3
  ga(i,j,k)=tmpga(k,k)
  enddo
  !F_ij = (C^T*dU*C)_ij / (V_j-V_i)
  do k=1,2
  do l=k+1,3
  f(i,j,k,l)=tmpga(k,l)/(a(l)-a(k))
  f(i,j,l,k)=-f(i,j,k,l)
  enddo
  enddo
  enddo
  enddo

!  write(13,*) "v1 analytical gradients:"
!  do i=1,4
!  write(13,*) ga(i,1,1), ga(i,2,1), ga(i,3,1)
!  end do
!
!  write(13,*) "v2 analytical gradients:"
!  do i=1,4
!  write(13,*) ga(i,1,2), ga(i,2,2), ga(i,3,2)
!  end do
!
!  write(13,*) "v3 analytical gradients:"
!  do i=1,4
!  write(13,*) ga(i,1,3), ga(i,2,3), ga(i,3,3)
!  end do
!
!  write(13,*) "f12 nonadiabatic coupling vectors:"
!  do i=1,4
!  write(13,*) f(i,1,1,2), f(i,2,1,2), f(i,3,1,2)
!  end do
!
!  write(13,*) "f13 nonadiabatic coupling vectors:"
!  do i=1,4
!  write(13,*) f(i,1,1,3), f(i,2,1,3), f(i,3,1,3)
!  end do
!
!  write(13,*) "f23 nonadiabatic coupling vectors:"
!  do i=1,4
!  write(13,*) f(i,1,2,3), f(i,2,2,3), f(i,3,2,3)
!  end do

end
