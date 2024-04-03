program evalsurf
  implicit none

  integer :: natoms, nstates
  double precision, dimension(:), allocatable :: e
  double precision, dimension(:,:), allocatable :: h
  double precision, dimension(:,:), allocatable :: geom
  double precision, dimension(:,:,:), allocatable :: cg, dcg
  integer :: i, j, ios
  
!  call get_command_argument(number=1,value=gfile,status=ios)
!  if (ios .ne. 0) then
!          print *, "Cannot get file name from command line. Exiting..."
!          stop
!  end if
  
  ! Initialize potential. Read in geometry. Evaluate surface. Return
  ! adiabatic energies and gradients. Return diabatic energies and
  ! gradients.
  call initPotential()
  call getinfo(natoms, nstates)

  allocate(e(nstates))
  allocate(h(nstates,nstates))
  allocate(geom(3,natoms))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))

  call read_colgeom(geom, natoms)

  call EvaluateSurfgen(geom, e, cg, h, dcg)

  print *, "Quasi-diabatic Hamiltonian"
  do i=1, nstates
          print "(2x,10f24.15)", h(i,:)
  end do
  print *, "Adiabatic energy(a.u.)"
  print "(2x,10f24.15)", e
  print *, ""
  do i=1, nstates
          do j=i, nstates
                  print "('Adiabatic gradients of states: ',2i3)", i, j
                  print "(3f25.15)", cg(:,i,j)
                  print *, ""
                  print "('Diabatic gradients of states: ',2i3)", i, j
                  print "(3f25.15)", dcg(:,i,j)
                  print *, ""
          end do
  end do

  deallocate(e)
  deallocate(h)
  deallocate(geom)
  deallocate(cg)
  deallocate(dcg)

contains

  !-------------------------------------------------------------------
  ! get_unit: get available file unit
  function get_unit () result(u)
    implicit none
    integer :: i
    integer :: u
    logical :: uexist, uopen
    u = 0
    do i = 15, 9999
            inquire(unit=i,exist=uexist,opened=uopen)
            if (uexist .and. .not. uopen) then
                    u=i
                    exit
            end if
    end do
    if (u .eq. 0) stop "No unit available."
  end function get_unit

  !-------------------------------------------------------------------
! read_colgeom: read geometry from standard input
subroutine read_colgeom(gm, na)
  implicit none
  integer, intent(in) :: na
  double precision, dimension(3,na), intent(out) :: gm
  integer :: i

  ! Read geometry from standard input (unit number 5)
  do i = 1, na
          read(5,fmt="(11x,3f14.8,14x)") gm(1:3,i)
  end do
  return
end subroutine read_colgeom
    
end program evalsurf