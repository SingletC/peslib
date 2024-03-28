subroutine evaluate_surface(geom, e, cg, h, dcg)
  implicit none

  ! Input parameters

  
  double precision, dimension(3,5), intent(in) :: geom
  double precision, dimension(2), intent(out) :: e
  double precision, dimension(2,2), intent(out) :: h
  double precision, dimension(3*5,2,2), intent(out) :: cg, dcg
  
  ! Local variables
  integer :: ios

  ! Initialize potential. Read in geometry. Evaluate surface. Return
  ! adiabatic energies and gradients. Return diabatic energies and
  ! gradients.

  call initPotential()
  call getinfo(5, 2)

  call EvaluateSurfgen(geom, e, cg, h, dcg)

  
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

    
end subroutine evaluate_surface
