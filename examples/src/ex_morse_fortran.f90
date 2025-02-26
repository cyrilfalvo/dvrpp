
module tools

  implicit none

integer, parameter :: mkind=kind(0.d0)
real(mkind), parameter :: alpha=0.3d0

contains

function morse(q)
  real(mkind), dimension(1), intent(in) :: q
  real(mkind) :: morse
  real(mkind) :: fact
  fact = 1.-exp(-alpha*q(1))
  morse = 0.5d0/alpha/alpha*fact*fact
end function morse

function analytic(n)
  integer, intent(in) :: n
  real(mkind) :: analytic
  real(mkind) :: D,xe,nn

  D = 0.5d0/alpha/alpha
  xe = 0.25d0/D
 
  nn = dble(n)+0.5d0
  analytic = nn*(1. - xe*nn)
end function analytic


end module tools


program ex_morse_fortran

  use tools
  use DVRpp_fortran
  implicit none

  integer :: np
  integer :: info
  integer :: size
  integer :: i
  integer :: dim
  type(C_PTR) :: dvr
  real(mkind) :: Edvr,Eana
  real(mkind), dimension(:), allocatable :: R

  !np = 20
  np = 31
  info = 0
  size = 0
  dim = 1

  dvr = DVRpp_create()
  call DVRpp_SetDimension(dvr, dim)
  call DVRpp_AddHermBasis(dvr,0,np,1.d0,1.2d0,2.5d0,info)
  !call DVRpp_AddSineBasis(dvr,0,np,1.d0,0.5d0,-2.5d0,info)
  call DVRpp_PrepareDvr(dvr,info)
  size = DVRpp_GetSize(dvr)

  allocate(R(dim))

  do i=0,size-1
    call DVRpp_GetNodeCoord(dvr,i, R,dim)    
    call Dvrpp_LoadPotential(dvr,i,morse(R))
    write(*,'(I3,2(F16.4))') i, R(1), morse(R)
  enddo
  
  call Dvrpp_SolveHamiltonianLapack(dvr,info)

  write(*,*)
  write(*,*)

  do i=0,6
   Edvr = Dvrpp_GetEnergy(dvr,i)
   Eana = analytic(i)
   write(*,'(I3,4(F16.4))') i, Eana,Edvr,Edvr-Eana
  enddo

  deallocate(R)
  call DVRpp_destroy(dvr)  


end program ex_morse_fortran
