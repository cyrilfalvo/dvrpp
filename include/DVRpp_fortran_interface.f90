
module DVRpp_fortran
  use, intrinsic :: iso_c_binding
  implicit none

  ! Interfaces for C functions
  interface
     function DVRpp_create() bind(C, name="DVRpp_create")
       import :: C_PTR
       type(C_PTR) :: DVRpp_create
     end function DVRpp_create

     subroutine DVRpp_destroy(dvr) bind(C, name="DVRpp_destroy")
       import :: C_PTR
       type(C_PTR), value :: dvr
     end subroutine DVRpp_destroy

     subroutine DVRpp_SetDimension(dvr, DimV) bind(C, name="DVRpp_SetDimension")
       import :: C_PTR, C_INT
       type(C_PTR), value :: dvr
       integer(C_INT), value :: DimV
     end subroutine DVRpp_SetDimension

     subroutine DVRpp_AddHermBasis(dvr, mode, N, Fact, Delta, Xeq, info) bind(C, name="DVRpp_AddHermBasis")
       import :: C_PTR, C_INT, C_DOUBLE
       type(C_PTR), value :: dvr
       integer(C_INT), intent(out) :: info
       integer(C_INT), value :: mode, N
       real(C_DOUBLE), value :: Fact, Delta, Xeq
     end subroutine DVRpp_AddHermBasis

     subroutine DVRpp_AddSineBasis(dvr, mode, N, Fact, Delta, Xeq, info) bind(C, name="DVRpp_AddSineBasis")
       import :: C_PTR, C_INT, C_DOUBLE
       type(C_PTR), value :: dvr
       integer(C_INT), intent(out) :: info
       integer(C_INT), value :: mode, N
       real(C_DOUBLE), value :: Fact, Delta, Xeq
     end subroutine DVRpp_AddSineBasis

     subroutine DVRpp_AddExpBasis(dvr, mode, N, Fact, Delta, Xeq, info) bind(C, name="DVRpp_AddExpBasis")
       import :: C_PTR, C_INT, C_DOUBLE
       type(C_PTR), value :: dvr
       integer(C_INT), intent(out) :: info
       integer(C_INT), value :: mode, N
       real(C_DOUBLE), value :: Fact, Delta, Xeq
     end subroutine DVRpp_AddExpBasis

     subroutine DVRpp_PrepareDvr(dvr, info) bind(C, name="DVRpp_PrepareDvr")
       import :: C_PTR, C_INT
       type(C_PTR), value :: dvr
       integer(C_INT), intent(out) :: info
     end subroutine DVRpp_PrepareDvr

     subroutine DVRpp_GetNodeCoord(dvr, index, R, dim) bind(C, name="DVRpp_GetNodeCoord")
       import :: C_PTR, C_INT, C_DOUBLE
       type(C_PTR), value :: dvr
       integer(C_INT), value :: index
       integer(C_INT), value :: dim
       real(C_DOUBLE), dimension(dim), intent(out)  :: R
     end subroutine DVRpp_GetNodeCoord

     subroutine DVRpp_LoadPotential_array(dvr, V, dim) bind(C, name="DVRpp_LoadPotential_array")
       import :: C_PTR, C_INT, C_DOUBLE
       type(C_PTR), value :: dvr
       integer(C_INT), value :: dim
       real(C_DOUBLE), dimension(dim), intent(in) :: V
     end subroutine DVRpp_LoadPotential_array

     subroutine DVRpp_LoadPotential(dvr, index, V) bind(C, name="DVRpp_LoadPotential")
       import :: C_PTR, C_INT, C_DOUBLE
       type(C_PTR), value :: dvr
       integer(C_INT), value :: index
       real(C_DOUBLE), value :: V
     end subroutine DVRpp_LoadPotential

#ifdef _ARPACK_
     subroutine DVRpp_SolveHamiltonianArpack(dvr, nev, ncv, info) bind(C, name="DVRpp_SolveHamiltonianArpack")
       import :: C_PTR, C_INT
       type(C_PTR), value :: dvr
       integer(C_INT), intent(out) :: info
       integer(C_INT), value :: nev, ncv
     end subroutine DVRpp_SolveHamiltonianArpack
#endif

     subroutine DVRpp_SolveHamiltonianLapack(dvr, info) bind(C, name="DVRpp_SolveHamiltonianLapack")
       import :: C_PTR, C_INT
       type(C_PTR), value :: dvr
       integer(C_INT), intent(out) :: info
     end subroutine DVRpp_SolveHamiltonianLapack

     function DVRpp_GetEnergy(dvr, i) bind(C, name="DVRpp_GetEnergy")
       import :: C_PTR, C_INT, C_DOUBLE
       type(C_PTR), value :: dvr
       integer(C_INT), value :: i
       real(C_DOUBLE) :: DVRpp_GetEnergy
     end function DVRpp_GetEnergy

     function DVRpp_GetElmDiagOperatorSingle(dvr, i, j, Ai) bind(C, name="DVRpp_GetElmDiagOperatorSingle")
       import :: C_PTR, C_INT, C_DOUBLE
       type(C_PTR), value :: dvr
       integer(C_INT), value :: i, j
       real(C_DOUBLE), dimension(:), intent(in) :: Ai
       real(C_DOUBLE) :: DVRpp_GetElmDiagOperatorSingle
     end function DVRpp_GetElmDiagOperatorSingle

     subroutine DVRpp_GetElmDiagOperator(dvr, i, j, Ai, rank, Op, info) bind(C, name="DVRpp_GetElmDiagOperator")
       import :: C_PTR, C_INT, C_DOUBLE
       type(C_PTR), value :: dvr
       integer(C_INT), intent(out) :: info
       integer(C_INT), value :: i, j, rank
       real(C_DOUBLE), dimension(:), intent(in) :: Ai
       real(C_DOUBLE), dimension(:), intent(out) :: Op
     end subroutine DVRpp_GetElmDiagOperator

     subroutine DVRpp_GetElmQOpt(dvr, PowMax, i, AvgOpt, info) bind(C, name="DVRpp_GetElmQOpt")
       import :: C_PTR, C_INT, C_DOUBLE
       type(C_PTR), value :: dvr
       integer(C_INT), intent(out) :: info
       integer(C_INT), value :: PowMax, i
       real(C_DOUBLE), dimension(:), intent(out) :: AvgOpt
     end subroutine DVRpp_GetElmQOpt

     subroutine DVRpp_GetPop(dvr, i, ni, info) bind(C, name="DVRpp_GetPop")
       import :: C_PTR, C_INT
       type(C_PTR), value :: dvr
       integer(C_INT), intent(out) :: info
       integer(C_INT), value :: i
       integer(C_INT), dimension(:), intent(out) :: ni
     end subroutine DVRpp_GetPop

     subroutine DVRpp_GetVector(dvr, i, Psi, info) bind(C, name="DVRpp_GetVector")
       import :: C_PTR, C_INT, C_DOUBLE
       type(C_PTR), value :: dvr
       integer(C_INT), intent(out) :: info
       integer(C_INT), value :: i
       real(C_DOUBLE), dimension(:), intent(out) :: Psi
     end subroutine DVRpp_GetVector

     function DVRpp_GetSize(dvr) bind(C, name="DVRpp_GetSize")
       import :: C_PTR, C_INT
       type(C_PTR), value :: dvr
       integer(C_INT) :: DVRpp_GetSize
     end function DVRpp_GetSize

  end interface


end module DVRpp_fortran

  

