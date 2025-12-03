module DQMC_Model_mod
    use NonInteract
    use Fields_mod
    use OperatorK_mod
    implicit none
    
    public
    type(SquareLattice), allocatable :: Latt
    type(SpinConf), allocatable :: NsigL_K
    type(OperatorKinetic), allocatable :: Op_T
    type(OperatorPhonon) :: Op_K
    
contains
    subroutine Model_init(iseed)
        integer, intent(out) :: iseed
! read in parameters
        call read_input()
        call Params_set()
        call write_info()
! initiate lattice lists
        allocate(Latt)
        call Lattice_make(Latt)
! initiate phonon and auxiliary field configuration
        allocate(NsigL_K)
        call NsigL_K%make()
        call conf_in(NsigL_K, iseed, Latt)
! set non-interacting exponential operator
        allocate(Op_T)
        call Op_T%make()
        call Op_T%set(Latt)
! set el-ph coupling exponential
        call Op_K%set()
        return
    end subroutine Model_init
    
    subroutine Model_clear(iseed)
        integer, intent(in) :: iseed
        call conf_out(NsigL_K, iseed)
        deallocate(Op_T)
        deallocate(NsigL_K)
        deallocate(Latt)
        return
    end subroutine Model_clear
end module DQMC_Model_mod