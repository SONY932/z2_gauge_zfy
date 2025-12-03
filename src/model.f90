module DQMC_Model_mod
    use MyLattice
    use calc_basic
    use NonInteract
    use Fields_mod
    use OperatorKin_mod
    implicit none

    public
    type(SquareLattice), allocatable :: Latt
    type(SigmaConf), allocatable :: NsigL_K
    type(OperatorMU), allocatable :: Op_MU
    type(OperatorKin) :: Op_K

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
! 初始化规范场构型
        allocate(NsigL_K)
        call NsigL_K%make()
        call conf_in(NsigL_K, iseed, Latt)
! 设置化学势
        allocate(Op_MU)
        call Op_MU%make()
        call Op_MU%set()
! 设置动能算符
        call Op_K%set()
        return
    end subroutine Model_init

    subroutine Model_clear(iseed)
        integer, intent(in) :: iseed
        call conf_out(NsigL_K, iseed)
        deallocate(Op_MU)
        deallocate(NsigL_K)
        deallocate(Latt)
        return
    end subroutine Model_clear

end module DQMC_Model_mod