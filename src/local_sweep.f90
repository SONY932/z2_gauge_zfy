module LocalSweep_mod
    use Dynamics_mod
    use LocalK_mod
    use GlobalK_mod
    use ObserEqual_mod
    use calc_basic
    use DQMC_Model_mod
    implicit none
    
    public
    private :: Dyn
    
    type :: LocalSweep
    contains
        procedure :: init => Local_sweep_init
        procedure :: clear => Local_sweep_clear
        procedure, private, nopass :: reset => Local_sweep_reset
        procedure :: therm => Local_sweep_therm
        procedure :: pre => Local_sweep_pre
        procedure, private, nopass :: sweep_L => Local_sweep_L
        procedure, private, nopass :: sweep_R => Local_sweep_R
        procedure :: sweep => Local_sweep
        procedure, nopass :: ctrl_print_l => Local_control_print
        procedure, nopass :: ctrl_print_t => Therm_control_print
    end type LocalSweep
    
    type(ObserTau),   allocatable :: Obs_tau
    type(ObserEqual), allocatable :: Obs_equal
    type(Dynamics) :: Dyn
    
contains
    subroutine Local_sweep_init(this)
        class(LocalSweep), intent(inout) :: this
        allocate(Obs_equal)
        call Obs_equal%make()
        if (is_tau) then
            allocate(Obs_tau)
            call Obs_tau%make()
            call Dyn%init()
        endif
        call LocalK_init()
        return
    end subroutine Local_sweep_init
    
    subroutine Local_sweep_clear(this)
        class(LocalSweep), intent(inout) :: this
        deallocate(Obs_equal)
        call LocalK_clear()
        if (is_tau) then
            deallocate(Obs_tau)
            call Dyn%clear()
        endif
        return
    end subroutine Local_sweep_clear
    
    subroutine Local_sweep_reset(toggle)
        logical, intent(in) :: toggle
        integer :: ii
        real(kind=8) :: sval
        call LocalK_reset()
        call Obs_equal%reset()
        if (toggle) call Obs_tau%reset()
! 当 h=0 时，沿虚时间同步每个格点的自旋，确保严格退化到 2D 经典极限
        if (abs(h) <= Zero) then
            do ii = 1, 2*Lq
                sval = NsigL_K%sigma(ii, 1)
                NsigL_K%sigma(ii, 1:Ltrot) = sval
            enddo
        endif
        return
    end subroutine Local_sweep_reset
    
    subroutine Local_sweep_therm(this, iseed)
        class(LocalSweep), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer :: ii, nt
        call this%reset(.false.)
        do nt = 1, Ltrot
            do ii = 1, 2*Lq
                call LocalK_therm(ii, nt, iseed)
            enddo
        enddo
        call Acc_Kt%ratio()
        return
    end subroutine Local_sweep_therm
    
    subroutine Local_sweep_pre(this, PropU, PropD, WrU, WrD)
        class(LocalSweep), intent(inout) :: this
        class(Propagator), intent(inout) :: PropU, PropD
        class(WrapList),   intent(inout) :: WrU, WrD
        integer :: nt
        call this%reset(.false.)
        call Wrap_pre(PropU, WrU, 0)
        call Wrap_pre(PropD, WrD, 0)
        do nt = 1, Ltrot
            if (RT > Zero) call propK_pre(PropU, PropD, NsigL_K%sigma, nt)
            call propT_pre(PropU, PropD, NsigL_K%lambda, nt)
            if (mod(nt, Nwrap) == 0) then
                call Wrap_pre(PropU, WrU, nt)
                call Wrap_pre(PropD, WrD, nt)
            endif
        enddo
        return
    end subroutine Local_sweep_pre
    
    subroutine Local_sweep_L(PropU, PropD, WrU, WrD, iseed, Nobs)
        class(Propagator), intent(inout) :: PropU, PropD
        class(WrapList),   intent(inout) :: WrU, WrD
        integer, intent(inout) :: iseed, Nobs
        integer :: nt
        do nt = Ltrot, 1, -1
            if (mod(nt, Nwrap) == 0) then
                call Wrap_L(PropU, WrU, nt, "S")
                call Wrap_L(PropD, WrD, nt, "S")
            endif
            call Obs_equal%calc(PropU, PropD, nt)
            Nobs = Nobs + 1
            call propT_L(PropU, PropD, NsigL_K%lambda, nt)
            call LocalK_prop_L(PropU, PropD, iseed, nt)
        enddo
        call Wrap_L(PropU, WrU, 0, "S")
        call Wrap_L(PropD, WrD, 0, "S")
        return
    end subroutine Local_sweep_L
    
    subroutine Local_sweep_R(PropU, PropD, WrU, WrD, iseed, toggle, Nobs, Nobst)
        class(Propagator), intent(inout) :: PropU, PropD
        class(WrapList),   intent(inout) :: WrU, WrD
        logical, intent(in) :: toggle
        integer, intent(inout) :: iseed, Nobs, Nobst
        integer :: nt
        if (toggle) then ! calculating time-sliced Green's function before sweep_right
            call Dyn%reset(PropU, PropD)
            call Dyn%sweep_R(Obs_tau, WrU, WrD)
            Nobst = Nobst + 1
        endif
        call Wrap_R(PropU, WrU, 0, "S")
        call Wrap_R(PropD, WrD, 0, "S")
        do nt = 1, Ltrot
            call LocalK_prop_R(PropU, PropD, iseed, nt)
            call propT_R(PropU, PropD, NsigL_K%lambda, nt)
            if (mod(nt, Nwrap) == 0) then
                call Wrap_R(PropU, WrU, nt, "S")
                call Wrap_R(PropD, WrD, nt, "S")
            endif
            call Obs_equal%calc(PropU, PropD, nt)
            Nobs = Nobs + 1
        enddo
        return
    end subroutine Local_sweep_R
    
    subroutine Local_sweep(this, PropU, PropD, WrU, WrD, iseed, is_beta, toggle)
        class(LocalSweep), intent(inout) :: this
        class(Propagator), intent(inout) :: PropU, PropD
        class(WrapList),   intent(inout) :: WrU, WrD
        integer, intent(inout) :: iseed
        logical, intent(inout) :: is_beta
        logical, intent(in) :: toggle
        integer :: Nobs, Nobst, nsw
        integer :: n_lambda_accept, n_lambda_total

        call this%reset(toggle)
        Nobs = 0; Nobst = 0
        n_lambda_accept = 0; n_lambda_total = 0
        do nsw = 1, Nsweep
            if(is_beta) then
                call this%sweep_L(PropU, PropD, WrU, WrD, iseed, Nobs)
                call this%sweep_R(PropU, PropD, WrU, WrD, iseed, toggle, Nobs, Nobst)
            else
                call this%sweep_R(PropU, PropD, WrU, WrD, iseed, toggle, Nobs, Nobst)
                call this%sweep_L(PropU, PropD, WrU, WrD, iseed, Nobs)
            endif
            ! 全局 λ 更新：暂时禁用
            ! 原因：当前的接受率公式使用 G_0，但正确的公式需要 G_λ
            ! 这会导致 λ 配置偏离正确分布
            ! TODO: 实现正确的 λ 更新需要在 wrap 时计算 det[1 + P[λ] B]
            ! call Global_lambda_update(PropU%Gr, PropD%Gr, iseed, n_lambda_accept, n_lambda_total)
        enddo
        ! 更新 λ 接受率统计（通过 count 方法）
        ! 注意：由于使用参数传递，这里需要手动调用 count
        ! Acc_lambda 的统计已在 Global_lambda_update 返回的计数器中
        ! 这里我们使用一个简化的方法：每次接受调用 count(.true.)，每次拒绝调用 count(.false.)
        ! 但由于 count 期望每次调用一次，我们需要循环
        do nsw = 1, n_lambda_accept
            call Acc_lambda%count(.true.)
        enddo
        do nsw = 1, n_lambda_total - n_lambda_accept
            call Acc_lambda%count(.false.)
        enddo
        call Obs_equal%ave(Nobs)
        if (toggle) call Obs_tau%ave(Nobst)
        call Acc_Kl%ratio()
        call Acc_lambda%ratio()
        return
    end subroutine Local_sweep
    
    subroutine Local_control_print(toggle)
        include 'mpif.h'
        logical, intent(in) :: toggle
        real(kind=8) :: collect
        collect = 0.d0
        call MPI_Reduce(Acc_Kl%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Acc_Kl%acc = collect / dble(ISIZE * Nbin)
        collect = 0.d0
        call MPI_Reduce(Acc_lambda%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) then
            Acc_lambda%acc = collect / dble(ISIZE * Nbin)
            write(50,*) 'Accept_Klocal_shift                            :', Acc_Kl%acc
            write(50,*) 'Accept_lambda_local                            :', Acc_lambda%acc
            call flush(50)
        endif
        if (toggle) call Dyn%ctrl_print()
    end subroutine Local_control_print
    
    subroutine Therm_control_print()
        include 'mpif.h'
        real(kind=8) :: collect
        if (Nwarm <= 0) return
        collect = 0.d0
        call MPI_Reduce(Acc_Kt%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) then
            Acc_Kt%acc = collect / dble(ISIZE * Nwarm)
            write(50,*) "Thermalize Accept Ratio                        :", Acc_Kt%acc
            call flush(50)
        endif
        return
    end subroutine Therm_control_print
end module LocalSweep_mod
