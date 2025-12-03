program SpinFluctuation
    use LocalSweep_mod
    use GlobalUpdate_mod
    use FourierTrans_mod
    use LocalK_mod,    only: Acc_Kl, Acc_lambda
    use ProcessMatrix
    use calc_basic
    implicit none
    include 'mpif.h'
    
    integer :: status(MPI_STATUS_SIZE)
    integer:: iseed, nth, nbc, N
    logical :: is_beta, istau_tmp
    real(kind=8) :: collect, CPUT
    integer(kind=8) :: ICPU_1, ICPU_2, N_P_SEC
    
    type(GlobalUpdate) :: Sweep_global
    type(LocalSweep) :: Sweep_local
    type(FourierTrans) :: Fourier
    type(Propagator), allocatable :: PropU, PropD
    type(WrapList),   allocatable :: WrU, WrD

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
    
    call SYSTEM_CLOCK(COUNT_RATE = N_P_SEC)
    call SYSTEM_CLOCK(COUNT = ICPU_1)
! initiate
    call Model_init(iseed)
    allocate(PropU); call PropU%make()
    allocate(PropD); call PropD%make()
    allocate(WrU);   call WrU%make()
    allocate(WrD);   call WrD%make()
    call Stabilize_init()
    call Sweep_local%init()
    if (is_global) call Sweep_global%init()
! boson warm-up
    if (is_warm) then
        ! 改进策略：热化阶段也使用Wolff更新来处理强时间耦合
        ! 混合使用 Wolff (时间方向) + local (空间方向)
        if (IRANK == 0) then
            write(50,*) "=============================================="
            write(50,*) "Starting thermalization (gauge-only):"
            if (is_global) then
                write(50,*) "  Strategy: Wolff (time) + Local (space)"
                write(50,*) "  Reason: Strong time coupling (gamma≈", &
                    -0.5d0*log(tanh(Dtau*h)), ") requires cluster updates"
            else
                write(50,*) "  Strategy: Local updates only"
                write(50,*) "  Warning: May be slow for strong time coupling!"
            endif
            write(50,*) "=============================================="
            call flush(50)
        endif
        
        do nth = 1, Nwarm
            call Sweep_local%therm(iseed)
            
            ! 每100步打印一次进度
            if (mod(nth, 100) == 0 .and. IRANK == 0) then
                write(50,'(A,I5,A,I5,A,F8.4,A)') '  Therm step ', nth, '/', Nwarm, &
                    ' - accept ratio: ', Acc_Kt%acc / dble(nth), ''
                call flush(50)
            endif
        enddo
        call Sweep_local%ctrl_print_t()
    else
        if (IRANK == 0) write(50,*) "Skipping Bosonic warm-up"
    endif
! Sweep
    call Sweep_local%pre(PropU, PropD, WrU, WrD)
    is_beta = .true.; istau_tmp = .false.
    do nbc = 1, Nbin
        if (nbc .gt. Nthermal) istau_tmp = is_tau ! 预热结束，开始进行动力学测量
        if (is_global) call Sweep_global%sweep(PropU, PropD, WrU, WrD, iseed, is_beta)
        call Sweep_local%sweep(PropU, PropD, WrU, WrD, iseed, is_beta, istau_tmp)
        call Fourier%preq(Obs_equal)
        if (istau_tmp) call Fourier%prtau(Obs_tau)
! 每 20 bin 打印一次接受率
        if (mod(nbc, 20) == 0) then
            if (RT > Zero) then
                collect = 0.d0
                call MPI_Reduce(Acc_Kl%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
                if (IRANK == 0) then
                    write(50,*) 'Progress ', nbc, ': Accept_Klocal_shift: ', collect / dble(ISIZE * nbc)
                endif
            endif
            collect = 0.d0
            call MPI_Reduce(Acc_lambda%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
            if (IRANK == 0) then
                write(50,*) 'Progress ', nbc, ': Accept_lambda_local: ', collect / dble(ISIZE * nbc)
            endif
            if (is_global) then
                collect = 0.d0
                call MPI_Reduce(Acc_Kg%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
                if (IRANK == 0) then
                    write(50,*) 'Progress ', nbc, ': Accept_Kglobal: ', collect / dble(ISIZE * nbc)
                endif
            endif
            if (IRANK == 0) call flush(50)
        endif
    enddo
! control print
    if (is_global) call Sweep_global%ctrl_print()
    call Sweep_local%ctrl_print_l(istau_tmp)
    collect = 0.d0
    call MPI_Reduce(PropU%Xmaxm, collect, 1, MPI_Real8, MPI_MAX, 0, MPI_COMM_WORLD, IERR)
    if (IRANK == 0) PropU%Xmaxm = collect
    N = 2* ISIZE * Nbin * Nst * Nsweep
    collect = 0.d0
    call MPI_Reduce(PropU%Xmeanm, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
    if (IRANK == 0) PropU%Xmeanm = collect / dble(N)
    call SYSTEM_CLOCK(COUNT = ICPU_2)
    CPUT = dble(ICPU_2 - ICPU_1) / dble(N_P_SEC)
    collect = 0.d0
    call MPI_Reduce(CPUT, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
    if (IRANK == 0) CPUT = collect / dble(ISIZE)
    if (IRANK == 0) then
        write(50,*) 'Max diff Matrix                                :', PropU%Xmaxm
        write(50,*) 'Mean diff Matrix                               :', PropU%Xmeanm
        write(50,*) 'Tot CPU time                                   :', CPUT
        call flush(50)
        close(50)
    endif
! deallocate
    if (is_global) call Sweep_global%clear()
    call Sweep_local%clear()
    call Stabilize_clear()
    deallocate(PropU)
    deallocate(PropD)
    deallocate(WrU)
    deallocate(WrD)
    call Model_clear(iseed) ! conf-out
    
    call MPI_FINALIZE(IERR)
end program SpinFluctuation
