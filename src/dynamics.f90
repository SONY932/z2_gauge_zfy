module Dynamics_mod
    use Stabilize_mod
    use Multiply_mod
    use ObserTau_mod
    use DQMC_Model_mod
    use calc_basic
    implicit none
    
    public
    private :: Prop_dU, Prop_dD, PropGrU, PropGrD
    
    type :: Dynamics
    contains
        procedure, nopass :: init      => Dyn_init
        procedure, nopass :: clear     => Dyn_clear
        procedure, nopass :: reset     => Dyn_reset
        procedure, nopass :: sweep_R   => Dyn_sweep_R
        procedure, nopass :: ctrl_print=> Dyn_control_print
    end type Dynamics
    
    type(Propagator), allocatable :: Prop_dU, Prop_dD
    type(PropGreen),  allocatable :: PropGrU, PropGrD
    
contains
    subroutine Dyn_init()
        allocate(Prop_dU); call Prop_dU%make()
        allocate(Prop_dD); call Prop_dD%make()
        allocate(PropGrU); call PropGrU%make()
        allocate(PropGrD); call PropGrD%make()
        return
    end subroutine Dyn_init
    
    subroutine Dyn_clear()
        deallocate(Prop_dU)
        deallocate(Prop_dD)
        deallocate(PropGrU)
        deallocate(PropGrD)
        return
    end subroutine Dyn_clear
    
    subroutine Dyn_reset(PropU, PropD)
        class(Propagator), intent(in) :: PropU, PropD
        call Prop_dU%asgn(PropU)
        call Prop_dD%asgn(PropD)
        call PropGrU%reset(PropU%Gr)
        call PropGrD%reset(PropD%Gr)
        return
    end subroutine Dyn_reset
    
    subroutine Dyn_sweep_R(Obs, WrU, WrD)
        class(ObserTau), intent(inout) :: Obs
        class(WrapList), intent(in) :: WrU, WrD
        integer :: nt
        do nt = 1, Ltrot
            call Obs%calc(PropGrU, PropGrD, nt)
            if (RT > Zero) call propgrK_R(Prop_dU, Prop_dD, PropGrU, PropGrD, NsigL_K%sigma, nt)
            call propgrT_R(Prop_dU, Prop_dD, PropGrU, PropGrD, NsigL_K%lambda, nt)
            if (mod(nt, Nwrap) == 0) call Wrap_tau(Prop_dU, PropGrU, WrU, nt)
            if (mod(nt, Nwrap) == 0) call Wrap_tau(Prop_dD, PropGrD, WrD, nt)
        enddo
        return
    end subroutine Dyn_sweep_R
    
    subroutine Dyn_control_print()
        include 'mpif.h'
        integer :: no, N
        real(kind=8) :: collect
        real(kind=8) :: tmpMax(3), tmpMean(3)

        do no = 1, 3
            tmpMax(no) = max(PropGrU%Xmaxm(no), PropGrD%Xmaxm(no))
            tmpMean(no) = (PropGrU%Xmeanm(no) + PropGrD%Xmeanm(no))
        enddo

        do no = 1, 3
            collect = 0.d0
            call MPI_Reduce(tmpMax(no), collect, 1, MPI_Real8, MPI_MAX, 0, MPI_COMM_WORLD, IERR)
            if (IRANK == 0) tmpMax(no) = collect
            N = ISIZE * Nbin * Nst * Nsweep
            collect = 0.d0
            call MPI_Reduce(tmpMean(no), collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
            if (IRANK == 0) tmpMean(no) = collect / dble(N)
        enddo
        if (IRANK == 0) then
            write(50,*) 'Max diff GRtt                                  :', tmpMax(1)
            write(50,*) 'Mean diff GRtt                                 :', tmpMean(1)
            write(50,*) 'Max diff GRt0                                  :', tmpMax(2)
            write(50,*) 'Mean diff GRt0                                 :', tmpMean(2)
            write(50,*) 'Max diff GR0t                                  :', tmpMax(3)
            write(50,*) 'Mean diff GR0t                                 :', tmpMean(3)
            call flush(50)
        endif
        return
    end subroutine Dyn_control_print
end module Dynamics_mod
