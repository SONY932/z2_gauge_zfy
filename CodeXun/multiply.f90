module Multiply_mod
    use DQMC_Model_mod
    use ProcessMatrix
    implicit none
    
contains
    subroutine propK_pre(Prop, phi, nt)
        class(Propagator), intent(inout) :: Prop
        real(kind=8), dimension(Nspin, Lq, Ltrot), intent(in) :: phi
        integer, intent(in) :: nt
        call Op_K%mmult_R(Prop%UUR, Latt, phi, nt, 1)
        return
    end subroutine propK_pre
    
    subroutine propK_L(Prop, phi, nt)
        class(Propagator), intent(inout) :: Prop
        real(kind=8), dimension(Nspin, Lq, Ltrot), intent(in) :: phi
        integer, intent(in) :: nt
        call Op_K%mmult_L(Prop%Gr, Latt, phi, nt, 1)
        call Op_K%mmult_R(Prop%Gr, Latt, phi, nt, -1)
        call Op_K%mmult_L(Prop%UUL, Latt, phi, nt, 1)
        return
    end subroutine propK_L
    
    subroutine propK_R(Prop, phi, nt)
        class(Propagator), intent(inout) :: Prop
        real(kind=8), dimension(Nspin, Lq, Ltrot), intent(in) :: phi
        integer, intent(in) :: nt
        call Op_K%mmult_R(Prop%Gr, Latt, phi, nt, 1)
        call Op_K%mmult_L(Prop%Gr, Latt, phi, nt, -1)
        call Op_K%mmult_R(Prop%UUR, Latt, phi, nt, 1)
        return
    end subroutine propK_R
    
    subroutine propgrK_R(Prop, Propgr, phi, nt)
        class(Propagator), intent(inout) :: Prop
        class(PropGreen), intent(inout) :: Propgr
        real(kind=8), dimension(Nspin, Lq, Ltrot), intent(in) :: phi
        integer, intent(in) :: nt
        call Op_K%mmult_R(Propgr%Grt0, Latt, phi, nt, 1)
        call Op_K%mmult_L(Propgr%Gr0t, Latt, phi, nt, -1)
        call Op_K%mmult_R(Propgr%Grtt, Latt, phi, nt, 1)
        call Op_K%mmult_L(Propgr%Grtt, Latt, phi, nt, -1)
        call Op_K%mmult_R(Prop%UUR, Latt, phi, nt, 1)
        return
    end subroutine propgrK_R
    
    subroutine propT_pre(Prop)
        class(Propagator), intent(inout) :: Prop
        call Op_T%mmult_R(Prop%UUR, 1)
        return
    end subroutine propT_pre
    
    subroutine propT_L(Prop)
        class(Propagator), intent(inout) :: Prop
        call Op_T%mmult_L(Prop%Gr, 1)
        call Op_T%mmult_R(Prop%Gr, -1)
        call Op_T%mmult_L(Prop%UUL, 1)
        return
    end subroutine propT_L
    
    subroutine propT_R(Prop)
        class(Propagator), intent(inout) :: Prop
        call Op_T%mmult_R(Prop%Gr, 1)
        call Op_T%mmult_L(Prop%Gr, -1)
        call Op_T%mmult_R(Prop%UUR, 1)
        return
    end subroutine propT_R
    
    subroutine propgrT_R(Prop, Propgr)
        class(Propagator), intent(inout) :: Prop
        class(PropGreen), intent(inout) :: Propgr
        call Op_T%mmult_R(Propgr%Grt0, 1)
        call Op_T%mmult_L(Propgr%Gr0t, -1)
        call Op_T%mmult_R(Propgr%Grtt, 1)
        call Op_T%mmult_L(Propgr%Grtt, -1)
        call Op_T%mmult_R(Prop%UUR, 1)
        return
    end subroutine propgrT_R
end module Multiply_mod