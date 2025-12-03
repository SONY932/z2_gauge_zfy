module Multiply_mod
    use DQMC_Model_mod
    use ProcessMatrix
    use calc_basic
    implicit none
    
contains
    subroutine propK_pre(Propu, Propd, sigma, nt)
        class(Propagator), intent(inout) :: Propu, Propd
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: nt
        call Op_K%mmult_R(Propu%UUR, Latt, sigma, nt, 1)
        call Op_K%mmult_R(Propd%UUR, Latt, sigma, nt, 1)
        return
    end subroutine propK_pre
    
    subroutine propK_L(Propu, Propd, sigma, nt)
        class(Propagator), intent(inout) :: Propu, Propd
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: nt
        call Op_K%mmult_L(Propu%Gr, Latt, sigma, nt,  1)
        call Op_K%mmult_L(Propd%Gr, Latt, sigma, nt,  1)
        call Op_K%mmult_R(Propu%Gr, Latt, sigma, nt, -1)
        call Op_K%mmult_R(Propd%Gr, Latt, sigma, nt, -1)
        call Op_K%mmult_L(Propu%UUL, Latt, sigma, nt,  1)
        call Op_K%mmult_L(Propd%UUL, Latt, sigma, nt,  1)
        return
    end subroutine propK_L
    
    subroutine propK_R(Propu, Propd, sigma, nt)
        class(Propagator), intent(inout) :: Propu, Propd
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: nt
        call Op_K%mmult_R(Propu%Gr, Latt, sigma, nt,  1)
        call Op_K%mmult_R(Propd%Gr, Latt, sigma, nt,  1)
        call Op_K%mmult_L(Propu%Gr, Latt, sigma, nt, -1)
        call Op_K%mmult_L(Propd%Gr, Latt, sigma, nt, -1)
        call Op_K%mmult_R(Propu%UUR, Latt, sigma, nt,  1)
        call Op_K%mmult_R(Propd%UUR, Latt, sigma, nt,  1)
        return
    end subroutine propK_R
    
    subroutine propgrK_R(PropU, PropD, PropgrU, PropgrD, sigma, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(PropGreen),  intent(inout) :: PropgrU, PropgrD
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: nt
        call Op_K%mmult_R(PropgrU%Grt0, Latt, sigma, nt,  1)
        call Op_K%mmult_R(PropgrD%Grt0, Latt, sigma, nt,  1)

        call Op_K%mmult_L(PropgrU%Gr0t, Latt, sigma, nt, -1)
        call Op_K%mmult_L(PropgrD%Gr0t, Latt, sigma, nt, -1)

        call Op_K%mmult_R(PropgrU%Grtt, Latt, sigma, nt,  1)
        call Op_K%mmult_R(PropgrD%Grtt, Latt, sigma, nt,  1)

        call Op_K%mmult_L(PropgrU%Grtt, Latt, sigma, nt, -1)
        call Op_K%mmult_L(PropgrD%Grtt, Latt, sigma, nt, -1)
        
        call Op_K%mmult_R(PropU%UUR,     Latt, sigma, nt,  1)
        call Op_K%mmult_R(PropD%UUR,     Latt, sigma, nt,  1)
        return
    end subroutine propgrK_R
    
    subroutine propT_pre(PropU, PropD, lambda, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        real(kind=8), dimension(Lq), intent(in) :: lambda
        integer, intent(in) :: nt
        call Op_MU%mmult_R(PropU%UUR, lambda, nt, 1)
        call Op_MU%mmult_R(PropD%UUR, lambda, nt, 1)
        return
    end subroutine propT_pre
    
    subroutine propT_L(PropU, PropD, lambda, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        real(kind=8), dimension(Lq), intent(in) :: lambda
        integer, intent(in) :: nt
        call Op_MU%mmult_L(PropU%Gr, lambda, nt, 1)
        call Op_MU%mmult_L(PropD%Gr, lambda, nt, 1)
        call Op_MU%mmult_R(PropU%Gr, lambda, nt, -1)
        call Op_MU%mmult_R(PropD%Gr, lambda, nt, -1)
        call Op_MU%mmult_L(PropU%UUL, lambda, nt, 1)
        call Op_MU%mmult_L(PropD%UUL, lambda, nt, 1)
        return
    end subroutine propT_L
    
    subroutine propT_R(PropU, PropD, lambda, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        real(kind=8), dimension(Lq), intent(in) :: lambda
        integer, intent(in) :: nt
        call Op_MU%mmult_R(PropU%Gr, lambda, nt, 1)
        call Op_MU%mmult_R(PropD%Gr, lambda, nt, 1)
        call Op_MU%mmult_L(PropU%Gr, lambda, nt, -1)
        call Op_MU%mmult_L(PropD%Gr, lambda, nt, -1)
        call Op_MU%mmult_R(PropU%UUR, lambda, nt, 1)
        call Op_MU%mmult_R(PropD%UUR, lambda, nt, 1)
        return
    end subroutine propT_R
    
    subroutine propgrT_R(PropU, PropD, PropgrU, PropgrD, lambda, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        class(PropGreen),  intent(inout) :: PropgrU, PropgrD
        real(kind=8), dimension(Lq), intent(in) :: lambda
        integer, intent(in) :: nt
        call Op_MU%mmult_R(PropgrU%Grt0, lambda, nt, 1)
        call Op_MU%mmult_R(PropgrD%Grt0, lambda, nt, 1)

        call Op_MU%mmult_L(PropgrU%Gr0t, lambda, nt, -1)
        call Op_MU%mmult_L(PropgrD%Gr0t, lambda, nt, -1)

        call Op_MU%mmult_R(PropgrU%Grtt, lambda, nt, 1)
        call Op_MU%mmult_R(PropgrD%Grtt, lambda, nt, 1)

        call Op_MU%mmult_L(PropgrU%Grtt, lambda, nt, -1)
        call Op_MU%mmult_L(PropgrD%Grtt, lambda, nt, -1)

        call Op_MU%mmult_R(PropU%UUR, lambda, nt, 1)
        call Op_MU%mmult_R(PropD%UUR, lambda, nt, 1)
        return
    end subroutine propgrT_R
end module Multiply_mod
