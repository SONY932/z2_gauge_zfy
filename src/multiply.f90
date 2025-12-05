module Multiply_mod
    use DQMC_Model_mod
    use ProcessMatrix
    use calc_basic
    implicit none
    
contains
    subroutine propK_pre(Propu, Propd, sigma, nt)
        ! 参考 CodeXun：只更新 UUR，不做 wrap
        class(Propagator), intent(inout) :: Propu, Propd
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: nt
        
        ! 使用完整的 mmult_R，不是分组版本
        call Op_K%mmult_R(Propu%UUR, Latt, sigma, nt, 1)
        call Op_K%mmult_R(Propd%UUR, Latt, sigma, nt, 1)
        return
    end subroutine propK_pre
    
    subroutine propK_L(Propu, Propd, sigma, nt)
        ! 参考 CodeXun：向左传播
        ! 1. wrap G（mmult_L + mmult_R）
        ! 2. 更新 UUL
        class(Propagator), intent(inout) :: Propu, Propd
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: nt
        
        ! 只做一次 wrap，不是 4 次！
        call Op_K%mmult_L(Propu%Gr, Latt, sigma, nt, 1)
        call Op_K%mmult_L(Propd%Gr, Latt, sigma, nt, 1)
        call Op_K%mmult_R(Propu%Gr, Latt, sigma, nt, -1)
        call Op_K%mmult_R(Propd%Gr, Latt, sigma, nt, -1)
        
        ! 更新 UUL
        call Op_K%mmult_L(Propu%UUL, Latt, sigma, nt, 1)
        call Op_K%mmult_L(Propd%UUL, Latt, sigma, nt, 1)
        return
    end subroutine propK_L
    
    subroutine propK_R(Propu, Propd, sigma, nt)
        ! 参考 CodeXun：向右传播
        ! 1. wrap G（mmult_R + mmult_L）
        ! 2. 更新 UUR
        class(Propagator), intent(inout) :: Propu, Propd
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: nt
        
        ! 只做一次 wrap，不是 4 次！
        call Op_K%mmult_R(Propu%Gr, Latt, sigma, nt, 1)
        call Op_K%mmult_R(Propd%Gr, Latt, sigma, nt, 1)
        call Op_K%mmult_L(Propu%Gr, Latt, sigma, nt, -1)
        call Op_K%mmult_L(Propd%Gr, Latt, sigma, nt, -1)
        
        ! 更新 UUR
        call Op_K%mmult_R(Propu%UUR, Latt, sigma, nt, 1)
        call Op_K%mmult_R(Propd%UUR, Latt, sigma, nt, 1)
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
        
        ! 标准的 B_μ 传播
        ! 注意：P[λ] 投影不在这里应用
        ! 由于 P[λ] 和 B 不对易，标准传播公式不能正确处理 P[λ]
        ! P[λ] 的效应需要在 λ 翻转时通过行列式比公式单独考虑
        call Op_MU%mmult_R(PropU%UUR, lambda, nt, 1)
        call Op_MU%mmult_R(PropD%UUR, lambda, nt, 1)
        return
    end subroutine propT_pre
    
    subroutine propT_L(PropU, PropD, lambda, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        real(kind=8), dimension(Lq), intent(in) :: lambda
        integer, intent(in) :: nt
        
        ! 标准 B_μ 传播
        ! P_λ 的乘法在 Local_lambda_flip 中同步处理
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
        
        ! 标准 B_μ 传播
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
