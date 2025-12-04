module Multiply_mod
    use DQMC_Model_mod
    use ProcessMatrix
    use calc_basic
    implicit none
    
contains
    subroutine propK_pre(Propu, Propd, sigma, nt)
        ! 使用分组方式更新 UUR，确保和 LocalK_prop_R 的顺序一致
        class(Propagator), intent(inout) :: Propu, Propd
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: nt
        ! Group 1 -> Group 2 -> Group 3 -> Group 4
        call Op_K%mmult_R_group(Propu%UUR, Latt%group_1, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%UUR, Latt%group_1, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propu%UUR, Latt%group_2, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%UUR, Latt%group_2, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propu%UUR, Latt%group_3, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%UUR, Latt%group_3, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propu%UUR, Latt%group_4, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%UUR, Latt%group_4, Latt, sigma, nt, 1)
        return
    end subroutine propK_pre
    
    subroutine propK_L(Propu, Propd, sigma, nt)
        ! 使用分组方式，确保和 LocalK_prop_L 的顺序一致
        ! Group 4 -> Group 3 -> Group 2 -> Group 1
        class(Propagator), intent(inout) :: Propu, Propd
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: nt
        ! Group 4: wrap and UUL
        call Op_K%mmult_L_group(Propu%Gr, Latt%group_4, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propd%Gr, Latt%group_4, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propu%Gr, Latt%group_4, Latt, sigma, nt, -1)
        call Op_K%mmult_R_group(Propd%Gr, Latt%group_4, Latt, sigma, nt, -1)
        call Op_K%mmult_L_group(Propu%UUL, Latt%group_4, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propd%UUL, Latt%group_4, Latt, sigma, nt, 1)
        ! Group 3
        call Op_K%mmult_L_group(Propu%Gr, Latt%group_3, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propd%Gr, Latt%group_3, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propu%Gr, Latt%group_3, Latt, sigma, nt, -1)
        call Op_K%mmult_R_group(Propd%Gr, Latt%group_3, Latt, sigma, nt, -1)
        call Op_K%mmult_L_group(Propu%UUL, Latt%group_3, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propd%UUL, Latt%group_3, Latt, sigma, nt, 1)
        ! Group 2
        call Op_K%mmult_L_group(Propu%Gr, Latt%group_2, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propd%Gr, Latt%group_2, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propu%Gr, Latt%group_2, Latt, sigma, nt, -1)
        call Op_K%mmult_R_group(Propd%Gr, Latt%group_2, Latt, sigma, nt, -1)
        call Op_K%mmult_L_group(Propu%UUL, Latt%group_2, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propd%UUL, Latt%group_2, Latt, sigma, nt, 1)
        ! Group 1
        call Op_K%mmult_L_group(Propu%Gr, Latt%group_1, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propd%Gr, Latt%group_1, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propu%Gr, Latt%group_1, Latt, sigma, nt, -1)
        call Op_K%mmult_R_group(Propd%Gr, Latt%group_1, Latt, sigma, nt, -1)
        call Op_K%mmult_L_group(Propu%UUL, Latt%group_1, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propd%UUL, Latt%group_1, Latt, sigma, nt, 1)
        return
    end subroutine propK_L
    
    subroutine propK_R(Propu, Propd, sigma, nt)
        ! 使用分组方式，确保和 LocalK_prop_R 的顺序一致
        ! Group 1 -> Group 2 -> Group 3 -> Group 4
        class(Propagator), intent(inout) :: Propu, Propd
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: nt
        ! Group 1: UUR, wrap
        call Op_K%mmult_R_group(Propu%UUR, Latt%group_1, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%UUR, Latt%group_1, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propu%Gr, Latt%group_1, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%Gr, Latt%group_1, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propu%Gr, Latt%group_1, Latt, sigma, nt, -1)
        call Op_K%mmult_L_group(Propd%Gr, Latt%group_1, Latt, sigma, nt, -1)
        ! Group 2
        call Op_K%mmult_R_group(Propu%UUR, Latt%group_2, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%UUR, Latt%group_2, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propu%Gr, Latt%group_2, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%Gr, Latt%group_2, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propu%Gr, Latt%group_2, Latt, sigma, nt, -1)
        call Op_K%mmult_L_group(Propd%Gr, Latt%group_2, Latt, sigma, nt, -1)
        ! Group 3
        call Op_K%mmult_R_group(Propu%UUR, Latt%group_3, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%UUR, Latt%group_3, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propu%Gr, Latt%group_3, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%Gr, Latt%group_3, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propu%Gr, Latt%group_3, Latt, sigma, nt, -1)
        call Op_K%mmult_L_group(Propd%Gr, Latt%group_3, Latt, sigma, nt, -1)
        ! Group 4
        call Op_K%mmult_R_group(Propu%UUR, Latt%group_4, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%UUR, Latt%group_4, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propu%Gr, Latt%group_4, Latt, sigma, nt, 1)
        call Op_K%mmult_R_group(Propd%Gr, Latt%group_4, Latt, sigma, nt, 1)
        call Op_K%mmult_L_group(Propu%Gr, Latt%group_4, Latt, sigma, nt, -1)
        call Op_K%mmult_L_group(Propd%Gr, Latt%group_4, Latt, sigma, nt, -1)
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
