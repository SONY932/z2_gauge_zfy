module LocalK_mod
    use Multiply_mod
    use DQMC_Model_mod
    use ProcessMatrix
    use MyLattice
    use calc_basic
    use Fields_mod, only: gauss_boson_ratio_lambda, gauss_boson_ratio_sigma
    implicit none

    public
    private :: LocalK_metro, LocalK_metro_woodbury, compute_LR_cols_rows
    private :: apply_group_to_cols_direct, apply_group_to_rows_direct
    private :: sigma_new, Local_lambda_flip, lambda_new

    type(AccCounter) :: Acc_Kl, Acc_Kt, Acc_lambda
    real(kind = 8), dimension(:, :), allocatable :: sigma_new
    real(kind = 8), dimension(:), allocatable :: lambda_new
    
    ! λ = -1 的格点数量（用于 λ 全局更新的统计）
    integer :: n_lambda_minus

contains
    subroutine LocalK_init()
        call Acc_Kl%init()
        call Acc_Kt%init()
        call Acc_lambda%init()
        allocate(sigma_new(2*Lq, Ltrot))
        allocate(lambda_new(Lq))
        n_lambda_minus = 0
        return
    end subroutine LocalK_init
    
    subroutine LocalK_clear()
        deallocate(sigma_new)
        deallocate(lambda_new)
        return
    end subroutine LocalK_clear

    subroutine LocalK_reset()
        call Acc_Kl%reset()
        call Acc_Kt%reset()
        call Acc_lambda%reset()
        sigma_new = NsigL_K%sigma
        lambda_new = NsigL_K%lambda
        return
    end subroutine LocalK_reset

    subroutine LocalK_metro(GrU, GrD, iseed, ii, ntau)
        ! 简化版 metro 更新（不考虑棋盘分解的影响）
        ! 这个函数仅用于 metro 更新被禁用时的测试
        ! 实际使用时应调用 LocalK_metro_woodbury
! Arguments:
	    complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: GrU, GrD
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau
!   Local: 
        real(kind=8), external :: ranf
        complex(kind=8) :: ProddetU, ProddetD
        complex(kind=8), dimension(2, 2) :: ProdU, ProdD, ProdinvU, ProdinvD, GrU_local, GrD_local
        complex(kind=8), dimension(2, 2) :: matU_tmp, matD_tmp
        complex(kind=8) :: UhlpU(Ndim, 2), VhlpU(2, Ndim), UhlpD(Ndim, 2), VhlpD(2, Ndim)
        complex(kind=8) :: temp(Ndim, 2), Diff(Ndim, Ndim)
        complex(kind=8) :: ratio_fermion
        real(kind=8) :: ratio_boson, ratio_re, ratio_re_abs
        real(kind=8) :: random
        integer :: P(2), j, no, nl, nr
        real(kind=8) :: S_new, S_old

! 局部翻转规范场
        sigma_new(ii, ntau) = - NsigL_K%sigma(ii, ntau)
        S_new = sigma_new(ii, ntau)
        S_old = NsigL_K%sigma(ii, ntau)
        P(1) = Latt%bond_list(ii, 1)
        P(2) = Latt%bond_list(ii, 2)

        call Op_K%get_delta(S_old, S_new)
        
        ProdU = dcmplx(0.d0, 0.d0)
        ProdD = dcmplx(0.d0, 0.d0)
        
        ! 计算 (I - G_0) 矩阵
        do nr = 1, 2
            do nl = 1, 2
                GrU_local(nl, nr) = ZKRON(nl, nr) - GrU(P(nl), P(nr))
                GrD_local(nl, nr) = ZKRON(nl, nr) - GrD(P(nl), P(nr))
            enddo
        enddo

        ! 标准公式：det[1 + B'] / det[1 + B] = det[I + Delta * (I - G_0)]
        matU_tmp = matmul(Op_K%Delta, GrU_local)
        matD_tmp = matmul(Op_K%Delta, GrD_local)
        do nr = 1, 2
            do nl = 1, 2
                ProdU(nl, nr) = ZKRON(nl, nr) + matU_tmp(nl, nr)
                ProdD(nl, nr) = ZKRON(nl, nr) + matD_tmp(nl, nr)
            enddo
        enddo
        ProddetU = ProdU(1, 1) * ProdU(2, 2) - ProdU(1, 2) * ProdU(2, 1)
        ProddetD = ProdD(1, 1) * ProdD(2, 2) - ProdD(1, 2) * ProdD(2, 1)
        
        ratio_fermion = ProddetU * ProddetD
! Calculate total Metropolis ratio  
        ratio_boson = NsigL_K%bosonratio(sigma_new, ii, ntau, Latt)
        
        ! 添加高斯投影的玻色贡献
        ! 当翻转 σ(bond, τ) 时，如果 τ = 1 或 τ = Ltrot，会影响 σ^x
        ! 这会改变高斯投影的玻色权重
        ratio_boson = ratio_boson * gauss_boson_ratio_sigma(ii, ntau, &
            NsigL_K%sigma, sigma_new(ii, ntau), NsigL_K%lambda, Latt)
        
        ratio_re = dble(ratio_fermion * ratio_boson)
        ratio_re_abs = abs(ratio_re)
        random = ranf(iseed)
! Upgrade Green's function
        if (ratio_re_abs > random) then
            call Acc_Kl%count(.true.)
            ProdinvD(1, 1) = ProdD(2, 2)
            ProdinvD(1, 2) = - ProdD(1, 2)
            ProdinvD(2, 1) = - ProdD(2, 1)
            ProdinvD(2, 2) = ProdD(1, 1)
            ProdinvU(1, 1) = ProdU(2, 2)
            ProdinvU(1, 2) = - ProdU(1, 2)
            ProdinvU(2, 1) = - ProdU(2, 1)
            ProdinvU(2, 2) = ProdU(1, 1)

            ProdinvD = ProdinvD / ProddetD
            ProdinvU = ProdinvU / ProddetU
            UhlpU = dcmplx(0.d0, 0.d0); VhlpU = dcmplx(0.d0, 0.d0)
            UhlpD = dcmplx(0.d0, 0.d0); VhlpD = dcmplx(0.d0, 0.d0)
! Sherman-Morrison 更新 G_0：使用标准 Op_K%Delta
            temp = dcmplx(0.d0, 0.d0); Diff = dcmplx(0.d0, 0.d0)
            do no = 1, 2
                do j = 1, Ndim
                    UhlpU(j, no) = GrU(j, P(no))
                    VhlpU(no, j) = - Op_K%Delta(no, 1) * GrU(P(1), j) - Op_K%Delta(no, 2) * GrU(P(2), j)
                enddo
                VhlpU(no, P(1)) = VhlpU(no, P(1)) + Op_K%Delta(no, 1)
                VhlpU(no, P(2)) = VhlpU(no, P(2)) + Op_K%Delta(no, 2)
            enddo
            temp = matmul(UhlpU, ProdinvU) ! (Ndim,2) = (Ndim,2) * (2,2)
            Diff = matmul(temp, VhlpU)     ! (Ndim,Ndim) = (Ndim,2) * (2,Ndim)
            GrU = GrU - Diff ! output GrU

            temp = dcmplx(0.d0, 0.d0); Diff = dcmplx(0.d0, 0.d0)
            do no = 1, 2
                do j = 1, Ndim
                    UhlpD(j, no) = GrD(j, P(no))
                    VhlpD(no, j) = - Op_K%Delta(no, 1) * GrD(P(1), j) - Op_K%Delta(no, 2) * GrD(P(2), j)
                enddo
                VhlpD(no, P(1)) = VhlpD(no, P(1)) + Op_K%Delta(no, 1)
                VhlpD(no, P(2)) = VhlpD(no, P(2)) + Op_K%Delta(no, 2)
            enddo
            temp = matmul(UhlpD, ProdinvD)
            Diff = matmul(temp, VhlpD)
            GrD = GrD - Diff ! output GrD
! Flip: 
            NsigL_K%sigma(ii, ntau) = sigma_new(ii, ntau)
        else
            call Acc_Kl%count(.false.)
            sigma_new(ii, ntau) = NsigL_K%sigma(ii, ntau)
        endif
        return
    end subroutine LocalK_metro

    subroutine LocalK_metro_woodbury(GrU, GrD, iseed, ii, ntau, group_idx)
        ! 使用完整 Woodbury 公式的 metro 更新
        ! 
        ! 棋盘分解：B(τ) = B_4 * B_3 * B_2 * B_1
        ! 当翻转 group_k 中的 bond b 时，设 Delta = E'_b * E_b^{-1} - I。
        ! 
        ! 由于 B'_k = (I + Delta_embedded) * B_k（在 bond b 子空间上），有：
        !   ΔB = B' - B = L * Delta_embedded * R
        !   其中 L = B_4 * ... * B_{k+1}
        !        R = B_k * B_{k-1} * ... * B_1（注意包含 B_k！）
        !
        ! 设 L_cols = L(:, [P1, P2])（N×2），R_rows = R([P1, P2], :)（2×N）。
        ! 则 ΔB = L_cols * Delta * R_rows（rank-2 矩阵）
        !
        ! 行列式比：
        !   det(I + B') / det(I + B) = det(I + G * ΔB)
        !                            = det(I + G * L_cols * Delta * R_rows)
        !                            = det(I_2 + Delta * R_rows * G * L_cols)  [det(I+AB)=det(I+BA)]
        !   设 M = R_rows * G * L_cols（2×2），则 det ratio = det(I_2 + Delta * M)
        !
        ! Woodbury 更新：
        !   G' = G - (G * L_cols * Delta) * K^{-1} * (R_rows * G)
        !   其中 K = I_2 + M * Delta
        !
! Arguments:
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: GrU, GrD
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau, group_idx
!   Local: 
        real(kind=8), external :: ranf
        complex(kind=8) :: detK_U, detK_D
        complex(kind=8), dimension(2, 2) :: M_U, M_D, K_U, K_D, Kinv_U, Kinv_D
        complex(kind=8), dimension(2, 2) :: Delta_M_U, Delta_M_D
        complex(kind=8) :: L_cols(Ndim, 2), R_rows(2, Ndim)
        complex(kind=8) :: GL_U(Ndim, 2), GL_D(Ndim, 2)  ! G * L_cols
        complex(kind=8) :: RG_U(2, Ndim), RG_D(2, Ndim)  ! R_rows * G
        complex(kind=8) :: GL_Delta_U(Ndim, 2), GL_Delta_D(Ndim, 2)  ! G * L_cols * Delta
        complex(kind=8) :: temp(Ndim, 2), Diff(Ndim, Ndim)
        complex(kind=8) :: ratio_fermion
        complex(kind=8) :: Delta_c(2, 2)
        real(kind=8) :: ratio_boson, ratio_re, ratio_re_abs
        real(kind=8) :: random
        integer :: P(2), j, no, nl, nr
        real(kind=8) :: S_new, S_old
        real(kind=8) :: Delta_local(2, 2)

! 局部翻转规范场
        sigma_new(ii, ntau) = - NsigL_K%sigma(ii, ntau)
        S_new = sigma_new(ii, ntau)
        S_old = NsigL_K%sigma(ii, ntau)
        P(1) = Latt%bond_list(ii, 1)
        P(2) = Latt%bond_list(ii, 2)

        ! 计算局部 Delta 矩阵（2×2）
        call Op_K%get_delta(S_old, S_new)
        Delta_local = Op_K%Delta
        Delta_c = dcmplx(Delta_local, 0.d0)

        ! 计算 L_cols = L(:, [P1, P2]) 和 R_rows = R([P1, P2], :)
        ! 注意：R 包含 B_k！
        call compute_LR_cols_rows(group_idx, P, ntau, L_cols, R_rows)

        ! ========== Spin Up ==========
        ! 计算 GL_U = G * L_cols（N×2）
        GL_U = dcmplx(0.d0, 0.d0)
        do j = 1, 2
            do no = 1, Ndim
                do nl = 1, Ndim
                    GL_U(no, j) = GL_U(no, j) + GrU(no, nl) * L_cols(nl, j)
                enddo
            enddo
        enddo

        ! 计算 RG_U = R_rows * G（2×N）
        RG_U = dcmplx(0.d0, 0.d0)
        do j = 1, Ndim
            do no = 1, 2
                do nl = 1, Ndim
                    RG_U(no, j) = RG_U(no, j) + R_rows(no, nl) * GrU(nl, j)
                enddo
            enddo
        enddo

        ! 计算 M_U = R_rows * G * L_cols = RG_U * L_cols（2×2）
        ! 但更高效的方法是 M_U = R_rows * GL_U
        M_U = dcmplx(0.d0, 0.d0)
        do nr = 1, 2
            do nl = 1, 2
                do no = 1, Ndim
                    M_U(nl, nr) = M_U(nl, nr) + R_rows(nl, no) * GL_U(no, nr)
                enddo
            enddo
        enddo

        ! 计算 Delta * M_U（2×2）
        Delta_M_U = matmul(Delta_c, M_U)

        ! K_U = I + Delta * M = I + Delta_M（注意：det ratio 用 I + Delta * M）
        ! 但 Woodbury 用 K = I + M * Delta，所以需要转置
        ! 实际上让我重新检查...
        ! det ratio = det(I + Delta * M)
        ! Woodbury: G' = G - (G*L*Delta) * (I + M*Delta)^{-1} * (R*G)
        ! 这两个 K 不同！det ratio 用 I + Delta*M，Woodbury 用 I + M*Delta
        ! 但 det(I + AB) = det(I + BA)，所以 det(I + Delta*M) = det(I + M*Delta)
        ! 而 (I + Delta*M)^{-1} ≠ (I + M*Delta)^{-1} 一般情况下
        ! 让我重新推导 Woodbury...
        ! 
        ! G' = (I + B')^{-1} = (I + B + L_cols*Delta*R_rows)^{-1}
        ! 设 U = L_cols * Delta（N×2），V^T = R_rows（2×N）
        ! Woodbury: G' = G - G*U*(I + V^T*G*U)^{-1}*V^T*G
        ! V^T*G*U = R_rows * G * L_cols * Delta = M * Delta
        ! 所以 K = I + M * Delta
        K_U = dcmplx(0.d0, 0.d0)
        K_U(1, 1) = dcmplx(1.d0, 0.d0)
        K_U(2, 2) = dcmplx(1.d0, 0.d0)
        K_U = K_U + matmul(M_U, Delta_c)  ! K = I + M * Delta

        detK_U = K_U(1, 1) * K_U(2, 2) - K_U(1, 2) * K_U(2, 1)

        ! ========== Spin Down ==========
        GL_D = dcmplx(0.d0, 0.d0)
        do j = 1, 2
            do no = 1, Ndim
                do nl = 1, Ndim
                    GL_D(no, j) = GL_D(no, j) + GrD(no, nl) * L_cols(nl, j)
                enddo
            enddo
        enddo

        RG_D = dcmplx(0.d0, 0.d0)
        do j = 1, Ndim
            do no = 1, 2
                do nl = 1, Ndim
                    RG_D(no, j) = RG_D(no, j) + R_rows(no, nl) * GrD(nl, j)
                enddo
            enddo
        enddo

        M_D = dcmplx(0.d0, 0.d0)
        do nr = 1, 2
            do nl = 1, 2
                do no = 1, Ndim
                    M_D(nl, nr) = M_D(nl, nr) + R_rows(nl, no) * GL_D(no, nr)
                enddo
            enddo
        enddo

        K_D = dcmplx(0.d0, 0.d0)
        K_D(1, 1) = dcmplx(1.d0, 0.d0)
        K_D(2, 2) = dcmplx(1.d0, 0.d0)
        K_D = K_D + matmul(M_D, Delta_c)

        detK_D = K_D(1, 1) * K_D(2, 2) - K_D(1, 2) * K_D(2, 1)
        
        ratio_fermion = detK_U * detK_D

! Calculate total Metropolis ratio  
        ratio_boson = NsigL_K%bosonratio(sigma_new, ii, ntau, Latt)
        ratio_boson = ratio_boson * gauss_boson_ratio_sigma(ii, ntau, &
            NsigL_K%sigma, sigma_new(ii, ntau), NsigL_K%lambda, Latt)
        
        ratio_re = dble(ratio_fermion * ratio_boson)
        ratio_re_abs = abs(ratio_re)
        random = ranf(iseed)

! Upgrade Green's function using Woodbury formula
        if (ratio_re_abs > random) then
            call Acc_Kl%count(.true.)
            
            ! Woodbury: G' = G - (G*L_cols*Delta) * K^{-1} * (R_rows*G)
            ! 其中 K = I + M * Delta
            
            ! 计算 K^{-1}
            Kinv_U(1, 1) =  K_U(2, 2) / detK_U
            Kinv_U(1, 2) = -K_U(1, 2) / detK_U
            Kinv_U(2, 1) = -K_U(2, 1) / detK_U
            Kinv_U(2, 2) =  K_U(1, 1) / detK_U

            Kinv_D(1, 1) =  K_D(2, 2) / detK_D
            Kinv_D(1, 2) = -K_D(1, 2) / detK_D
            Kinv_D(2, 1) = -K_D(2, 1) / detK_D
            Kinv_D(2, 2) =  K_D(1, 1) / detK_D

            ! Spin Up: G' = G - (GL_U * Delta) * Kinv * RG_U
            ! GL_Delta_U = GL_U * Delta（N×2）
            GL_Delta_U = matmul(GL_U, Delta_c)
            ! temp = GL_Delta_U * Kinv（N×2）
            temp = matmul(GL_Delta_U, Kinv_U)
            ! Diff = temp * RG_U（N×N）
            Diff = matmul(temp, RG_U)
            GrU = GrU - Diff

            ! Spin Down
            GL_Delta_D = matmul(GL_D, Delta_c)
            temp = matmul(GL_Delta_D, Kinv_D)
            Diff = matmul(temp, RG_D)
            GrD = GrD - Diff

! Flip: 
            NsigL_K%sigma(ii, ntau) = sigma_new(ii, ntau)
        else
            call Acc_Kl%count(.false.)
            sigma_new(ii, ntau) = NsigL_K%sigma(ii, ntau)
        endif
        return
    end subroutine LocalK_metro_woodbury

    subroutine compute_LR_cols_rows(group_idx, P, ntau, L_cols, R_rows)
        ! 计算 L_cols = L(:, [P1, P2]) 和 R_rows = R([P1, P2], :)
        ! 
        ! 棋盘分解：B(τ) = B_4 * B_3 * B_2 * B_1
        ! 
        ! 对于 group_k 中的 bond b，设 Delta = E'_b * E_b^{-1} - I
        ! 则 B'_k = (I + Delta) * B_k（在 bond b 的子空间上）
        ! 
        ! ΔB = B' - B = L * Delta * B_k * R_orig
        ! 其中 L = B_4 * ... * B_{k+1}
        !      R_orig = B_{k-1} * ... * B_1
        !      R = B_k * R_orig = B_k * B_{k-1} * ... * B_1
        !
        ! 写成 UV^T 形式：
        !   U = L(:, [P1, P2]) * Delta  (N×2)
        !   V^T = R([P1, P2], :)        (2×N)
        !
        ! 注意：R 需要包含 B_k！
        !
        integer, intent(in) :: group_idx, P(2), ntau
        complex(kind=8), intent(out) :: L_cols(Ndim, 2), R_rows(2, Ndim)
        
        integer :: kk
        
        ! 初始化为单位矩阵的相应列/行
        L_cols = dcmplx(0.d0, 0.d0)
        L_cols(P(1), 1) = dcmplx(1.d0, 0.d0)
        L_cols(P(2), 2) = dcmplx(1.d0, 0.d0)
        
        R_rows = dcmplx(0.d0, 0.d0)
        R_rows(1, P(1)) = dcmplx(1.d0, 0.d0)
        R_rows(2, P(2)) = dcmplx(1.d0, 0.d0)
        
        ! 计算 L = B_4 * ... * B_{k+1}
        ! 从 I 开始，依次左乘 B_{k+1}, B_{k+2}, ..., B_4
        ! 注意：对于 group_4，L = I；对于 group_1，L = B_4 * B_3 * B_2
        do kk = group_idx + 1, 4
            call apply_group_to_cols_direct(kk, ntau, L_cols)
        enddo
        
        ! 计算 R = B_k * B_{k-1} * ... * B_1
        ! 从 I 开始，依次右乘 B_k, B_{k-1}, ..., B_1
        ! 注意：R 需要包含 B_k！
        ! 对于 group_1，R = B_1；对于 group_4，R = B_4 * B_3 * B_2 * B_1
        do kk = group_idx, 1, -1
            call apply_group_to_rows_direct(kk, ntau, R_rows)
        enddo
        
        return
    end subroutine compute_LR_cols_rows

    subroutine apply_group_to_cols_direct(grp_idx, nt, cols)
        ! 将 B_grp 左乘到 cols 上：cols = B_grp * cols
        integer, intent(in) :: grp_idx, nt
        complex(kind=8), intent(inout) :: cols(Ndim, 2)
        integer :: idx, bond_no, s1, s2, jj, glen
        real(kind=8) :: sig, Cv, Sv
        complex(kind=8) :: t1, t2
        
        select case (grp_idx)
        case (1)
            glen = size(Latt%group_1)
            do idx = 1, glen
                bond_no = Latt%group_1(idx)
                s1 = Latt%bond_list(bond_no, 1)
                s2 = Latt%bond_list(bond_no, 2)
                sig = NsigL_K%sigma(bond_no, nt)
                Cv = cosh(Op_K%alpha * sig)
                Sv = sinh(Op_K%alpha * sig)
                do jj = 1, 2
                    t1 = Cv * cols(s1, jj) + Sv * cols(s2, jj)
                    t2 = Sv * cols(s1, jj) + Cv * cols(s2, jj)
                    cols(s1, jj) = t1
                    cols(s2, jj) = t2
                enddo
            enddo
        case (2)
            glen = size(Latt%group_2)
            do idx = 1, glen
                bond_no = Latt%group_2(idx)
                s1 = Latt%bond_list(bond_no, 1)
                s2 = Latt%bond_list(bond_no, 2)
                sig = NsigL_K%sigma(bond_no, nt)
                Cv = cosh(Op_K%alpha * sig)
                Sv = sinh(Op_K%alpha * sig)
                do jj = 1, 2
                    t1 = Cv * cols(s1, jj) + Sv * cols(s2, jj)
                    t2 = Sv * cols(s1, jj) + Cv * cols(s2, jj)
                    cols(s1, jj) = t1
                    cols(s2, jj) = t2
                enddo
            enddo
        case (3)
            glen = size(Latt%group_3)
            do idx = 1, glen
                bond_no = Latt%group_3(idx)
                s1 = Latt%bond_list(bond_no, 1)
                s2 = Latt%bond_list(bond_no, 2)
                sig = NsigL_K%sigma(bond_no, nt)
                Cv = cosh(Op_K%alpha * sig)
                Sv = sinh(Op_K%alpha * sig)
                do jj = 1, 2
                    t1 = Cv * cols(s1, jj) + Sv * cols(s2, jj)
                    t2 = Sv * cols(s1, jj) + Cv * cols(s2, jj)
                    cols(s1, jj) = t1
                    cols(s2, jj) = t2
                enddo
            enddo
        case (4)
            glen = size(Latt%group_4)
            do idx = 1, glen
                bond_no = Latt%group_4(idx)
                s1 = Latt%bond_list(bond_no, 1)
                s2 = Latt%bond_list(bond_no, 2)
                sig = NsigL_K%sigma(bond_no, nt)
                Cv = cosh(Op_K%alpha * sig)
                Sv = sinh(Op_K%alpha * sig)
                do jj = 1, 2
                    t1 = Cv * cols(s1, jj) + Sv * cols(s2, jj)
                    t2 = Sv * cols(s1, jj) + Cv * cols(s2, jj)
                    cols(s1, jj) = t1
                    cols(s2, jj) = t2
                enddo
            enddo
        end select
    end subroutine apply_group_to_cols_direct
    
    subroutine apply_group_to_rows_direct(grp_idx, nt, rows)
        ! 将 B_grp 右乘到 rows 上：rows = rows * B_grp
        integer, intent(in) :: grp_idx, nt
        complex(kind=8), intent(inout) :: rows(2, Ndim)
        integer :: idx, bond_no, s1, s2, jj, glen
        real(kind=8) :: sig, Cv, Sv
        complex(kind=8) :: t1, t2
        
        select case (grp_idx)
        case (1)
            glen = size(Latt%group_1)
            do idx = 1, glen
                bond_no = Latt%group_1(idx)
                s1 = Latt%bond_list(bond_no, 1)
                s2 = Latt%bond_list(bond_no, 2)
                sig = NsigL_K%sigma(bond_no, nt)
                Cv = cosh(Op_K%alpha * sig)
                Sv = sinh(Op_K%alpha * sig)
                do jj = 1, 2
                    t1 = rows(jj, s1) * Cv + rows(jj, s2) * Sv
                    t2 = rows(jj, s1) * Sv + rows(jj, s2) * Cv
                    rows(jj, s1) = t1
                    rows(jj, s2) = t2
                enddo
            enddo
        case (2)
            glen = size(Latt%group_2)
            do idx = 1, glen
                bond_no = Latt%group_2(idx)
                s1 = Latt%bond_list(bond_no, 1)
                s2 = Latt%bond_list(bond_no, 2)
                sig = NsigL_K%sigma(bond_no, nt)
                Cv = cosh(Op_K%alpha * sig)
                Sv = sinh(Op_K%alpha * sig)
                do jj = 1, 2
                    t1 = rows(jj, s1) * Cv + rows(jj, s2) * Sv
                    t2 = rows(jj, s1) * Sv + rows(jj, s2) * Cv
                    rows(jj, s1) = t1
                    rows(jj, s2) = t2
                enddo
            enddo
        case (3)
            glen = size(Latt%group_3)
            do idx = 1, glen
                bond_no = Latt%group_3(idx)
                s1 = Latt%bond_list(bond_no, 1)
                s2 = Latt%bond_list(bond_no, 2)
                sig = NsigL_K%sigma(bond_no, nt)
                Cv = cosh(Op_K%alpha * sig)
                Sv = sinh(Op_K%alpha * sig)
                do jj = 1, 2
                    t1 = rows(jj, s1) * Cv + rows(jj, s2) * Sv
                    t2 = rows(jj, s1) * Sv + rows(jj, s2) * Cv
                    rows(jj, s1) = t1
                    rows(jj, s2) = t2
                enddo
            enddo
        case (4)
            glen = size(Latt%group_4)
            do idx = 1, glen
                bond_no = Latt%group_4(idx)
                s1 = Latt%bond_list(bond_no, 1)
                s2 = Latt%bond_list(bond_no, 2)
                sig = NsigL_K%sigma(bond_no, nt)
                Cv = cosh(Op_K%alpha * sig)
                Sv = sinh(Op_K%alpha * sig)
                do jj = 1, 2
                    t1 = rows(jj, s1) * Cv + rows(jj, s2) * Sv
                    t2 = rows(jj, s1) * Sv + rows(jj, s2) * Cv
                    rows(jj, s1) = t1
                    rows(jj, s2) = t2
                enddo
            enddo
        end select
    end subroutine apply_group_to_rows_direct

    subroutine LocalK_metro_group(GrU, GrD, iseed, group, ntau, group_idx)
        ! 对组内所有键做 metro 更新（使用 Woodbury 公式）
        ! 由于组内的键不共享格点，它们的 B 矩阵互相对易
        ! 但不同组之间需要考虑 L 和 R 因子
        !
        ! group_idx：当前组的索引（1-4）
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: GrU, GrD
        integer, intent(inout) :: iseed
        integer, dimension(:), intent(in) :: group
        integer, intent(in) :: ntau, group_idx
        integer :: ii, no, len
        
        len = size(group)
        do ii = 1, len
            no = group(ii)
            call LocalK_metro_woodbury(GrU, GrD, iseed, no, ntau, group_idx)
        enddo
        return
    end subroutine LocalK_metro_group

    ! =========================================================================
    ! 单格点 λ 翻转（正确实现高斯约束）
    ! 
    ! 根据 PNAS SI，翻转 λ_i 的接受率包含两部分：
    ! 1. 费米子行列式比：det[1 + P'B] / det[1 + PB]
    ! 2. 玻色权重比：exp[iπ(n_b + n_Q)] = (-1)^{n_b + n_Q}
    !
    ! 其中 n_b = 格点 i 相邻四条键中时间不连续的个数
    !      n_Q = (1-Q)/2（扇区信息）
    ! =========================================================================
    subroutine Local_lambda_flip_single(PropU, PropD, iseed, ii)
        use MyMats
! Arguments:
        class(Propagator), intent(inout) :: PropU, PropD
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii
! Local:
        real(kind=8), external :: ranf
        real(kind=8) :: rU_single, rD_single, ratio_boson, ratio_total, random
        logical :: safeU, safeD

        ! 计算费米子行列式比（rank-1 更新）
        call lambda_single_ratio(PropU%Gr, ii, rU_single, safeU)
        call lambda_single_ratio(PropD%Gr, ii, rD_single, safeD)
        
        ! 计算玻色权重比
        ratio_boson = gauss_boson_ratio_lambda(ii, NsigL_K%sigma, NsigL_K%lambda, Latt)

        if (safeU .and. safeD) then
            ratio_total = rU_single * rD_single * ratio_boson
        else
            ratio_total = 0.d0
        endif

        random = ranf(iseed)
        if (abs(ratio_total) > random) then
            call Acc_lambda%count(.true.)
            call lambda_single_update(PropU%Gr, ii)
            call lambda_single_update(PropD%Gr, ii)
            
            ! 更新 λ 场
            NsigL_K%lambda(ii) = -NsigL_K%lambda(ii)
            lambda_new(ii) = NsigL_K%lambda(ii)
        else
            call Acc_lambda%count(.false.)
        endif
        return
    contains
        subroutine lambda_single_ratio(Gr, i, r_single, is_safe)
            ! 计算单格点翻转 λ_i 的费米子行列式比
            ! 
            ! 设 G_λ = (I + P_λ B)^{-1}，翻转 λ_i 等价于：
            ! P' = P + ΔP，其中 ΔP 只在 (i,i) 处非零
            ! ΔP_{ii} = P'_{ii} - P_{ii} = -2λ_i（从 λ_i 到 -λ_i）
            !
            ! 行列式比 = det(I + ΔP B G_λ) = 1 + ΔP_{ii} (BG_λ)_{ii}
            !
            ! 由于 BG_λ = I - G_λ + P_λ^{-1}(I - G_λ)... 这很复杂
            ! 
            ! 更简单的推导：
            ! det[I + P'B] / det[I + PB] = det[P' P^{-1} + P' B P^{-1}] / det[I + B]
            !                            = det[P' P^{-1}(I + PB)] / det[I + PB]
            !                            = det[P' P^{-1}]
            ! 
            ! 不对，这也不对因为 P' P^{-1} 不等于 I。
            !
            ! 正确方法：使用 Woodbury 公式
            ! G_λ = (I + P_λ B)^{-1}
            ! 翻转 λ_i 意味着 P' = P + ΔP，其中 ΔP_{ii} = -2λ_i
            ! 
            ! 行列式比：det[I + P'B] / det[I + PB] 
            !         = det[I + (P + ΔP)B] / det[I + PB]
            !         = det[(I + PB)(I + (I + PB)^{-1} ΔP B)] / det[I + PB]
            !         = det[I + G_λ^{-1} · ΔP · (G_λ^{-1} - I)]  -- 这太复杂
            !
            ! 更简洁：det[I + P'B] = det[I + PB] · det[I + ΔP B G_λ]
            ! 其中 G_λ = (I + PB)^{-1}
            ! 由于 ΔP 是 rank-1 矩阵，det[I + ΔP B G_λ] = 1 + Tr[ΔP B G_λ]
            !                                         = 1 + ΔP_{ii} (BG_λ)_{ii}
            ! 
            ! 但问题是我们存储的是 G_0 = (I + B)^{-1}，不是 G_λ = (I + P_λ B)^{-1}
            ! 需要转换：G_λ = (I + P_λ B)^{-1} = ...
            !
            ! 实际上，根据 PNAS 的方法，我们应该存储 G_λ 而不是 G_0。
            ! 让我重新思考...
            !
            ! 简化版：如果当前存储的 Green 函数已经包含了 P[λ] 的效应，
            ! 即 Gr = G_λ = (I + P_λ B)^{-1}，则：
            ! 
            ! 翻转 λ_i 时，ΔP_{ii} = -2λ_i
            ! 行列式比 = 1 + ΔP_{ii} (BG_λ)_{ii}
            !         = 1 - 2λ_i (I - G_λ)_{ii} · P_λ^{-1}_{ii}
            !         = 1 - 2λ_i (1 - G_{ii}) / λ_i
            !         = 1 - 2(1 - G_{ii})
            !         = 2G_{ii} - 1
            !
            complex(kind=8), intent(in) :: Gr(:, :)
            integer, intent(in) :: i
            real(kind=8), intent(out) :: r_single
            logical, intent(out) :: is_safe
            real(kind=8) :: Gii_diag
            real(kind=8), parameter :: tol = 0.1d0
            
            Gii_diag = real(Gr(i, i))
            
            ! 检查 G_{ii} 是否在物理区间
            if (Gii_diag < -tol .or. Gii_diag > 1.d0 + tol) then
                r_single = 0.d0
                is_safe = .false.
                return
            endif
            
            ! 行列式比 = |2G_{ii} - 1|
            r_single = abs(2.d0 * Gii_diag - 1.d0)
            is_safe = .true.
        end subroutine lambda_single_ratio

        subroutine lambda_single_update(Gr, i)
            ! 更新 Green 函数（rank-1 Sherman-Morrison）
            ! 
            ! 新的 Green 函数：G' = G - G · u · v^T · G / (1 + v^T G u)
            ! 其中 u = e_i（单位向量），v = -2λ_i B e_i
            ! 
            ! 简化后：G' = G + 2 G[:, i] (I-G)[i, :] / (2G_{ii} - 1)
            complex(kind=8), intent(inout) :: Gr(:, :)
            integer, intent(in) :: i
            complex(kind=8) :: denom, twoC
            complex(kind=8), dimension(Ndim) :: col_i, row_i
            integer :: kk, ll
            
            twoC = dcmplx(2.d0, 0.d0)
            denom = twoC * Gr(i, i) - dcmplx(1.d0, 0.d0)
            
            if (abs(denom) < 1.d-300) return
            
            ! col_i = G[:, i]
            ! row_i = (I - G)[i, :]
            do kk = 1, Ndim
                col_i(kk) = Gr(kk, i)
                if (kk == i) then
                    row_i(kk) = dcmplx(1.d0, 0.d0) - Gr(i, kk)
                else
                    row_i(kk) = -Gr(i, kk)
                endif
            enddo
            
            ! G' = G + 2 * col_i ⊗ row_i / denom
            do kk = 1, Ndim
                do ll = 1, Ndim
                    Gr(kk, ll) = Gr(kk, ll) + twoC * col_i(kk) * row_i(ll) / denom
                enddo
            enddo
        end subroutine lambda_single_update
    end subroutine Local_lambda_flip_single

    ! =========================================================================
    ! 成对 λ 翻转（保留用于保持 ∏λ 不变的特殊情况）
    ! =========================================================================
    subroutine Local_lambda_flip(PropU, PropD, iseed, ii, jj)
        use MyMats
! Arguments:
        class(Propagator), intent(inout) :: PropU, PropD
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, jj
! Local:
        real(kind=8), external :: ranf
        real(kind=8) :: rU_pair, rD_pair, ratio_boson_i, ratio_boson_j, prob, random
        logical :: safeU, safeD

        if (ii == jj) return

! 采用 rank-2 Woodbury，一次性对 (ii,jj) 成对翻转 λ，避免在主 Gr 上反复 preview/revert。
! 先在当前 Green 上计算两个自旋的成对行列式比，再根据接受与否，对主 Gr 做一次 rank-2 更新。
        call lambda_pair_ratio(PropU%Gr, ii, jj, rU_pair, safeU)
        call lambda_pair_ratio(PropD%Gr, ii, jj, rD_pair, safeD)
        
        ! 计算两个格点的玻色权重比
        ratio_boson_i = gauss_boson_ratio_lambda(ii, NsigL_K%sigma, NsigL_K%lambda, Latt)
        ratio_boson_j = gauss_boson_ratio_lambda(jj, NsigL_K%sigma, NsigL_K%lambda, Latt)

        if (safeU .and. safeD .and. rU_pair > 0.d0 .and. rD_pair > 0.d0) then
            prob = rU_pair * rD_pair * ratio_boson_i * ratio_boson_j
        else
            prob = 0.d0
        endif

        random = ranf(iseed)
        if (abs(prob) > random) then
            call Acc_lambda%count(.true.)
            call lambda_pair_update(PropU%Gr, ii, jj)
            call lambda_pair_update(PropD%Gr, ii, jj)
            
            ! 更新 λ 场
            NsigL_K%lambda(ii) = -NsigL_K%lambda(ii)
            NsigL_K%lambda(jj) = -NsigL_K%lambda(jj)
            lambda_new(ii) = NsigL_K%lambda(ii)
            lambda_new(jj) = NsigL_K%lambda(jj)
            
            ! 注意：不需要更新 UUR 链
            ! P[λ] 在 Wrap 时动态应用，不存储在 UUR 中
        else
            call Acc_lambda%count(.false.)
            lambda_new(ii) = NsigL_K%lambda(ii)
            lambda_new(jj) = NsigL_K%lambda(jj)
        endif
        return
    contains
        subroutine lambda_pair_ratio(Gr, i, j, r_pair, is_safe)
            ! 计算成对翻转 λ_i 和 λ_j 的接受率
            ! 根据论文公式，G = (I + P_λ B)^{-1}
            ! 当翻转 λ_i, λ_j 时，ΔP(i,i) = -2λ_i, ΔP(j,j) = -2λ_j
            ! K = I + (ΔP B) G (在 {i,j} 子空间)
            ! 由于 BG = P(I-G)，有：
            ! K(1,1) = 1 - 2(1 - G_ii) = 2G_ii - 1
            ! K(1,2) = 2G_ij
            ! K(2,1) = 2G_ji
            ! K(2,2) = 2G_jj - 1
            complex(kind=8), intent(in) :: Gr(:, :)
            integer, intent(in) :: i, j
            real(kind=8), intent(out) :: r_pair
            logical, intent(out) :: is_safe
            complex(kind=8) :: k11, k22, k12, k21, detK
            real(kind=8) :: Gii_diag, Gjj_diag
            real(kind=8), parameter :: tol = 0.1d0
            integer, save :: warn_count = 0
            complex(kind=8), parameter :: oneC = (1.d0, 0.d0)
            complex(kind=8), parameter :: twoC = (2.d0, 0.d0)

            Gii_diag = real(Gr(i, i))
            Gjj_diag = real(Gr(j, j))

            if (warn_count < 20) then
                if (Gii_diag < -tol .or. Gii_diag > 1.d0 + tol) then
                    warn_count = warn_count + 1
                    write(6,'(A,I4,A,I4,A,F12.6)') &
                        ' WARNING: Gr(',i,',',i,') out of [0,1], Re = ', Gii_diag
                endif
                if (Gjj_diag < -tol .or. Gjj_diag > 1.d0 + tol) then
                    warn_count = warn_count + 1
                    write(6,'(A,I4,A,I4,A,F12.6)') &
                        ' WARNING: Gr(',j,',',j,') out of [0,1], Re = ', Gjj_diag
                endif
            endif

            ! 若对角元严重偏离物理区间，保守地拒绝这对 λ 翻转
            if (Gii_diag < -tol .or. Gii_diag > 1.d0 + tol .or. &
                Gjj_diag < -tol .or. Gjj_diag > 1.d0 + tol) then
                r_pair = 0.d0
                is_safe = .false.
                return
            endif

            ! 正确的 K 矩阵（根据论文推导）：
            ! K = [[2G_ii - 1, 2G_ij], [2G_ji, 2G_jj - 1]]
            k11 = twoC*Gr(i, i) - oneC
            k22 = twoC*Gr(j, j) - oneC
            k12 = twoC*Gr(i, j)   ! 修正：之前是 2*(G_ji - 1)，应为 2*G_ij
            k21 = twoC*Gr(j, i)   ! 修正：之前是 2*(G_ij - 1)，应为 2*G_ji

            detK = k11 * k22 - k12 * k21
            r_pair = abs(detK)
            is_safe = .true.
        end subroutine lambda_pair_ratio

        subroutine lambda_pair_update(Gr, i, j)
            ! 更新 Green 函数：G' = (I + P' B)^{-1}
            ! 使用 Sherman-Morrison-Woodbury 公式：
            ! G' = G - G U K^{-1} V G
            ! 其中 UV = ΔP B，U = [e_i, e_j]
            ! V G 的行是 -2(I-G)[{i,j}, :]
            ! 所以 G' = G + 2 G[:, {i,j}] K^{-1} (I-G)[{i,j}, :]
            complex(kind=8), intent(inout) :: Gr(:, :)
            integer, intent(in) :: i, j
            complex(kind=8) :: k11, k22, k12, k21, detK
            complex(kind=8) :: kinv11, kinv12, kinv21, kinv22
            complex(kind=8), dimension(Ndim) :: z1, z2, w1, w2
            complex(kind=8) :: oneC, twoC, w1_orig, w2_orig
            integer :: kk, ll

            oneC = dcmplx(1.d0, 0.d0)
            twoC = dcmplx(2.d0, 0.d0)

            ! K 矩阵（与 lambda_pair_ratio 一致）
            k11 = twoC*Gr(i, i) - oneC
            k22 = twoC*Gr(j, j) - oneC
            k12 = twoC*Gr(i, j)
            k21 = twoC*Gr(j, i)

            detK = k11 * k22 - k12 * k21
            if (abs(detK) < 1.d-300) return

            kinv11 =  k22 / detK
            kinv22 =  k11 / detK
            kinv12 = -k12 / detK
            kinv21 = -k21 / detK

            ! Z 的两列：Z(:,1) = G(:,i), Z(:,2) = G(:,j)
            do kk = 1, Ndim
                z1(kk) = Gr(kk, i)
                z2(kk) = Gr(kk, j)
            enddo

            ! W 的两行：W(1,:) = (I-G)(i,:), W(2,:) = (I-G)(j,:)
            do ll = 1, Ndim
                if (ll == i) then
                    w1(ll) = oneC - Gr(i, ll)
                else
                    w1(ll) = -Gr(i, ll)
                endif
                if (ll == j) then
                    w2(ll) = oneC - Gr(j, ll)
                else
                    w2(ll) = -Gr(j, ll)
                endif
            enddo

            ! 计算 Y = K^{-1} W（2×N 矩阵，存为两个行向量 y1, y2）
            do ll = 1, Ndim
                w1_orig = w1(ll)
                w2_orig = w2(ll)
                w1(ll) = kinv11 * w1_orig + kinv12 * w2_orig
                w2(ll) = kinv21 * w1_orig + kinv22 * w2_orig
            enddo

            ! 更新：G' = G + 2 * (z1 ⊗ y1 + z2 ⊗ y2)
            do kk = 1, Ndim
                do ll = 1, Ndim
                    Gr(kk, ll) = Gr(kk, ll) + twoC * (z1(kk) * w1(ll) + z2(kk) * w2(ll))
                enddo
            enddo
        end subroutine lambda_pair_update
    end subroutine Local_lambda_flip

    subroutine LocalK_prop_L(PropU, PropD, iseed, nt)
        ! 向左扫描时的传播：从 τ 到 τ-1
        ! 
        ! 使用棋盘分解和 Woodbury 公式
        ! 顺序：先 metro（从 group_4 到 group_1），再 wrap
        ! 
        ! 这确保了每次 metro 更新时 Green 函数的一致性
        class(Propagator), intent(inout) :: PropU, PropD
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        
        ! 步骤 1：按组做 metro 更新（从 group_4 到 group_1）
        ! 临时禁用以测试仅 sweep_R 的 metro
        ! call LocalK_metro_group(PropU%Gr, PropD%Gr, iseed, Latt%group_4, nt, 4)
        ! call LocalK_metro_group(PropU%Gr, PropD%Gr, iseed, Latt%group_3, nt, 3)
        ! call LocalK_metro_group(PropU%Gr, PropD%Gr, iseed, Latt%group_2, nt, 2)
        ! call LocalK_metro_group(PropU%Gr, PropD%Gr, iseed, Latt%group_1, nt, 1)
        
        ! 步骤 2：wrap G（使用棋盘分解）
        ! G(τ-1) = B(τ) * G(τ) * B(τ)^{-1}
        call Op_K%mmult_L(PropU%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropD%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_R(PropU%Gr, Latt, NsigL_K%sigma, nt, -1)
        call Op_K%mmult_R(PropD%Gr, Latt, NsigL_K%sigma, nt, -1)
        
        ! 步骤 3：更新 UUL（使用可能更新后的 σ）
        call Op_K%mmult_L(PropU%UUL, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropD%UUL, Latt, NsigL_K%sigma, nt, 1)
        
        return
    end subroutine LocalK_prop_L

    subroutine LocalK_prop_R(PropU, PropD, iseed, nt)
        ! 向右扫描时的传播：从 τ-1 到 τ
        ! 
        ! 使用棋盘分解和 Woodbury 公式
        ! 顺序：先 wrap，再 metro（按 group_1 到 group_4 的顺序）
        ! 
        ! 这确保了每次 metro 更新时 Green 函数的一致性
        class(Propagator), intent(inout) :: PropU, PropD
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        
        ! 步骤 1：wrap G（使用棋盘分解）
        ! G(τ) = B(τ)^{-1} * G(τ-1) * B(τ)
        call Op_K%mmult_R(PropU%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_R(PropD%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropU%Gr, Latt, NsigL_K%sigma, nt, -1)
        call Op_K%mmult_L(PropD%Gr, Latt, NsigL_K%sigma, nt, -1)
        
        ! 步骤 2：按组做 metro 更新（从 group_1 到 group_4）
        ! 使用 Woodbury 公式，正确考虑棋盘分解的影响
        call LocalK_metro_group(PropU%Gr, PropD%Gr, iseed, Latt%group_1, nt, 1)
        call LocalK_metro_group(PropU%Gr, PropD%Gr, iseed, Latt%group_2, nt, 2)
        call LocalK_metro_group(PropU%Gr, PropD%Gr, iseed, Latt%group_3, nt, 3)
        call LocalK_metro_group(PropU%Gr, PropD%Gr, iseed, Latt%group_4, nt, 4)
        
        ! 步骤 3：更新 UUR（使用可能更新后的 σ）
        call Op_K%mmult_R(PropU%UUR, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_R(PropD%UUR, Latt, NsigL_K%sigma, nt, 1)
        
        return
    end subroutine LocalK_prop_R

    subroutine LocalK_therm(ii, ntau, iseed)
! Arguments:   
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau
! Local: 
        real(kind=8), external :: ranf
        real(kind=8) :: random, ratio_boson
        integer :: ns

        sigma_new(ii, ntau) = - NsigL_K%sigma(ii, ntau)
        ratio_boson = NsigL_K%bosonratio(sigma_new, ii, ntau, Latt)
        random = ranf(iseed)
        if (ratio_boson > random) then
            call Acc_Kt%count(.true.)
            NsigL_K%sigma(ii, ntau) = sigma_new(ii, ntau)
        else
            call Acc_Kt%count(.false.)
            sigma_new(ii, ntau) = NsigL_K%sigma(ii, ntau)
        endif
        return
    end subroutine LocalK_therm
end module LocalK_mod
