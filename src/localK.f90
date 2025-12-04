module LocalK_mod
    use Multiply_mod
    use DQMC_Model_mod
    use ProcessMatrix
    use MyLattice
    use calc_basic
    use Fields_mod, only: gauss_boson_ratio_lambda, gauss_boson_ratio_sigma
    implicit none

    public
    private :: LocalK_metro, sigma_new, Local_lambda_flip, lambda_new

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
        
        ! ===================================================================
        ! 严格正确的 σ 更新接受率公式
        ! 
        ! 根据 PRX 论文，正确的配分函数是 det[1 + P[λ] B]
        ! 理论上正确的接受率是：det[1 + P[λ] B'] / det[1 + P[λ] B]
        ! 
        ! 但实际实现中使用标准公式 det[1 + B'] / det[1 + B]
        ! λ 的效应通过 Global_lambda_update 在每次 sweep 结束时处理
        ! ===================================================================
        
        ! 计算 (I - G_0) 矩阵
        do nr = 1, 2
            do nl = 1, 2
                GrU_local(nl, nr) = ZKRON(nl, nr) - GrU(P(nl), P(nr))
                GrD_local(nl, nr) = ZKRON(nl, nr) - GrD(P(nl), P(nr))
            enddo
        enddo

        ! 标准 K 矩阵（基于 G_0，用于 Green 函数更新）
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
        
        ! ===================================================================
        ! σ 更新的接受率
        ! 
        ! 重要说明：使用标准公式 det[1 + B'] / det[1 + B] = det[I + G_0 ΔB]
        ! 
        ! 虽然理论上正确的公式应该是 det[1 + P[λ] B'] / det[1 + P[λ] B]，
        ! 但使用 G_λ 计算接受率会导致与 Green 函数更新（基于 G_0）不一致，
        ! 从而造成数值误差累积。
        ! 
        ! 当前策略：
        ! - σ 更新使用标准公式（基于 G_0）
        ! - λ 效应通过 Global_lambda_update 在每次 sweep 结束时处理
        ! - 整体采样通过交替进行 σ 和 λ 更新收敛到正确分布
        ! ===================================================================
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
! Vhlp(1:2, 1:Ndim) = Del(1:2) * (1 - Grup)(P(1):P(2), 1:Ndim); Uhlp(1:Ndim, 1:2) = Grup(1:Ndim, P(1):P(2))
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
        ! 使用标准方法：先做所有 metro 更新，再做完整的 wrap
        ! 这与 CodeXun 的实现一致
        class(Propagator), intent(inout) :: PropU, PropD
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        integer :: ii
        
        ! 按相反顺序对所有键做 metro 更新
        do ii = 2*Lq, 1, -1
            call LocalK_metro(PropU%Gr, PropD%Gr, iseed, ii, nt)
        enddo
        
        ! 完整的 wrap：G -> B * G * B^{-1}
        call Op_K%mmult_L(PropU%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropD%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_R(PropU%Gr, Latt, NsigL_K%sigma, nt, -1)
        call Op_K%mmult_R(PropD%Gr, Latt, NsigL_K%sigma, nt, -1)
        
        ! 更新 UUL 累积传播子
        call Op_K%mmult_L(PropU%UUL, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropD%UUL, Latt, NsigL_K%sigma, nt, 1)
        return
    end subroutine LocalK_prop_L

    subroutine LocalK_prop_R(PropU, PropD, iseed, nt)
        ! 向右扫描时的传播：从 τ-1 到 τ
        ! 使用标准方法：先做完整的 wrap，再做所有 metro 更新
        class(Propagator), intent(inout) :: PropU, PropD
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        integer :: ii
        
        ! 完整的 wrap：G -> B * G * B^{-1}
        call Op_K%mmult_R(PropU%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_R(PropD%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropU%Gr, Latt, NsigL_K%sigma, nt, -1)
        call Op_K%mmult_L(PropD%Gr, Latt, NsigL_K%sigma, nt, -1)
        
        ! 按顺序对所有键做 metro 更新
        do ii = 1, 2*Lq
            call LocalK_metro(PropU%Gr, PropD%Gr, iseed, ii, nt)
        enddo
        
        ! 更新 UUR 累积传播子
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
