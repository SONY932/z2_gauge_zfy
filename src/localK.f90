module LocalK_mod
    use Multiply_mod
    use DQMC_Model_mod
    use ProcessMatrix
    use MyLattice
    use calc_basic
    implicit none

    public
    private :: LocalK_metro, sigma_new, Local_lambda_flip, lambda_new

    type(AccCounter) :: Acc_Kl, Acc_Kt, Acc_lambda
    real(kind = 8), dimension(:, :), allocatable :: sigma_new
    real(kind = 8), dimension(:), allocatable :: lambda_new

contains
    subroutine LocalK_init()
        call Acc_Kl%init()
        call Acc_Kt%init()
        call Acc_lambda%init()
        allocate(sigma_new(2*Lq, Ltrot))
        allocate(lambda_new(Lq))
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
        real(kind=8) :: lambda_P1, lambda_P2  ! λ 投影因子

! 局部翻转规范场
        sigma_new(ii, ntau) = - NsigL_K%sigma(ii, ntau)
        S_new = sigma_new(ii, ntau)
        S_old = NsigL_K%sigma(ii, ntau)
        P(1) = Latt%bond_list(ii, 1)
        P(2) = Latt%bond_list(ii, 2)
        
        ! 获取 λ 投影因子：λ_{P(1)} 和 λ_{P(2)}
        ! 费米子矩阵 M = I + P_λ B_tot，其中 P_λ = diag(λ_r)
        lambda_P1 = NsigL_K%lambda(P(1))
        lambda_P2 = NsigL_K%lambda(P(2))

        call Op_K%get_delta(S_old, S_new)
        ProdU = dcmplx(0.d0, 0.d0)
        ProdD = dcmplx(0.d0, 0.d0)
        ! 标准 Sherman-Morrison 公式：K = I + Δ (I - G)_{P,P}
        ! 注意：当 λ 投影被启用时，G 已经是 (I + P_λ B)^{-1} 的形式
        ! 但 Sherman-Morrison 更新仍然使用 (I - G) 的因子，因为更新来自
        ! B' = B(I + Δ_local)，即 ΔB = B Δ_local，而不是直接的 ΔB
        ! 在这种情况下，K 矩阵不需要额外的 λ 因子
        do nr = 1, 2
            do nl = 1, 2
                GrU_local(nl, nr) = ZKRON(nl, nr) - GrU(P(nl), P(nr))
                GrD_local(nl, nr) = ZKRON(nl, nr) - GrD(P(nl), P(nr))
            enddo
        enddo

        matU_tmp = matmul(Op_K%Delta, GrU_local) ! 2x2 * 2x2
        matD_tmp = matmul(Op_K%Delta, GrD_local) ! 2x2 * 2x2
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

    subroutine Local_lambda_flip(GrU, GrD, iseed, ii, jj)
        use MyMats
! Arguments:
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: GrU, GrD
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, jj
! Local:
        real(kind=8), external :: ranf
        real(kind=8) :: rU_pair, rD_pair, prob, random
        logical :: safeU, safeD

        if (ii == jj) return

! 采用 rank-2 Woodbury，一次性对 (ii,jj) 成对翻转 λ，避免在主 Gr 上反复 preview/revert。
! 先在当前 Green 上计算两个自旋的成对行列式比，再根据接受与否，对主 Gr 做一次 rank-2 更新。
        call lambda_pair_ratio(GrU, ii, jj, rU_pair, safeU)
        call lambda_pair_ratio(GrD, ii, jj, rD_pair, safeD)

        if (safeU .and. safeD .and. rU_pair > 0.d0 .and. rD_pair > 0.d0) then
            prob = exp(min(0.d0, log(rU_pair) + log(rD_pair)))
        else
            prob = 0.d0
        endif
        prob = max(0.d0, min(1.d0, prob))

        random = ranf(iseed)
        if (prob > random) then
            call Acc_lambda%count(.true.)
            call lambda_pair_update(GrU, ii, jj)
            call lambda_pair_update(GrD, ii, jj)
            NsigL_K%lambda(ii) = -NsigL_K%lambda(ii)
            NsigL_K%lambda(jj) = -NsigL_K%lambda(jj)
            lambda_new(ii) = NsigL_K%lambda(ii)
            lambda_new(jj) = NsigL_K%lambda(jj)
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
        class(Propagator), intent(inout) :: PropU, PropD
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        integer :: ii, jj
        do ii = 2*Lq, 1, -1
            call LocalK_metro(PropU%Gr, PropD%Gr, iseed, ii, nt)
        enddo
        ! λ 场更新：在最后一个时间片 (nt == Ltrot) 进行
        ! σ 更新的 Sherman-Morrison 公式已修改为支持 G = (I + P_λ B)^{-1}
        ! λ 成对翻转保持 Q = Π_r λ_r 不变
        if (nt == Ltrot) then
            do ii = 1, Lq - 1
                do jj = ii + 1, Lq
                    call Local_lambda_flip(PropU%Gr, PropD%Gr, iseed, ii, jj)
                enddo
            enddo
        endif
        call Op_K%mmult_L(PropU%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropD%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_R(PropU%Gr, Latt, NsigL_K%sigma, nt, -1)
        call Op_K%mmult_R(PropD%Gr, Latt, NsigL_K%sigma, nt, -1)
        call Op_K%mmult_L(PropU%UUL, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropD%UUL, Latt, NsigL_K%sigma, nt, 1)
        return
    end subroutine LocalK_prop_L

    subroutine LocalK_prop_R(PropU, PropD, iseed, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        integer :: ii, jj
        call Op_K%mmult_R(PropU%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_R(PropD%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropU%Gr, Latt, NsigL_K%sigma, nt, -1)
        call Op_K%mmult_L(PropD%Gr, Latt, NsigL_K%sigma, nt, -1)
        do ii = 1, 2*Lq
            call LocalK_metro(PropU%Gr, PropD%Gr, iseed, ii, nt)
        enddo
        ! λ 场更新：在最后一个时间片 (nt == Ltrot) 进行
        if (nt == Ltrot) then
            do ii = 1, Lq - 1
                do jj = ii + 1, Lq
                    call Local_lambda_flip(PropU%Gr, PropD%Gr, iseed, ii, jj)
                enddo
            enddo
        endif
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
