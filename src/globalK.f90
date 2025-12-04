module GlobalK_mod
    use Multiply_mod
    use DQMC_Model_mod
    use ProcessMatrix
    use MyLattice
    use calc_basic
    implicit none

    real(kind=8), parameter :: RATIO_EPS = 1.d-300

    public
    private :: ratioK_fermion

contains
    real(kind=8) function ratioK_fermion(GrU, GrD, sigma_new, ii, ntau) result(ratio_fermion)
! Arguments:
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: GrU, GrD
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_new
        integer, intent(in) :: ii, ntau
! Local:
        complex(kind=8), dimension(2, 2) :: ProdU, ProdD, ProdinvU, ProdinvD
        complex(kind=8), dimension(2, 2) :: GrU_local, GrD_local, matU_tmp, matD_tmp
        complex(kind=8), dimension(Ndim, 2) :: UhlpU, UhlpD
        complex(kind=8), dimension(2, Ndim) :: VhlpU, VhlpD
        complex(kind=8), dimension(Ndim, 2) :: temp
        complex(kind=8), dimension(Ndim, Ndim) :: Diff
        complex(kind=8) :: ProddetU, ProddetD
        real(kind=8) :: S_old, S_new
        integer :: P(2), nl, nr, j

        S_old = NsigL_K%sigma(ii, ntau)
        S_new = sigma_new(ii, ntau)
        if (S_new == S_old) then
            ratio_fermion = 1.d0
            return
        endif

        P(1) = Latt%bond_list(ii, 1)
        P(2) = Latt%bond_list(ii, 2)
        call Op_K%get_delta(S_old, S_new)

        ProdU = dcmplx(0.d0, 0.d0)
        ProdD = dcmplx(0.d0, 0.d0)
        do nr = 1, 2
            do nl = 1, 2
                GrU_local(nl, nr) = ZKRON(nl, nr) - GrU(P(nl), P(nr))
                GrD_local(nl, nr) = ZKRON(nl, nr) - GrD(P(nl), P(nr))
            enddo
        enddo

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
        ratio_fermion = max(abs(ProddetU * ProddetD), RATIO_EPS)

! Update spin-up Green's function
        ProdinvU(1, 1) = ProdU(2, 2)
        ProdinvU(2, 2) = ProdU(1, 1)
        ProdinvU(1, 2) = -ProdU(1, 2)
        ProdinvU(2, 1) = -ProdU(2, 1)
        ProdinvU = ProdinvU / ProddetU
        UhlpU = dcmplx(0.d0, 0.d0)
        VhlpU = dcmplx(0.d0, 0.d0)
        temp = dcmplx(0.d0, 0.d0)
        Diff = dcmplx(0.d0, 0.d0)
        do nl = 1, 2
            do j = 1, Ndim
                UhlpU(j, nl) = GrU(j, P(nl))
                VhlpU(nl, j) = - Op_K%Delta(nl, 1) * GrU(P(1), j) - Op_K%Delta(nl, 2) * GrU(P(2), j)
            enddo
            VhlpU(nl, P(1)) = VhlpU(nl, P(1)) + Op_K%Delta(nl, 1)
            VhlpU(nl, P(2)) = VhlpU(nl, P(2)) + Op_K%Delta(nl, 2)
        enddo
        temp = matmul(UhlpU, ProdinvU)
        Diff = matmul(temp, VhlpU)
        GrU = GrU - Diff

! Update spin-down Green's function
        ProdinvD(1, 1) = ProdD(2, 2)
        ProdinvD(2, 2) = ProdD(1, 1)
        ProdinvD(1, 2) = -ProdD(1, 2)
        ProdinvD(2, 1) = -ProdD(2, 1)
        ProdinvD = ProdinvD / ProddetD
        UhlpD = dcmplx(0.d0, 0.d0)
        VhlpD = dcmplx(0.d0, 0.d0)
        temp = dcmplx(0.d0, 0.d0)
        Diff = dcmplx(0.d0, 0.d0)
        do nl = 1, 2
            do j = 1, Ndim
                UhlpD(j, nl) = GrD(j, P(nl))
                VhlpD(nl, j) = - Op_K%Delta(nl, 1) * GrD(P(1), j) - Op_K%Delta(nl, 2) * GrD(P(2), j)
            enddo
            VhlpD(nl, P(1)) = VhlpD(nl, P(1)) + Op_K%Delta(nl, 1)
            VhlpD(nl, P(2)) = VhlpD(nl, P(2)) + Op_K%Delta(nl, 2)
        enddo
        temp = matmul(UhlpD, ProdinvD)
        Diff = matmul(temp, VhlpD)
        GrD = GrD - Diff

        return
    end function ratioK_fermion

    subroutine GlobalK_prop_L(PropU, PropD, log_ratio_fermion, sigma_new, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        real(kind=8), intent(inout) :: log_ratio_fermion
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_new
        integer, intent(in) :: nt
        integer :: ii
        real(kind=8) :: ratio

        do ii = 2*Lq, 1, -1
            ratio = ratioK_fermion(PropU%Gr, PropD%Gr, sigma_new, ii, nt)
            log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
        enddo
        call Op_K%mmult_L(PropU%Gr, Latt, sigma_new, nt, 1)
        call Op_K%mmult_L(PropD%Gr, Latt, sigma_new, nt, 1)
        call Op_K%mmult_R(PropU%Gr, Latt, sigma_new, nt, -1)
        call Op_K%mmult_R(PropD%Gr, Latt, sigma_new, nt, -1)
        call Op_K%mmult_L(PropU%UUL, Latt, sigma_new, nt, 1)
        call Op_K%mmult_L(PropD%UUL, Latt, sigma_new, nt, 1)
        return
    end subroutine GlobalK_prop_L

    subroutine GlobalK_prop_R(PropU, PropD, log_ratio_fermion, sigma_new, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        real(kind=8), intent(inout) :: log_ratio_fermion
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_new
        integer, intent(in) :: nt
        integer :: ii
        real(kind=8) :: ratio

        call Op_K%mmult_R(PropU%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_R(PropD%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropU%Gr, Latt, NsigL_K%sigma, nt, -1)
        call Op_K%mmult_L(PropD%Gr, Latt, NsigL_K%sigma, nt, -1)
        do ii = 1, 2*Lq
            ratio = ratioK_fermion(PropU%Gr, PropD%Gr, sigma_new, ii, nt)
            log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
        enddo
        call Op_K%mmult_R(PropU%UUR, Latt, sigma_new, nt, 1)
        call Op_K%mmult_R(PropD%UUR, Latt, sigma_new, nt, 1)
        return
    end subroutine GlobalK_prop_R

    !=========================================================================
    ! 全局 λ 更新函数
    ! 根据 PRX Appendix A，λ 更新需要在全局进行
    ! 因为 P[λ] 和传播子 B 不对易，不能使用标准传播公式
    !=========================================================================
    
    subroutine Global_lambda_update(GrU, GrD, iseed, n_accept, n_total)
        ! 全局 λ 更新：成对翻转 λ_i 和 λ_j
        ! 行列式比：det[1 + P[λ'] B] / det[1 + P[λ] B]
        ! 使用 Woodbury 公式计算
        complex(kind=8), dimension(Ndim, Ndim), intent(in) :: GrU, GrD
        integer, intent(inout) :: iseed
        integer, intent(inout) :: n_accept, n_total
        
        real(kind=8), external :: ranf
        integer :: ii, jj, attempt
        integer :: n_attempts
        real(kind=8) :: rU, rD, prob, random
        logical :: safeU, safeD
        
        ! 每次全局更新尝试多次成对翻转
        n_attempts = Lq / 2  ! 尝试 Lq/2 次成对翻转
        
        do attempt = 1, n_attempts
            ! 随机选择两个不同的格点进行成对翻转
            ii = int(ranf(iseed) * Lq) + 1
            if (ii > Lq) ii = Lq
            jj = int(ranf(iseed) * (Lq - 1)) + 1
            if (jj >= ii) jj = jj + 1
            if (jj > Lq) jj = Lq
            
            ! 计算行列式比
            ! 使用与 lambda_pair_ratio 相同的公式
            ! K = | 2G_ii - 1    2G_ij     |
            !     | 2G_ji        2G_jj - 1 |
            ! det K = (2G_ii - 1)(2G_jj - 1) - 4 G_ij G_ji
            call compute_lambda_ratio(GrU, ii, jj, rU, safeU)
            call compute_lambda_ratio(GrD, ii, jj, rD, safeD)
            
            if (safeU .and. safeD .and. rU > 0.d0 .and. rD > 0.d0) then
                prob = min(1.d0, rU * rD)
            else
                prob = 0.d0
            endif
            
            random = ranf(iseed)
            n_total = n_total + 1
            if (prob > random) then
                ! 接受翻转
                n_accept = n_accept + 1
                NsigL_K%lambda(ii) = -NsigL_K%lambda(ii)
                NsigL_K%lambda(jj) = -NsigL_K%lambda(jj)
            endif
        enddo
        
        return
    contains
        subroutine compute_lambda_ratio(Gr, i, j, r, is_safe)
            ! 计算成对翻转 λ_i 和 λ_j 的行列式比
            ! 公式基于 G_0 = (1 + B)^{-1}
            ! 当 λ_i, λ_j 当前都是 +1 或都是 -1 时，翻转后 S 的大小变化
            complex(kind=8), dimension(:, :), intent(in) :: Gr
            integer, intent(in) :: i, j
            real(kind=8), intent(out) :: r
            logical, intent(out) :: is_safe
            
            complex(kind=8) :: k11, k12, k21, k22, detK
            complex(kind=8), parameter :: oneC = dcmplx(1.d0, 0.d0)
            complex(kind=8), parameter :: twoC = dcmplx(2.d0, 0.d0)
            real(kind=8), parameter :: SAFE_THRESHOLD = 1.d-10
            
            ! K 矩阵：当 λ 翻转时，det[1 + P[λ'] B] / det[1 + P[λ] B]
            ! 对于 G_0 = (1 + B)^{-1}，(I - G_0) = B G_0
            ! K_ij = δ_ij - 2 λ_i (B G_0)_ij = δ_ij - 2 λ_i (I - G_0)_ij
            ! 当 λ_i = λ_j = 1（翻转到 -1）：
            !   K_ii = 1 - 2(1 - G_ii) = 2G_ii - 1
            !   K_ij = -2(-G_ij) = 2G_ij
            ! 当 λ_i = λ_j = -1（翻转到 1）：
            !   K_ii = 1 + 2(1 - G_ii) = 3 - 2G_ii
            !   K_ij = 2(-G_ij) = -2G_ij
            
            if (NsigL_K%lambda(i) > 0.d0 .and. NsigL_K%lambda(j) > 0.d0) then
                ! 当前 λ_i = λ_j = 1，翻转到 -1
                k11 = twoC * Gr(i, i) - oneC
                k22 = twoC * Gr(j, j) - oneC
                k12 = twoC * Gr(i, j)
                k21 = twoC * Gr(j, i)
            else if (NsigL_K%lambda(i) < 0.d0 .and. NsigL_K%lambda(j) < 0.d0) then
                ! 当前 λ_i = λ_j = -1，翻转到 1
                k11 = dcmplx(3.d0, 0.d0) - twoC * Gr(i, i)
                k22 = dcmplx(3.d0, 0.d0) - twoC * Gr(j, j)
                k12 = -twoC * Gr(i, j)
                k21 = -twoC * Gr(j, i)
            else
                ! 一个是 +1，一个是 -1：不能成对翻转（会改变 Q = Π λ）
                r = 0.d0
                is_safe = .false.
                return
            endif
            
            detK = k11 * k22 - k12 * k21
            r = dble(detK)
            is_safe = abs(detK) > SAFE_THRESHOLD
            
            return
        end subroutine compute_lambda_ratio
    end subroutine Global_lambda_update
    
end module GlobalK_mod