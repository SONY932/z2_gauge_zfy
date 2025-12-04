module GlobalK_mod
    use Multiply_mod
    use DQMC_Model_mod
    use ProcessMatrix
    use MyLattice
    use calc_basic
    use Fields_mod, only: gauss_boson_ratio_lambda
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
    ! 单格点 λ 更新（正确实现高斯约束）
    !=========================================================================
    !
    ! 根据 PNAS SI，翻转单个 λ_i 的接受率包含两部分：
    ! 1. 费米子行列式比：det[1 + P'B] / det[1 + PB]
    ! 2. 玻色权重比：(-1)^{n_b + n_Q}
    !
    ! 重要说明：
    ! - λ 不需要成对翻转，单格点翻转是合法的
    ! - 扇区由 Q 参数决定，不是由 ∏λ 决定
    ! - 在 global update 中进行 λ 更新可以避免 local update 中的数值不稳定问题
    !   因为此时 Green 函数刚刚被稳定化
    !=========================================================================
    
    subroutine Global_lambda_update(GrU, GrD, iseed, n_accept, n_total)
! ===================================================================
! λ 场更新
!
! 重要说明：
! - 传入的 GrU, GrD 是 G_0 = (1 + B)^{-1}，不是 G_λ
! - λ 更新需要 G_λ = (1 + P[λ]B)^{-1}
! - 但 λ 翻转不改变 B，所以 G_0 不需要更新
! - 我们临时计算 G_λ 进行接受率计算
!
! 方案：
! 1. 一次性从 G_0 计算 G_λ
! 2. 遍历所有格点，使用 G_λ 计算接受率
! 3. 如果接受，用 rank-1 公式更新 G_λ，翻转 λ
! 4. 最后 G_0 不变（因为 B 没变）
! ===================================================================
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: GrU, GrD
        integer, intent(inout) :: iseed
        integer, intent(out) :: n_accept, n_total
        
        real(kind=8), external :: ranf
        real(kind=8) :: ratio_fermion_U, ratio_fermion_D, ratio_boson, ratio_total, random
        real(kind=8) :: Gii_U, Gii_D
        integer :: ii
        real(kind=8), parameter :: tol = 0.1d0
        logical :: safe
        
        ! 临时存储 G_λ
        complex(kind=8), dimension(Ndim, Ndim) :: Glambda_U, Glambda_D
        
        n_accept = 0
        n_total = 0
        
        ! 从 G_0 计算 G_λ
        Glambda_U = GrU
        Glambda_D = GrD
        call compute_G_lambda_from_G0(Glambda_U)
        call compute_G_lambda_from_G0(Glambda_D)
        
        ! 遍历所有格点，对每个格点尝试翻转 λ
        do ii = 1, Lq
            n_total = n_total + 1
            
            Gii_U = real(Glambda_U(ii, ii))
            Gii_D = real(Glambda_D(ii, ii))
            
            ! 检查 Green 函数对角元是否在物理区间
            safe = .true.
            if (Gii_U < -tol .or. Gii_U > 1.d0 + tol) safe = .false.
            if (Gii_D < -tol .or. Gii_D > 1.d0 + tol) safe = .false.
            
            if (.not. safe) then
                cycle
            endif
            
            ! 费米子行列式比：|2G_λ(i, i) - 1|
            ratio_fermion_U = abs(2.d0 * Gii_U - 1.d0)
            ratio_fermion_D = abs(2.d0 * Gii_D - 1.d0)
            
            ! 玻色权重比
            ratio_boson = gauss_boson_ratio_lambda(ii, NsigL_K%sigma, NsigL_K%lambda, Latt)
            
            ! 总接受率
            ratio_total = ratio_fermion_U * ratio_fermion_D * ratio_boson
            
            random = ranf(iseed)
            if (abs(ratio_total) > random) then
                n_accept = n_accept + 1
                
                ! 更新 G_λ（用于后续的 λ 翻转）
                call lambda_update_Green(Glambda_U, ii)
                call lambda_update_Green(Glambda_D, ii)
                
                ! 翻转 λ
                NsigL_K%lambda(ii) = -NsigL_K%lambda(ii)
            endif
        enddo
        
        ! 注意：G_0 (GrU, GrD) 不需要更新，因为 λ 翻转不改变 B
        return
    contains
        subroutine compute_G_lambda_from_G0(Gr)
            ! 从 G_0 = (1 + B)^{-1} 计算 G_λ = (1 + P[λ]B)^{-1}
            ! 
            ! M_0 = G_0^{-1} = 1 + B
            ! M_λ = 1 + P[λ]B = P[λ]M_0 + (I - P[λ])
            ! G_λ = M_λ^{-1}
            complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gr
            complex(kind=8), allocatable :: M0(:,:), Mlambda(:,:)
            real(kind=8) :: lam_i
            integer :: i, j, info
            integer, allocatable :: ipiv(:)
            complex(kind=8), allocatable :: work(:)
            integer :: lwork
            logical :: all_one
            
            ! 如果 λ 全为 1，G_λ = G_0，直接返回
            all_one = .true.
            do i = 1, Ndim
                if (abs(NsigL_K%lambda(i) - 1.d0) > 1.d-12) then
                    all_one = .false.
                    exit
                endif
            enddo
            if (all_one) return
            
            allocate(M0(Ndim, Ndim), Mlambda(Ndim, Ndim))
            allocate(ipiv(Ndim))
            lwork = Ndim * Ndim
            allocate(work(lwork))
            
            ! M0 = Gr^{-1}
            M0 = Gr
            call ZGETRF(Ndim, Ndim, M0, Ndim, ipiv, info)
            if (info /= 0) then
                deallocate(M0, Mlambda, ipiv, work)
                return
            endif
            call ZGETRI(Ndim, M0, Ndim, ipiv, work, lwork, info)
            if (info /= 0) then
                deallocate(M0, Mlambda, ipiv, work)
                return
            endif
            
            ! M_λ = P[λ]*M0 + (I - P[λ])
            do i = 1, Ndim
                lam_i = NsigL_K%lambda(i)
                do j = 1, Ndim
                    if (i == j) then
                        Mlambda(i, j) = dcmplx(lam_i, 0.d0) * M0(i, j) + dcmplx(1.d0 - lam_i, 0.d0)
                    else
                        Mlambda(i, j) = dcmplx(lam_i, 0.d0) * M0(i, j)
                    endif
                enddo
            enddo
            
            ! G_λ = M_λ^{-1}
            call ZGETRF(Ndim, Ndim, Mlambda, Ndim, ipiv, info)
            if (info /= 0) then
                deallocate(M0, Mlambda, ipiv, work)
                return
            endif
            call ZGETRI(Ndim, Mlambda, Ndim, ipiv, work, lwork, info)
            if (info == 0) then
                Gr = Mlambda
            endif
            
            deallocate(M0, Mlambda, ipiv, work)
        end subroutine compute_G_lambda_from_G0
        
        subroutine lambda_update_Green(Gr, i)
            ! 更新 G_λ（rank-1 Sherman-Morrison）
            ! G' = G + 2 * G[:, i] * (I - G)[i, :] / (2G(i,i) - 1)
            complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gr
            integer, intent(in) :: i
            complex(kind=8) :: denom
            complex(kind=8), dimension(Ndim) :: col_i, row_i
            integer :: kk, ll
            
            denom = 2.d0 * Gr(i, i) - dcmplx(1.d0, 0.d0)
            
            if (abs(denom) < 1.d-300) return
            
            do kk = 1, Ndim
                col_i(kk) = Gr(kk, i)
                if (kk == i) then
                    row_i(kk) = dcmplx(1.d0, 0.d0) - Gr(i, kk)
                else
                    row_i(kk) = -Gr(i, kk)
                endif
            enddo
            
            do kk = 1, Ndim
                do ll = 1, Ndim
                    Gr(kk, ll) = Gr(kk, ll) + 2.d0 * col_i(kk) * row_i(ll) / denom
                enddo
            enddo
        end subroutine lambda_update_Green
    end subroutine Global_lambda_update
    
end module GlobalK_mod