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
    private :: ratioK_fermion_simple, ratioK_fermion_woodbury
    private :: compute_LR_cols_rows_global
    private :: apply_group_to_cols_global, apply_group_to_rows_global, apply_group_inv_to_rows_global

contains
    subroutine ratioK_fermion_woodbury(GrU, GrD, sigma_base, sigma_new, ii, ntau, group_idx, ratio_fermion)
        ! =========================================================================
        ! 使用 Woodbury 公式计算费米子行列式比并更新 Green 函数
        ! 
        ! 与 LocalK_metro_woodbury 使用相同的逻辑
        ! 
        ! 重要说明：传入的 G 是 G_0 = (1 + B)^{-1}
        ! 当所有 λ = 1 时，公式完全正确
        ! =========================================================================
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: GrU, GrD
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_base, sigma_new
        integer, intent(in) :: ii, ntau, group_idx
        real(kind=8), intent(out) :: ratio_fermion
        
        ! Local
        complex(kind=8) :: detR_U, detR_D
        complex(kind=8), dimension(2, 2) :: R_small_U, R_small_D, Rinv_U, Rinv_D
        complex(kind=8) :: L_cols(Ndim, 2)
        complex(kind=8) :: Linv_rows(2, Ndim)
        complex(kind=8) :: R_Binv_rows(2, Ndim)
        complex(kind=8) :: VG_U(2, Ndim), VG_D(2, Ndim)
        complex(kind=8) :: U(Ndim, 2)
        complex(kind=8) :: GU_U(Ndim, 2), GU_D(Ndim, 2)
        complex(kind=8) :: temp(Ndim, 2), Diff(Ndim, Ndim)
        complex(kind=8) :: Delta_c(2, 2)
        complex(kind=8) :: Bk_inv_2x2(2, 2)
        real(kind=8) :: Delta_local(2, 2)
        real(kind=8) :: S_old, S_new_val
        real(kind=8) :: Cv_inv, Sv_inv
        integer :: P(2), j, no, nl, nr, kk

        S_old = sigma_base(ii, ntau)
        S_new_val = sigma_new(ii, ntau)
        if (abs(S_new_val - S_old) < 1.d-12) then
            ratio_fermion = 1.d0
            return
        endif

        P(1) = Latt%bond_list(ii, 1)
        P(2) = Latt%bond_list(ii, 2)
        
        ! 计算 Delta 矩阵
        ! 注意：Woodbury 公式需要 ΔE = E_new - E_old，不是 E_new E_old^{-1} - I
        block
            real(kind=8) :: C_new, S_new_tmp, C_old, S_old_tmp
            C_new = cosh(Op_K%alpha * S_new_val)
            S_new_tmp = sinh(Op_K%alpha * S_new_val)
            C_old = cosh(Op_K%alpha * S_old)
            S_old_tmp = sinh(Op_K%alpha * S_old)
            Delta_local(1, 1) = C_new - C_old
            Delta_local(1, 2) = S_new_tmp - S_old_tmp
            Delta_local(2, 1) = S_new_tmp - S_old_tmp
            Delta_local(2, 2) = C_new - C_old
        end block
        Delta_c = dcmplx(Delta_local, 0.d0)

        ! ========== 步骤 1：计算 L_cols = L_ℓ(:, [i, j]) ==========
        L_cols = dcmplx(0.d0, 0.d0)
        L_cols(P(1), 1) = dcmplx(1.d0, 0.d0)
        L_cols(P(2), 2) = dcmplx(1.d0, 0.d0)
        do kk = group_idx + 1, 4
            call apply_group_to_cols_global(kk, ntau, sigma_base, L_cols)
        enddo

        ! ========== 步骤 2：计算 R_Binv_rows = (R_ℓ × B(τ)^{-1})([i,j], :) ==========
        Linv_rows = dcmplx(0.d0, 0.d0)
        Linv_rows(1, P(1)) = dcmplx(1.d0, 0.d0)
        Linv_rows(2, P(2)) = dcmplx(1.d0, 0.d0)
        do kk = 4, group_idx + 1, -1
            call apply_group_inv_to_rows_global(kk, ntau, sigma_base, Linv_rows)
        enddo

        ! 左乘 B_k^{-1}([i,j], [i,j])
        Cv_inv = cosh(-Op_K%alpha * S_old)
        Sv_inv = sinh(-Op_K%alpha * S_old)
        Bk_inv_2x2(1, 1) = dcmplx(Cv_inv, 0.d0)
        Bk_inv_2x2(1, 2) = dcmplx(Sv_inv, 0.d0)
        Bk_inv_2x2(2, 1) = dcmplx(Sv_inv, 0.d0)
        Bk_inv_2x2(2, 2) = dcmplx(Cv_inv, 0.d0)
        
        R_Binv_rows = dcmplx(0.d0, 0.d0)
        do j = 1, Ndim
            R_Binv_rows(1, j) = Bk_inv_2x2(1, 1) * Linv_rows(1, j) + Bk_inv_2x2(1, 2) * Linv_rows(2, j)
            R_Binv_rows(2, j) = Bk_inv_2x2(2, 1) * Linv_rows(1, j) + Bk_inv_2x2(2, 2) * Linv_rows(2, j)
        enddo

        ! ========== 步骤 3：计算 U = L_cols × Δ ==========
        U = matmul(L_cols, Delta_c)

        ! ========== Spin Up ==========
        VG_U = R_Binv_rows
        do j = 1, Ndim
            do no = 1, 2
                do nl = 1, Ndim
                    VG_U(no, j) = VG_U(no, j) - R_Binv_rows(no, nl) * GrU(nl, j)
                enddo
            enddo
        enddo

        GU_U = dcmplx(0.d0, 0.d0)
        do j = 1, 2
            do no = 1, Ndim
                do nl = 1, Ndim
                    GU_U(no, j) = GU_U(no, j) + GrU(no, nl) * U(nl, j)
                enddo
            enddo
        enddo

        R_small_U = dcmplx(0.d0, 0.d0)
        R_small_U(1, 1) = dcmplx(1.d0, 0.d0)
        R_small_U(2, 2) = dcmplx(1.d0, 0.d0)
        do nr = 1, 2
            do nl = 1, 2
                do no = 1, Ndim
                    R_small_U(nl, nr) = R_small_U(nl, nr) + VG_U(nl, no) * U(no, nr)
                enddo
            enddo
        enddo
        detR_U = R_small_U(1, 1) * R_small_U(2, 2) - R_small_U(1, 2) * R_small_U(2, 1)

        ! ========== Spin Down ==========
        VG_D = R_Binv_rows
        do j = 1, Ndim
            do no = 1, 2
                do nl = 1, Ndim
                    VG_D(no, j) = VG_D(no, j) - R_Binv_rows(no, nl) * GrD(nl, j)
                enddo
            enddo
        enddo

        GU_D = dcmplx(0.d0, 0.d0)
        do j = 1, 2
            do no = 1, Ndim
                do nl = 1, Ndim
                    GU_D(no, j) = GU_D(no, j) + GrD(no, nl) * U(nl, j)
                enddo
            enddo
        enddo

        R_small_D = dcmplx(0.d0, 0.d0)
        R_small_D(1, 1) = dcmplx(1.d0, 0.d0)
        R_small_D(2, 2) = dcmplx(1.d0, 0.d0)
        do nr = 1, 2
            do nl = 1, 2
                do no = 1, Ndim
                    R_small_D(nl, nr) = R_small_D(nl, nr) + VG_D(nl, no) * U(no, nr)
                enddo
            enddo
        enddo
        detR_D = R_small_D(1, 1) * R_small_D(2, 2) - R_small_D(1, 2) * R_small_D(2, 1)
        
        ratio_fermion = max(abs(detR_U * detR_D), RATIO_EPS)

        ! ========== 更新 Green 函数 ==========
        Rinv_U(1, 1) =  R_small_U(2, 2) / detR_U
        Rinv_U(1, 2) = -R_small_U(1, 2) / detR_U
        Rinv_U(2, 1) = -R_small_U(2, 1) / detR_U
        Rinv_U(2, 2) =  R_small_U(1, 1) / detR_U

        Rinv_D(1, 1) =  R_small_D(2, 2) / detR_D
        Rinv_D(1, 2) = -R_small_D(1, 2) / detR_D
        Rinv_D(2, 1) = -R_small_D(2, 1) / detR_D
        Rinv_D(2, 2) =  R_small_D(1, 1) / detR_D

        ! Spin Up: G' = G - GU × Rinv × VG
        temp = matmul(GU_U, Rinv_U)
        Diff = matmul(temp, VG_U)
        GrU = GrU - Diff

        ! Spin Down
        temp = matmul(GU_D, Rinv_D)
        Diff = matmul(temp, VG_D)
        GrD = GrD - Diff

        return
    end subroutine ratioK_fermion_woodbury

    subroutine compute_LR_cols_rows_global(group_idx, P, ntau, sigma, L_cols, R_rows)
        ! 计算 L_cols 和 R_rows（用于 Woodbury 更新）
        ! 
        ! L = B_4 * ... * B_{k+1}
        ! R = B_{k-1} * ... * B_1（不包含 B_k！）
        !
        integer, intent(in) :: group_idx, P(2), ntau
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        complex(kind=8), intent(out) :: L_cols(Ndim, 2), R_rows(2, Ndim)
        
        integer :: kk
        
        L_cols = dcmplx(0.d0, 0.d0)
        L_cols(P(1), 1) = dcmplx(1.d0, 0.d0)
        L_cols(P(2), 2) = dcmplx(1.d0, 0.d0)
        
        R_rows = dcmplx(0.d0, 0.d0)
        R_rows(1, P(1)) = dcmplx(1.d0, 0.d0)
        R_rows(2, P(2)) = dcmplx(1.d0, 0.d0)
        
        ! L = B_4 * ... * B_{k+1}
        do kk = group_idx + 1, 4
            call apply_group_to_cols_global(kk, ntau, sigma, L_cols)
        enddo
        
        ! R = B_{k-1} * ... * B_1（不包含 B_k！）
        do kk = group_idx - 1, 1, -1
            call apply_group_to_rows_global(kk, ntau, sigma, R_rows)
        enddo
        
        return
    end subroutine compute_LR_cols_rows_global

    subroutine apply_group_to_cols_global(grp_idx, nt, sigma, cols)
        integer, intent(in) :: grp_idx, nt
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
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
                sig = sigma(bond_no, nt)
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
                sig = sigma(bond_no, nt)
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
                sig = sigma(bond_no, nt)
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
                sig = sigma(bond_no, nt)
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
    end subroutine apply_group_to_cols_global
    
    subroutine apply_group_to_rows_global(grp_idx, nt, sigma, rows)
        integer, intent(in) :: grp_idx, nt
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
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
                sig = sigma(bond_no, nt)
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
                sig = sigma(bond_no, nt)
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
                sig = sigma(bond_no, nt)
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
                sig = sigma(bond_no, nt)
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
    end subroutine apply_group_to_rows_global
    
    subroutine apply_group_inv_to_rows_global(grp_idx, nt, sigma, rows)
        ! 将 B_grp^{-1} 右乘到 rows 上：rows = rows * B_grp^{-1}
        integer, intent(in) :: grp_idx, nt
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma
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
                sig = sigma(bond_no, nt)
                Cv = cosh(-Op_K%alpha * sig)
                Sv = sinh(-Op_K%alpha * sig)
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
                sig = sigma(bond_no, nt)
                Cv = cosh(-Op_K%alpha * sig)
                Sv = sinh(-Op_K%alpha * sig)
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
                sig = sigma(bond_no, nt)
                Cv = cosh(-Op_K%alpha * sig)
                Sv = sinh(-Op_K%alpha * sig)
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
                sig = sigma(bond_no, nt)
                Cv = cosh(-Op_K%alpha * sig)
                Sv = sinh(-Op_K%alpha * sig)
                do jj = 1, 2
                    t1 = rows(jj, s1) * Cv + rows(jj, s2) * Sv
                    t2 = rows(jj, s1) * Sv + rows(jj, s2) * Cv
                    rows(jj, s1) = t1
                    rows(jj, s2) = t2
                enddo
            enddo
        end select
    end subroutine apply_group_inv_to_rows_global

    subroutine ratioK_fermion_simple(GrU, GrD, sigma_new, ii, ntau, ratio_fermion)
        ! 使用简单的 Sherman-Morrison 公式计算费米子行列式比并更新 Green 函数
        ! det ratio = det(I + Delta * (I - G))
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: GrU, GrD
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_new
        integer, intent(in) :: ii, ntau
        real(kind=8), intent(out) :: ratio_fermion
        
        complex(kind=8), dimension(2, 2) :: Prod_U, Prod_D, Prodinv_U, Prodinv_D
        complex(kind=8), dimension(2, 2) :: GrU_local, GrD_local, mat_tmp
        complex(kind=8), dimension(2, Ndim) :: Vhlp_U, Vhlp_D
        complex(kind=8), dimension(Ndim, 2) :: Uhlp_U, Uhlp_D, temp
        complex(kind=8), dimension(Ndim, Ndim) :: Diff
        complex(kind=8) :: Proddet_U, Proddet_D
        real(kind=8) :: S_old, S_new
        integer :: P(2), nl, nr, j

        S_old = NsigL_K%sigma(ii, ntau)
        S_new = sigma_new(ii, ntau)
        if (abs(S_new - S_old) < 1.d-12) then
            ratio_fermion = 1.d0
            return
        endif

        P(1) = Latt%bond_list(ii, 1)
        P(2) = Latt%bond_list(ii, 2)
        call Op_K%get_delta(S_old, S_new)

        ! det ratio = det(I + Delta * (I - G))
        ! 使用与 CodeXun 相同的形式
        do nr = 1, 2
            do nl = 1, 2
                GrU_local(nl, nr) = ZKRON(nl, nr) - GrU(P(nl), P(nr))
                GrD_local(nl, nr) = ZKRON(nl, nr) - GrD(P(nl), P(nr))
            enddo
        enddo

        mat_tmp = matmul(dcmplx(Op_K%Delta, 0.d0), GrU_local)
        do nr = 1, 2
            do nl = 1, 2
                Prod_U(nl, nr) = ZKRON(nl, nr) + mat_tmp(nl, nr)
            enddo
        enddo
        
        mat_tmp = matmul(dcmplx(Op_K%Delta, 0.d0), GrD_local)
        do nr = 1, 2
            do nl = 1, 2
                Prod_D(nl, nr) = ZKRON(nl, nr) + mat_tmp(nl, nr)
            enddo
        enddo

        Proddet_U = Prod_U(1,1) * Prod_U(2,2) - Prod_U(1,2) * Prod_U(2,1)
        Proddet_D = Prod_D(1,1) * Prod_D(2,2) - Prod_D(1,2) * Prod_D(2,1)
        ratio_fermion = max(abs(Proddet_U * Proddet_D), RATIO_EPS)

        ! 更新 Green 函数
        Prodinv_U(1,1) =  Prod_U(2,2) / Proddet_U
        Prodinv_U(2,2) =  Prod_U(1,1) / Proddet_U
        Prodinv_U(1,2) = -Prod_U(1,2) / Proddet_U
        Prodinv_U(2,1) = -Prod_U(2,1) / Proddet_U
        
        Prodinv_D(1,1) =  Prod_D(2,2) / Proddet_D
        Prodinv_D(2,2) =  Prod_D(1,1) / Proddet_D
        Prodinv_D(1,2) = -Prod_D(1,2) / Proddet_D
        Prodinv_D(2,1) = -Prod_D(2,1) / Proddet_D

        ! Sherman-Morrison 更新（与 CodeXun 相同）
        ! Spin Up
        Uhlp_U = dcmplx(0.d0, 0.d0)
        Vhlp_U = dcmplx(0.d0, 0.d0)
        do nl = 1, 2
            do j = 1, Ndim
                Uhlp_U(j, nl) = GrU(j, P(nl))
                Vhlp_U(nl, j) = -Op_K%Delta(nl, 1) * GrU(P(1), j) &
                               -Op_K%Delta(nl, 2) * GrU(P(2), j)
            enddo
            Vhlp_U(nl, P(1)) = Vhlp_U(nl, P(1)) + Op_K%Delta(nl, 1)
            Vhlp_U(nl, P(2)) = Vhlp_U(nl, P(2)) + Op_K%Delta(nl, 2)
        enddo
        temp = matmul(Uhlp_U, Prodinv_U)
        Diff = matmul(temp, Vhlp_U)
        GrU = GrU - Diff

        ! Spin Down
        Uhlp_D = dcmplx(0.d0, 0.d0)
        Vhlp_D = dcmplx(0.d0, 0.d0)
        do nl = 1, 2
            do j = 1, Ndim
                Uhlp_D(j, nl) = GrD(j, P(nl))
                Vhlp_D(nl, j) = -Op_K%Delta(nl, 1) * GrD(P(1), j) &
                               -Op_K%Delta(nl, 2) * GrD(P(2), j)
            enddo
            Vhlp_D(nl, P(1)) = Vhlp_D(nl, P(1)) + Op_K%Delta(nl, 1)
            Vhlp_D(nl, P(2)) = Vhlp_D(nl, P(2)) + Op_K%Delta(nl, 2)
        enddo
        temp = matmul(Uhlp_D, Prodinv_D)
        Diff = matmul(temp, Vhlp_D)
        GrD = GrD - Diff

        return
    end subroutine ratioK_fermion_simple
    
    function get_bond_group(bond_idx) result(group_idx)
        ! 确定 bond 属于哪个 group
        integer, intent(in) :: bond_idx
        integer :: group_idx
        integer :: ii
        
        ! 搜索 group_1
        do ii = 1, size(Latt%group_1)
            if (Latt%group_1(ii) == bond_idx) then
                group_idx = 1
                return
            endif
        enddo
        ! 搜索 group_2
        do ii = 1, size(Latt%group_2)
            if (Latt%group_2(ii) == bond_idx) then
                group_idx = 2
                return
            endif
        enddo
        ! 搜索 group_3
        do ii = 1, size(Latt%group_3)
            if (Latt%group_3(ii) == bond_idx) then
                group_idx = 3
                return
            endif
        enddo
        ! 搜索 group_4
        do ii = 1, size(Latt%group_4)
            if (Latt%group_4(ii) == bond_idx) then
                group_idx = 4
                return
            endif
        enddo
        
        ! 如果找不到，默认返回 1
        group_idx = 1
    end function get_bond_group

    subroutine GlobalK_prop_L(PropU, PropD, log_ratio_fermion, sigma_new, nt)
        ! 向左传播时计算费米子行列式比并更新 Green 函数
        ! 使用正确的 Woodbury 公式（包含 P[λ]）
        class(Propagator), intent(inout) :: PropU, PropD
        real(kind=8), intent(inout) :: log_ratio_fermion
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_new
        integer, intent(in) :: nt
        integer :: ii, grp, idx
        real(kind=8) :: ratio

        ! 按组顺序处理（从 group_4 到 group_1）
        ! 对于每个组内的 bonds，使用 Woodbury 公式
        do grp = 4, 1, -1
            select case (grp)
            case (4)
                do idx = size(Latt%group_4), 1, -1
                    ii = Latt%group_4(idx)
                    call ratioK_fermion_woodbury(PropU%Gr, PropD%Gr, NsigL_K%sigma, &
                        sigma_new, ii, nt, grp, ratio)
                    log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
                enddo
            case (3)
                do idx = size(Latt%group_3), 1, -1
                    ii = Latt%group_3(idx)
                    call ratioK_fermion_woodbury(PropU%Gr, PropD%Gr, NsigL_K%sigma, &
                        sigma_new, ii, nt, grp, ratio)
                    log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
                enddo
            case (2)
                do idx = size(Latt%group_2), 1, -1
                    ii = Latt%group_2(idx)
                    call ratioK_fermion_woodbury(PropU%Gr, PropD%Gr, NsigL_K%sigma, &
                        sigma_new, ii, nt, grp, ratio)
                    log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
                enddo
            case (1)
                do idx = size(Latt%group_1), 1, -1
                    ii = Latt%group_1(idx)
                    call ratioK_fermion_woodbury(PropU%Gr, PropD%Gr, NsigL_K%sigma, &
                        sigma_new, ii, nt, grp, ratio)
                    log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
                enddo
            end select
        enddo
        
        ! wrap G
        call Op_K%mmult_L(PropU%Gr, Latt, sigma_new, nt, 1)
        call Op_K%mmult_L(PropD%Gr, Latt, sigma_new, nt, 1)
        call Op_K%mmult_R(PropU%Gr, Latt, sigma_new, nt, -1)
        call Op_K%mmult_R(PropD%Gr, Latt, sigma_new, nt, -1)
        call Op_K%mmult_L(PropU%UUL, Latt, sigma_new, nt, 1)
        call Op_K%mmult_L(PropD%UUL, Latt, sigma_new, nt, 1)
        return
    end subroutine GlobalK_prop_L

    subroutine GlobalK_prop_R(PropU, PropD, log_ratio_fermion, sigma_new, nt)
        ! 向右传播时计算费米子行列式比并更新 Green 函数
        ! 使用正确的 Woodbury 公式（包含 P[λ]）
        class(Propagator), intent(inout) :: PropU, PropD
        real(kind=8), intent(inout) :: log_ratio_fermion
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_new
        integer, intent(in) :: nt
        integer :: ii, grp, idx
        real(kind=8) :: ratio

        ! 先 wrap G（使用旧的 sigma）
        call Op_K%mmult_R(PropU%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_R(PropD%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropU%Gr, Latt, NsigL_K%sigma, nt, -1)
        call Op_K%mmult_L(PropD%Gr, Latt, NsigL_K%sigma, nt, -1)
        
        ! 按组顺序处理（从 group_1 到 group_4）
        do grp = 1, 4
            select case (grp)
            case (1)
                do idx = 1, size(Latt%group_1)
                    ii = Latt%group_1(idx)
                    call ratioK_fermion_woodbury(PropU%Gr, PropD%Gr, NsigL_K%sigma, &
                        sigma_new, ii, nt, grp, ratio)
                    log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
                enddo
            case (2)
                do idx = 1, size(Latt%group_2)
                    ii = Latt%group_2(idx)
                    call ratioK_fermion_woodbury(PropU%Gr, PropD%Gr, NsigL_K%sigma, &
                        sigma_new, ii, nt, grp, ratio)
                    log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
                enddo
            case (3)
                do idx = 1, size(Latt%group_3)
                    ii = Latt%group_3(idx)
                    call ratioK_fermion_woodbury(PropU%Gr, PropD%Gr, NsigL_K%sigma, &
                        sigma_new, ii, nt, grp, ratio)
                    log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
                enddo
            case (4)
                do idx = 1, size(Latt%group_4)
                    ii = Latt%group_4(idx)
                    call ratioK_fermion_woodbury(PropU%Gr, PropD%Gr, NsigL_K%sigma, &
                        sigma_new, ii, nt, grp, ratio)
                    log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
                enddo
            end select
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