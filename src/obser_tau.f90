module ObserTau_mod
    use ProcessMatrix
    use DQMC_Model_mod
    use calc_basic
    implicit none
    
    type, public :: ObserTau
        complex(kind=8), dimension(:,:), allocatable :: single_tau, spin_tau, charge_d_tau, charge_s_tau
        complex(kind=8), dimension(:,:), allocatable :: pair_d_tau, pair_s_tau, curxx_tau
    contains
        procedure :: make  => Obs_tau_make
        final     :: Obs_tau_clear
        procedure :: reset => Obs_tau_reset
        procedure :: ave   => Obs_tau_ave
        procedure :: calc  => Obs_tau_calc
    end type ObserTau
    
contains
    subroutine Obs_tau_make(this)
        class(ObserTau), intent(inout) :: this
        allocate(this%spin_tau(Lq, Ltrot), this%single_tau(Lq, Ltrot))
        allocate(this%pair_d_tau(Lq, Ltrot), this%pair_s_tau(Lq, Ltrot))
        allocate(this%charge_d_tau(Lq, Ltrot), this%charge_s_tau(Lq, Ltrot))
        allocate(this%curxx_tau(Lq, Ltrot))
        return
    end subroutine Obs_tau_make
    
    subroutine Obs_tau_reset(this)
        class(ObserTau), intent(inout) :: this
        this%spin_tau = dcmplx(0.d0, 0.d0); this%single_tau = dcmplx(0.d0, 0.d0)
        this%charge_d_tau = dcmplx(0.d0, 0.d0); this%charge_s_tau = dcmplx(0.d0, 0.d0)
        this%pair_d_tau = dcmplx(0.d0, 0.d0); this%pair_s_tau = dcmplx(0.d0, 0.d0)
        this%curxx_tau = dcmplx(0.d0, 0.d0)
        return
    end subroutine Obs_tau_reset
    
    subroutine Obs_tau_ave(this, Nobst)
        class(ObserTau), intent(inout) :: this
        integer, intent(in) :: Nobst
        complex(kind=8) :: znorm
        znorm = dcmplx(1.d0 / dble(Nobst), 0.d0)
        this%spin_tau = znorm * this%spin_tau; this%single_tau = znorm * this%single_tau
        this%pair_d_tau = znorm * this%pair_d_tau; this%pair_s_tau = znorm * this%pair_s_tau
        this%charge_d_tau = znorm * this%charge_d_tau; this%charge_s_tau = znorm * this%charge_s_tau
        this%curxx_tau = znorm * this%curxx_tau
        return
    end subroutine Obs_tau_ave
    
    subroutine Obs_tau_clear(this)
        type(ObserTau), intent(inout) :: this
        deallocate(this%single_tau, this%spin_tau, this%pair_d_tau, this%pair_s_tau)
        deallocate(this%charge_d_tau, this%charge_s_tau, this%curxx_tau)
        return
    end subroutine Obs_tau_clear
    
    subroutine Obs_tau_calc(this, PropGrU, PropGrD, ntau) ! ntau: 1..Ltrot, 1 表示等时
        class(ObserTau), intent(inout) :: this
        class(PropGreen), intent(in) :: PropGrU, PropGrD
        integer, intent(in) :: ntau
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Gr00u, Gr00uc, Grttu, Grttuc, Gr0tu, Grt0u
        complex(kind=8), dimension(Ndim, Ndim) :: Gr00d, Gr00dc, Grttd, Grttdc, Gr0td, Grt0d
        complex(kind=8), dimension(Ndim) :: tmp_t_u, tmp_0_u, tmp_t_d, tmp_0_d
        integer :: ii, jj, imj, ip, jp, iiy, jjy
        complex(kind=8) :: diag1, diag2, term
        complex(kind=8) :: phase_i_p, phase_i_m, phase_j_p, phase_j_m
        real(kind=8) :: Ai, Aj

        ! 从 G_0 获取 Green 函数
        ! 注意：P[λ] 投影对含时 Green 函数的影响通过简单的对角缩放处理
        ! G(i,j) → λ_i × G(i,j) × λ_j，这是一个 O(N²) 的高效操作
        Gr00u  = PropGrU%Gr00; Grttu  = PropGrU%Grtt; Grt0u = PropGrU%Grt0; Gr0tu = PropGrU%Gr0t
        Gr00d  = PropGrD%Gr00; Grttd  = PropGrD%Grtt; Grt0d = PropGrD%Grt0; Gr0td = PropGrD%Gr0t
        
        ! 应用 P[λ] 投影（简单对角缩放，O(N²)）
        call apply_lambda_projection(Gr00u)
        call apply_lambda_projection(Grttu)
        call apply_lambda_projection(Grt0u)
        call apply_lambda_projection(Gr0tu)
        call apply_lambda_projection(Gr00d)
        call apply_lambda_projection(Grttd)
        call apply_lambda_projection(Grt0d)
        call apply_lambda_projection(Gr0td)
        
        Gr00uc = ZKRON - transpose(Gr00u); Grttuc = ZKRON - transpose(Grttu)
        Gr00dc = ZKRON - transpose(Gr00d); Grttdc = ZKRON - transpose(Grttd)

        do ii = 1, Ndim
            tmp_t_u(ii) = Grttuc(ii, ii) + dconjg(Grttuc(ii, ii))
            tmp_0_u(ii) = Gr00uc(ii, ii) + dconjg(Gr00uc(ii, ii))
            tmp_t_d(ii) = Grttdc(ii, ii) + dconjg(Grttdc(ii, ii))
            tmp_0_d(ii) = Gr00dc(ii, ii) + dconjg(Gr00dc(ii, ii))
        enddo

        do jj = 1, Lq
            do ii = 1, Lq
                imj = Latt%imj(ii, jj)
! SDW susceptibility (bosonic, scalar Ising)
                this%spin_tau(imj, ntau) = this%spin_tau(imj, ntau) + dcmplx( NsigL_K%sigma(ii, ntau) * NsigL_K%sigma(jj, 1), 0.d0 )
! pairing susceptibility（自旋分离：无自旋翻转项，cross=0）
                diag1 = Gr0tu(jj, ii) * dconjg(Gr0tu(jj, ii))
                diag2 = Gr0td(jj, ii) * dconjg(Gr0td(jj, ii))
                this%pair_d_tau(imj, ntau) = this%pair_d_tau(imj, ntau) + (diag1 + diag2) * 4.0
                this%pair_s_tau(imj, ntau) = this%pair_s_tau(imj, ntau) + (diag1 + diag2) * 4.0
! charge susceptibility（cross=0）
                this%charge_d_tau(imj, ntau) = this%charge_d_tau(imj, ntau) &
                    &   - Gr0tu(jj, ii) * Grt0u(ii, jj) - dconjg(Gr0tu(jj, ii) * Grt0u(ii, jj)) &
                    &   - Gr0td(jj, ii) * Grt0d(ii, jj) - dconjg(Gr0td(jj, ii) * Grt0d(ii, jj)) &
                    &   + (tmp_t_u(ii) - tmp_t_d(ii)) * (tmp_0_u(jj) - tmp_0_d(jj))
                this%charge_s_tau(imj, ntau) = this%charge_s_tau(imj, ntau) &
                    &   - Gr0tu(jj, ii) * Grt0u(ii, jj) - dconjg(Gr0tu(jj, ii) * Grt0u(ii, jj)) &
                    &   - Gr0td(jj, ii) * Grt0d(ii, jj) - dconjg(Gr0td(jj, ii) * Grt0d(ii, jj)) &
                    &   + (tmp_t_u(ii) + tmp_t_d(ii)) * (tmp_0_u(jj) + tmp_0_d(jj))
! single-particle correlation
                term = Grt0u(ii, jj) + Grt0d(ii, jj)
                this%single_tau(imj, ntau) = this%single_tau(imj, ntau) + term + dconjg(term)
            enddo
        enddo

! superfluid density（按旧代码公式分别对两自旋求和）
        do jj = 1, Lq
            do ii = 1, Lq
                imj = Latt%imj(ii, jj)
                ip = Latt%L_bonds(ii, 1)
                jp = Latt%L_bonds(jj, 1)
                
                iiy = Latt%n_list(ii, 2)
                jjy = Latt%n_list(jj, 2)
                Ai = -2.d0 * Pi * NB_field * dble(iiy)/dble(Lq)
                Aj = -2.d0 * Pi * NB_field * dble(jjy)/dble(Lq)
                phase_i_p = exp( dcmplx(0.d0, 1.d0) * Ai)
                phase_i_m = exp(-dcmplx(0.d0, 1.d0) * Ai)
                phase_j_p = exp( dcmplx(0.d0, 1.d0) * Aj)
                phase_j_m = exp(-dcmplx(0.d0, 1.d0) * Aj)
                
                ! ↑ 自旋贡献
                term = - Gr0tu(jj, ip) * Grt0u(ii, jp) * phase_i_m * phase_j_m - Gr0tu(jp, ii) * Grt0u(ip, jj) * phase_i_p * phase_j_p &
                    &   + Gr0tu(jj, ii) * Grt0u(ip, jp) * phase_i_p * phase_j_m + Gr0tu(jp, ip) * Grt0u(ii, jj) * phase_i_m * phase_j_p
                this%curxx_tau(imj, ntau) = this%curxx_tau(imj, ntau) - RT * RT * term - RT * RT * dconjg(term)
                ! ↓ 自旋贡献
                term = - Gr0td(jj, ip) * Grt0d(ii, jp) * phase_i_m * phase_j_m - Gr0td(jp, ii) * Grt0d(ip, jj) * phase_i_p * phase_j_p &
                    &   + Gr0td(jj, ii) * Grt0d(ip, jp) * phase_i_p * phase_j_m + Gr0td(jp, ip) * Grt0d(ii, jj) * phase_i_m * phase_j_p
                this%curxx_tau(imj, ntau) = this%curxx_tau(imj, ntau) - RT * RT * term - RT * RT * dconjg(term)
            enddo
        enddo
        return
    end subroutine Obs_tau_calc

    subroutine apply_lambda_projection(Gr)
        ! 对 Green 函数应用简单的 P[λ] 投影（对角缩放）
        ! 这是一个 O(N²) 的高效操作
        ! G(i, j) → λ_i × G(i, j) × λ_j
        ! 这等价于 P[λ] × G × P[λ]
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gr
        integer :: i, j
        real(kind=8) :: lam_i, lam_j
        logical :: all_one
        
        ! 快速检查：如果 λ 全为 1，直接返回
        all_one = .true.
        do i = 1, Lq
            if (abs(NsigL_K%lambda(i) - 1.d0) > 1.d-12) then
                all_one = .false.
                exit
            endif
        enddo
        if (all_one) return
        
        ! 应用 P[λ] × G × P[λ]
        do j = 1, Ndim
            lam_j = NsigL_K%lambda(j)
            do i = 1, Ndim
                lam_i = NsigL_K%lambda(i)
                Gr(i, j) = Gr(i, j) * dcmplx(lam_i * lam_j, 0.d0)
            enddo
        enddo
        return
    end subroutine apply_lambda_projection

end module ObserTau_mod
