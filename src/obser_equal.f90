module ObserEqual_mod
    use ProcessMatrix
    use DQMC_Model_mod
    use calc_basic
    implicit none
    
    type, public :: ObserEqual
        real(kind=8) :: flux_sum, flux_sq
        real(kind=8) :: flux_avg, chi_flux
        real(kind=8) :: pair_sc_acc, pair_sc_q0
        real(kind=8) :: cdw_pi_acc, cdw_shift_acc
        real(kind=8) :: cdw_pi, cdw_shift, cdw_ratio
        real(kind=8) :: density
    contains
        procedure :: make  => Obs_equal_make
        final     :: Obs_equal_clear
        procedure :: reset => Obs_equal_reset
        procedure :: ave   => Obs_equal_ave
        procedure :: calc  => Obs_equal_calc
    end type ObserEqual
    
contains
    subroutine Obs_equal_make(this)
        class(ObserEqual), intent(inout) :: this
        return
    end subroutine Obs_equal_make
    
    subroutine Obs_equal_clear(this)
        type(ObserEqual), intent(inout) :: this
        return
    end subroutine Obs_equal_clear
    
    subroutine Obs_equal_reset(this)
        class(ObserEqual), intent(inout) :: this
        this%flux_sum = 0.d0
        this%flux_sq = 0.d0
        this%flux_avg = 0.d0
        this%chi_flux = 0.d0
        this%pair_sc_acc = 0.d0
        this%pair_sc_q0 = 0.d0
        this%cdw_pi_acc = 0.d0
        this%cdw_shift_acc = 0.d0
        this%cdw_pi = 0.d0
        this%cdw_shift = 0.d0
        this%cdw_ratio = 0.d0
        this%density = 0.d0
        return
    end subroutine Obs_equal_reset
    
    subroutine Obs_equal_ave(this, Nobs)
        class(ObserEqual), intent(inout) :: this
        integer, intent(in) :: Nobs
        real(kind=8) :: znorm
        real(kind=8) :: flux_mean, flux_sq_mean
        if (Nobs <= 0) return
        znorm = 1.d0 / dble(Nobs)
        flux_mean = this%flux_sum * znorm
        flux_sq_mean = this%flux_sq * znorm
        this%flux_avg = flux_mean
        this%chi_flux = Beta * dble(Lq) * max(0.d0, flux_sq_mean - flux_mean*flux_mean)
        this%pair_sc_q0 = (this%pair_sc_acc * znorm) / dble(Lq)
        this%cdw_pi = (this%cdw_pi_acc * znorm) / dble(Lq)
        this%cdw_shift = (this%cdw_shift_acc * znorm) / dble(Lq)
        if (this%cdw_pi > 1.d-12) then
            this%cdw_ratio = 1.d0 - this%cdw_shift / this%cdw_pi
        else
            this%cdw_ratio = 0.d0
        endif
        this%density = this%density * znorm
        return
    end subroutine Obs_equal_ave
    
    subroutine Obs_equal_calc(this, PropU, PropD, ntau)
        class(ObserEqual), intent(inout) :: this
        class(Propagator), intent(in)    :: PropU, PropD
        integer, intent(in) :: ntau
        complex(kind=8), dimension(Ndim, Ndim) :: Gru, Grd
        real(kind=8) :: flux_tau, pair_sum, cdw_sum, cdw_shift_sum
        real(kind=8) :: n_up_i, n_dn_i, n_up_j, n_dn_j
        real(kind=8) :: density_corr, phase_pi, phase_shift
        complex(kind=8) :: Gup_ij, Gup_ji, Gdn_ij, Gdn_ji
        real(kind=8) :: delta_ij, qshift_x
        integer :: ii, jj
        integer :: site, right_site, up_site, diag_site
        integer :: bond_bottom, bond_right, bond_top, bond_left
        real(kind=8) :: density_slice
        real(kind=8) :: stag_i, stag_j
        real(kind=8) :: dx, dy

        Gru = PropU%Gr
        Grd = PropD%Gr
        qshift_x = PI - 2.d0 * PI / dble(Nlx)

! 规范场 Z₂ 磁通：B(\tau) = (1/N_\square) \sum_{\square} \prod_{b\in\square} \sigma_b^z(\tau)
        flux_tau = 0.d0
        do site = 1, Lq
            right_site = Latt%L_bonds(site, 1)
            up_site    = Latt%L_bonds(site, 2)
            diag_site  = Latt%L_bonds(up_site, 1)
            bond_bottom = Latt%inv_bond_list(site, right_site)
            bond_left   = Latt%inv_bond_list(site, up_site)
            bond_top    = Latt%inv_bond_list(up_site, diag_site)
            bond_right  = Latt%inv_bond_list(right_site, Latt%L_bonds(right_site, 2))
            flux_tau = flux_tau + NsigL_K%sigma(bond_bottom, ntau) &
                & * NsigL_K%sigma(bond_right, ntau) * NsigL_K%sigma(bond_top, ntau) * NsigL_K%sigma(bond_left, ntau)
        enddo
        flux_tau = flux_tau / dble(Lq)
        this%flux_sum = this%flux_sum + flux_tau
        this%flux_sq = this%flux_sq + flux_tau * flux_tau

! s-wave 对结构因子：P_SC(q=0) = (1/N) \sum_{ij} \langle b_i^\dagger b_j + h.c. \rangle
        pair_sum = 0.d0
        cdw_sum = 0.d0
        cdw_shift_sum = 0.d0
        density_slice = 0.d0
        do ii = 1, Lq
            n_up_i = real(1.d0 - Gru(ii, ii))
            n_dn_i = real(1.d0 - Grd(ii, ii))
            density_slice = density_slice + n_up_i + n_dn_i
            do jj = 1, Lq
                n_up_j = real(1.d0 - Gru(jj, jj))
                n_dn_j = real(1.d0 - Grd(jj, jj))
                delta_ij = 0.d0
                if (ii == jj) delta_ij = 1.d0
                Gup_ij = Gru(ii, jj)
                Gup_ji = Gru(jj, ii)
                Gdn_ij = Grd(ii, jj)
                Gdn_ji = Grd(jj, ii)
                pair_sum = pair_sum + real((delta_ij - Gup_ij) * (delta_ij - Gdn_ij))
! CDW 结构因子：P_D(q) = (1/N) \sum_{ij} e^{i q\cdot(r_i - r_j)} \langle n_i n_j \rangle
                density_corr = n_up_i * n_up_j - real((delta_ij - Gup_ji) * Gup_ij)
                density_corr = density_corr + n_dn_i * n_dn_j - real((delta_ij - Gdn_ji) * Gdn_ij)
                density_corr = density_corr + n_up_i * n_dn_j + n_dn_i * n_up_j
                stag_i = merge(1.d0, -1.d0, mod(Latt%n_list(ii,1) + Latt%n_list(ii,2), 2) == 0)
                stag_j = merge(1.d0, -1.d0, mod(Latt%n_list(jj,1) + Latt%n_list(jj,2), 2) == 0)
                phase_pi = stag_i * stag_j
                dx = dble(Latt%n_list(ii,1) - Latt%n_list(jj,1))
                dy = dble(Latt%n_list(ii,2) - Latt%n_list(jj,2))
                phase_shift = cos(qshift_x * dx + PI * dy)
                cdw_sum = cdw_sum + phase_pi * density_corr
                cdw_shift_sum = cdw_shift_sum + phase_shift * density_corr
            enddo
        enddo
        this%pair_sc_acc = this%pair_sc_acc + pair_sum
        this%cdw_pi_acc = this%cdw_pi_acc + cdw_sum
        this%cdw_shift_acc = this%cdw_shift_acc + cdw_shift_sum
        this%density = this%density + density_slice / dble(Lq)
        return
    end subroutine Obs_equal_calc
end module ObserEqual_mod
