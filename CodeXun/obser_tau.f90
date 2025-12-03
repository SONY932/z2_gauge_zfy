module ObserTau_mod
    use ProcessMatrix
    use DQMC_Model_mod
    implicit none
    
    type, public :: ObserTau
        complex(kind=8), dimension(:,:), allocatable :: single_tau, spin_tau, charge_d_tau, charge_s_tau
        complex(kind=8), dimension(:,:), allocatable :: pair_d_tau, pair_s_tau, curxx_tau
    contains
        procedure :: make => Obs_tau_make
        final :: Obs_tau_clear
        procedure :: reset => Obs_tau_reset
        procedure :: ave => Obs_tau_ave
        procedure :: calc => Obs_tau_calc
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
    
    subroutine Obs_tau_calc(this, PropGr, ntau) ! here ntau ranges from 1 to Ltrot, with ntau=1 denotes equal time correlation
        class(ObserTau), intent(inout) :: this
        class(PropGreen), intent(in) :: PropGr
        integer, intent(in) :: ntau
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Gr00up, Gr00upc, Grttup, Grttupc, Gr0tup, Grt0up
        complex(kind=8), dimension(Ndim) :: tmp_t, tmp_0
        integer :: i, ii, jj, imj, ip, i1, i2, jp, j1, j2, ip1, ip2, jp1, jp2, iiy, jjy
        complex(kind=8) :: diag1, diag2, cross1, cross2, term11, term12, term21, term22
        complex(kind=8) :: phase_i_p, phase_i_m, phase_j_p, phase_j_m
        real(kind=8) :: Ai, Aj
        
        Gr00up = PropGr%Gr00 ! Gr00(i, j) = <c_i(0) * c^+_j(0)>
        Gr00upc = ZKRON - transpose(Gr00up) ! Gr00c(i, j) = <c^+_i(0) * c_j(0)>
        Grttup = PropGr%Grtt ! Grtt(i, j) = <c_i(\tau) * c^+_j(\tau)>
        Grttupc = ZKRON - transpose(Grttup) ! Grttc(i, j) = <c^+_i(\tau) * c_j(\tau)>
        Grt0up = PropGr%Grt0 ! Grt0(i, j) = <c_i(\tau) * c^+_j(0)>
        Gr0tup = PropGr%Gr0t ! Gr0t(i, j) = - <c^+_j(\tau) * c_i(0)>

        do i = 1, Ndim
            tmp_t(i) = Grttupc(i, i) + dconjg(Grttupc(i, i))
            tmp_0(i) = Gr00upc(i, i) + dconjg(Gr00upc(i, i))
        enddo
        
        do jj = 1, Lq
            do ii = 1, Lq
                imj = Latt%imj(ii, jj)
                i1 = Latt%inv_o_list(ii, 1)
                i2 = Latt%inv_o_list(ii, 2)
                j1 = Latt%inv_o_list(jj, 1)
                j2 = Latt%inv_o_list(jj, 2)
! SDW susceptibility (bosonic)
                this%spin_tau(imj, ntau) = this%spin_tau(imj, ntau) &
                    &   + NsigL_K%phi(1, ii, ntau) * NsigL_K%phi(1, jj, 1) + NsigL_K%phi(2, ii, ntau) * NsigL_K%phi(2, jj, 1)
! pairing susceptibility
                diag1 = Gr0tup(j1, i1) * dconjg(Gr0tup(j1, i1))
                diag2 = Gr0tup(j2, i2) * dconjg(Gr0tup(j2, i2))
                cross1 = Gr0tup(j1, i2) * dconjg(Gr0tup(j1, i2))
                cross2 = Gr0tup(j2, i1) * dconjg(Gr0tup(j2, i1))
                this%pair_d_tau(imj, ntau) = this%pair_d_tau(imj, ntau) + (diag1 + diag2 + cross1 + cross2) * 4.0
                this%pair_s_tau(imj, ntau) = this%pair_s_tau(imj, ntau) + (diag1 + diag2 - cross1 - cross2) * 4.0
! charge susceptibility
                diag1 = - Gr0tup(j1, i1) * Grt0up(i1, j1) - dconjg(Gr0tup(j1, i1) * Grt0up(i1, j1))
                diag2 = - Gr0tup(j2, i2) * Grt0up(i2, j2) - dconjg(Gr0tup(j2, i2) * Grt0up(i2, j2))
                cross1 = - Gr0tup(j2, i1) * Grt0up(i1, j2) - dconjg(Gr0tup(j2, i1) * Grt0up(i1, j2))
                cross2 = - Gr0tup(j1, i2) * Grt0up(i2, j1) - dconjg(Gr0tup(j1, i2) * Grt0up(i2, j1))
                this%charge_d_tau(imj, ntau) = this%charge_d_tau(imj, ntau) &
                    &   + diag1 + diag2 - cross1 - cross2 + (tmp_t(i1) - tmp_t(i2)) * (tmp_0(j1) - tmp_0(j2))
                this%charge_s_tau(imj, ntau) = this%charge_s_tau(imj, ntau) &
                    &   + diag1 + diag2 + cross1 + cross2 + (tmp_t(i1) + tmp_t(i2)) * (tmp_0(j1) + tmp_0(j2))
! single-particle correlation
                term11 = Grt0up(i1, j1) + Grt0up(i2, j2)
                this%single_tau(imj, ntau) = this%single_tau(imj, ntau) + term11 + dconjg(term11)
            enddo
        enddo
! superfluid density
        do jj = 1, Lq
            do ii = 1, Lq
                imj = Latt%imj(ii, jj)
                ip = Latt%L_bonds(ii, 1)
                jp = Latt%L_bonds(jj, 1)
                i1 = Latt%inv_o_list(ii, 1)
                ip1 = Latt%inv_o_list(ip, 1)
                i2 = Latt%inv_o_list(ii, 2)
                ip2 = Latt%inv_o_list(ip, 2)
                j1 = Latt%inv_o_list(jj, 1)
                jp1 = Latt%inv_o_list(jp, 1)
                j2 = Latt%inv_o_list(jj, 2)
                jp2 = Latt%inv_o_list(jp, 2)
                
                iiy = Latt%n_list(ii, 2)
                jjy = Latt%n_list(jj, 2)
                Ai = -2.d0 * Pi * NB_field * dble(iiy)/dble(Lq)
                Aj = -2.d0 * Pi * NB_field * dble(jjy)/dble(Lq)
                phase_i_p = exp( dcmplx(0.d0, 1.d0) * Ai)
                phase_i_m = exp(-dcmplx(0.d0, 1.d0) * Ai)
                phase_j_p = exp( dcmplx(0.d0, 1.d0) * Aj)
                phase_j_m = exp(-dcmplx(0.d0, 1.d0) * Aj)
                
                term11 = - Gr0tup(j1, ip1) * Grt0up(i1, jp1) * phase_i_m * phase_j_m - Gr0tup(jp1, i1) * Grt0up(ip1, j1) * phase_i_p * phase_j_p &
                    &   + Gr0tup(j1, i1) * Grt0up(ip1, jp1) * phase_i_p * phase_j_m + Gr0tup(jp1, ip1) * Grt0up(i1, j1) * phase_i_m * phase_j_p
                term12 = - Gr0tup(j2, ip1) * Grt0up(i1, jp2) * phase_i_m * phase_j_m - Gr0tup(jp2, i1) * Grt0up(ip1, j2) * phase_i_p * phase_j_p &
                    &   + Gr0tup(j2, i1) * Grt0up(ip1, jp2) * phase_i_p * phase_j_m + Gr0tup(jp2, ip1) * Grt0up(i1, j2) * phase_i_m * phase_j_p
                term21 = - Gr0tup(j1, ip2) * Grt0up(i2, jp1) * phase_i_m * phase_j_m - Gr0tup(jp1, i2) * Grt0up(ip2, j1) * phase_i_p * phase_j_p &
                    &   + Gr0tup(j1, i2) * Grt0up(ip2, jp1) * phase_i_p * phase_j_m + Gr0tup(jp1, ip2) * Grt0up(i2, j1) * phase_i_m * phase_j_p
                term22 = - Gr0tup(j2, ip2) * Grt0up(i2, jp2) * phase_i_m * phase_j_m - Gr0tup(jp2, i2) * Grt0up(ip2, j2) * phase_i_p * phase_j_p &
                    &   + Gr0tup(j2, i2) * Grt0up(ip2, jp2) * phase_i_p * phase_j_m + Gr0tup(jp2, ip2) * Grt0up(i2, j2) * phase_i_m * phase_j_p
                this%curxx_tau(imj, ntau) = this%curxx_tau(imj, ntau) &
                    &   - RT(1,1) * RT(1,1) * term11 - RT(1,1) * RT(2,1) * (term12+term21) - RT(2,1) * RT(2,1) * term22
                term11 = dconjg(term11); term12 = dconjg(term12)
                term21 = dconjg(term21); term22 = dconjg(term22)
                this%curxx_tau(imj, ntau) = this%curxx_tau(imj, ntau) &
                    &   - RT(1,1) * RT(1,1) * term11 - RT(1,1) * RT(2,1) * (term12+term21) - RT(2,1) * RT(2,1) * term22
            enddo
        enddo
        return
    end subroutine Obs_tau_calc
end module ObserTau_mod