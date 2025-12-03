module ObserEqual_mod
    use ProcessMatrix
    use DQMC_Model_mod
    implicit none
    
    type, public :: ObserEqual
        real(kind=8), dimension(:), allocatable :: den_occ
        real(kind=8), dimension(:,:), allocatable :: bond_occ, spin_occ
        real(kind=8) :: el_ke(Nbond), density, spin_avg(Nspin), spin_order(Nspin), diam
    contains
        procedure :: make => Obs_equal_make
        final :: Obs_equal_clear
        procedure :: reset => Obs_equal_reset
        procedure :: ave => Obs_equal_ave
        procedure :: calc => Obs_equal_calc
    end type ObserEqual
    
contains
    subroutine Obs_equal_make(this)
        class(ObserEqual), intent(inout) :: this
        allocate( this%spin_occ(Nspin, Lq), this%den_occ(Lq), this%bond_occ(Lq, Nbond) )
        return
    end subroutine Obs_equal_make
    
    subroutine Obs_equal_clear(this)
        type(ObserEqual), intent(inout) :: this
        deallocate(this%den_occ, this%spin_occ, this%bond_occ)
        return
    end subroutine Obs_equal_clear
    
    subroutine Obs_equal_reset(this)
        class(ObserEqual), intent(inout) :: this
        this%bond_occ = dcmplx(0.d0, 0.d0)
        this%spin_occ = 0.d0
        this%den_occ = 0.d0
        this%density = 0.d0
        this%el_ke = 0.d0
        this%diam = 0.d0
        this%spin_avg = 0.d0
        this%spin_order = 0.d0
        return
    end subroutine Obs_equal_reset
    
    subroutine Obs_equal_ave(this, Nobs)
        class(ObserEqual), intent(inout) :: this
        integer, intent(in) :: Nobs
        real(kind=8) :: znorm
        znorm = 1.d0 / dble(Nobs)
        this%bond_occ = znorm * this%bond_occ
        this%spin_occ = znorm * this%spin_occ
        this%den_occ = znorm * this%den_occ
        this%density = znorm * this%density
        this%el_ke = znorm * this%el_ke
        this%diam = znorm * this%diam
        this%spin_avg = znorm * this%spin_avg
        this%spin_order = znorm * this%spin_order
        return
    end subroutine Obs_equal_ave
    
    subroutine Obs_equal_calc(this, Prop, ntau)
!   Arguments: 
        class(ObserEqual), intent(inout) :: this
        class(Propagator), intent(in) :: Prop
        integer, intent(in) :: ntau
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Grupc, Grup
        real(kind=8), dimension(Ndim) :: tmpd
        integer :: i, ii, no, ns, i1, nf, sign, iy
        real(kind=8) :: tmp, Ai
        complex(kind=8) :: phase_p, phase_m
        
        Grup = Prop%Gr !   Gr(i, j) = <c_i c^+_j >
        Grupc = ZKRON - transpose(Grup) !   Grc(i, j) = <c^+_i c_j >
        
        do i = 1, Ndim
            tmpd(i) = real(Grupc(i, i) + dconjg(Grupc(i, i)))
        enddo
!   Density, spin and bond density
		do i = 1, Ndim
            ii = Latt%o_list(i, 1)
            no = Latt%o_list(i, 2)
            if (Latt%b_list(ii, 2) == 1) sign = 1
            if (Latt%b_list(ii, 2) == 2) sign = -1
            do ns = 1, Nspin
                this%spin_occ(ns, ii) = this%spin_occ(ns, ii) + NsigL_K%phi(ns, ii, ntau) * sign
                this%spin_avg(ns) = this%spin_avg(ns) + NsigL_K%phi(ns, ii, ntau) / dble(Lq)
                this%spin_order(ns) = this%spin_order(ns) + NsigL_K%phi(ns, ii, ntau) * sign / dble(Lq)
            enddo
            this%den_occ(ii) = this%den_occ(ii) + tmpd(i)
            this%density = this%density + tmpd(i) / dble(Lq)
            do nf = 1, Nbond
                i1 = Latt%inv_o_list(Latt%L_bonds(ii, nf), no)
                tmp = real(Grupc(i, i1) + Grupc(i1, i) + dconjg(Grupc(i, i1) + Grupc(i1, i)))
                this%bond_occ(ii, nf) = this%bond_occ(ii, nf) + tmp
                this%el_ke(nf) = this%el_ke(nf) - RT(no, nf) * tmp / dble(Lq)
                if (nf == 1) then
                    iy = Latt%n_list(ii, 2)
                    Ai = -2.d0 * Pi * NB_field * dble(iy)/dble(Lq)
                    phase_p = exp( dcmplx(0.d0, 1.d0) * Ai)
                    phase_m = exp(-dcmplx(0.d0, 1.d0) * Ai)
                    this%diam = this%diam + 2.0 * RT(no, 1) * real(Grupc(i, i1) * phase_p+ Grupc(i1, i) * phase_m) / dble(Lq)
                endif
            enddo
        enddo
        return
    end subroutine Obs_equal_calc
end module ObserEqual_mod