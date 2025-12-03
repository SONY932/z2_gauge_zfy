module FourierTrans_mod
    use DQMC_Model_mod
    use ObserEqual_mod
    use ObserTau_mod
    implicit none
    
    type, public :: FourierTrans
    contains
        procedure, private, nopass :: m_write_real_1
        procedure, private, nopass :: m_write_real_2
        generic :: write_real => m_write_real_1
        generic :: write_real => m_write_real_2
        
        procedure, private, nopass :: m_write_reciprocal_1
        procedure, private, nopass :: m_write_reciprocal_2
        procedure, private, nopass :: m_write_reciprocal_3
        generic :: write_reciprocal => m_write_reciprocal_1
        generic :: write_reciprocal => m_write_reciprocal_2
        generic :: write_reciprocal => m_write_reciprocal_3
        
        procedure, private, nopass :: m_write_k_1
        procedure, private, nopass :: m_write_k_2
        procedure, private, nopass :: m_write_k_3
        generic :: write_k => m_write_k_1
        generic :: write_k => m_write_k_2
        generic :: write_k => m_write_k_3
        
        procedure, private, nopass :: integrate_susc => m_integrate_susc_2_mom
        procedure, private, nopass :: integrate_susc_freq => m_integrate_susc_2_freq
        
        procedure, private, nopass :: m_write_k_tau_2
        procedure, private, nopass :: m_write_k_tau_4
        generic :: write_k_tau => m_write_k_tau_2
        generic :: write_k_tau => m_write_k_tau_4
        
        procedure, private, nopass :: write_r_tau => m_write_r_tau_2
        
        procedure, private, nopass :: write_w => m_write_w
        
        procedure, private :: write_obs_equal => m_write_obs_equal
        procedure, private :: write_obs_tau => m_write_obs_tau
        
        procedure, public :: preq => m_process_obs_equal
        procedure, public :: prtau => m_process_obs_tau
    end type FourierTrans

contains
    subroutine m_write_real_1(gr, filek) ! overloading routine for other correlations
! Arguments:
        real(kind=8), dimension(Lq), intent(in) :: gr
        character(len=*), intent(in) :: filek
! Local: 
        integer :: nr
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nr = 1, Lq
            write(20,*) Latt%aimj_v(nr, 1), Latt%aimj_v(nr, 2)
            write(20,*) gr(nr)
        enddo
        close(20)
        return
    end subroutine m_write_real_1
    
    subroutine m_write_real_2(gr, filek) ! overloading routine
! Arguments:
        real(kind=8), dimension(:,:), intent(in) :: gr
        character(len=*), intent(in) :: filek
! Local: 
        integer :: nr, nf
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nr = 1, Lq
            write(20,*) Latt%aimj_v(nr, 1), Latt%aimj_v(nr, 2)
            if (size(gr, 1) == Lq) then ! (Lq, Nbond)
                do nf = 1, Nbond
                    write(20,*) gr(nr, nf)
                enddo
            elseif (size(gr, 2) == Lq) then ! (Nspin, Lq)
                do nf = 1, Nspin
                    write(20,*) gr(nf, nr)
                enddo
            else
                write(6,*) "ERROR: incorrect input size in write_real_2"; stop
            endif
        enddo
        close(20)
        return
    end subroutine m_write_real_2
   
    subroutine m_write_reciprocal_1(gk, filek)
        complex(kind=8), dimension(Lq), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer :: nk
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nk = 1, Lq
            write(20,*) Latt%xk_v(nk, 1), Latt%xk_v(nk, 2)
            write(20,*) gk(nk)
        enddo
        close(20)
        return
    end subroutine m_write_reciprocal_1
    
    subroutine m_write_reciprocal_2(gk, filek)
        complex(kind=8), dimension(Lq, Nbond), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer :: nk, no
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nk = 1, Lq
            write(20,*) Latt%xk_v(nk, 1), Latt%xk_v(nk, 2)
            do no = 1, Nbond
                write(20,*) gk(nk, no)
            enddo
        enddo
        close(20)
        return
    end subroutine m_write_reciprocal_2
    
    subroutine m_write_reciprocal_3(gk, filek)
        complex(kind=8), dimension(Lq, Nbond, Nbond), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer :: nk, no1, no2
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nk = 1, Lq
            write(20,*) Latt%xk_v(nk, 1), Latt%xk_v(nk, 2)
            do no2 = 1, Nbond
                do no1 = 1, Nbond
                    write(20,*) gk(nk, no1, no2)
                enddo
            enddo
        enddo
        close(20)
        return
    end subroutine m_write_reciprocal_3

    subroutine m_write_k_1(gk, filek, momindex)
        complex(kind=8), dimension(Lq), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex
        open(unit=30, file=filek, status='unknown', action="write", position="append")
        write(30, *) real(gk(momindex)), imag(gk(momindex))
        close(30)
        return
    end subroutine m_write_k_1
    
    subroutine m_write_k_2(gk, filek, momindex, nf)
        complex(kind=8), dimension(Lq, Nbond), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex, nf
        open(unit=30, file=filek, status='unknown', action="write", position="append")
        write(30,*) real(gk(momindex, nf)), imag(gk(momindex, nf))
        close(30)
        return
    end subroutine m_write_k_2
    
    subroutine m_write_k_3(gk, filek, momindex, no1, no2)
        complex(kind=8), dimension(Lq, Nbond, Nbond), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex, no1, no2
        open(unit=30, file=filek, status='unknown', action="write", position="append")
        write(30,*) real(gk(momindex, no1, no2)), imag(gk(momindex, no1, no2))
        close(30)
        return
    end subroutine m_write_k_3

    subroutine m_integrate_susc_2_mom(gr, gk)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gr
        complex(kind=8), dimension(Lq), intent(out) :: gk
        integer :: nt, nk
        gk = dcmplx(0.d0, 0.d0)
        do nt = 1, Ltrot
            do nk = 1, Lq
                gk(nk) = gk(nk) + gr(nk, nt)
            enddo
        enddo
        gk = gk * dcmplx(Dtau, 0.d0)
        return
    end subroutine m_integrate_susc_2_mom
    
    subroutine m_integrate_susc_2_freq(gr, gk)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gr
        complex(kind=8), dimension(Ltrot), intent(out) :: gk
        integer :: nt, nw
        gk = dcmplx(0.d0, 0.d0)
        do nw = 1, Ltrot
            do nt = 1, Ltrot
                gk(nw) = gk(nw) + exp( dcmplx(0.d0, 2.d0*Pi*dble((nt-1)*(nw-1))/dble(Ltrot)) ) * gr(1, nt)
            enddo
        enddo
        gk = gk * dcmplx(Dtau, 0.d0)
        return
    end subroutine m_integrate_susc_2_freq
    
    subroutine m_write_r_tau_2(gr, filek, momindex)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gr
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex
        integer :: nt
        complex(kind=8) :: tmp
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        write(20,*) Latt%aimj_v(momindex, 1), Latt%aimj_v(momindex, 2)
        do nt = 1, Ltrot
            tmp = gr(momindex, nt) / dble(Lq)
            write(20,*) real(tmp), imag(tmp)
        enddo
        close(20)
        return
    end subroutine m_write_r_tau_2
    
    subroutine m_write_w(gk, filek)
        complex(kind=8), dimension(Ltrot), intent(in) :: gk
        character(len=*), intent(in) :: filek
        real(kind=8) :: omega
        integer :: nw
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nw = 1, Ltrot
            omega = 2.d0 * Pi * dble(nw-1) / beta ! Matsubara frequency
            write(20,*) omega
            write(20,*) real(gk(nw)), imag(gk(nw))
        enddo
        close(20)
        return
    end subroutine m_write_w
    
    subroutine m_write_k_tau_2(gk, filek, momindex)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex
        integer :: nt
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        write(20,*) Latt%xk_v(momindex, 1), Latt%xk_v(momindex, 2)
        do nt = 1, Ltrot
            write(20,*) real(gk(momindex, nt)), imag(gk(momindex, nt))
        enddo
        close(20)
        return
    end subroutine m_write_k_tau_2
    
    subroutine m_write_k_tau_4(gk, filek, momindex, no1, no2)
        complex(kind=8), dimension(Lq, Nbond, Nbond, Ltrot), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex, no1, no2
        integer :: nt
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        write(20,*) Latt%xk_v(momindex, 1), Latt%xk_v(momindex, 2)
        do nt = 1, Ltrot
            write(20,*) real(gk(momindex, no1, no2, nt))
        enddo
        close(20)
        return
    end subroutine m_write_k_tau_4
    
    subroutine m_write_obs_equal(this, Obs)
        class(FourierTrans), intent(inout) :: this
        class(ObserEqual), intent(in) :: Obs
!        complex(kind=8) :: gk1(Lq), gk2(Lq, Nbond), gk3(Lq, Nbond, Nbond)
        character(len=25) :: filek
        
        filek = "den_occ"
        call this%write_real(Obs%den_occ, filek)
        
        filek = "spin_occ"
        call this%write_real(Obs%spin_occ, filek)
            
        filek = "bond_occ"
        call this%write_real(Obs%bond_occ, filek)
            
        open(unit=80, file='filling', status='unknown', action="write", position="append")
        write(80,*) Obs%density
        close(80)
        
        open(unit=80, file='diamagnetic', status='unknown', action="write", position="append")
        write(80,*) Obs%diam
        close(80)
            
        open(unit=80, file='kinetic', status='unknown', action="write", position="append")
        write(80,*) Obs%el_ke(1),  Obs%el_ke(2)
        close(80)
        
        open(unit=80, file='spin_average', status='unknown', action="write", position="append")
        write(80,*) Obs%spin_avg(1),  Obs%spin_avg(2)
        close(80)
        
        open(unit=80, file='spin_order', status='unknown', action="write", position="append")
        write(80,*) Obs%spin_order(1),  Obs%spin_order(2)
        close(80)
        return
    end subroutine m_write_obs_equal

    subroutine m_process_obs_equal(this, Obs)
!#define DEC
        include 'mpif.h'
! Arguments:
        class(FourierTrans), intent(inout) :: this
        class(ObserEqual), intent(inout) :: Obs
! Local:
!        complex(kind=8), dimension(Lq, Nbond, Nbond) :: Collect3
        real(kind=8), dimension(Lq) :: Collect1
        real(kind=8), dimension(Lq, Nbond) :: Collect2
        real(kind=8), dimension(Nspin, Lq) :: Collect2prime
        real(kind=8) :: Collect0, Collect1prime(Nbond)
        integer :: N
! test:
!        real(kind=8) :: start, finish
! collect data from different nodes
!        if (IRANK == 0 .or. IRANK == 1) call CPU_TIME(start)
        N = Lq
        Collect1 = 0.d0
        call MPI_REDUCE(Obs%den_occ, Collect1, N, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%den_occ = Collect1/dble(ISIZE)
        
        N = Nspin * Lq
        Collect2prime = 0.d0
        call MPI_REDUCE(Obs%spin_occ, Collect2prime, N, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%spin_occ = Collect2prime/dble(ISIZE)
        
        Collect1prime = 0.d0
        call MPI_REDUCE(Obs%spin_avg, Collect1prime, Nspin, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%spin_avg = Collect1prime/dble(ISIZE)
        
        Collect1prime = 0.d0
        call MPI_REDUCE(Obs%spin_order, Collect1prime, Nspin, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%spin_order = Collect1prime/dble(ISIZE)
        
        Collect0 = 0.d0
        call MPI_REDUCE(Obs%density, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%density = Collect0/dble(ISIZE)
        
        Collect0 = 0.d0
        call MPI_REDUCE(Obs%diam, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%diam = Collect0/dble(ISIZE)
        
        Collect1prime = 0.d0
        call MPI_REDUCE(Obs%el_ke, Collect1prime, Nbond, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%el_ke = Collect1prime/dble(ISIZE)

        N = Lq * Nbond
        Collect2 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%bond_occ, Collect2, N, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%bond_occ = Collect2/dcmplx(dble(ISIZE),0.d0)

!        if (IRANK == 0 .or. IRANK == 1) then
!            call CPU_TIME(finish)
!            print '("Time MPI_Reduce on IRANK ", i3.1, " = ", f10.5, " seconds.")', IRANK, finish - start
!        endif
        if (IRANK == 0) then
!            call CPU_TIME(start)
            call this%write_obs_equal(Obs)
!            call CPU_TIME(finish)
!            print '("Time fourier_trans on IRANK ", i3.1, " = ", f10.5, " seconds.")', IRANK, finish - start
        endif
        return
    end subroutine m_process_obs_equal
    
    subroutine m_write_obs_tau(this, Obs)
        class(FourierTrans), intent(inout) :: this
        class(ObserTau), intent(in) :: Obs
        complex(kind=8) :: gk1(Lq), gk1w(Ltrot), gk2(Lq, Ltrot)
        character(len=25) :: filek
        integer :: indexzero, indexpi, indexqx, indexqy
        
        indexpi = int((Lq + Nlx) / 2) + 1
        indexzero = 1
        indexqx = 2
        indexqy = Nlx + 1
        
        call Fourier_R_to_K(Obs%spin_tau, gk2, Latt)
        call this%integrate_susc(gk2, gk1)
        filek = "sdw_susc"
        call this%write_reciprocal(gk1, filek)
        filek = "sdw_pipi_susc"
        call this%write_k(gk1, filek, indexpi)
        filek = "sdw_zero_susc"
        call this%write_k(gk1, filek, indexzero)
        
        call Fourier_R_to_K(Obs%charge_d_tau, gk2, Latt)
        filek = "cdw_d_pipi_tau"
        call this%write_k_tau(gk2, filek, indexpi)
        call this%integrate_susc(gk2, gk1)
        filek = "cdw_d_susc"
        call this%write_reciprocal(gk1, filek)
        filek = "cdw_d_pipi_susc"
        call this%write_k(gk1, filek, indexpi)
        filek = "cdw_d_zero_susc"
        call this%write_k(gk1, filek, indexzero)
        
        call Fourier_R_to_K(Obs%charge_s_tau, gk2, Latt)
        filek = "cdw_s_pipi_tau"
        call this%write_k_tau(gk2, filek, indexpi)
        call this%integrate_susc(gk2, gk1)
        filek = "cdw_s_susc"
        call this%write_reciprocal(gk1, filek)
        filek = "cdw_s_pipi_susc"
        call this%write_k(gk1, filek, indexpi)
        filek = "cdw_s_zero_susc"
        call this%write_k(gk1, filek, indexzero)
        
        filek = "g_r0_tau"
        call this%write_r_tau(Obs%single_tau, filek, Lq)
        call Fourier_R_to_K(Obs%single_tau, gk2, Latt)
        call this%integrate_susc(gk2, gk1)
        filek = "g_susc"
        call this%write_reciprocal(gk1, filek)
        filek = "g_pipi_susc"
        call this%write_k(gk1, filek, indexpi)
        filek = "g_zero_susc"
        call this%write_k(gk1, filek, indexzero)
        
        call Fourier_R_to_K(Obs%curxx_tau, gk2, Latt)
        filek = "curxx_L_tau"
        call this%write_k_tau(gk2, filek, indexqx)
        filek = "curxx_T_tau"
        call this%write_k_tau(gk2, filek, indexqy)
        filek = "curxx_zero_tau"
        call this%write_k_tau(gk2, filek, indexzero)
        call this%integrate_susc(gk2, gk1)
        call this%integrate_susc_freq(gk2, gk1w)
        filek = "curxx_freq"
        call this%write_w(gk1w, filek)
        filek = "curxx_susc"
        call this%write_reciprocal(gk1, filek)
        filek = "curxx_L_susc"
        call this%write_k(gk1, filek, indexqx)
        filek = "curxx_T_susc"
        call this%write_k(gk1, filek, indexqy)
        filek = "curxx_zero_susc"
        call this%write_k(gk1, filek, indexzero)
        
        call Fourier_R_to_K(Obs%pair_s_tau, gk2, Latt)
        filek = "pair_s_zero_tau"
        call this%write_k_tau(gk2, filek, indexzero)
        call this%integrate_susc(gk2, gk1)
        filek = "pair_s_susc"
        call this%write_reciprocal(gk1, filek)
        filek = "pair_s_zero_susc"
        call this%write_k(gk1, filek, indexzero)
        
        call Fourier_R_to_K(Obs%pair_d_tau, gk2, Latt)
        filek = "pair_d_zero_tau"
        call this%write_k_tau(gk2, filek, indexzero)
        call this%integrate_susc(gk2, gk1)
        filek = "pair_d_susc"
        call this%write_reciprocal(gk1, filek)
        filek = "pair_d_zero_susc"
        call this%write_k(gk1, filek, indexzero)
        return
    end subroutine m_write_obs_tau
    
    subroutine m_process_obs_tau(this, Obs)
!#define DEC
        include 'mpif.h'
! Arguments:
        class(FourierTrans), intent(inout) :: this
        class(ObserTau), intent(inout) :: Obs
! Local:
        complex(kind=8), dimension(Lq, Ltrot) :: Collect2
!        complex(kind=8), dimension(Lq, Nbond, Nbond, Ltrot) :: Collect4
        integer :: N

        N = Lq * Ltrot
        Collect2 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%spin_tau, Collect2, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%spin_tau = Collect2/dcmplx(dble(ISIZE),0.d0)
        
        Collect2 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%charge_d_tau, Collect2, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%charge_d_tau = Collect2/dcmplx(dble(ISIZE),0.d0)
        
        Collect2 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%charge_s_tau, Collect2, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%charge_s_tau = Collect2/dcmplx(dble(ISIZE),0.d0)
        
        Collect2 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%single_tau, Collect2, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%single_tau = Collect2/dcmplx(dble(ISIZE),0.d0)
        
        Collect2 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%pair_s_tau, Collect2, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%pair_s_tau = Collect2/dcmplx(dble(ISIZE),0.d0)
        
        Collect2 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%pair_d_tau, Collect2, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%pair_d_tau = Collect2/dcmplx(dble(ISIZE),0.d0)
        
        Collect2 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%curxx_tau, Collect2, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%curxx_tau = Collect2/dcmplx(dble(ISIZE),0.d0)

        if (IRANK == 0) call this%write_obs_tau(Obs)
        return
    end subroutine m_process_obs_tau
end module FourierTrans_mod