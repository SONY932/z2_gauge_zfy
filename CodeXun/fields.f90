module Fields_mod
    use MyLattice
    implicit none
    
    public :: SpinConf, conf_in, conf_out
    private
    
    type :: SpinConf
        real(kind=8), dimension(:,:,:), public, allocatable :: phi
    contains
        procedure :: make => SpinConf_make
        final :: SpinConf_clear
        procedure :: m_bosonratio_local
!        procedure :: m_bosonratio_global_time
        procedure :: m_bosonratio_global_spacetime
        generic :: bosonratio => m_bosonratio_local
!        generic :: bosonratio => m_bosonratio_global_time
        generic :: bosonratio => m_bosonratio_global_spacetime
    end type SpinConf

contains
    subroutine SpinConf_make(this)
        class(SpinConf), intent(inout) :: this
        allocate(this%phi(Nspin, Lq, Ltrot))
        return
    end subroutine SpinConf_make
    
    subroutine SpinConf_clear(this)
        type(SpinConf), intent(inout) :: this
        deallocate(this%phi)
        return
    end subroutine SpinConf_clear
    
    real(kind=8) function m_bosonratio_local(this, phi_new, ii, ntau, Latt) result(ratio_boson)
! Arguments: 
        class(SpinConf), intent(in) :: this
        real(kind=8), dimension(Nspin, Lq, Ltrot), intent(in) :: phi_new
        integer, intent(in) :: ii, ntau
        class(SquareLattice), intent(in) :: Latt
! Local: 
        real(kind=8) :: new1, old1, new2, old2, new3, old3, new4, old4, action_dif
        real(kind=8) :: vec1(Nspin), vec2(Nspin)
        integer :: jj, nf
! phonon kinetics
        vec1(:) = this%phi(:, ii, npbc(ntau+1, Ltrot)) - phi_new(:, ii, ntau)
        vec2(:) = this%phi(:, ii, npbc(ntau-1, Ltrot)) - phi_new(:, ii, ntau)
        vec1 = vec1 / Dtau
        vec2 = vec2 / Dtau
        new1 = sqr_vec(vec1) + sqr_vec(vec2)
        vec1(:) = this%phi(:, ii, npbc(ntau+1, Ltrot)) - this%phi(:, ii, ntau)
        vec2(:) = this%phi(:, ii, npbc(ntau-1, Ltrot)) - this%phi(:, ii, ntau)
        vec1 = vec1 / Dtau
        vec2 = vec2 / Dtau
        old1 = sqr_vec(vec1) + sqr_vec(vec2)
! phonon-phonon interactions
        new2 = 0.d0; old2 = 0.d0
        do nf = 1, 2*Nbond
            jj = Latt%L_bonds(ii, nf)
            vec1(:) = phi_new(:, ii, ntau) - this%phi(:, jj, ntau)
            new2 = new2 + sqr_vec(vec1)
            vec2(:) = this%phi(:, ii, ntau) - this%phi(:, jj, ntau)
            old2 = old2 + sqr_vec(vec2)
        enddo
! phonon mass term and quartic terms
        vec1(:) = phi_new(:, ii, ntau)
        new3 = sqr_vec(vec1)
        new4 = new3 * new3
        vec2(:) = this%phi(:, ii, ntau)
        old3 = sqr_vec(vec2)
        old4 = old3 * old3
        
        action_dif = (new1 - old1) * Dtau/(2.0*velocity*velocity) + (new2 - old2) * (Dtau/2.0) &
            &   + (new3 - old3) * (Dtau*mass/2.0) + (new4 - old4) * (Dtau*quartic/4.0)
        ratio_boson = exp(-action_dif)
        return
    end function m_bosonratio_local
    
!    real(kind=8) function m_bosonratio_global_time(this, phi_new, ii, Latt) result(ratio_boson)
! Arguments: 
!        class(SpinConf), intent(in) :: this
!        real(kind=8), dimension(Nspin, Lq, Ltrot), intent(in) :: phi_new
!        integer, intent(in) :: ii
!        class(SquareLattice), intent(in) :: Latt
! Local: 
!        real(kind=8) :: new1, old1, new2, old2, new3, old3, new4, old4, action_dif, tmp
!        integer :: nt, jj, nf
!        real(kind=8) :: vec(Nspin)
        
!        new1 = 0.d0; old1 = 0.d0
!        new2 = 0.d0; old2 = 0.d0
!        new3 = 0.d0; old3 = 0.d0
!        new4 = 0.d0; old4 = 0.d0
!        do nt = 1, Ltrot
! phonon kinetics
!            vec(:) = phi_new(:, ii, npbc(nt+1, Ltrot)) - phi_new(:, ii, nt)
!            vec = vec/Dtau
!            new1 = new1 + sqr_vec(vec)
!            vec(:) = this%phi(:, ii, npbc(nt+1, Ltrot)) - this%phi(:, ii, nt)
!            vec = vec/Dtau
!            old1 = old1 + sqr_vec(vec)
! phonon-phonon interactions
!            do nf = 1, 2*Nbond
!                jj = Latt%L_bonds(ii, nf)
!                vec(:) = phi_new(:, ii, nt) - this%phi(:, jj, nt)
!                new2 = new2 + sqr_vec(vec)
!                vec(:) = this%phi(:, ii, nt) - this%phi(:, jj, nt)
!                old2 = old2 + sqr_vec(vec)
!            enddo
! phonon mass term and quartic terms
!            vec(:) = phi_new(:, ii, nt)
!            tmp = sqr_vec(vec)
!            new3 = new3 + tmp
!            new4 = new4 + tmp * tmp
!            vec(:) = this%phi(:, ii, nt)
!            tmp = sqr_vec(vec)
!            old3 = old3 + tmp
!            old4 = old4 + tmp * tmp
!        enddo
!        action_dif = (new1 - old1) * Dtau/(2.0*velocity*velocity) + (new2 - old2) * (Dtau/2.0) &
!            &   + (new3 - old3) * (Dtau*mass/2.0) + (new4 - old4) * (Dtau*quartic/4.0)
!        ratio_boson = exp(-action_dif)
!        return
!    end function m_bosonratio_global_time
    
! without phonon kinetics and mutual phonon interactions
    real(kind=8) function m_bosonratio_global_spacetime(this, phi_new) result(ratio_boson)
! Arguments: 
        class(SpinConf), intent(in) :: this
        real(kind=8), dimension(Nspin, Lq, Ltrot), intent(in) :: phi_new
! Local: 
        real(kind=8) :: new3, old3, new4, old4, action_dif, tmp
        integer :: nt, ii
        real(kind=8) :: vec(Nspin)
        
        new3 = 0.d0; old3 = 0.d0
        new4 = 0.d0; old4 = 0.d0
        do nt = 1, Ltrot
            do ii = 1, Lq
! phonon mass term and quartic terms
                vec(:) = phi_new(:, ii, nt)
                tmp = sqr_vec(vec)
                new3 = new3 + tmp
                new4 = new4 + tmp * tmp
                vec(:) = this%phi(:, ii, nt)
                tmp = sqr_vec(vec)
                old3 = old3 + tmp
                old4 = old4 + tmp * tmp
            enddo
        enddo
        action_dif = (new3 - old3) * (Dtau*mass/2.0) + (new4 - old4) * (Dtau*quartic/4.0)
        ratio_boson = exp(-action_dif)
        return
    end function m_bosonratio_global_spacetime
    
    subroutine conf_in(Conf, iseed, Latt)
        class(SpinConf), intent(inout) :: Conf
        integer, intent(out) :: iseed
        class(SquareLattice), intent(in) :: Latt
! Local: 
        real(kind=8), dimension(Nspin, LqTherm, LtrotTherm) :: phi_therm
        
        call conf_therm_in(phi_therm, iseed)
        call conf_transfer(Conf%phi, phi_therm, Latt)
        return
    end subroutine conf_in
    
    subroutine conf_therm_in(phi_therm, iseed)
!#define DEC
        include 'mpif.h'
! Arguments: 
        real(kind=8), dimension(Nspin, LqTherm, LtrotTherm), intent(inout) :: phi_therm
        integer, intent(out) :: iseed
! Local: 
        integer :: status(MPI_STATUS_SIZE)
        integer :: iseed0, itmp
        integer :: ii, ns, nt, N
        real(kind=8), external :: ranf
        real(kind=8), dimension(Nspin, LqTherm, LtrotTherm) :: phi_itmp
        
        if (IRANK == 0 ) then
            open(unit=30, file='confin.txt', status='unknown')
            open(unit=10, file='seeds.txt', status='unknown')
        endif
        if ( IRANK == 0 ) then
            write(6,*) 'Number of process', ISIZE
            read(30,*) iseed ! read seeds from "confin"
            if (iseed == 0) then
                read(10,*) iseed0 ! if confin=0, read the first row of "seeds" as the seed of the master process.
                do N = 1, ISIZE - 1
!   Setup node I and send data.
                    read(10,*) itmp ! read the following rows of "seeds" as the seeds for other nodes. Each node has only one seed to generate the space-time auxiliary fields
                    call conf_set(phi_itmp, itmp)
                    call MPI_SEND(itmp, 1, MPI_Integer, N, N, MPI_COMM_WORLD, IERR)
					call MPI_SEND(phi_itmp, LqTherm*Nspin*LtrotTherm, MPI_Real8, N, N+1024, MPI_COMM_WORLD, IERR)
                    print '("after sending message from master ", i3.1, " to Rank ", i3.1)', IRANK, N
                enddo
! Set node zero.
                iseed = iseed0
                call conf_set(phi_therm, iseed)
                print '("initiating master rank ", i3.1)', IRANK
            else
! read all confins from NODE 0.
! Setup Node 0
                do nt = 1, LtrotTherm
                    do ii = 1, LqTherm
                        do ns = 1, Nspin
                            read(30,*) phi_therm(ns, ii, nt)
                        enddo
                    enddo
                enddo
                do N = 1, ISIZE - 1
                    read(30,*) itmp
                    do nt = 1, LtrotTherm
                        do ii = 1, LqTherm
                            do ns = 1, Nspin
                                read(30,*) phi_itmp(ns, ii, nt)
                            enddo
                        enddo
                    enddo
                    call MPI_SEND(itmp, 1, MPI_Integer, N, N, MPI_COMM_WORLD, IERR)
					call MPI_SEND(phi_itmp, Nspin*LqTherm*LtrotTherm, MPI_Real8, N, N+1024, MPI_COMM_WORLD, IERR)
                enddo
            endif
        else
            call MPI_RECV(iseed, 1, MPI_Integer, 0, IRANK, MPI_COMM_WORLD, STATUS, IERR)
			call MPI_RECV(phi_therm, Nspin*LqTherm*LtrotTherm, MPI_Real8, 0, IRANK + 1024, MPI_COMM_WORLD, STATUS, IERR)
            print '("Rank ", i3.1, " after receiving message. ")', IRANK
        endif
        if (IRANK == 0 ) then
            close(30); close(10)
        endif
        return
    end subroutine conf_therm_in

    subroutine conf_set(phi, itmp)
        real(kind=8), dimension(Nspin, LqTherm, LtrotTherm), intent(out) :: phi
        integer, intent(inout) :: itmp
! Local: 
        integer :: ii, ns, nt
        real(kind=8) :: X
        real(kind=8), external :: ranf
        
        if (absolute == 1) then
            do nt = 1, LtrotTherm
                do ii = 1, LqTherm
                    do ns = 1, Nspin
                        X = ranf(itmp)
                        phi(ns, ii, nt) = init(ns) + inconf * (X - 0.5)
                    enddo
                enddo
            enddo
        elseif (absolute == 2) then
            do nt = 1, LtrotTherm
                do ii = 1, LqTherm
                    do ns = 1, Nspin
                        X = rng_Gaussian(itmp)
                        phi(ns, ii, nt) = init(ns) + inconf * X
                    enddo
                enddo
            enddo
        elseif (absolute == 3) then
            do nt = 1, LtrotTherm
                do ii = 1, LqTherm
                    do ns = 1, Nspin
                        X = rng_Gaussian(itmp)
                        phi(ns, ii, nt) = init(ns) + inconf * abs(X)
                    enddo
                enddo
            enddo
        else
            write(6,*) "incorrect absolute input", absolute, " on rank ", IRANK
        endif
        return
    end subroutine conf_set
    
    subroutine conf_transfer(phi, phi_therm, Latt)
        real(kind=8), dimension(Nspin, Lq, Ltrot), intent(inout) :: phi
        real(kind=8), dimension(Nspin, LqTherm, LtrotTherm), intent(in) :: phi_therm
        class(SquareLattice), intent(in) :: Latt
        integer :: nt, ntt, nx, ny, ii, iit, nc
        integer :: n_list_therm(LqTherm, 1:2), inv_n_list_therm(NlxTherm, NlyTherm)

        nc = 0
        do ny = 1, NlyTherm
            do nx = 1, NlxTherm
                nc = nc + 1
                n_list_therm(nc, 1) = nx
                n_list_therm(nc, 2) = ny
                inv_n_list_therm(nx, ny) = nc
            enddo
        enddo
! copy NsigL_K
        do nt = 1, LtrotTherm
            do ny = 1, NlyTherm
                do nx = 1, NlxTherm
                    iit = inv_n_list_therm(nx, ny)
                    ii = Latt%inv_n_list(nx, ny)
                    phi(:, ii, nt) = phi_therm(:, iit, nt)
                enddo
            enddo
        enddo
        do nt = 1, LtrotTherm
            do ny = 1, NlyTherm
                do nx = NlxTherm+1, Nlx
                    iit = inv_n_list_therm(nx - NlxTherm, ny)
                    ii = Latt%inv_n_list(nx, ny)
                    phi(:, ii, nt) = phi_therm(:, iit, nt)
                enddo
            enddo
        enddo
        do nt = 1, LtrotTherm
            do ny = NlyTherm+1, Nly
                do nx = 1, NlxTherm
                    iit = inv_n_list_therm(nx, ny - NlyTherm)
                    ii = Latt%inv_n_list(nx, ny)
                    phi(:, ii, nt) = phi_therm(:, iit, nt)
                enddo
            enddo
        enddo
        do nt = 1, LtrotTherm
            do ny = NlyTherm+1, Nly
                do nx = NlxTherm+1, Nlx
                    iit = inv_n_list_therm(nx - NlxTherm, ny - NlyTherm)
                    ii = Latt%inv_n_list(nx, ny)
                    phi(:, ii, nt) = phi_therm(:, iit, nt)
                enddo
            enddo
        enddo
        do nt = LtrotTherm+1, Ltrot
            do ny = 1, Nly
                do nx = 1, Nlx
                    ii = Latt%inv_n_list(nx, ny)
                    ntt = nt - LtrotTherm
                    phi(:, ii, nt) = phi(:, ii, ntt)
                enddo
            enddo
        enddo
        return
    end subroutine conf_transfer

    subroutine conf_out(Conf, iseed)
!#define DEC
        include 'mpif.h'
! Arguments: 
        class(SpinConf), intent(in) :: Conf
        integer, intent(in) :: iseed
! Local:
        integer :: ii, ns, nt, N
        integer :: status(MPI_STATUS_SIZE)

        if (IRANK == 0) open(unit=35, file='confout.txt', status='unknown')
        if (IRANK .NE. 0) then
            call MPI_SEND(iseed, 1, MPI_Integer, 0, IRANK, MPI_COMM_WORLD, IERR)
			call MPI_SEND(Conf%phi, Lq*Nspin*Ltrot, MPI_Real8, 0, IRANK+1024, MPI_COMM_WORLD, IERR)
            print '("Rank", i3.1, " after sending message. ")', IRANK
        endif
        if (IRANK == 0) then
            write(35,*) iseed
            do nt = 1, Ltrot
                do ii = 1, Lq
                    do ns = 1, Nspin
                        write(35,*) Conf%phi(ns, ii, nt)
                    enddo
                enddo
            enddo
            print '("Rank ", i3.1, " after writing down own auxfields. ")', IRANK
            do N = 1, ISIZE - 1
                call MPI_RECV(iseed, 1, MPI_Integer, N, N, MPI_COMM_WORLD, STATUS, IERR)
				call MPI_RECV(Conf%phi, Lq*Nspin*Ltrot, MPI_Real8, N, N+1024, MPI_COMM_WORLD, STATUS, IERR)
                write(35,*) iseed
                do nt = 1, Ltrot
                    do ii = 1, Lq
                        do ns = 1, Nspin
                            write(35,*) Conf%phi(ns, ii, nt)
                        enddo
                    enddo
                enddo
                print '("Master rank ", i3.1, " after writing down auxfields on RANK ",i3.1 ,".")', IRANK, N
            enddo
        endif
        if (IRANK == 0) close(35)
        return
    end subroutine conf_out
    
    real(kind=8) function rng_Gaussian(iseed) result(X)
        integer, intent(inout) :: iseed
        real(kind=8), external :: ranf
        real(kind=8) :: X1, X2
        X1 = ranf(iseed)
        X2 = ranf(iseed)
        X = sqrt(-2.0*log(X1))*cos(2.0*Pi*X2) ! Gaussian distribution
        return
    end function rng_Gaussian
end module Fields_mod