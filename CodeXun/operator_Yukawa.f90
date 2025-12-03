module OperatorK_mod
    use MyLattice
    implicit none
    public
    
    type :: OperatorPhonon
        real(kind=8), private :: alpha ! = Dtau*phg
        real(kind=8), private :: entryC
        complex(kind=8), private :: entryS
        
        complex(kind=8), public :: Delta(2,2)
    contains
        procedure :: set => opK_set
        procedure, private :: get_exp => opK_get_exp
        procedure :: get_delta => opK_get_delta
        procedure :: mmult_R => opK_mmult_R
        procedure :: mmult_L => opK_mmult_L
    end type OperatorPhonon
    
contains
    subroutine opK_set(this)
        class(OperatorPhonon), intent(inout) :: this
        this%alpha = Dtau * lambda
        return
    end subroutine opK_set
    
    subroutine opK_get_exp(this, vec, nflag)
        class(OperatorPhonon), intent(inout) :: this
        integer, intent(in) :: nflag ! +1 or -1; propagating direction
        real(kind=8), dimension(Nspin), intent(in) :: vec ! two-component order parameter with space-time coordinate (ii, ntau)
        real(kind=8) :: magnitude
        magnitude = sqrt(sqr_vec(vec))
        this%entryC = cosh( nflag * this%alpha * magnitude )
        this%entryS = sinh( nflag * this%alpha * magnitude ) * dcmplx(-vec(1), vec(2)) / magnitude
        return
    end subroutine opK_get_exp
    
    subroutine opK_get_delta(this, vec_old, vec_new, sign)
        class(OperatorPhonon), intent(inout) :: this
        real(kind=8), dimension(Nspin), intent(in) :: vec_old, vec_new
        integer, intent(in) :: sign ! = \pm 1
        real(kind=8) :: C_new, C_old, mag_new, mag_old
        complex(kind=8) :: S_new, S_old
        
        mag_new = sqrt(sqr_vec(vec_new))
        C_new = cosh( this%alpha * mag_new )
        S_new = sinh( this%alpha * mag_new ) * dcmplx(- vec_new(1), vec_new(2)) / mag_new
        mag_old = sqrt(sqr_vec(vec_old))
        C_old = cosh( this%alpha * mag_old )
        S_old = sinh( -this%alpha * mag_old ) * dcmplx(- vec_old(1), vec_old(2)) / mag_old
        this%Delta(1,1) = C_new * C_old + S_new * dconjg(S_old) - 1.d0
        this%Delta(1,2) = (C_new * S_old + S_new * C_old) * dble(sign)
        this%Delta(2,1) = (C_new * dconjg(S_old) + dconjg(S_new) * C_old) * dble(sign)
        this%Delta(2,2) = C_new * C_old + dconjg(S_new) * S_old - 1.d0
        return
    end subroutine opK_get_delta
    
    subroutine opK_mmult_R(this, Mat, Latt, phi, ntau, nflag)
! In Mat Out U(NF) * EXP(D(NF)) * Mat
! Arguments: 
        class(OperatorPhonon), intent(inout) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(SquareLattice), intent(in) :: Latt
        real(kind=8), dimension(Nspin, Lq, Ltrot), intent(in) :: phi
        integer, intent(in) :: ntau, nflag
! Local: 
        integer :: P(Norb), ii, no, j, sign
        real(kind=8), dimension(Nspin) :: vec
        complex(kind=8), dimension(2, Ndim) :: Vhlp

        do ii = 1, Lq
            vec(:) = phi(:, ii, ntau)
            call this%get_exp(vec, nflag) ! output entryC and entryS
            if (Latt%b_list(ii, 2) == 1) sign = 1
            if (Latt%b_list(ii, 2) == 2) sign = -1
            do no = 1, Norb
                P(no) = Latt%inv_o_list(ii, no)
            enddo
            Vhlp = dcmplx(0.d0, 0.d0)
            do j = 1, Ndim
                Vhlp(1, j) = this%entryC * Mat(P(1), j) + sign * this%entryS * Mat(P(2), j)
                Vhlp(2, j) = sign * dconjg(this%entryS) * Mat(P(1), j) + this%entryC * Mat(P(2), j)
            enddo
            do j = 1, Ndim
                Mat(P(1), j) = Vhlp(1, j)
                Mat(P(2), j) = Vhlp(2, j)
            enddo
        enddo
        return
    end subroutine opK_mmult_R
    
    subroutine opK_mmult_L(this, Mat, Latt, phi, ntau, nflag) 
! In Mat Out Mat* EXP(D(NF)) * UT(NF)
! Arguments:
        class(OperatorPhonon), intent(inout) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(SquareLattice), intent(in) :: Latt
        real(kind=8), dimension(Nspin, Lq, Ltrot), intent(in) :: phi
        integer, intent(in) :: ntau, nflag
! Local: 
        integer :: P(Norb), j, ii, no, sign
        real(kind=8), dimension(Nspin) :: vec
        complex(kind=8), dimension(Ndim, 2) :: Uhlp

        do ii = 1, Lq
            vec(:) = phi(:, ii, ntau)
            call this%get_exp(vec, nflag) ! output entryC and entryS
            if (Latt%b_list(ii, 2) == 1) sign = 1
            if (Latt%b_list(ii, 2) == 2) sign = -1
            do no = 1, Norb
                P(no) = Latt%inv_o_list(ii, no)
            enddo
            Uhlp = dcmplx(0.d0, 0.d0)
            do j = 1, Ndim
                Uhlp(j, 1) = Mat(j, P(1)) * this%entryC + sign * Mat(j, P(2)) * dconjg(this%entryS)
                Uhlp(j, 2) = sign * Mat(j, P(1)) * this%entryS + Mat(j, P(2)) * this%entryC
            enddo
            do j = 1, Ndim
                Mat(j, P(1)) = Uhlp(j, 1)
                Mat(j, P(2)) = Uhlp(j, 2)
            enddo
        enddo
        return
    end subroutine opK_mmult_L
end module OperatorK_mod