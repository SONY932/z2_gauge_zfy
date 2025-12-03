module OperatorKin_mod
    use MyLattice
    use calc_basic
    implicit none
    public

    type :: OperatorKin
        real(kind = 8), private :: alpha, alpha_2
        real(kind = 8), private :: entryC
        real(kind = 8), private :: entryS

        real(kind = 8), public :: Delta(2, 2)
    contains
        procedure :: set => opK_set
        procedure, private :: get_exp => opK_get_exp
        procedure :: get_delta => opK_get_delta
        procedure :: mmult_R => opK_mmult_R
        procedure :: mmult_L => opK_mmult_L
    end type OperatorKin

contains
    subroutine opK_set(this)
        class(OperatorKin), intent(inout) :: this
        this%alpha = 0.5d0 * Dtau * RT
        this%alpha_2 = Dtau * RT
        return
    end subroutine opK_set

    subroutine opK_get_exp(this, sigma, nflag)
        class(OperatorKin), intent(inout) :: this
        integer, intent(in) :: nflag
        real(kind = 8), intent(in) :: sigma
        this%entryC = cosh( nflag * this%alpha * sigma )
        this%entryS = sinh( nflag * this%alpha * sigma )
        return
    end subroutine opK_get_exp

    subroutine opK_get_delta(this, sigma_old, sigma_new)
! Calculate Delta = exp(alpha*sigma_new) * exp(-alpha*sigma_old) - I
! Used for efficient Green's function update via Sherman-Morrison formula
        class(OperatorKin), intent(inout) :: this
        real(kind = 8), intent(in) :: sigma_old, sigma_new
        real(kind = 8) :: C_new, C_old, S_new, S_old

        C_new = cosh( this%alpha_2 * sigma_new )
        C_old = cosh( this%alpha_2 * sigma_old )
        S_new = sinh( this%alpha_2 * sigma_new )
        S_old = sinh( this%alpha_2 * sigma_old )
        
        ! Delta = exp(+alpha*sigma_new) * exp(-alpha*sigma_old) - I
        ! Note: sinh is odd function, so exp(-alpha*sigma_old) has -S_old entries
        this%Delta(1, 1) = C_new * C_old - S_new * S_old - 1.d0
        this%Delta(1, 2) = S_new * C_old - C_new * S_old
        this%Delta(2, 1) = S_new * C_old - C_new * S_old
        this%Delta(2, 2) = C_new * C_old - S_new * S_old - 1.d0
        return
    end subroutine opK_get_delta

    subroutine opK_mmult_R(this, Mat, Latt, sigma, ntau, nflag)
! In Mat Out U(NF) * EXP(D(NF)) * Mat
! Symmetrized checkerboard decomposition: forward (1->2->3->4) then backward (4->3->2->1)
! Arguments: 
        class(OperatorKin), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(SquareLattice), intent(in) :: Latt
        real(kind = 8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: ntau, nflag
! Local:
        integer :: ig

        ! Forward sweep: group 1 -> 2 -> 3 -> 4
        call process_group(this, Mat, Latt%group_1, Latt, sigma, ntau, nflag)
        call process_group(this, Mat, Latt%group_2, Latt, sigma, ntau, nflag)
        call process_group(this, Mat, Latt%group_3, Latt, sigma, ntau, nflag)
        call process_group(this, Mat, Latt%group_4, Latt, sigma, ntau, nflag)

        ! Backward sweep: group 4 -> 3 -> 2 -> 1 (for symmetrization)
        call process_group(this, Mat, Latt%group_4, Latt, sigma, ntau, nflag)
        call process_group(this, Mat, Latt%group_3, Latt, sigma, ntau, nflag)
        call process_group(this, Mat, Latt%group_2, Latt, sigma, ntau, nflag)
        call process_group(this, Mat, Latt%group_1, Latt, sigma, ntau, nflag)

        return
    end subroutine opK_mmult_R

    subroutine process_group(this, Mat, group, Latt, sigma, ntau, nflag)
! Helper subroutine to process one checkerboard group
        class(OperatorKin), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        integer, dimension(:), intent(in) :: group
        class(SquareLattice), intent(in) :: Latt
        real(kind = 8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: ntau, nflag
! Local:
        integer :: P(2), ii, jj, len, no
        real(kind = 8) :: vec
        complex(kind=8), dimension(2, Ndim) :: Vhlp

        len = size(group)
        do ii = 1, len
            no = group(ii)
            vec = sigma(no, ntau)
            call this%get_exp(vec, nflag)  ! output entryC and entryS
            P(1) = Latt%bond_list(no, 1)
            P(2) = Latt%bond_list(no, 2)
            
            ! Apply 2x2 Givens rotation: [C S; S C] to rows P(1) and P(2)
            do jj = 1, Ndim
                Vhlp(1, jj) = this%entryC * Mat(P(1), jj) + this%entryS * Mat(P(2), jj)
                Vhlp(2, jj) = this%entryS * Mat(P(1), jj) + this%entryC * Mat(P(2), jj)
            enddo
            do jj = 1, Ndim
                Mat(P(1), jj) = Vhlp(1, jj)
                Mat(P(2), jj) = Vhlp(2, jj)
            enddo
        enddo
        return
    end subroutine process_group

    subroutine opK_mmult_L(this, Mat, Latt, sigma, ntau, nflag)
! In Mat Out Mat * EXP(D(NF)) * UT(NF)
! Symmetrized checkerboard decomposition: forward (1->2->3->4) then backward (4->3->2->1)
! Arguments: 
        class(OperatorKin), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(SquareLattice), intent(in) :: Latt
        real(kind = 8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: ntau, nflag

        ! Forward sweep: group 1 -> 2 -> 3 -> 4
        call process_group_L(this, Mat, Latt%group_1, Latt, sigma, ntau, nflag)
        call process_group_L(this, Mat, Latt%group_2, Latt, sigma, ntau, nflag)
        call process_group_L(this, Mat, Latt%group_3, Latt, sigma, ntau, nflag)
        call process_group_L(this, Mat, Latt%group_4, Latt, sigma, ntau, nflag)

        ! Backward sweep: group 4 -> 3 -> 2 -> 1 (for symmetrization)
        call process_group_L(this, Mat, Latt%group_4, Latt, sigma, ntau, nflag)
        call process_group_L(this, Mat, Latt%group_3, Latt, sigma, ntau, nflag)
        call process_group_L(this, Mat, Latt%group_2, Latt, sigma, ntau, nflag)
        call process_group_L(this, Mat, Latt%group_1, Latt, sigma, ntau, nflag)

        return
    end subroutine opK_mmult_L

    subroutine process_group_L(this, Mat, group, Latt, sigma, ntau, nflag)
! Helper subroutine to process one checkerboard group (left multiplication)
        class(OperatorKin), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        integer, dimension(:), intent(in) :: group
        class(SquareLattice), intent(in) :: Latt
        real(kind = 8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: ntau, nflag
! Local:
        integer :: P(2), ii, jj, len, no
        real(kind = 8) :: vec
        complex(kind=8), dimension(Ndim, 2) :: Uhlp

        len = size(group)
        do ii = 1, len
            no = group(ii)
            vec = sigma(no, ntau)
            call this%get_exp(vec, nflag)  ! output entryC and entryS
            P(1) = Latt%bond_list(no, 1)
            P(2) = Latt%bond_list(no, 2)
            
            ! Apply 2x2 Givens rotation: [C S; S C] to columns P(1) and P(2)
            do jj = 1, Ndim
                Uhlp(jj, 1) = Mat(jj, P(1)) * this%entryC + Mat(jj, P(2)) * this%entryS
                Uhlp(jj, 2) = Mat(jj, P(1)) * this%entryS + Mat(jj, P(2)) * this%entryC
            enddo
            do jj = 1, Ndim
                Mat(jj, P(1)) = Uhlp(jj, 1)
                Mat(jj, P(2)) = Uhlp(jj, 2)
            enddo
        enddo
        return
    end subroutine process_group_L

end module OperatorKin_mod
