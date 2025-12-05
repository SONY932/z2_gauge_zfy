module OperatorKin_mod
    use MyLattice
    use calc_basic
    use, intrinsic :: ieee_arithmetic
    implicit none
    public

    type :: OperatorKin
        real(kind = 8), public :: alpha, alpha_2  ! 改为 public 以便在 localK 中使用
        real(kind = 8), private :: entryC
        real(kind = 8), private :: entryS

        real(kind = 8), public :: Delta(2, 2)
    contains
        procedure :: set => opK_set
        procedure, private :: get_exp => opK_get_exp
        procedure :: get_delta => opK_get_delta
        procedure :: mmult_R => opK_mmult_R
        procedure :: mmult_L => opK_mmult_L
        procedure :: mmult_R_group => opK_mmult_R_group
        procedure :: mmult_L_group => opK_mmult_L_group
    end type OperatorKin

    real(kind=8), parameter :: HYPER_LIMIT = 40.d0
    logical :: hyper_warned = .false.

contains
    subroutine safe_cosh_sinh(arg, Cv, Sv)
        real(kind=8), intent(in)  :: arg
        real(kind=8), intent(out) :: Cv, Sv
        real(kind=8) :: arg_clamped

        arg_clamped = arg
        if (abs(arg_clamped) > HYPER_LIMIT) then
            arg_clamped = sign(HYPER_LIMIT, arg_clamped)
            if (.not. hyper_warned) then
                write(6,*) "WARNING: opK cosh/sinh argument clipped:", arg, "->", arg_clamped
                hyper_warned = .true.
            endif
        endif

        Cv = cosh(arg_clamped)
        Sv = sinh(arg_clamped)
    end subroutine safe_cosh_sinh

    subroutine opK_set(this)
        class(OperatorKin), intent(inout) :: this
        ! 标准棋盘分解：每个 group 操作一次，使用完整步长
        ! 这与 get_delta 中的计算一致
        this%alpha = Dtau * RT
        this%alpha_2 = Dtau * RT
        return
    end subroutine opK_set

    subroutine opK_get_exp(this, sigma, nflag)
        class(OperatorKin), intent(inout) :: this
        integer, intent(in) :: nflag
        real(kind = 8), intent(in) :: sigma
        call safe_cosh_sinh(nflag * this%alpha * sigma, this%entryC, this%entryS)
        return
    end subroutine opK_get_exp

    subroutine opK_get_delta(this, sigma_old, sigma_new)
! Calculate Delta = exp(alpha*sigma_new) * exp(-alpha*sigma_old) - I
! Used for efficient Green's function update via Sherman-Morrison formula
        class(OperatorKin), intent(inout) :: this
        real(kind = 8), intent(in) :: sigma_old, sigma_new
        real(kind = 8) :: C_new, C_old, S_new, S_old

        call safe_cosh_sinh(this%alpha_2 * sigma_new, C_new, S_new)
        call safe_cosh_sinh(this%alpha_2 * sigma_old, C_old, S_old)
        
        ! Delta = exp(+alpha*sigma_new) * exp(-alpha*sigma_old) - I
        ! Note: sinh is odd function, so exp(-alpha*sigma_old) has -S_old entries
        this%Delta(1, 1) = C_new * C_old - S_new * S_old - 1.d0
        this%Delta(1, 2) = S_new * C_old - C_new * S_old
        this%Delta(2, 1) = S_new * C_old - C_new * S_old
        this%Delta(2, 2) = C_new * C_old - S_new * S_old - 1.d0
        return
    end subroutine opK_get_delta

    subroutine opK_mmult_R(this, Mat, Latt, sigma, ntau, nflag)
! In Mat Out B * Mat（nflag=1）或 B^{-1} * Mat（nflag=-1）
! B = B_4 * B_3 * B_2 * B_1（棋盘分解）
! 
! 对于 nflag=1 (B * Mat)：
!   B * Mat = B_4 * B_3 * B_2 * B_1 * Mat
!   计算顺序：group_1 -> group_2 -> group_3 -> group_4
!
! 对于 nflag=-1 (B^{-1} * Mat)：
!   B^{-1} = B_1^{-1} * B_2^{-1} * B_3^{-1} * B_4^{-1}
!   计算顺序：group_4 -> group_3 -> group_2 -> group_1（反序）
!
! Arguments: 
        class(OperatorKin), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(SquareLattice), intent(in) :: Latt
        real(kind = 8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: ntau, nflag

        if (nflag > 0) then
            ! B * Mat：group_1 -> group_2 -> group_3 -> group_4
            call process_group(this, Mat, Latt%group_1, Latt, sigma, ntau, nflag)
            call process_group(this, Mat, Latt%group_2, Latt, sigma, ntau, nflag)
            call process_group(this, Mat, Latt%group_3, Latt, sigma, ntau, nflag)
            call process_group(this, Mat, Latt%group_4, Latt, sigma, ntau, nflag)
        else
            ! B^{-1} * Mat：group_4 -> group_3 -> group_2 -> group_1（反序）
            call process_group(this, Mat, Latt%group_4, Latt, sigma, ntau, nflag)
            call process_group(this, Mat, Latt%group_3, Latt, sigma, ntau, nflag)
            call process_group(this, Mat, Latt%group_2, Latt, sigma, ntau, nflag)
            call process_group(this, Mat, Latt%group_1, Latt, sigma, ntau, nflag)
        endif

        return
    end subroutine opK_mmult_R
    
    ! 按组进行右乘（用于 LocalK_prop_R 中的分组更新）
    subroutine opK_mmult_R_group(this, Mat, group, Latt, sigma, ntau, nflag)
        class(OperatorKin), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        integer, dimension(:), intent(in) :: group
        class(SquareLattice), intent(in) :: Latt
        real(kind = 8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: ntau, nflag

        call process_group(this, Mat, group, Latt, sigma, ntau, nflag)
        return
    end subroutine opK_mmult_R_group

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
! In Mat Out Mat * B（nflag=1）或 Mat * B^{-1}（nflag=-1）
! B = B_4 * B_3 * B_2 * B_1（棋盘分解）
! 
! 对于 nflag=1 (Mat * B)：
!   Mat * B = Mat * B_4 * B_3 * B_2 * B_1
!   计算顺序：group_4 -> group_3 -> group_2 -> group_1（从右往左乘）
!
! 对于 nflag=-1 (Mat * B^{-1})：
!   Mat * B^{-1} = Mat * B_1^{-1} * B_2^{-1} * B_3^{-1} * B_4^{-1}
!   计算顺序：group_1 -> group_2 -> group_3 -> group_4（从右往左乘逆）
!
! Arguments: 
        class(OperatorKin), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(SquareLattice), intent(in) :: Latt
        real(kind = 8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: ntau, nflag

        if (nflag > 0) then
            ! Mat * B：group_4 -> group_3 -> group_2 -> group_1
            call process_group_L(this, Mat, Latt%group_4, Latt, sigma, ntau, nflag)
            call process_group_L(this, Mat, Latt%group_3, Latt, sigma, ntau, nflag)
            call process_group_L(this, Mat, Latt%group_2, Latt, sigma, ntau, nflag)
            call process_group_L(this, Mat, Latt%group_1, Latt, sigma, ntau, nflag)
        else
            ! Mat * B^{-1}：group_1 -> group_2 -> group_3 -> group_4（反序）
            call process_group_L(this, Mat, Latt%group_1, Latt, sigma, ntau, nflag)
            call process_group_L(this, Mat, Latt%group_2, Latt, sigma, ntau, nflag)
            call process_group_L(this, Mat, Latt%group_3, Latt, sigma, ntau, nflag)
            call process_group_L(this, Mat, Latt%group_4, Latt, sigma, ntau, nflag)
        endif

        return
    end subroutine opK_mmult_L
    
    ! 按组进行左乘（用于 LocalK_prop_L 中的分组更新）
    subroutine opK_mmult_L_group(this, Mat, group, Latt, sigma, ntau, nflag)
        class(OperatorKin), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        integer, dimension(:), intent(in) :: group
        class(SquareLattice), intent(in) :: Latt
        real(kind = 8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: ntau, nflag

        call process_group_L(this, Mat, group, Latt, sigma, ntau, nflag)
        return
    end subroutine opK_mmult_L_group

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
