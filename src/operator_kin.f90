module OperatorKin_mod
    use MyLattice
    use calc_basic
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

contains
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
! 按顺序处理所有键（不使用棋盘分解）
! Arguments: 
        class(OperatorKin), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(SquareLattice), intent(in) :: Latt
        real(kind = 8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: ntau, nflag
! Local:
        integer :: ii, P(2), jj
        real(kind = 8) :: vec
        complex(kind=8), dimension(2, Ndim) :: Vhlp

        ! 按顺序处理所有键
        do ii = 1, 2*Lq
            vec = sigma(ii, ntau)
            call this%get_exp(vec, nflag)
            P(1) = Latt%bond_list(ii, 1)
            P(2) = Latt%bond_list(ii, 2)
            
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
! In Mat Out Mat * EXP(D(NF)) * UT(NF)
! 按相反顺序处理所有键，确保和 mmult_R 给出相同的 B 矩阵
! mmult_R 给出 B * G = B_n * ... * B_1 * G（按正向顺序处理）
! mmult_L 给出 G * B = G * B_n * ... * B_1（按相反顺序处理）
! 这样 B 矩阵定义一致：B = B_n * ... * B_1
! 按相反顺序处理时：先乘 B_n，再乘 B_{n-1}，...，最后乘 B_1
! 结果是 Mat * B_n * B_{n-1} * ... * B_1 = Mat * B
! Arguments: 
        class(OperatorKin), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(SquareLattice), intent(in) :: Latt
        real(kind = 8), dimension(2*Lq, Ltrot), intent(in) :: sigma
        integer, intent(in) :: ntau, nflag
! Local:
        integer :: ii, P(2), jj
        real(kind = 8) :: vec
        complex(kind=8), dimension(Ndim, 2) :: Uhlp

        ! 按相反顺序处理所有键：先乘 B_n，最后乘 B_1
        ! Mat -> Mat * B_n -> Mat * B_n * B_{n-1} -> ... -> Mat * B_n * ... * B_1
        do ii = 2*Lq, 1, -1
            vec = sigma(ii, ntau)
            call this%get_exp(vec, nflag)
            P(1) = Latt%bond_list(ii, 1)
            P(2) = Latt%bond_list(ii, 2)
            
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
