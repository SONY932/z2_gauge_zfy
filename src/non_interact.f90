module NonInteract
! 该模块目前只支持化学势 mu 的乘法处理
    use MyLattice
    use calc_basic
    implicit none

    public

    type :: OperatorMU
        real(kind=8) :: expMU_P, expMU_M
    contains
        procedure :: make => opMU_make
        procedure :: set => opMU_set
        procedure :: mmult_R => opMU_mmult_R
        procedure :: mmult_L => opMU_mmult_L
    end type OperatorMU

    contains
    subroutine opMU_make(this)
        class(OperatorMU), intent(inout) :: this
        ! 标量不需要allocate，直接初始化
        this%expMU_P = 0.d0
        this%expMU_M = 0.d0
        return
    end subroutine opMU_make

    subroutine opMU_set(this)
        class(OperatorMU), intent(inout) :: this
        this%expMU_P = exp( Dtau * mu)
        this%expMU_M = exp(- Dtau * mu)
        return
    end subroutine opMU_set

    subroutine opMU_mmult_R(this, Mat, lambda, nt, nflag)
!   In Mat Out exp(+Dtau*mu) * Mat for nflag = 1 (B_mu part of propagator)
!   In Mat Out exp(-Dtau*mu) * Mat for nflag = -1 (B_mu^{-1} part)
!   传播阶段只处理化学势 μ，λ 的投影只在构造 M=I+P_lambda·B_tot 时出现。
        class(OperatorMU), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        real(kind=8), dimension(Lq), intent(in) :: lambda  ! 当前不在此处使用
        integer, intent(in) :: nt
        integer, intent(in) :: nflag
        real(kind=8) :: factor

        ! nflag=1: multiply by B_mu = exp(+Dtau*mu)
        ! nflag=-1: multiply by B_mu^{-1} = exp(-Dtau*mu)
        if (nflag == 1) then
            factor = this%expMU_P  ! exp(+Dtau*mu)
        elseif (nflag == -1) then
            factor = this%expMU_M  ! exp(-Dtau*mu)
        else
            write(6,*) "incorrect nflag in opMU_mmult"; stop
        endif

        Mat = Mat * factor
        return
    end subroutine opMU_mmult_R

    subroutine opMU_mmult_L(this, Mat, lambda, nt, nflag)
!   In Mat Out Mat * exp(+Dtau*mu) for nflag = 1 (B_mu part of propagator)
!   In Mat Out Mat * exp(-Dtau*mu) for nflag = -1 (B_mu^{-1} part)
!   同上：传播阶段只处理 μ，不在此处插入 λ。
        class(OperatorMU), intent(inout) :: this
        complex(kind = 8), dimension(Ndim, Ndim), intent(inout) :: Mat
        real(kind=8), dimension(Lq), intent(in) :: lambda  ! 当前不在此处使用
        integer, intent(in) :: nt
        integer, intent(in) :: nflag
        real(kind=8) :: factor

        ! nflag=1: multiply by B_mu = exp(+Dtau*mu)
        ! nflag=-1: multiply by B_mu^{-1} = exp(-Dtau*mu)
        if (nflag == 1) then
            factor = this%expMU_P  ! exp(+Dtau*mu)
        elseif (nflag == -1) then
            factor = this%expMU_M  ! exp(-Dtau*mu)
        else
            write(6,*) "incorrect nflag in opMU_mmult_L"; stop
        endif

        Mat = Mat * factor
        return
    end subroutine opMU_mmult_L

end module NonInteract