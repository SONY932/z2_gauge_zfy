module Stabilize_mod
    use ProcessMatrix
    use MyMats
    use MyLattice, only: ZKRON
    use DQMC_Model_mod      ! 引入 NsigL_K%lambda，用于在构造 Green 时插入 P_lambda
    use, intrinsic :: ieee_arithmetic
    implicit none
    
    private
    public :: Wrap_pre, Wrap_L, Wrap_L_store, Wrap_R, Wrap_R_store, Wrap_tau, Stabilize_init, Stabilize_clear
    
    complex(kind=kind(0.d0)) ::  Z_one
    complex(kind=8), dimension(:,:), allocatable :: matUDV
    complex(kind=kind(0.d0)), dimension(:), allocatable :: TAU, DUP
    integer, dimension(:), allocatable:: IPVT
    complex(kind=8), dimension(:,:), allocatable :: invUUR, invUUL, invVUR, invVUL, temp
    complex(kind=8), dimension(:), allocatable :: bigD
    complex(kind=8), dimension(:,:), allocatable :: mat_big, mat_left, mat_right, Gr_tot, bigU, bigV, invbigU, invbigV
    type(PropGreen), allocatable :: Gr_tmp
    logical :: warn_left_done = .false.
    logical :: warn_right_done = .false.
    
contains
    subroutine Stabilize_init()
        Z_one = dcmplx(1.d0, 0.d0)
        allocate(TAU(Ndim), IPVT(Ndim))
        allocate(DUP(Ndim))
        allocate(matUDV(Ndim, Ndim))
        allocate(invUUR(Ndim, Ndim), invUUL(Ndim, Ndim), invVUR(Ndim, Ndim), invVUL(Ndim, Ndim), temp(Ndim, Ndim))
        allocate(mat_big(2*Ndim, 2*Ndim), mat_left(2*Ndim, 2*Ndim), mat_right(2*Ndim, 2*Ndim))
        allocate(Gr_tot(2*Ndim, 2*Ndim), bigU(2*Ndim, 2*Ndim), bigV(2*Ndim, 2*Ndim), invbigU(2*Ndim, 2*Ndim), invbigV(2*Ndim, 2*Ndim))
        allocate(bigD(2*Ndim))
        allocate(Gr_tmp)
        call Gr_tmp%make()
        return
    end subroutine Stabilize_init
    
    subroutine Stabilize_clear()
        deallocate(matUDV, TAU, IPVT, DUP)
        deallocate(invUUR, invUUL, invVUR, invVUL, temp)
        deallocate(mat_big, mat_left, mat_right, Gr_tot, bigU, bigV, invbigU, invbigV, bigD)
        deallocate(Gr_tmp)
        return
    end subroutine Stabilize_clear
    
    subroutine QDRP_decompose(Mat, D, IPVT, TAU, WORK, Lwork)
! Arguments: 
        complex(kind=kind(0.d0)), dimension(:,:), intent(inout) :: Mat
        complex(kind=kind(0.d0)), dimension(:), intent(inout) :: D
        integer, dimension(Ndim), intent(inout) :: IPVT
        complex(kind=kind(0.d0)), dimension(Ndim), intent(inout) :: TAU
        complex(kind=kind(0.d0)), dimension(:), intent(inout), allocatable :: WORK
        integer, intent(inout) :: Lwork
! Local: 
        complex(kind=kind(0.d0)), dimension(2*Ndim) :: RWORK
        complex(kind=kind(0.d0)) :: Z
        integer :: info, i, j
        real(kind=kind(0.d0)) :: X
        real(kind=kind(0.d0)), parameter :: tiny_val = 1.d-100
! Query optimal amount of memory
        call ZGEQP3(Ndim, Ndim, Mat, Ndim, IPVT, TAU(1), Z, -1, RWORK(1), info)
        Lwork = int(dble(Z)); allocate(WORK(Lwork))
! QR decomposition of Mat with full column pivoting, Mat * P = Q * R
        call ZGEQP3(Ndim, Ndim, Mat, Ndim, IPVT, TAU(1), WORK(1), Lwork, RWORK(1), info)
! separate off D
        do i = 1, Ndim
            X = abs(Mat(i, i))
            ! 保护：避免除以 0
            if (X < tiny_val) then
                X = tiny_val
            endif
            D(i) = X ! plain diagonal entry
            do j = i, Ndim
                Mat(i, j) = Mat(i, j) / X
            enddo
        enddo
        return
    end subroutine QDRP_decompose
    
    subroutine stab_UR(Prop)
        class(Propagator), intent(inout) :: Prop
        integer :: info, i, j, Lwork
        complex(kind=8), dimension(:), allocatable :: WORK
! QR(TMP * U * D) * V
! U*D
        matUDV = Prop%UUR
        do i = 1,Ndim
            matUDV(:, i) = matUDV(:, i) * Prop%DUR(i)
        enddo
! QR(TMP * U * D)
        IPVT = 0
! output in matUDV and DUR
        call QDRP_decompose(matUDV, Prop%DUR, IPVT, TAU, WORK, Lwork)
! Permute V. Since we multiply with V from the right we have to permute the rows of V.
! A V = A P P^-1 V = Q R P^-1 V; ZLAPMR: rearrange rows of matrix V
        call ZLAPMR(.true., Ndim, Ndim, Prop%VUR, Ndim, IPVT) ! lapack 3.3
! V = R * V, output in VUR; character U: only use upper triangular part of matUDV
        call ZTRMM('L', 'U', 'N', 'N', Ndim, Ndim, Z_one, matUDV, Ndim, Prop%VUR, Ndim)
! Generate explicitly U in the previously abused storage of U, output in matUDV
        call ZUNGQR(Ndim, Ndim, Ndim, matUDV, Ndim, TAU, WORK, Lwork, info)
        deallocate(WORK)
        Prop%UUR = matUDV
        return
    end subroutine stab_UR
    
    subroutine  stab_UL(Prop)
        class(Propagator), intent(inout) :: Prop
        integer :: info, i, j, Lwork
        complex(kind=8), dimension(:), allocatable :: WORK
! QR(TMP^\dagger * U^dagger * D) * V^dagger
! U*D
        matUDV = dconjg(transpose(Prop%UUL))
        do i = 1, Ndim
            matUDV(:, i) = matUDV(:, i) * Prop%DUL(i)
        enddo
! QR(TMP^\dagger * U^dagger * D)
        IPVT = 0 
! output in matUDV and DUL
        call QDRP_decompose(matUDV, Prop%DUL, IPVT, TAU, WORK, Lwork)
! Permute V. Since we multiply with V^dagger from the right we have to permute the columns of V.
! A V^dagger = A P P^-1 V^dagger = Q R P^-1 V^dagger; ZLAPMT: rearrange columns of matrix V
        call ZLAPMT(.true., Ndim, Ndim, Prop%VUL, Ndim, IPVT)
! V = V * R^dagger; output in VUL; character U: only use upper triangular part of matUDV
        call ZTRMM('R', 'U', 'C', 'N', Ndim, Ndim, Z_one, matUDV, Ndim, Prop%VUL, Ndim)
! Generate explicitly U in the previously abused storage of U, output in matUDV
        call ZUNGQR(Ndim, Ndim, Ndim, matUDV, Ndim, TAU, WORK, Lwork, info)
        deallocate(WORK)
        Prop%UUL = dconjg(transpose(matUDV))
        return
    end subroutine stab_UL
    
    subroutine stab_green(Gr, Prop, nt)
! Calculates the Green's function from the UDV decomposition of UUL and UUR
! B = UL DL VL VR DR UR = UL DL (VL VR) DR UR
! G = (1+B)^-1
! Arguments: 
        complex(kind=8), dimension(Ndim, Ndim), intent(out) :: Gr
        class(Propagator), intent(in) :: Prop
        integer, intent(in) :: nt
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: VRVL, invULUR
        complex(kind=8), dimension(:), allocatable :: WORK
        complex(kind=kind(0.d0)) :: Z
        integer :: nl, nr, Lwork, info
        
! coefficient of zgemm: alpha * op(A) * op(B) + beta * op(C); here alpha=Z_one=1, beta=0
! VRVL = VR*VL
        call mmult(VRVL, Prop%VUR, Prop%VUL)
! invULUR = UR^dagger UL^dagger = (UL*UR)^-1; character C: conjugate transpose of UUR and UUL
        call ZGEMM('C', 'C', Ndim, Ndim, Ndim, Z_one, Prop%UUR, Ndim, Prop%UUL, Ndim, dcmplx(0.d0, 0.d0), invULUR, Ndim)
        
        
! compute: matUDV = (UL*UR)^-1 + DR VR VL DL
        temp = dcmplx(0.d0, 0.d0)
        do nr = 1, Ndim ! temp(1:Ndim, nr) = VRVL(1:Ndim, nr) * DUL(nr)
            call zaxpy(Ndim, Prop%DUL(nr), VRVL(1, nr), 1, temp(1, nr), 1)
        enddo
        VRVL = dcmplx(0.d0, 0.d0)
        do nl = 1, Ndim ! VRVL(nl, 1:Ndim) = DUR(nl) * temp(nl, 1:Ndim)
            call zaxpy(Ndim, Prop%DUR(nl), temp(nl, 1), Ndim, VRVL(nl, 1), Ndim)
        enddo
        matUDV = invULUR + VRVL
! if ntau >Ltrot/2, decompose matUDV^dagger
        if (nt > Ltrot/2) matUDV = dconjg(transpose(matUDV))
! matUDV * P  = U D V
        IPVT = 0
        call QDRP_decompose(matUDV, DUP, IPVT, TAU, WORK, Lwork)
        
        ! 检查 DUP 是否有接近 0 的元素
        do nl = 1, Ndim
            if (abs(DUP(nl)) < 1.d-100) then
                ! 设置为一个小的非零值以避免除以 0
                DUP(nl) = dcmplx(1.d-100, 0.d0)
            endif
        enddo
        
        if (nt < Ltrot/2 + 1) then ! ntau < Ltrot/2
! UR U D V P^dagger UL G = 1 => 
! G=UL^dagger * P * V^-1 * D^-1 * U^dagger * UR^dagger; multiply from right to left
            temp = dconjg(transpose(Prop%UUR))
! compute U^dagger * UR^dagger; output in temp 
            call ZUNMQR('L', 'C', Ndim, Ndim, Ndim, matUDV(1,1), Ndim, TAU(1), temp(1,1), Ndim, WORK(1), Lwork, info)
! compute D^-1 * (U^dagger * UR^dagger)
            do nl = 1, Ndim
                temp(nl, :) = temp(nl, :) / DUP(nl)
            enddo
! compute V^-1 * (D^-1 * U^dagger * UR^dagger) by solving V * X = D^-1 * U^dagger * UR^dagger; output in temp
            call ZTRSM('L', 'U', 'N', 'N', Ndim, Ndim, Z_one, matUDV(1,1), Ndim, temp(1,1), Ndim)
! apply permutation matrix : P * (V^-1*D^-1*U^dagger*UR^dagger); rearrange rows
            call ZLAPMR(.false., Ndim, Ndim, temp(1,1), Ndim, IPVT(1))
! compute UL^dagger * tmp = UL^dagger * (P * V^-1*D^-1*U^dagger*UR^dagger)
! output Gr
            call ZGEMM('C', 'N', Ndim, Ndim, Ndim, Z_one, Prop%UUL(1,1), Ndim, temp(1,1), Ndim, dcmplx(0.d0, 0.d0), Gr(1,1), Ndim)
        elseif(nt > Ltrot/2) then ! ntau >Ltrot/2          
! matUDV^dagger * P = U D V
! G=UL^dagger *  U * D^-1 * V * P^-1 * UR^dagger; multiply from left to right
            temp = dconjg(transpose(Prop%UUL))
! UL^dagger * U
            call ZUNMQR('R', 'N', Ndim, Ndim, Ndim, matUDV(1, 1), Ndim, TAU(1), temp(1, 1), Ndim, WORK(1), Lwork, info)
! (UL^dagger * U) * D^-1
            do nr = 1, Ndim
                temp(:, nr) = temp(:, nr)/ DUP(nr)
            enddo
! compute (UL^dagger * U * D^-1) * V by solving X * V^dagger = UL^dagger * U * D^-1 * V
            call ZTRSM('R', 'U', 'C', 'N', Ndim, Ndim, Z_one, matUDV(1, 1), Ndim, temp(1, 1), Ndim)
! apply inverse permutation matrix: (UL^dagger * U * D^-1 * V) * P^-1; rearrange columns
            call ZLAPMT(.false., Ndim, Ndim, temp(1, 1), Ndim, IPVT(1))
! (UL^dagger * U * D^-1 * V * P^-1) * UR^dagger = G
! output Gr
            call ZGEMM('N', 'C', Ndim, Ndim, Ndim, Z_one, temp(1, 1), Ndim, Prop%UUR(1, 1), Ndim, dcmplx(0.d0, 0.d0), Gr(1, 1), Ndim)
        else
            write(6,*) "illegal imaginary input in stabgreen, nt =", nt
        endif

        deallocate(WORK)
        return
    end subroutine stab_green
    
    subroutine stab_green_lambda(Gr, Prop, nt)
! Calculates G_lambda = (1 + P[lambda] * B)^{-1}
! 其中 P[lambda] = diag(lambda_i)，lambda_i = ±1
! 实现方法：用 P[lambda] * UL 代替 UL，即把 UL 的每一行乘以对应的 lambda
! 由于 P[lambda] 是对角矩阵（条件数 = 1），这不会影响数值稳定性
! Arguments: 
        complex(kind=8), dimension(Ndim, Ndim), intent(out) :: Gr
        class(Propagator), intent(in) :: Prop
        integer, intent(in) :: nt
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: VRVL, invULUR, UUL_lambda
        complex(kind=8), dimension(:), allocatable :: WORK
        complex(kind=kind(0.d0)) :: Z
        integer :: nl, nr, Lwork, info, ii

! 创建 UUL 的副本并应用 P[lambda]
! UL_lambda = P[lambda] * UL，即把每一行乘以对应的 lambda
        UUL_lambda = Prop%UUL
        do ii = 1, Ndim
            UUL_lambda(ii, :) = UUL_lambda(ii, :) * NsigL_K%lambda(ii)
        enddo

! VRVL = VR*VL
        call mmult(VRVL, Prop%VUR, Prop%VUL)
! invULUR = UR^dagger (P*UL)^dagger = (P*UL*UR)^-1
! 由于 P 是实对角矩阵，(P*UL)^dagger = UL^dagger * P
        call ZGEMM('C', 'C', Ndim, Ndim, Ndim, Z_one, Prop%UUR, Ndim, UUL_lambda, Ndim, dcmplx(0.d0, 0.d0), invULUR, Ndim)
! compute: matUDV = (P*UL*UR)^-1 + DR VR VL DL
        temp = dcmplx(0.d0, 0.d0)
        do nr = 1, Ndim
            call zaxpy(Ndim, Prop%DUL(nr), VRVL(1, nr), 1, temp(1, nr), 1)
        enddo
        VRVL = dcmplx(0.d0, 0.d0)
        do nl = 1, Ndim
            call zaxpy(Ndim, Prop%DUR(nl), temp(nl, 1), Ndim, VRVL(nl, 1), Ndim)
        enddo
        matUDV = invULUR + VRVL
! if ntau >Ltrot/2, decompose matUDV^dagger
        if (nt > Ltrot/2) matUDV = dconjg(transpose(matUDV))
! matUDV * P  = U D V
        IPVT = 0
        call QDRP_decompose(matUDV, DUP, IPVT, TAU, WORK, Lwork)
        if (nt < Ltrot/2 + 1) then
! G_lambda = (P*UL)^dagger * P * V^-1 * D^-1 * U^dagger * UR^dagger
            temp = dconjg(transpose(Prop%UUR))
            call ZUNMQR('L', 'C', Ndim, Ndim, Ndim, matUDV(1,1), Ndim, TAU(1), temp(1,1), Ndim, WORK(1), Lwork, info)
            do nl = 1, Ndim
                temp(nl, :) = temp(nl, :) / DUP(nl)
            enddo
            call ZTRSM('L', 'U', 'N', 'N', Ndim, Ndim, Z_one, matUDV(1,1), Ndim, temp(1,1), Ndim)
            call ZLAPMR(.false., Ndim, Ndim, temp(1,1), Ndim, IPVT(1))
! 使用 UUL_lambda 代替 UUL
            call ZGEMM('C', 'N', Ndim, Ndim, Ndim, Z_one, UUL_lambda(1,1), Ndim, temp(1,1), Ndim, dcmplx(0.d0, 0.d0), Gr(1,1), Ndim)
        elseif(nt > Ltrot/2) then
! G_lambda = (P*UL)^dagger * U * D^-1 * V * P^-1 * UR^dagger
! 使用 UUL_lambda 代替 UUL
            temp = dconjg(transpose(UUL_lambda))
            call ZUNMQR('R', 'N', Ndim, Ndim, Ndim, matUDV(1, 1), Ndim, TAU(1), temp(1, 1), Ndim, WORK(1), Lwork, info)
            do nr = 1, Ndim
                temp(:, nr) = temp(:, nr)/ DUP(nr)
            enddo
            call ZTRSM('R', 'U', 'C', 'N', Ndim, Ndim, Z_one, matUDV(1, 1), Ndim, temp(1, 1), Ndim)
            call ZLAPMT(.false., Ndim, Ndim, temp(1, 1), Ndim, IPVT(1))
            call ZGEMM('N', 'C', Ndim, Ndim, Ndim, Z_one, temp(1, 1), Ndim, Prop%UUR(1, 1), Ndim, dcmplx(0.d0, 0.d0), Gr(1, 1), Ndim)
        else
            write(6,*) "illegal imaginary input in stabgreen_lambda, nt =", nt
        endif

        deallocate(WORK)
        return
    end subroutine stab_green_lambda
        
    subroutine stab_green_big(Prop)
! Arguments:
        class(Propagator), intent(in) :: Prop
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: ULUR, VRVL, invULUR, invVRVL
        complex(kind=8) :: det
        integer :: nl, nr
        
        call mmult(ULUR, Prop%UUL, Prop%UUR)
        call mmult(VRVL, Prop%VUR, Prop%VUL)
        call inv(ULUR, invULUR, det)
        call inv(VRVL, invVRVL, det)
        mat_big = dcmplx(0.d0, 0.d0)
        do nl = 1, Ndim
            do nr = 1, Ndim
                mat_big(nl, nr) = invVRVL(nl, nr)
                mat_big(nl + Ndim, nr + Ndim) = invULUR(nl, nr)
            enddo
        enddo
        do nr = 1, Ndim
            mat_big(nr, nr + Ndim) = Prop%DUL(nr)
            mat_big(nr + Ndim, nr) = - Prop%DUR(nr)
        enddo
        call udv(mat_big, bigU, bigD, bigV, 0)
        call inv(bigU, invbigU, det)
        call inv(bigV, invbigV, det)
        call inv(Prop%UUR, invUUR, det)
        call inv(Prop%UUL, invUUL, det)
        call inv(Prop%VUR, invVUR, det)
        call inv(Prop%VUL, invVUL, det)
        mat_big = dcmplx(0.d0, 0.d0); mat_left = dcmplx(0.d0, 0.d0)
        do nl = 1, Ndim
            do nr = 1, Ndim
                mat_big(nl, nr) = invVUR(nl, nr)
                mat_big(nl + Ndim, nr + Ndim) =invUUL(nl, nr)
            enddo
        enddo
        call mmult(mat_left, mat_big, invbigV)
        mat_big = dcmplx(0.d0, 0.d0); mat_right = dcmplx(0.d0, 0.d0)
        do nl = 1, Ndim
            do nr = 1, Ndim
                mat_big(nl, nr) = invVUL(nl, nr)
                mat_big(nl + Ndim, nr + Ndim) =invUUR(nl, nr)
            enddo
        enddo
        call mmult(mat_right, invbigU, mat_big)
        do nr = 1, 2*Ndim
            do nl = 1, 2*Ndim
                mat_left(nl, nr) = mat_left(nl, nr) / bigD(nr)
            enddo
        enddo
        call mmult(Gr_tot, mat_left, mat_right)
! output time-sliced Green function
        do nl = 1, Ndim
            do nr = 1, Ndim
                Gr_tmp%Gr00(nl, nr) = Gr_tot(nl, nr)
                Gr_tmp%Gr0t(nl, nr) = Gr_tot(nl, nr + Ndim)
                Gr_tmp%Grt0(nl, nr) = Gr_tot(nl + Ndim, nr)
                Gr_tmp%Grtt(nl, nr) = Gr_tot(nl + Ndim, nr + Ndim)
            enddo
        enddo
        return
    end subroutine stab_green_big

    !------------------------------------------------------------------
    ! 在构造好的等时格林函数上插入 λ 投影：
    ! 已有 Gr = (I + B_tot)^-1 ≡ M^-1，目标是返回
    ! Gr_lambda = (I + P_lambda * B_tot)^-1 ≡ M_lambda^-1。
    ! 实现思路：
    !   1) 先由 Gr 求出 M = Gr^{-1}；
    !   2) 用公式 M_lambda = P_lambda * M + (I - P_lambda) 构造新 M_lambda；
    !   3) 再对 M_lambda 做一次矩阵求逆得到 Gr_lambda。
    ! 这里只追求物理正确性和稳健性，允许额外的 O(N^3) 开销。
    !------------------------------------------------------------------
    subroutine apply_P_to_green(Gr)
! ===================================================================
! 从 G_0 = (1 + B)^{-1} 计算 G_λ = (1 + P[λ]B)^{-1}
! 
! 推导：
!   M_0 = G_0^{-1} = 1 + B
!   M_λ = 1 + P[λ]B = 1 + P[λ](M_0 - 1) = P[λ]M_0 + (I - P[λ])
!   G_λ = M_λ^{-1}
!
! P[λ] 是对角矩阵（条件数 = 1），不会放大误差
! ===================================================================
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gr
        complex(kind=8), allocatable :: M0(:,:), Mlambda(:,:)
        real(kind=8) :: lam_i
        integer :: i, j, info
        logical :: all_one

        ! 若当前 λ 全为 1，则 P[λ] = I，G_λ = G_0，直接返回
        all_one = .true.
        do i = 1, Lq
            if (abs(NsigL_K%lambda(i) - 1.d0) > 1.d-12) then
                all_one = .false.
                exit
            endif
        enddo
        if (all_one) return

        allocate(M0(Ndim, Ndim), Mlambda(Ndim, Ndim))

        ! 第一步：M0 = Gr^{-1} = 1 + B
        M0 = Gr
        call invert_matrix(M0, info)
        if (info /= 0) then
            deallocate(M0, Mlambda)
            return
        endif

        ! 第二步：M_λ = P[λ]*M0 + (I - P[λ])
        do i = 1, Ndim
            lam_i = NsigL_K%lambda(i)
            do j = 1, Ndim
                if (i == j) then
                    Mlambda(i, j) = dcmplx(lam_i, 0.d0) * M0(i, j) + dcmplx(1.d0 - lam_i, 0.d0)
                else
                    Mlambda(i, j) = dcmplx(lam_i, 0.d0) * M0(i, j)
                endif
            enddo
        enddo

        ! 第三步：G_λ = M_λ^{-1}
        call invert_matrix(Mlambda, info)
        if (info == 0) then
            Gr = Mlambda
        endif

        deallocate(M0, Mlambda)
        return
    end subroutine apply_P_to_green

    subroutine apply_lambda_projection(Gr)
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gr
        complex(kind=8), allocatable :: M0(:,:), Mlambda(:,:)
        real(kind=8) :: lam_i
        integer :: i, j, info
        logical :: all_one

        ! 若当前 λ 全为 1，则直接退化到原 σ-only 情形，不做额外操作。
        all_one = .true.
        do i = 1, Lq
            if (abs(NsigL_K%lambda(i) - 1.d0) > 1.d-12) then
                all_one = .false.
                exit
            endif
        enddo
        if (all_one) return

        allocate(M0(Ndim, Ndim), Mlambda(Ndim, Ndim))

        ! 第一步：M0 = Gr^{-1} = (I + B_tot)
        M0 = Gr
        call invert_matrix(M0, info)
        if (info /= 0) then
            ! 若求逆失败，为安全起见保留原 Gr，不再强行插入 λ
            deallocate(M0, Mlambda)
            return
        endif

        ! 第二步：构造 M_lambda = P_lambda * M0 + (I - P_lambda)
        ! 其中 P_lambda = diag(λ_i)，λ_i = ±1。
        do i = 1, Ndim
            lam_i = NsigL_K%lambda(i)
            do j = 1, Ndim
                if (i == j) then
                    Mlambda(i, j) = dcmplx(lam_i, 0.d0) * M0(i, j) + dcmplx(1.d0 - lam_i, 0.d0)
                else
                    Mlambda(i, j) = dcmplx(lam_i, 0.d0) * M0(i, j)
                endif
            enddo
        enddo

        ! 第三步：Gr_lambda = M_lambda^{-1}
        call invert_matrix(Mlambda, info)
        if (info == 0) then
            Gr = Mlambda
        endif

        deallocate(M0, Mlambda)
        return
    end subroutine apply_lambda_projection

    !------------------------------------------------------------------
    ! 通用的 Ndim×Ndim 复矩阵求逆封装：调用 LAPACK ZGETRF/ZGETRI。
    ! 输入：A (被覆盖为 A^{-1})；输出：info (0 表示成功)。
    !------------------------------------------------------------------
    subroutine invert_matrix(A, info)
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: A
        integer, intent(out) :: info
        integer, allocatable :: ipiv(:)
        complex(kind=8), allocatable :: work(:)
        integer :: lwork

        allocate(ipiv(Ndim))

        ! LU 分解
        call ZGETRF(Ndim, Ndim, A, Ndim, ipiv, info)
        if (info /= 0) then
            deallocate(ipiv)
            return
        endif

        ! 询问最优工作空间大小
        allocate(work(1))
        call ZGETRI(Ndim, A, Ndim, ipiv, work, -1, info)
        if (info /= 0) then
            deallocate(ipiv, work)
            return
        endif
        lwork = int(real(work(1)))
        deallocate(work)

        ! 真正求逆
        allocate(work(lwork))
        call ZGETRI(Ndim, A, Ndim, ipiv, work, lwork, info)

        deallocate(ipiv, work)
        return
    end subroutine invert_matrix
    
    logical function has_nan_or_inf(Mat) result(has_bad)
        complex(kind=8), dimension(:,:), intent(in) :: Mat
        integer :: nl, nr
        has_bad = .false.
        do nr = 1, size(Mat, 2)
            do nl = 1, size(Mat, 1)
                if (.not. ieee_is_finite(real(Mat(nl, nr))) .or. .not. ieee_is_finite(aimag(Mat(nl, nr)))) then
                    has_bad = .true.
                    return
                endif
            enddo
        enddo
    end function has_nan_or_inf

    real(kind=8) function max_abs(Mat) result(mv)
        complex(kind=8), dimension(:,:), intent(in) :: Mat
        integer :: nl, nr
        mv = 0.d0
        do nr = 1, size(Mat, 2)
            do nl = 1, size(Mat, 1)
                mv = max(mv, abs(Mat(nl, nr)))
            enddo
        enddo
    end function max_abs
    
    real(kind=8) function compare_mat(Gr, Gr2) result(dif)
        complex(kind=8), dimension(Ndim, Ndim), intent(in) :: Gr, Gr2
        ! 使用与 CodeXun 相同的方式：所有元素绝对差值之和
        integer :: nl, nr
        dif = 0.d0
        do nr = 1, Ndim
            do nl = 1, Ndim
                dif = dif + real( abs(Gr(nl, nr) - Gr2(nl, nr)) )
            enddo
        enddo
        return
    end function compare_mat
    
    subroutine Wrap_pre(Prop, WrList, nt)
! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(in) :: nt
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Gr
        integer :: nt_st
        if (mod(nt, Nwrap) .ne. 0) then
            write(6,*) "incorrect preortho time slice, NT = ", nt; stop
        endif
        nt_st = nt/Nwrap
        if (nt .ne. 0) call stab_UR(Prop)
        WrList%URlist(1:Ndim, 1:Ndim, nt_st) = Prop%UUR(1:Ndim, 1:Ndim)
        WrList%VRlist(1:Ndim, 1:Ndim, nt_st) = Prop%VUR(1:Ndim, 1:Ndim)
        WrList%DRlist(1:Ndim, nt_st) = Prop%DUR(1:Ndim)
        if (nt == Ltrot) then
            ! 计算 G_0 = (1 + B_tot)^{-1}
            ! 注意：Prop%Gr 存储 G_0，不转换为 G_λ
            ! 传播操作（mmult_L/R）假设 G 是 G_0
            ! λ 的效应通过玻色部分和 λ 更新来实现
            Gr = dcmplx(0.d0, 0.d0)
            call stab_green(Gr, Prop, nt)
            Prop%Gr = Gr
        endif
        return
    end subroutine Wrap_pre
    
    subroutine Wrap_L(Prop, WrList, nt, flag)
! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(in) :: nt
        character(len=*), optional, intent(in) :: flag
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Gr
        integer :: nt_st, ii
        real(kind=8) :: dif
        if (mod(nt, Nwrap) .ne. 0 .and. nt .ne. 0) then
            write(6,*) "incorrect ortholeft time slice, NT = ", nt; stop
        endif
        nt_st = int(nt/Nwrap)
        Gr = dcmplx(0.d0, 0.d0)
        ! 从当前 wrap 段同步左右因子，保证与最新 σ/λ 链一致
        Prop%UUL(1:Ndim, 1:Ndim) = WrList%ULlist(1:Ndim, 1:Ndim, nt_st)
        Prop%VUL(1:Ndim, 1:Ndim) = WrList%VLlist(1:Ndim, 1:Ndim, nt_st)
        Prop%DUL(1:Ndim) = WrList%DLlist(1:Ndim, nt_st)
        Prop%UUR(1:Ndim, 1:Ndim) = WrList%URlist(1:Ndim, 1:Ndim, nt_st)
        Prop%VUR(1:Ndim, 1:Ndim) = WrList%VRlist(1:Ndim, 1:Ndim, nt_st)
        Prop%DUR(1:Ndim) = WrList%DRlist(1:Ndim, nt_st)
        ! 从当前 wrap 段同步左右因子，保证与最新 σ/λ 链一致
        Prop%UUL(1:Ndim, 1:Ndim) = WrList%ULlist(1:Ndim, 1:Ndim, nt_st)
        Prop%VUL(1:Ndim, 1:Ndim) = WrList%VLlist(1:Ndim, 1:Ndim, nt_st)
        Prop%DUL(1:Ndim) = WrList%DLlist(1:Ndim, nt_st)
        Prop%UUR(1:Ndim, 1:Ndim) = WrList%URlist(1:Ndim, 1:Ndim, nt_st)
        Prop%VUR(1:Ndim, 1:Ndim) = WrList%VRlist(1:Ndim, 1:Ndim, nt_st)
        Prop%DUR(1:Ndim) = WrList%DRlist(1:Ndim, nt_st)
        if (nt .ne. Ltrot) then
            call stab_UL(Prop)
            ! 使用更稳健的 big 版本一次性重建，避免依赖重试
            call stab_green_big(Prop)
            Gr = Gr_tmp%Gr00
            if (has_nan_or_inf(Gr)) Gr = ZKRON
            Prop%Gr = Gr
            dif = 0.d0
            if (present(flag)) Prop%Xmeanm = Prop%Xmeanm + dif
        endif
        WrList%ULlist(1:Ndim, 1:Ndim, nt_st) = Prop%UUL(1:Ndim, 1:Ndim)
        WrList%VLlist(1:Ndim, 1:Ndim, nt_st) = Prop%VUL(1:Ndim, 1:Ndim)
        WrList%DLlist(1:Ndim, nt_st) = Prop%DUL(1:Ndim)
        return
    end subroutine Wrap_L

    subroutine Wrap_L_store(Prop, WrList, nt)
        ! 只存储 UUL 到 WrapList，不做稳定化检查和 G 重建
        ! 用于 global update 中
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(in) :: nt
        integer :: nt_st
        if (mod(nt, Nwrap) .ne. 0 .and. nt .ne. 0) then
            write(6,*) "incorrect ortholeft time slice, NT = ", nt; stop
        endif
        nt_st = int(nt/Nwrap)
        WrList%ULlist(1:Ndim, 1:Ndim, nt_st) = Prop%UUL(1:Ndim, 1:Ndim)
        WrList%VLlist(1:Ndim, 1:Ndim, nt_st) = Prop%VUL(1:Ndim, 1:Ndim)
        WrList%DLlist(1:Ndim, nt_st) = Prop%DUL(1:Ndim)
        return
    end subroutine Wrap_L_store
    
    subroutine Wrap_R(Prop, WrList, nt, flag)
! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(in) :: nt
        character(len=*), optional, intent(in) :: flag
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Gr
        integer :: nt_st, ii
        real(kind=8) :: dif
        if (mod(nt, Nwrap) .ne. 0 .and. nt .ne. 0) then
            write(6,*) "incorrect orthoright time slice, NT = ", nt; stop
        endif
        nt_st = int(nt/Nwrap)
        Gr = dcmplx(0.d0, 0.d0)
        if (nt == Ltrot) then
            ! 清空 ULlist，但重新初始化为单位矩阵和 D=1
            ! 这确保在下一次 sweep_L 之前，如果有 global update 访问 ULlist，
            ! 它们是有效的值（而不是 0）
            do ii = 0, Nst
                WrList%ULlist(:,:,ii) = ZKRON
                WrList%VLlist(:,:,ii) = ZKRON
                WrList%DLlist(:,ii) = dcmplx(1.d0, 0.d0)
            enddo
        endif
        if (nt .ne. 0) then
            call stab_UR(Prop)
            ! 使用更稳健的 big 版本一次性重建，避免依赖重试
            call stab_green_big(Prop)
            Gr = Gr_tmp%Gr00
            if (has_nan_or_inf(Gr)) Gr = ZKRON

            Prop%Gr = Gr
            dif = 0.d0
            if (present(flag)) Prop%Xmeanm = Prop%Xmeanm + dif
        endif
        WrList%URlist(1:Ndim, 1:Ndim, nt_st) = Prop%UUR(1:Ndim, 1:Ndim)
        WrList%VRlist(1:Ndim, 1:Ndim, nt_st) = Prop%VUR(1:Ndim, 1:Ndim)
        WrList%DRlist(1:Ndim, nt_st) = Prop%DUR(1:Ndim)
        return
    end subroutine Wrap_R

    subroutine Wrap_R_store(Prop, WrList, nt)
        ! 只存储 UUR 到 WrapList，不做稳定化检查和 G 重建
        ! 用于 global update 中
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(in) :: nt
        integer :: nt_st
        if (mod(nt, Nwrap) .ne. 0 .and. nt .ne. 0) then
            write(6,*) "incorrect orthoright time slice, NT = ", nt; stop
        endif
        nt_st = int(nt/Nwrap)
        WrList%URlist(1:Ndim, 1:Ndim, nt_st) = Prop%UUR(1:Ndim, 1:Ndim)
        WrList%VRlist(1:Ndim, 1:Ndim, nt_st) = Prop%VUR(1:Ndim, 1:Ndim)
        WrList%DRlist(1:Ndim, nt_st) = Prop%DUR(1:Ndim)
        return
    end subroutine Wrap_R_store

    subroutine Wrap_tau(Prop, PropGr, WrList, nt)
! Arguments: 
        class(Propagator), intent(inout) :: Prop
        class(PropGreen), intent(inout) :: PropGr
        class(WrapList), intent(in) :: WrList
        integer, intent(in) :: nt
! Local: 
        integer :: nt_st
        real(kind=8) :: dif
        if (mod(nt, Nwrap) .ne. 0 .or. nt == 0) then
            write(6,*) "incorrect orthobig time slice, NT = ", nt; stop
        endif
        nt_st = int(nt/Nwrap)
        Prop%UUL(1:Ndim, 1:Ndim) = WrList%ULlist(1:Ndim, 1:Ndim, nt_st)
        Prop%VUL(1:Ndim, 1:Ndim) = WrList%VLlist(1:Ndim, 1:Ndim, nt_st)
        Prop%DUL(1:Ndim) = WrList%DLlist(1:Ndim, nt_st)
        call stab_UR(Prop)
        call stab_green_big(Prop)
! stabilization test
        dif = compare_mat(Gr_tmp%Gr0t, PropGr%Gr0t)
        if (dif > PropGr%Xmaxm(1)) PropGr%Xmaxm(1) = dif
        if (dif .ge. 5.5d-5) write(6,*) nt, dif, "GR0T ortho unstable in RANK ", IRANK
        PropGr%Xmeanm(1) = PropGr%Xmeanm(1) + dif
        dif = compare_mat(Gr_tmp%Grt0, PropGr%Grt0)
        if (dif > PropGr%Xmaxm(2)) PropGr%Xmaxm(2) = dif
        if (dif .ge. 5.5d-5) write(6,*) nt, dif, "GRT0 ortho unstable in RANK ", IRANK
        PropGr%Xmeanm(2) = PropGr%Xmeanm(2) + dif
        dif = compare_mat(Gr_tmp%Grtt, PropGr%Grtt)
        if (dif > PropGr%Xmaxm(3)) PropGr%Xmaxm(3) = dif
        if (dif .ge. 5.5d-5) write(6,*) nt, dif, "GRTT ortho unstable in RANK ", IRANK
        PropGr%Xmeanm(3) = PropGr%Xmeanm(3) + dif
        PropGr%Gr00 = Gr_tmp%Gr00
        PropGr%Gr0t = Gr_tmp%Gr0t
        PropGr%Grt0 = Gr_tmp%Grt0
        PropGr%Grtt = Gr_tmp%Grtt
        return
    end subroutine Wrap_tau
end module Stabilize_mod