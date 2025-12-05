module GlobalUpdate_mod ! a combination of shift&Wolff update
    use GlobalK_mod
    use Stabilize_mod, only: Wrap_pre, Wrap_L, Wrap_L_store, Wrap_R, Wrap_R_store, Wrap_tau, Stabilize_init, Stabilize_clear, stab_UR, stab_UL, stab_green_big_out
    use Multiply_mod, only: propK_pre, propT_pre
    use MyLattice, only: ZKRON
    use Fields_mod, only: gauss_boson_ratio_sigma
    use DQMC_Model_mod
    use calc_basic
    implicit none
    logical, parameter :: DEBUG_WOLFF = .false.
    
    public
    private :: Wolff, sigma_new
    
    type WolffContainer
        logical, dimension(:), allocatable :: tau_mask
        integer, dimension(:), allocatable :: queue
        real(kind=8) :: size_cluster
    contains
        procedure :: init => Wolff_init
        procedure :: reset => Wolff_reset
        procedure :: flip  => Wolff_flip
        final :: Wolff_clear
    end type WolffContainer
    
    type :: GlobalUpdate
        type(Propagator), allocatable, private :: propU, propD
        type(WrapList),   allocatable, private :: wrU, wrD
    contains
        procedure :: init  => Global_init
        procedure :: clear => Global_clear
        procedure, private :: reset => Global_reset
        procedure, private :: flip  => Global_flip
        procedure, private :: sweep_L => Global_sweep_L
        procedure, private :: sweep_R => Global_sweep_R
        procedure :: sweep => Global_sweep
        procedure, nopass :: ctrl_print => Global_control_print
    end type GlobalUpdate
    
    type(AccCounter) :: Acc_Kg
    type(WolffContainer), allocatable :: Wolff
    real(kind=8), dimension(:,:), allocatable :: sigma_new
    
contains
    subroutine Wolff_init(this)
        class(WolffContainer), intent(inout) :: this
        allocate(this%tau_mask(Ltrot))
        this%tau_mask = .false.
        allocate(this%queue(Ltrot))
        this%queue = 0
        this%size_cluster = 0.d0
        return
    end subroutine Wolff_init
    
    subroutine Wolff_reset(this)
        class(WolffContainer), intent(inout) :: this
        if (allocated(this%tau_mask)) this%tau_mask = .false.
        return
    end subroutine Wolff_reset
    
    subroutine Wolff_clear(this)
        type(WolffContainer), intent(inout) :: this
        if (allocated(this%tau_mask)) deallocate(this%tau_mask)
        if (allocated(this%queue)) deallocate(this%queue)
        return
    end subroutine Wolff_clear
    
    subroutine Wolff_flip(this, iseed, size_cluster, log_ratio_space)
        class(WolffContainer), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer, intent(out) :: size_cluster
        real(kind=8), intent(out) :: log_ratio_space
        real(kind=8) :: p_time
        integer :: bond_idx
        bond_idx = nranf(iseed, 2*Lq)
        if (bond_idx < 1) bond_idx = 1

        p_time = compute_time_bond_prob()
        call build_time_cluster(bond_idx, this%tau_mask, this%queue, iseed, p_time)
        size_cluster = count(this%tau_mask)
        if (size_cluster == 0) then
            log_ratio_space = 0.d0
            return
        endif

        log_ratio_space = compute_space_log_ratio(bond_idx, this%tau_mask)
        call apply_cluster_flip(bond_idx, this%tau_mask)
        return
    end subroutine Wolff_flip

    real(kind=8) function compute_time_bond_prob()
        real(kind=8) :: tanh_arg, gamma
        if (abs(h) > Zero) then
            tanh_arg = tanh(Dtau * abs(h))
            tanh_arg = max(tanh_arg, RATIO_EPS)
            gamma = -0.5d0 * log(tanh_arg)
            compute_time_bond_prob = 1.d0 - exp(-2.d0 * gamma)
        else
            compute_time_bond_prob = 0.d0
        endif
        compute_time_bond_prob = max(0.d0, min(1.d0, compute_time_bond_prob))
        return
    end function compute_time_bond_prob

    subroutine build_time_cluster(bond_idx, mask, queue, iseed, p_time)
        integer, intent(in) :: bond_idx
        logical, intent(inout) :: mask(:)
        integer, intent(inout) :: queue(:)
        integer, intent(inout) :: iseed
        real(kind=8), intent(in) :: p_time
        real(kind=8), external :: ranf
        integer :: head, tail, current, neighbor, offset
        mask = .false.
        if (Ltrot <= 0) return
        queue = 0
        queue(1) = nranf(iseed, Ltrot)
        if (queue(1) < 1) queue(1) = 1
        head = 1; tail = 1
        mask(queue(1)) = .true.
        do while (head <= tail)
            current = queue(head)
            head = head + 1
            do offset = -1, 1, 2
                neighbor = npbc(current + offset, Ltrot)
                if (mask(neighbor)) cycle
                if (sigma_new(bond_idx, neighbor) == sigma_new(bond_idx, current)) then
                    if (p_time > 0.d0 .and. ranf(iseed) < p_time) then
                        tail = tail + 1
                        queue(tail) = neighbor
                        mask(neighbor) = .true.
                    endif
                endif
            enddo
        enddo
        return
    end subroutine build_time_cluster

    real(kind=8) function compute_space_log_ratio(bond_idx, mask)
        integer, intent(in) :: bond_idx
        logical, intent(in) :: mask(:)
        real(kind=8) :: alpha, gauss_ratio
        integer :: tau
        compute_space_log_ratio = 0.d0
        alpha = Dtau * J
        
        do tau = 1, Ltrot
            if (.not. mask(tau)) cycle
            ! Plaquette 贡献
            if (alpha /= 0.d0) then
                compute_space_log_ratio = compute_space_log_ratio + plaquette_log_ratio(bond_idx, tau, 1, alpha)
                compute_space_log_ratio = compute_space_log_ratio + plaquette_log_ratio(bond_idx, tau, 4, alpha)
            endif
            ! 高斯投影玻色贡献
            ! 注意：gauss_boson_ratio_sigma 返回 1 或 -1
            ! 如果是 -1，log 未定义，我们使用 abs() 处理（与 local update 一致）
            gauss_ratio = gauss_boson_ratio_sigma(bond_idx, tau, sigma_new, &
                -sigma_new(bond_idx, tau), NsigL_K%lambda, Latt)
            ! 使用绝对值，符号问题通过 reweighting 处理
            compute_space_log_ratio = compute_space_log_ratio + log(abs(gauss_ratio) + 1.d-300)
        enddo
        return
    end function compute_space_log_ratio

    real(kind=8) function plaquette_log_ratio(bond_idx, tau, start_idx, alpha)
        integer, intent(in) :: bond_idx, tau, start_idx
        real(kind=8), intent(in) :: alpha
        real(kind=8) :: prod
        integer :: k, neighbor
        ! 计算格点格子的旧构型值（sigma_new此时还等于旧构型）
        prod = sigma_new(bond_idx, tau)
        do k = 0, 2
            neighbor = Latt%bond_bonds(bond_idx, start_idx + k)
            if (neighbor <= 0) then
                plaquette_log_ratio = 0.d0
                return
            endif
            prod = prod * sigma_new(neighbor, tau)
        enddo
        ! Z2规范场作用量：S = -J*Dtau * sum(plaquette)
        ! 翻转一个键后：plaquette_new = -plaquette_old
        ! Delta_S = S_new - S_old = -alpha * (-prod - prod) = 2*alpha*prod
        ! 但Metropolis接受率是 exp(-Delta_S)，所以log_ratio = -Delta_S
        plaquette_log_ratio = -2.d0 * alpha * prod
        return
    end function plaquette_log_ratio

    subroutine apply_cluster_flip(bond_idx, mask)
        integer, intent(in) :: bond_idx
        logical, intent(in) :: mask(:)
        integer :: tau
        do tau = 1, Ltrot
            if (mask(tau)) sigma_new(bond_idx, tau) = -sigma_new(bond_idx, tau)
        enddo
    end subroutine apply_cluster_flip

    subroutine Global_init(this)
        class(GlobalUpdate), intent(inout) :: this
        allocate(sigma_new(2*Lq, Ltrot))
        sigma_new = 0.d0
        allocate(Wolff)
        call Wolff%init()
        call Acc_Kg%init()
        allocate(this%propU); call this%propU%make()
        allocate(this%propD); call this%propD%make()
        allocate(this%wrU);   call this%wrU%make()
        allocate(this%wrD);   call this%wrD%make()
        return
    end subroutine Global_init

    subroutine Global_clear(this)
        class(GlobalUpdate), intent(inout) :: this
        deallocate(this%propU); deallocate(this%propD)
        deallocate(this%wrU);   deallocate(this%wrD)
        deallocate(Wolff)
        deallocate(sigma_new)
        return
    end subroutine Global_clear
    
    subroutine Global_reset(this, PropU, PropD, WrU, WrD)
        class(GlobalUpdate), intent(inout) :: this
        class(Propagator),  intent(in)    :: PropU, PropD
        class(WrapList),    intent(in)    :: WrU, WrD
        call this%propU%asgn(PropU)
        call this%propD%asgn(PropD)
        call this%wrU%asgn(WrU)
        call this%wrD%asgn(WrD)
        call Wolff%reset()
        sigma_new = NsigL_K%sigma
        return
    end subroutine Global_reset
    
    subroutine Global_flip(this, iseed, size_cluster, log_ratio_space)
        class(GlobalUpdate), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer, intent(out) :: size_cluster
        real(kind=8), intent(out) :: log_ratio_space
        call Wolff%reset()
        call Wolff%flip(iseed, size_cluster, log_ratio_space) ! update sigma_new and size_cluster
        return
    end subroutine Global_flip
    
    subroutine Global_sweep_L(this, PropU, PropD, WrU, WrD, iseed, size_cluster, is_beta)
        class(GlobalUpdate), intent(inout) :: this
        class(Propagator), intent(inout) :: PropU, PropD
        class(WrapList),   intent(inout) :: WrU, WrD
        integer, intent(inout) :: iseed
        integer, intent(out) :: size_cluster
        logical, intent(inout) :: is_beta
! Local:
        real(kind=8), external :: ranf
        real(kind=8) :: log_ratio_fermion, log_ratio_space, log_ratio_total, random
        real(kind=8), dimension(2*Lq, Ltrot) :: sigma_curr
        integer :: nt
        
        log_ratio_fermion = 0.d0
        sigma_curr = NsigL_K%sigma
        call this%flip(iseed, size_cluster, log_ratio_space)
        ! 用提出的 sigma_new 重建内部稳定化链，避免旧 ULlist/URlist 失配
        call rebuild_stabilization_chain(this%propU, this%propD, this%wrU, this%wrD, sigma_new)
        sigma_curr = sigma_new
        
        ! 按照 CodeXun 的顺序：先 Wrap_L，再 propT_L，再 GlobalK_prop_L
        ! 这样 Wrap_L 会在每个 Nwrap 间隔重建 Green 函数，消除累积误差
        do nt = Ltrot, 1, -1
            if (mod(nt, Nwrap) == 0) then
                call Wrap_L(this%propU, this%wrU, nt)
                call Wrap_L(this%propD, this%wrD, nt)
            endif
            call propT_L(this%propU, this%propD, NsigL_K%lambda, nt)
            call GlobalK_prop_L(this%propU, this%propD, log_ratio_fermion, sigma_new, sigma_curr, nt)
            if (mod(nt, Nwrap) == 0) then
                call Wrap_L_store(this%propU, this%wrU, nt)
                call Wrap_L_store(this%propD, this%wrD, nt)
                call rebuild_stabilization_chain_factors_only(this%propU, this%propD, this%wrU, this%wrD, sigma_curr)
                call stab_UL(this%propU)
                call stab_green_big_out(this%propU, this%propU%Gr)
                call stab_UL(this%propD)
                call stab_green_big_out(this%propD, this%propD%Gr)
            endif
        enddo
        call Wrap_L(this%propU, this%wrU, 0)
        call Wrap_L(this%propD, this%wrD, 0)
        log_ratio_total = log_ratio_fermion + log_ratio_space
        random = ranf(iseed)
        if (min(0.d0, log_ratio_total) .ge. log(max(random, RATIO_EPS))) then
            call Acc_Kg%count(.true.)
            NsigL_K%sigma = sigma_new
            ! 接受：复制 prop 和 wrlist
            call PropU%asgn(this%propU); call PropD%asgn(this%propD)
            call WrU%asgn(this%wrU);     call WrD%asgn(this%wrD)
            is_beta = .false.
        else
            call Acc_Kg%count(.false.)
            sigma_new = NsigL_K%sigma
            ! 拒绝：恢复原始 prop 和 wrlist
            call this%propU%asgn(PropU); call this%propD%asgn(PropD)
            call this%wrU%asgn(WrU);     call this%wrD%asgn(WrD)
            is_beta = .true. 
        endif
        return
    end subroutine Global_sweep_L
    
    subroutine Global_sweep_R(this, PropU, PropD, WrU, WrD, iseed, size_cluster, is_beta)
        class(GlobalUpdate), intent(inout) :: this
        class(Propagator), intent(inout) :: PropU, PropD
        class(WrapList),   intent(inout) :: WrU, WrD
        integer, intent(inout) :: iseed
        integer, intent(out) :: size_cluster
        logical, intent(inout) :: is_beta
! Local:
        real(kind=8), external :: ranf
        real(kind=8) :: log_ratio_fermion, log_ratio_space, log_ratio_total, random
        real(kind=8), dimension(2*Lq, Ltrot) :: sigma_curr
        integer :: nt
        
        log_ratio_fermion = 0.d0
        sigma_curr = NsigL_K%sigma
        call this%flip(iseed, size_cluster, log_ratio_space)
        ! 用提出的 sigma_new 重建内部稳定化链，避免旧 ULlist/URlist 失配
        call rebuild_stabilization_chain(this%propU, this%propD, this%wrU, this%wrD, sigma_new)
        sigma_curr = sigma_new
        call Wrap_R(this%propU, this%wrU, 0)
        call Wrap_R(this%propD, this%wrD, 0)
        
        ! 按照 CodeXun 的顺序：先 GlobalK_prop_R，再 propT_R，最后 Wrap_R
        ! 这样 Wrap_R 会在每个 Nwrap 间隔重建 Green 函数，消除累积误差
        do nt = 1, Ltrot
            call GlobalK_prop_R(this%propU, this%propD, log_ratio_fermion, sigma_new, sigma_curr, nt)
            call propT_R(this%propU, this%propD, NsigL_K%lambda, nt)
            if (mod(nt, Nwrap) == 0) then
                call Wrap_R(this%propU, this%wrU, nt)
                call Wrap_R(this%propD, this%wrD, nt)
                call Wrap_R_store(this%propU, this%wrU, nt)
                call Wrap_R_store(this%propD, this%wrD, nt)
                ! 同步左向因子到当前 σ：确保后续 stab_green 使用的 ULlist 与 sigma_curr 一致
                call Wrap_L_store(this%propU, this%wrU, nt)
                call Wrap_L_store(this%propD, this%wrD, nt)
                call rebuild_stabilization_chain_factors_only(this%propU, this%propD, this%wrU, this%wrD, sigma_curr)
                call stab_UR(this%propU)
                call stab_green_big_out(this%propU, this%propU%Gr)
                call stab_UR(this%propD)
                call stab_green_big_out(this%propD, this%propD%Gr)
            endif
        enddo
        log_ratio_total = log_ratio_fermion + log_ratio_space
        random = ranf(iseed)
        if (min(0.d0, log_ratio_total) .ge. log(max(random, RATIO_EPS))) then
            call Acc_Kg%count(.true.)
            NsigL_K%sigma = sigma_new
            ! 接受：复制 prop 和 wrlist
            call PropU%asgn(this%propU); call PropD%asgn(this%propD)
            call WrU%asgn(this%wrU);     call WrD%asgn(this%wrD)
            is_beta = .true.
        else
            call Acc_Kg%count(.false.)
            sigma_new = NsigL_K%sigma
            ! 拒绝：恢复原始 prop 和 wrlist
            call this%propU%asgn(PropU); call this%propD%asgn(PropD)
            call this%wrU%asgn(WrU);     call this%wrD%asgn(WrD)
            is_beta = .false. 
        endif
        return
    end subroutine Global_sweep_R
    
    subroutine Global_sweep(this, PropU, PropD, WrU, WrD, iseed, is_beta)
        class(GlobalUpdate), intent(inout) :: this
        class(Propagator), intent(inout) :: PropU, PropD
        class(WrapList),   intent(inout) :: WrU, WrD
        integer, intent(inout) :: iseed
        logical, intent(inout) :: is_beta
        integer :: ngl, size_count
        real(kind=8) :: size_cluster
        
        call this%reset(PropU, PropD, WrU, WrD)
        call Acc_Kg%reset()
        size_cluster = 0.d0
        do ngl = 1, Nglobal
            if (is_beta) then
                call this%sweep_L(PropU, PropD, WrU, WrD, iseed, size_count, is_beta)
            else
                call this%sweep_R(PropU, PropD, WrU, WrD, iseed, size_count, is_beta)
            endif
            size_cluster = size_cluster + dble(size_count)
        enddo
        call Acc_Kg%ratio()
        size_cluster = size_cluster / dble(Nglobal)
        Wolff%size_cluster = Wolff%size_cluster + size_cluster
        return
    end subroutine Global_sweep
    
    subroutine rebuild_stabilization_chain(PropU, PropD, WrU, WrD, sigma_in)
        ! 重新准备稳定化链
        ! 在 global update 被接受后调用，确保 WrList 中的 UDV 分解与当前 sigma 一致
        ! 完全按照 Local_sweep_pre 的方式进行初始化
        class(Propagator), intent(inout) :: PropU, PropD
        class(WrapList),   intent(inout) :: WrU, WrD
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_in
        integer :: nt
        
        ! 重置 Propagator 的 UDV 分解（与 Prop_make 一致：UUR/VUR 为单位矩阵，DUR 为 1）
        PropU%UUR = ZKRON; PropU%VUR = ZKRON; PropU%DUR = dcmplx(1.d0, 0.d0)
        PropU%UUL = ZKRON; PropU%VUL = ZKRON; PropU%DUL = dcmplx(1.d0, 0.d0)
        PropU%Gr = ZKRON  ! 重置 Green 函数
        PropD%UUR = ZKRON; PropD%VUR = ZKRON; PropD%DUR = dcmplx(1.d0, 0.d0)
        PropD%UUL = ZKRON; PropD%VUL = ZKRON; PropD%DUL = dcmplx(1.d0, 0.d0)
        PropD%Gr = ZKRON  ! 重置 Green 函数
        
        ! 清除 WrapList（与 Wrlist_make 一致）
        WrU%URlist = dcmplx(0.d0, 0.d0)
        WrU%VRlist = dcmplx(0.d0, 0.d0)
        WrU%DRlist = dcmplx(0.d0, 0.d0)
        WrD%URlist = dcmplx(0.d0, 0.d0)
        WrD%VRlist = dcmplx(0.d0, 0.d0)
        WrD%DRlist = dcmplx(0.d0, 0.d0)
        ! 初始化 ULlist 为单位矩阵，DLlist 为 1
        ! 这确保在 stab_green 中 matUDV 不会是奇异的
        do nt = 0, Ltrot/Nwrap
            WrU%ULlist(:,:,nt) = ZKRON
            WrU%VLlist(:,:,nt) = ZKRON
            WrU%DLlist(:,nt) = dcmplx(1.d0, 0.d0)
            WrD%ULlist(:,:,nt) = ZKRON
            WrD%VLlist(:,:,nt) = ZKRON
            WrD%DLlist(:,nt) = dcmplx(1.d0, 0.d0)
        enddo
        
        ! 重新准备稳定化链（与 Local_sweep_pre 完全一致）
        call Wrap_pre(PropU, WrU, 0)
        call Wrap_pre(PropD, WrD, 0)
        WrU%ULlist(:,:,0) = PropU%UUL
        WrU%VLlist(:,:,0) = PropU%VUL
        WrU%DLlist(:,0)   = PropU%DUL
        WrD%ULlist(:,:,0) = PropD%UUL
        WrD%VLlist(:,:,0) = PropD%VUL
        WrD%DLlist(:,0)   = PropD%DUL
        do nt = 1, Ltrot
            if (RT > Zero) call propK_pre(PropU, PropD, sigma_in, nt)
            call propT_pre(PropU, PropD, NsigL_K%lambda, nt)
            if (mod(nt, Nwrap) == 0) then
                call Wrap_pre(PropU, WrU, nt)
                call Wrap_pre(PropD, WrD, nt)
                WrU%ULlist(:,:,nt/Nwrap) = PropU%UUL
                WrU%VLlist(:,:,nt/Nwrap) = PropU%VUL
                WrU%DLlist(:,nt/Nwrap)   = PropU%DUL
                WrD%ULlist(:,:,nt/Nwrap) = PropD%UUL
                WrD%VLlist(:,:,nt/Nwrap) = PropD%VUL
                WrD%DLlist(:,nt/Nwrap)   = PropD%DUL
            endif
        enddo
        return
    end subroutine rebuild_stabilization_chain

    subroutine rebuild_stabilization_chain_factors_only(PropU, PropD, WrU, WrD, sigma_in)
        ! 重建因子和 WrList，但保持当前 G 不变
        class(Propagator), intent(inout) :: PropU, PropD
        class(WrapList),   intent(inout) :: WrU, WrD
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_in
        complex(kind=8), dimension(Ndim, Ndim) :: GrU_keep, GrD_keep
        GrU_keep = PropU%Gr
        GrD_keep = PropD%Gr
        call rebuild_stabilization_chain(PropU, PropD, WrU, WrD, sigma_in)
        PropU%Gr = GrU_keep
        PropD%Gr = GrD_keep
        return
    end subroutine rebuild_stabilization_chain_factors_only

    subroutine Global_control_print()
        include 'mpif.h'
        real(kind=8) :: collect
        collect = 0.d0
        call MPI_Reduce(Acc_Kg%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Acc_Kg%acc = collect / dble(ISIZE * Nbin)
        collect = 0
        call MPI_Reduce(Wolff%size_cluster, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Wolff%size_cluster = collect / dble(ISIZE * Nbin)
        if (IRANK == 0) then
            write(50,*) 'Average size of Wolff cluster                  :', Wolff%size_cluster
            write(50,*) 'Accept_Kglobal                                 :', Acc_Kg%acc
            call flush(50)
        endif
        return
    end subroutine Global_control_print

end module GlobalUpdate_mod
