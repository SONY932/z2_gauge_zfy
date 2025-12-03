module Fields_mod
    use Mylattice
    use calc_basic
    implicit none

    public :: SigmaConf, conf_in, conf_out
    private

    type :: SigmaConf
        real(kind=8), dimension(:,:), public, allocatable :: sigma ! 这里是规范场，sigma(ii,nt)
        real(kind=8), dimension(:),   public, allocatable :: lambda ! 这是高斯规范(Lq), 末片轴规 λ_i ∈ {±1}
    contains
        procedure :: make => SigmaConf_make
        final :: SigmaConf_final
        procedure :: m_bosonratio_local
        procedure :: m_bosonratio_global_spacetime
        generic :: bosonratio => m_bosonratio_local
        generic :: bosonratio => m_bosonratio_global_spacetime
    end type SigmaConf

contains
! 初始化规范场
    subroutine SigmaConf_make(this)
        class(SigmaConf), intent(inout) :: this
        allocate( this%sigma(2*Lq, Ltrot) )
        allocate( this%lambda(Lq) )
        this%lambda = 1.d0
        return
    end subroutine SigmaConf_make

! 释放规范场
    subroutine SigmaConf_final(this)
        type(SigmaConf), intent(inout) :: this
        if (allocated(this%sigma))  deallocate(this%sigma)
        if (allocated(this%lambda)) deallocate(this%lambda)
        return
    end subroutine SigmaConf_final
    
! 计算规范场的更新 ratio
    real(kind=8) function m_bosonratio_local(this, sigma_new, ii, ntau, Latt) result(ratio_boson)
! Arguments: 
        class(SigmaConf), intent(in) :: this
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_new
        integer, intent(in) :: ii, ntau
        class(SquareLattice), intent(in) :: Latt
! Local: 
        real(kind=8) :: alpha, gamma
        real(kind=8) :: tmp_new, tmp_old, tmp
        real(kind=8) :: tanh_arg, log_upbound, expo, expo_space, expo_time
        logical :: use_temporal
        integer :: jj, nf
! 首先定义常数系数（h -> 0 的鲁棒极限处理）
        alpha = Dtau * J
        if (abs(h) > Zero) then
            tanh_arg = tanh(Dtau * h)
            if (tanh_arg < 1.d-300) tanh_arg = 1.d-300
            gamma = - 0.5d0 * log(tanh_arg)
            use_temporal = .true.
        else
            gamma = 0.d0
            use_temporal = .false.
        endif
        log_upbound = log(upbound)

! 首先计算 Pauli Z方向上的权重比值
        tmp_new = 0.d0; tmp_old = 0.d0

        tmp = sigma_new(ii, ntau)
        do nf = 1, 3
            jj = Latt%bond_bonds(ii, nf)
            tmp = tmp  * this%sigma(jj, ntau)
        enddo
        tmp_new = tmp_new + tmp
        tmp = sigma_new(ii, ntau)
        do nf = 4, 6
            jj = Latt%bond_bonds(ii, nf)
            tmp = tmp  * this%sigma(jj, ntau)
        enddo
        tmp_new = tmp_new + tmp

        tmp = this%sigma(ii, ntau)
        do nf = 1, 3
            jj = Latt%bond_bonds(ii, nf)
            tmp = tmp  * this%sigma(jj, ntau)
        enddo
        tmp_old = tmp_old + tmp
        tmp = this%sigma(ii, ntau)
        do nf = 4, 6
            jj = Latt%bond_bonds(ii, nf)
            tmp = tmp  * this%sigma(jj, ntau)
        enddo
        tmp_old = tmp_old + tmp

        ! Z2规范场作用量：S = -J * sum(plaquette)
        ! Delta_S = S_new - S_old = -alpha * (tmp_new - tmp_old)
        ! 接受率：ratio = exp(-Delta_S) = exp(alpha * (tmp_new - tmp_old))
        expo = alpha * ( tmp_new - tmp_old )
        if (expo >  log_upbound) expo =  log_upbound
        if (expo < -log_upbound) expo = -log_upbound
        expo_space = expo
        ratio_boson = exp(expo)

! 然后计算 Pauli X方向上的权重比值（仅当 h 非零时考虑时间耦合）
        tmp_new = 0.d0; tmp_old = 0.d0
        if (use_temporal) then
            do nf = -1, 1, 2 ! 遍历ii键对应的虚时间最近邻时间片
                tmp_new = tmp_new + sigma_new(ii, ntau) * this%sigma(ii, npbc(ntau + nf, Ltrot))
                tmp_old = tmp_old + this%sigma(ii, ntau) * this%sigma(ii, npbc(ntau + nf, Ltrot))
            enddo
            expo = gamma * ( tmp_new - tmp_old )
            if (expo >  log_upbound) expo =  log_upbound
            if (expo < -log_upbound) expo = -log_upbound
            expo_time = expo
            ratio_boson = ratio_boson * exp(expo)
        else
            expo_time = 0.d0
        endif
        return
    end function m_bosonratio_local

    real(kind=8) function m_bosonratio_global_spacetime(this, sigma_new, Latt) result(ratio_boson)
! Arguments: 
        class(SigmaConf), intent(in) :: this
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_new
        class(SquareLattice), intent(in) :: Latt
! Local: 
        real(kind=8) :: alpha, gamma
        real(kind=8) :: new_spatial, old_spatial, tmp_new, tmp_old, tmp
        real(kind=8) :: new_temporal, old_temporal
        real(kind=8) :: tanh_arg, log_upbound, expo
        logical :: use_temporal
        integer :: nt, ii, nf, jj
        
        alpha = Dtau * J
        if (abs(h) > Zero) then
            tanh_arg = tanh(Dtau * h)
            if (tanh_arg < 1.d-300) tanh_arg = 1.d-300
            gamma = - 0.5d0 * log(tanh_arg)
            use_temporal = .true.
        else
            gamma = 0.d0
            use_temporal = .false.
        endif
        log_upbound = log(upbound)

! 空间项：遍历所有键，每个键通过6个近邻键计算plaquette贡献（每个plaquette被计数4次，需要除以4）
        new_spatial = 0.d0; old_spatial = 0.d0
        do nt = 1, Ltrot
            do ii = 1, 2*Lq
                tmp_new = 0.d0; tmp_old = 0.d0
                
                ! 计算上半部分的三个plaquette
                tmp = sigma_new(ii, nt)
                do nf = 1, 3
                    jj = Latt%bond_bonds(ii, nf)
                    tmp = tmp * sigma_new(jj, nt)
                enddo
                tmp_new = tmp_new + tmp
                
                tmp = this%sigma(ii, nt)
                do nf = 1, 3
                    jj = Latt%bond_bonds(ii, nf)
                    tmp = tmp * this%sigma(jj, nt)
                enddo
                tmp_old = tmp_old + tmp
                
                ! 计算下半部分的三个plaquette
                tmp = sigma_new(ii, nt)
                do nf = 4, 6
                    jj = Latt%bond_bonds(ii, nf)
                    tmp = tmp * sigma_new(jj, nt)
                enddo
                tmp_new = tmp_new + tmp
                
                tmp = this%sigma(ii, nt)
                do nf = 4, 6
                    jj = Latt%bond_bonds(ii, nf)
                    tmp = tmp * this%sigma(jj, nt)
                enddo
                tmp_old = tmp_old + tmp
                
                new_spatial = new_spatial + tmp_new
                old_spatial = old_spatial + tmp_old
            enddo
        enddo
        ! 每个plaquette被4个键重复计数，需要除以4
        new_spatial = new_spatial / 4.d0
        old_spatial = old_spatial / 4.d0

! 时间项（h=0 时跳过以得到正确极限）
        new_temporal = 0.d0; old_temporal = 0.d0
        if (use_temporal) then
            do nt = 1, Ltrot
                do ii = 1, 2*Lq
                    new_temporal = new_temporal + sigma_new(ii, nt) * sigma_new(ii, npbc(nt + 1, Ltrot))
                    old_temporal = old_temporal + this%sigma(ii, nt) * this%sigma(ii, npbc(nt + 1, Ltrot))
                enddo
            enddo
        endif

        expo = alpha * ( new_spatial - old_spatial )
        if (expo >  log_upbound) expo =  log_upbound
        if (expo < -log_upbound) expo = -log_upbound
        ratio_boson = exp(expo)
        if (use_temporal) then
            expo = gamma * ( new_temporal - old_temporal )
            if (expo >  log_upbound) expo =  log_upbound
            if (expo < -log_upbound) expo = -log_upbound
            ratio_boson = ratio_boson * exp(expo)
        endif
        return
    end function m_bosonratio_global_spacetime

    subroutine conf_in(Conf, iseed, Latt)
        class(SigmaConf), intent(inout) :: Conf
        integer, intent(out) :: iseed
        class(SquareLattice), intent(in) :: Latt
! Local: 
        real(kind=8), dimension(2*LqTherm, LtrotTherm) :: sigma_therm
        real(kind=8), dimension(LqTherm) :: lambda_therm
        
        call conf_therm_in(sigma_therm, lambda_therm, iseed)
        call conf_transfer(Conf%sigma, Conf%lambda, sigma_therm, lambda_therm, Latt)
        return
    end subroutine conf_in

    subroutine conf_therm_in(sigma_therm, lambda_therm, iseed)
!#define DEC
        include 'mpif.h'
! Arguments: 
        real(kind=8), dimension(2*LqTherm, LtrotTherm), intent(inout) :: sigma_therm
        real(kind=8), dimension(LqTherm), intent(inout) :: lambda_therm
        integer, intent(out) :: iseed
! Local: 
        integer :: status(MPI_STATUS_SIZE)
        integer :: iseed0, itmp
        integer :: ii, nt, N
        real(kind=8), external :: ranf
        real(kind=8), dimension(2*LqTherm, LtrotTherm) :: sigma_itmp
        real(kind=8), dimension(LqTherm) :: lambda_itmp
        integer :: ios_lambda
        
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
                    call conf_set(sigma_itmp, lambda_itmp, itmp)
                    call MPI_SEND(itmp, 1, MPI_Integer, N, N, MPI_COMM_WORLD, IERR)
                    call MPI_SEND(sigma_itmp, 2*LqTherm*LtrotTherm, MPI_Real8, N, N+1024, MPI_COMM_WORLD, IERR)
                    call MPI_SEND(lambda_itmp, LqTherm, MPI_Real8, N, N+2048, MPI_COMM_WORLD, IERR)
                    print '("after sending message from master ", i3.1, " to Rank ", i3.1)', IRANK, N
                enddo
! Set node zero.
                iseed = iseed0
                call conf_set(sigma_therm, lambda_therm, iseed)
                print '("initiating master rank ", i3.1)', IRANK
            else
! read all confins from NODE 0.
! Setup Node 0
                do nt = 1, LtrotTherm
                    do ii = 1, 2*LqTherm
                        read(30,*) sigma_therm(ii, nt)
                    enddo
                enddo
                do ii = 1, LqTherm
                    read(30,*, iostat=ios_lambda) lambda_therm(ii)
                    if (ios_lambda /= 0) then
                        call lambda_random_fill(lambda_therm, iseed)
                        exit
                    endif
                enddo
                do N = 1, ISIZE - 1
                    read(30,*) itmp
                    do nt = 1, LtrotTherm
                        do ii = 1, 2*LqTherm
                            read(30,*) sigma_itmp(ii, nt)
                        enddo
                    enddo
                    do ii = 1, LqTherm
                        read(30,*, iostat=ios_lambda) lambda_itmp(ii)
                        if (ios_lambda /= 0) then
                            call lambda_random_fill(lambda_itmp, itmp)
                            exit
                        endif
                    enddo
                    call MPI_SEND(itmp, 1, MPI_Integer, N, N, MPI_COMM_WORLD, IERR)
                    call MPI_SEND(sigma_itmp, 2*LqTherm*LtrotTherm, MPI_Real8, N, N+1024, MPI_COMM_WORLD, IERR)
                    call MPI_SEND(lambda_itmp, LqTherm, MPI_Real8, N, N+2048, MPI_COMM_WORLD, IERR)
                enddo
            endif
        else
            call MPI_RECV(iseed, 1, MPI_Integer, 0, IRANK, MPI_COMM_WORLD, STATUS, IERR)
            call MPI_RECV(sigma_therm, 2*LqTherm*LtrotTherm, MPI_Real8, 0, IRANK + 1024, MPI_COMM_WORLD, STATUS, IERR)
            call MPI_RECV(lambda_therm, LqTherm, MPI_Real8, 0, IRANK + 2048, MPI_COMM_WORLD, STATUS, IERR)
            print '("Rank ", i3.1, " after receiving message. ")', IRANK
        endif
        if (IRANK == 0 ) then
            close(30); close(10)
        endif
        return
    end subroutine conf_therm_in

    subroutine conf_set(sigma, lambda, itmp)
        real(kind=8), dimension(2*LqTherm, LtrotTherm), intent(out) :: sigma
        real(kind=8), dimension(LqTherm), intent(out) :: lambda
        integer, intent(inout) :: itmp
! Local: 
        integer :: ii, nt, nc, nx, ny
        integer :: n_list_therm(LqTherm, 1:2), bond_list_therm(2*LqTherm, 1:2)
        integer :: jx, jy, jj
        real(kind=8) :: X
        real(kind=8), external :: ranf
        
        if (absolute == 1) then
! 随机初态：每个键上的规范场随机为 +1 或 -1
            do nt = 1, LtrotTherm
                do ii = 1, 2*LqTherm
                    X = ranf(itmp)
                    if (X >= 0.5d0) then
                        sigma(ii, nt) = 1.d0
                    else
                        sigma(ii, nt) = -1.d0
                    endif
                enddo
            enddo
        elseif (absolute == 2) then
! 高斯随机初态
            do nt = 1, LtrotTherm
                do ii = 1, 2*LqTherm
                    X = rng_Gaussian(itmp)
                    if (X >= 0.d0) then
                        sigma(ii, nt) = 1.d0
                    else
                        sigma(ii, nt) = -1.d0
                    endif
                enddo
            enddo
        elseif (absolute == 3) then
! 铁磁初态：所有键上的规范场都为 +1
            do nt = 1, LtrotTherm
                do ii = 1, 2*LqTherm
                    sigma(ii, nt) = 1.d0
                enddo
            enddo
        elseif (absolute == 4) then
! 棋盘初态：对于规范场，在每个plaquette上的4个键的乘积为 -1
! 这对应一个π通量的初始配置
! 建立格点编号到坐标的映射
            nc = 0
            do ny = 1, NlyTherm
                do nx = 1, NlxTherm
                    nc = nc + 1
                    n_list_therm(nc, 1) = nx
                    n_list_therm(nc, 2) = ny
                enddo
            enddo
! 建立键列表（与lattice.f90中的逻辑一致）
            nc = 0
            do ny = 1, NlyTherm
                do nx = 1, NlxTherm
                    ii = (ny - 1) * NlxTherm + nx
                    ! 横向键
                    nc = nc + 1
                    jx = mod(nx, NlxTherm) + 1
                    jj = (ny - 1) * NlxTherm + jx
                    bond_list_therm(nc, 1) = ii
                    bond_list_therm(nc, 2) = jj
                    ! 纵向键
                    nc = nc + 1
                    jy = mod(ny, NlyTherm) + 1
                    jj = (jy - 1) * NlxTherm + nx
                    bond_list_therm(nc, 1) = ii
                    bond_list_therm(nc, 2) = jj
                enddo
            enddo
! 设置棋盘规范场配置：在 (nx+ny) 为奇数的plaquette上设置π通量
            do nt = 1, LtrotTherm
                do nc = 1, 2*LqTherm
                    ii = bond_list_therm(nc, 1)
                    nx = n_list_therm(ii, 1)
                    ny = n_list_therm(ii, 2)
                    ! 键编号是交替的：奇数是横向键，偶数是纵向键
                    if (mod(nc, 2) == 1) then  ! 奇数编号是横向键
                        if (mod(nx + ny, 2) == 0) then
                            sigma(nc, nt) = 1.d0
                        else
                            sigma(nc, nt) = -1.d0
                        endif
                    else  ! 偶数编号是纵向键
                        if (mod(nx + ny, 2) == 0) then
                            sigma(nc, nt) = 1.d0
                        else
                            sigma(nc, nt) = -1.d0
                        endif
                    endif
                enddo
            enddo
        elseif (absolute == 5) then
! 时间铁磁初态：每个键在所有虚时间片上取相同的值，但不同键之间随机
! 这适用于强时间耦合 (h 很大) 的情况，在这种情况下平衡态自然倾向于时间方向的铁磁序
! 从时间铁磁初态开始可以大大加速热化过程
            do ii = 1, 2*LqTherm
                X = ranf(itmp)
                if (X >= 0.5d0) then
                    do nt = 1, LtrotTherm
                        sigma(ii, nt) = 1.d0
                    enddo
                else
                    do nt = 1, LtrotTherm
                        sigma(ii, nt) = -1.d0
                    enddo
                endif
            enddo
        else
            write(6,*) "incorrect absolute input", absolute, " on rank ", IRANK
        endif
        call lambda_random_fill(lambda, itmp)
        return
    end subroutine conf_set

    subroutine conf_transfer(sigma, lambda, sigma_therm, lambda_therm, Latt)
        real(kind=8), dimension(2*Lq, Ltrot), intent(inout) :: sigma
        real(kind=8), dimension(2*LqTherm, LtrotTherm), intent(in) :: sigma_therm
        real(kind=8), dimension(Lq), intent(inout) :: lambda
        real(kind=8), dimension(LqTherm), intent(in) :: lambda_therm
        class(SquareLattice), intent(in) :: Latt
        integer :: nt, ntt, nx, ny, ii, jj, iit, jjt, nc, bond_idx, bond_idx_therm
        integer :: n_list_therm(LqTherm, 1:2), inv_n_list_therm(NlxTherm, NlyTherm)
        integer :: nxt, nyt, jx, jy, jxt, jyt
        integer :: site_idx

! 建立热化尺寸的格点映射（用于确定键的端点）
        nc = 0
        do ny = 1, NlyTherm
            do nx = 1, NlxTherm
                nc = nc + 1
                n_list_therm(nc, 1) = nx
                n_list_therm(nc, 2) = ny
                inv_n_list_therm(nx, ny) = nc
            enddo
        enddo

! 复制所有键的配置，通过映射键的两个端点坐标
! 对于每个实际尺寸的键，找到对应的热化尺寸的键
        do nt = 1, LtrotTherm
            do bond_idx = 1, 2*Lq
                ! 获取实际尺寸中这个键的两个端点
                ii = Latt%bond_list(bond_idx, 1)
                jj = Latt%bond_list(bond_idx, 2)
                nx = Latt%n_list(ii, 1)
                ny = Latt%n_list(ii, 2)
                jx = Latt%n_list(jj, 1)
                jy = Latt%n_list(jj, 2)
                
                ! 将坐标映射回热化尺寸（周期性）
                nxt = mod(nx - 1, NlxTherm) + 1
                nyt = mod(ny - 1, NlyTherm) + 1
                jxt = mod(jx - 1, NlxTherm) + 1
                jyt = mod(jy - 1, NlyTherm) + 1
                
                ! 找到热化尺寸中对应的格点编号
                iit = inv_n_list_therm(nxt, nyt)
                jjt = inv_n_list_therm(jxt, jyt)
                
                ! 判断是横向键还是纵向键，找到对应的键编号
                if (npbc(nx + 1, Nlx) == jx) then
                    ! 横向键：在热化尺寸中找对应的横向键
                    bond_idx_therm = (nyt - 1) * NlxTherm * 2 + (nxt - 1) * 2 + 1
                elseif (npbc(ny + 1, Nly) == jy) then
                    ! 纵向键：在热化尺寸中找对应的纵向键
                    bond_idx_therm = (nyt - 1) * NlxTherm * 2 + (nxt - 1) * 2 + 2
                endif
                
                ! 复制规范场值
                sigma(bond_idx, nt) = sigma_therm(bond_idx_therm, nt)
            enddo
        enddo

! 时间方向的周期复制
        do nt = LtrotTherm+1, Ltrot
            ntt = nt - LtrotTherm
            do bond_idx = 1, 2*Lq
                sigma(bond_idx, nt) = sigma(bond_idx, ntt)
            enddo
        enddo

! 映射 lambda：每个空间格点对应热化体系中的一个格点
        do ii = 1, Lq
            nx = Latt%n_list(ii, 1)
            ny = Latt%n_list(ii, 2)
            nxt = mod(nx - 1, NlxTherm) + 1
            nyt = mod(ny - 1, NlyTherm) + 1
            site_idx = inv_n_list_therm(nxt, nyt)
            lambda(ii) = lambda_therm(site_idx)
        enddo
        call lambda_enforce_sector(lambda)
        
        return
    end subroutine conf_transfer

    subroutine conf_out(Conf, iseed)
!#define DEC
        include 'mpif.h'
! Arguments: 
        class(SigmaConf), intent(in) :: Conf
        integer, intent(in) :: iseed
! Local:
        integer :: ii, nt, N
        integer :: status(MPI_STATUS_SIZE)
        real(kind=8), dimension(2*Lq, Ltrot) :: sigma_tmp
        real(kind=8), dimension(Lq) :: lambda_tmp

        if (IRANK == 0) open(unit=35, file='confout.txt', status='unknown')
        if (IRANK .NE. 0) then
            call MPI_SEND(iseed, 1, MPI_Integer, 0, IRANK, MPI_COMM_WORLD, IERR)
            call MPI_SEND(Conf%sigma, 2*Lq*Ltrot, MPI_Real8, 0, IRANK+1024, MPI_COMM_WORLD, IERR)
            call MPI_SEND(Conf%lambda, Lq, MPI_Real8, 0, IRANK+2048, MPI_COMM_WORLD, IERR)
            print '("Rank", i3.1, " after sending message. ")', IRANK
        endif
        if (IRANK == 0) then
            write(35,*) iseed
            do nt = 1, Ltrot
                do ii = 1, 2*Lq
                    write(35,*) Conf%sigma(ii, nt)
                enddo
            enddo
            do ii = 1, Lq
                write(35,*) Conf%lambda(ii)
            enddo
            print '("Rank ", i3.1, " after writing down own gauge fields. ")', IRANK
            do N = 1, ISIZE - 1
                call MPI_RECV(iseed, 1, MPI_Integer, N, N, MPI_COMM_WORLD, STATUS, IERR)
                call MPI_RECV(sigma_tmp, 2*Lq*Ltrot, MPI_Real8, N, N+1024, MPI_COMM_WORLD, STATUS, IERR)
                call MPI_RECV(lambda_tmp, Lq, MPI_Real8, N, N+2048, MPI_COMM_WORLD, STATUS, IERR)
                write(35,*) iseed
                do nt = 1, Ltrot
                    do ii = 1, 2*Lq
                        write(35,*) sigma_tmp(ii, nt)
                    enddo
                enddo
                do ii = 1, Lq
                    write(35,*) lambda_tmp(ii)
                enddo
                print '("Master rank ", i3.1, " after writing down gauge fields on RANK ",i3.1 ,".")', IRANK, N
            enddo
        endif
        if (IRANK == 0) close(35)
        return
    end subroutine conf_out
    
    subroutine lambda_random_fill(lambda_vec, seed)
        real(kind=8), dimension(:), intent(inout) :: lambda_vec
        integer, intent(inout) :: seed
        real(kind=8), external :: ranf
        real(kind=8) :: X
        integer :: ii

        ! 暂时禁用 λ 随机化，保持 λ = 1，以测试基本框架
        ! 当 λ = 1 时，P_λ = I，程序应该等价于不带 λ 投影的版本
        do ii = 1, size(lambda_vec)
            lambda_vec(ii) = 1.d0
        enddo
        ! TODO: 当基本框架验证通过后，可以启用以下随机化代码
        ! do ii = 1, size(lambda_vec)
        !     X = ranf(seed)
        !     if (X >= 0.5d0) then
        !         lambda_vec(ii) = 1.d0
        !     else
        !         lambda_vec(ii) = -1.d0
        !     endif
        ! enddo
        ! call lambda_enforce_sector(lambda_vec)
        return
    end subroutine lambda_random_fill

    subroutine lambda_enforce_sector(lambda_vec)
        real(kind=8), dimension(:), intent(inout) :: lambda_vec
        real(kind=8) :: prod, target
        integer :: ii

        prod = 1.d0
        do ii = 1, size(lambda_vec)
            if (lambda_vec(ii) >= 0.d0) then
                lambda_vec(ii) = 1.d0
            else
                lambda_vec(ii) = -1.d0
            endif
            prod = prod * lambda_vec(ii)
        enddo
        if (Q >= 0.d0) then
            target = 1.d0
        else
            target = -1.d0
        endif
        if (prod * target < 0.d0) then
            lambda_vec(size(lambda_vec)) = -lambda_vec(size(lambda_vec))
        endif
        return
    end subroutine lambda_enforce_sector

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