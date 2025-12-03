module calc_basic
    implicit none

! 常数
    real(kind = 8), parameter :: Zero = 1.0d-10
    real(kind = 8), parameter :: PI = acos(-1.d0)
    real(kind = 8), parameter :: upbound = 1.0d+200
! 晶格参数
    integer, parameter :: Norb = 1 ! 该模型只有 1 个轨道
    integer, parameter :: Nbond = 2 ! 晶格有两个方向，水平和竖直
    integer, parameter :: Nspin = 2 ! 自旋σ = 1、2
    integer, public :: Nlx, Nly, NlxTherm, NlyTherm ! 晶格大小和热化时的晶格大小
    integer, public :: Lq, LqTherm ! 总格点数目
    integer, public :: Ndim ! 格林函数维度
    real(kind = 8), public :: Dtau ! 虚时间隔
    real(kind = 8), public :: Beta ! 逆温度
    integer, public :: Ltrot, LtrotTherm ! 虚时长度
! 哈密顿量参数
    real(kind = 8), public :: RT ! 跃迁项的系数 t
    real(kind = 8), public :: J ! 规范场耦合系数
    real(kind = 8), public :: h ! 规范场横场强度
    real(kind = 8), public :: mu ! 化学势
    real(kind = 8), public :: Q ! 所选取的规范扇区（+1 为 even 扇区； -1 为 odd 扇区）
! 更新参数
    logical, public :: is_global ! 是否全局更新
    integer, public :: Nglobal ! 全局更新频率（多少步更新一次）
! 初始化状态参数
    integer, public :: absolute ! 选择初始化Ising场的方式
    ! 说明：Ising场 phi 严格取值为 ±1。初始化模式如下：
    ! absolute = 1: 均匀随机 ±1
    ! absolute = 2: 取高斯随机数的符号作为 ±1
    ! absolute = 3: 全部设为 +1（铁磁初态）
    ! absolute = 4: 反铁磁初态（棋盘模式，相邻格点相反磁化）
    ! absolute = 5: 时间铁磁初态（每个bond在所有τ上同值，空间随机）
    ! 实现见 Fields_mod::conf_set
! 过程控制参数
    logical, public :: is_tau ! 是否计算含时观测量
    integer, public :: Nthermal ! 计算含时观测量前热化步数（自定义）
    logical, public :: is_warm ! 是否进行热化（只使用Ising作用量更新场）
    integer, public :: Nwarm ! 热化步数（自定义）
    integer, public :: Nst ! 存储稳定矩阵的格式，Ltrot/Nwrap
    integer, public :: Nwrap ! 稳定化频率，隔 Nwrap 步存储一次
    integer, public :: Nbin ! 数据存储频率，隔 Nbin 步存储一次
    integer, public :: Nsweep ! 来回扫的更新次数（蒙卡步）
    integer, public :: ISIZE, IRANK, IERR ! 并行参数
! 额外控制参数
    real(kind = 8), public :: NB_field = 0.d0 ! 默认为 0，用于兼容含时观测量中的磁通相位

contains
    subroutine read_input()
        include 'mpif.h'
        if (IRANK == 0) then
            open(unit = 20, file = 'paramC_sets.txt', status = 'unknown')
            read(20, *) RT, J, h, mu, Q
            read(20, *) Nlx, Nly, Ltrot, Beta
            read(20, *) NlxTherm, NlyTherm, LtrotTherm
            read(20, *) Nwrap, Nbin, Nsweep
            read(20, *) is_tau, Nthermal
            read(20, *) is_warm, Nwarm
            read(20, *) is_global, Nglobal
            read(20, *) absolute
            close(20)
        endif

! 广播参数
        call MPI_BCAST(Beta, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(RT, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(J, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(h, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(mu, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Q, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(absolute, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nlx, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nly, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Ltrot, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NlxTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NlyTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(LtrotTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nwrap, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nbin, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nsweep, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_tau, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nthermal, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_warm, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nwarm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_global, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nglobal, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        return
    end subroutine read_input

    subroutine Params_set()
        Lq = Nlx * Nly
        LqTherm = NlxTherm * NlyTherm
        Dtau = Beta / dble(Ltrot)
        Ndim = Lq
        if (mod(Ltrot, Nwrap) == 0) then
            Nst = Ltrot / Nwrap
        else
            write(6, *) "Ltrot必须为Nwrap的整数倍"; stop
        endif
        return
    end subroutine Params_set

    ! 该函数用于生成1到N范围内的随机整数
    ! 输入参数:
    !   iseed - 随机数种子(输入/输出)
    !   N     - 随机数上限(输入)
    ! 返回值: 1到N之间的随机整数
    ! 实现细节:
    ! 1. 调用ranf(iseed)生成[0,1)区间的随机实数
    ! 2. 缩放并偏移到[0.5,N+0.5]区间
    ! 3. 取整后得到[1,N]区间的随机整数
    integer function nranf(iseed, N)
        integer, intent( inout ) :: iseed
        integer, intent( in ) :: N
        real(kind = 8), external :: ranf
        nranf = nint( ranf(iseed) * dble(N) + 0.5 )
        if (nranf < 1) nranf = 1
        if (nranf > N) nranf = N
        return
    end function nranf

    ! 定义周期边界条件
    integer function npbc(nr, L)
        integer, intent( in ) :: nr, L
        npbc = nr
        if (nr < 1) npbc = nr + L
        if (nr > L) npbc = nr - L
        return
    end function npbc

    subroutine write_info()
        if (IRANK == 0) then
            open (unit = 50, file = 'info.txt', status = 'unknown', action = 'write')
            write(50, *) '=============================================='
            write(50, *) 'z2规范场模型'
            write(50, *) '晶格 x 方向长度 Lx          :', Nlx
            write(50, *) '晶格 y 方向长度 Ly          :', Nly
            write(50, *) '跃迁系数 t                  :', RT
            write(50, *) '规范场耦合系数 J            :', J
            write(50, *) '规范场横场强度 h            :', h
            write(50, *) '化学势 mu                   :', mu
            write(50, *) '所选取的规范扇区 Q          :', Q
            if (is_global) then
                write(50, *) '全局更新频率 Nglobal          :', Nglobal
            endif
            if (is_warm) then
                write(50, *) '热化步数 Nwarm               :', Nwarm
            endif
            write(50, *) '选择的初始化场构型方式 absolute :', absolute
            write(50, *) '逆温度 Beta                   :', Beta
            write(50, *) '虚时长度 Ltrot                :', Ltrot
            write(50, *) '=> Dtau = Beta / Ltrot      :', Dtau
            write(50, *) '稳定化频率 Nwrap              :', Nwrap
            write(50, *) '# Bins                       :', Nbin
            write(50, *) '来回扫的更新次数 Nsweep        :', Nsweep
            write(50, *) '# No. of MPI processes       :', ISIZE
            if (mod(Ltrot, Nwrap) /= 0) then
                write(50, *) 'Ltrot 必须为 Nwrap 的整数倍'
            endif
            write(50, *) '=============================================='
            call flush(50)
        endif
        return
    end subroutine write_info

end module calc_basic