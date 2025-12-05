module Mylattice ! definition on space geometry
    use calc_basic
    implicit none

    type, public :: SquareLattice
        integer, dimension(:,:), allocatable :: n_list, inv_n_list 
        integer, dimension(:,:), allocatable :: bond_list, inv_bond_list, t_list, inv_t_list ! 专门为规范场定义的数组
        integer, dimension(:,:), allocatable :: L_bonds, nn_bonds, LT_bonds, imj
        integer, dimension(:),   allocatable :: group_1, group_2, group_3, group_4 ! 棋盘分解用的组
        integer, dimension(:,:), allocatable :: bond_bonds ! 近邻键列表
        integer, dimension(:,:), allocatable :: site_bonds ! 格点到相邻四条键的映射（用于高斯约束）
        real(kind=8), dimension(:,:), allocatable :: xk_v, aimj_v, k_dot_r
        real(kind=8) :: a1_v(2), a2_v(2), b1_v(2), b2_v(2)
    contains
        final :: Lattice_clear
    end type SquareLattice

    complex(kind=8), dimension(:,:), public, allocatable, save :: ZKRON


    interface Fourier_R_to_K
        module procedure FFT_RtoK_1, FFT_RtoK_2, FFT_RtoK_3,  FFT_RtoK_4
    end interface

contains
    subroutine Lattice_make(Latt)
        class(SquareLattice), intent(inout) :: Latt
        integer :: nc, nx, ny
        integer :: ii, jj, ix, jx, iy, jy, nt, iit, imjx, imjy
        integer :: cnt1, cnt2, cnt3, cnt4  ! 棋盘分解的计数器
! 定义晶格坐标到格点的映射
        allocate(Latt%n_list(Lq, 1:2), Latt%inv_n_list(Nlx, Nly))
        nc = 0
        do ny = 1, Nly
            do nx = 1, Nlx
                nc = nc + 1
                Latt%n_list(nc, 1) = nx
                Latt%n_list(nc, 2) = ny
                Latt%inv_n_list(nx, ny) = nc
            enddo
        enddo
! 定义键列表到格点的映射
        allocate(Latt%bond_list(2*Lq, 1:2), Latt%inv_bond_list(Lq, Lq))
        nc = 0
        do ny = 1, Nly
            do nx = 1, Nlx
                ! ii 的左边一个键
                nc = nc + 1
                ii = Latt%inv_n_list(nx, ny)
                ix = Latt%inv_n_list( npbc(nx+1, Nlx), ny )
                Latt%bond_list(nc, 1) = ii
                Latt%bond_list(nc, 2) = ix
                Latt%inv_bond_list(ii, ix) = nc
                ! ii 的上边一个键
                nc = nc + 1
                iy = Latt%inv_n_list( nx, npbc(ny+1, Nly) )
                Latt%bond_list(nc, 1) = ii
                Latt%bond_list(nc, 2) = iy
                Latt%inv_bond_list(ii, iy) = nc
            enddo
        enddo
! 定义虚时间和键列表到总规范场采样空间的编号映射
        allocate(Latt%t_list(2*Lq*Ltrot, 2), Latt%inv_t_list(2*Lq, Ltrot))
        nc = 0
        do nt = 1, Ltrot
            do ii = 1, 2*Lq
                nc = nc + 1
                Latt%t_list(nc, 1) = ii
                Latt%t_list(nc, 2) = nt
                Latt%inv_t_list(ii, nt) = nc
            enddo
        enddo
! 定义棋盘分解用的组
        allocate(Latt%group_1( Nlx / 2 * Nly ), Latt%group_2( (Nlx + 1) / 2 * Nly ) )
        allocate(Latt%group_3( Nly / 2 * Nlx ), Latt%group_4( (Nly + 1) / 2 * Nlx ) )
        cnt1 = 0; cnt2 = 0; cnt3 = 0; cnt4 = 0  ! 初始化计数器
        do nc = 1, 2 * Lq
            ii = Latt%bond_list(nc, 1)
            jj = Latt%bond_list(nc, 2)
            ix = Latt%n_list(ii, 1)
            iy = Latt%n_list(ii, 2)
            jx = Latt%n_list(jj, 1)
            jy = Latt%n_list(jj, 2)
            if (npbc(ix + 1, Nlx) == jx) then
            ! 说明是横向的键
                if (mod(ix, 2) == 0) then
                    cnt1 = cnt1 + 1
                    Latt%group_1(cnt1) = nc
                else
                    cnt2 = cnt2 + 1
                    Latt%group_2(cnt2) = nc
                endif
            elseif (npbc(iy + 1, Nly) == jy) then
            ! 说明是纵向的键
                if (mod(iy, 2) == 0) then
                    cnt3 = cnt3 + 1
                    Latt%group_3(cnt3) = nc
                else
                    cnt4 = cnt4 + 1
                    Latt%group_4(cnt4) = nc
                endif
            endif
        enddo 
! 定义格点到相邻键的映射（高斯约束用）
! site_bonds(ii, 1:4) 存储与格点 ii 相邻的四条键的编号
! 对于格点 ii，相邻的四条键是：
!   1. ii -> right neighbor (横向键，从 ii 出发)
!   2. ii -> up neighbor (纵向键，从 ii 出发)
!   3. left neighbor -> ii (横向键，终点是 ii)
!   4. down neighbor -> ii (纵向键，终点是 ii)
        allocate(Latt%site_bonds(Lq, 1:4))
        do ii = 1, Lq
            ix = Latt%n_list(ii, 1)
            iy = Latt%n_list(ii, 2)
            ! 键 1: ii -> right (横向键从 ii 出发)
            jx = npbc(ix + 1, Nlx)
            jj = Latt%inv_n_list(jx, iy)
            Latt%site_bonds(ii, 1) = Latt%inv_bond_list(ii, jj)
            ! 键 2: ii -> up (纵向键从 ii 出发)
            jy = npbc(iy + 1, Nly)
            jj = Latt%inv_n_list(ix, jy)
            Latt%site_bonds(ii, 2) = Latt%inv_bond_list(ii, jj)
            ! 键 3: left -> ii (横向键终点是 ii)
            jx = npbc(ix - 1, Nlx)
            jj = Latt%inv_n_list(jx, iy)
            Latt%site_bonds(ii, 3) = Latt%inv_bond_list(jj, ii)
            ! 键 4: down -> ii (纵向键终点是 ii)
            jy = npbc(iy - 1, Nly)
            jj = Latt%inv_n_list(ix, jy)
            Latt%site_bonds(ii, 4) = Latt%inv_bond_list(jj, ii)
        enddo

! 定义近邻键列表(规范场专用)
        allocate(Latt%bond_bonds(2*Lq, 0:6))
        do nc = 1, 2*Lq
            ii = Latt%bond_list(nc, 1)
            jj = Latt%bond_list(nc, 2)
            ix = Latt%n_list(ii, 1)
            iy = Latt%n_list(ii, 2)
            jx = Latt%n_list(jj, 1)
            jy = Latt%n_list(jj, 2)

            Latt%bond_bonds(nc, 0) = nc
            if (npbc(ix + 1, Nlx) == jx) then
            ! 说明是横向的键
                nx = Latt%inv_n_list( ix, npbc(iy+1, Nly) )
                Latt%bond_bonds(nc, 1) = Latt%inv_bond_list(ii, nx)
                ny = Latt%inv_n_list( jx, npbc(jy+1, Nly) )
                Latt%bond_bonds(nc, 2) = Latt%inv_bond_list(nx, ny)
                Latt%bond_bonds(nc, 3) = Latt%inv_bond_list(jj, ny)

                nx = Latt%inv_n_list( ix, npbc(iy-1, Nly) )
                Latt%bond_bonds(nc, 4) = Latt%inv_bond_list(nx, ii)
                ny = Latt%inv_n_list( jx, npbc(jy-1, Nly) )
                Latt%bond_bonds(nc, 5) = Latt%inv_bond_list(nx, ny)
                Latt%bond_bonds(nc, 6) = Latt%inv_bond_list(ny, jj)
            elseif (npbc(iy + 1, Nly) == jy) then
            ! 说明是纵向的键
                nx = Latt%inv_n_list( npbc(jx-1, Nlx), jy)
                Latt%bond_bonds(nc, 1) = Latt%inv_bond_list(nx, jj)
                ny = Latt%inv_n_list( npbc(ix-1, Nlx), iy)
                Latt%bond_bonds(nc, 2) = Latt%inv_bond_list(ny, nx)
                Latt%bond_bonds(nc, 3) = Latt%inv_bond_list(ny, ii)

                nx = Latt%inv_n_list( npbc(jx+1, Nlx), jy)
                Latt%bond_bonds(nc, 4) = Latt%inv_bond_list(jj, nx)
                ny = Latt%inv_n_list( npbc(ix+1, Nlx), iy)
                Latt%bond_bonds(nc, 5) = Latt%inv_bond_list(ny, nx)
                Latt%bond_bonds(nc, 6) = Latt%inv_bond_list(ii, ny)
            endif
        enddo

! 定义两个格点差到格点编号的映射（利用平移对称性）
        allocate(Latt%imj(Lq, Lq))
        do jj = 1, Lq
            do ii =1, Lq
                ix = Latt%n_list(ii, 1)
                iy = Latt%n_list(ii, 2)
                jx = Latt%n_list(jj, 1)
                jy = Latt%n_list(jj, 2)
                imjx = npbc(ix - jx, Nlx)
                imjy = npbc(iy - jy, Nly)
                Latt%imj(ii, jj) = Latt%inv_n_list(imjx, imjy)
            enddo
        enddo

! 定义晶格基矢和倒格矢
        Latt%a1_v(1) = 1.d0; Latt%a1_v(2) = 0.d0
        Latt%a2_v(1) = 0.d0; Latt%a2_v(2) = 1.d0
        Latt%b1_v(1) = 2.d0 * PI; Latt%b1_v(2) = 0.d0
        Latt%b2_v(1) = 0.d0; Latt%b2_v(2) = 2.d0 * PI
! 定义傅里叶变换的数组（预计算）
        allocate(Latt%xk_v(Lq, 2), Latt%aimj_v(Lq, 2), Latt%k_dot_r(Lq, Lq))
        do ii = 1, Lq
            ix = Latt%n_list(ii, 1)
            iy = Latt%n_list(ii, 2)
            Latt%aimj_v(ii, 1) = dble(ix) * Latt%a1_v(1) + dble(iy) * Latt%a2_v(1)
            Latt%aimj_v(ii, 2) = dble(ix) * Latt%a1_v(2) + dble(iy) * Latt%a2_v(2)
            Latt%xk_v(ii, 1) = dble(ix - 1) * Latt%b1_v(1)/dble(Nlx) + dble(iy - 1) * Latt%b2_v(1)/dble(Nly)
            Latt%xk_v(ii, 2) = dble(ix - 1) * Latt%b1_v(2)/dble(Nlx) + dble(iy - 1) * Latt%b2_v(2)/dble(Nly)
        enddo
        do jj = 1, Lq
            do ii = 1, Lq
                Latt%k_dot_r(ii, jj) = Latt%xk_v(ii, 1)*Latt%aimj_v(jj, 1) + Latt%xk_v(ii, 2)*Latt%aimj_v(jj, 2)
            enddo
        enddo
! 定义单位矩阵，用来计算 sm 更新中的 Delta函数
        allocate(ZKRON(Ndim, Ndim))
        ZKRON = dcmplx(0.d0, 0.d0)
        do ii = 1, Ndim
            ZKRON(ii, ii) = dcmplx(1.d0, 0.d0)
        enddo !define delta function

        allocate(Latt%L_Bonds(Lq, 0:4), Latt%nn_bonds(Lq, 0:4), Latt%LT_bonds(2*Lq*Ltrot, 0:6))
!define the nearest neighbor bonds
        do iy = 1, Nly
            do ix = 1, Nlx
                ii = Latt%inv_n_list(ix, iy)
                Latt%L_bonds(ii, 0) = ii
                Latt%L_bonds(ii, 1) = Latt%inv_n_list( npbc(ix+1, Nlx), iy )
                Latt%L_bonds(ii, 2) = Latt%inv_n_list( ix, npbc(iy+1, Nly) )
                Latt%L_bonds(ii, 3) = Latt%inv_n_list( npbc(ix-1, Nlx), iy )
                Latt%L_bonds(ii, 4) = Latt%inv_n_list( ix, npbc(iy-1, Nly) )
            enddo
        enddo
! define the second nearest neighbor bonds
        do iy = 1, Nly
            do ix = 1, Nlx
                ii = Latt%inv_n_list(ix, iy)
                Latt%nn_bonds(ii, 0) = ii
                Latt%nn_bonds(ii, 1) = Latt%inv_n_list( npbc(ix+1, Nlx), npbc(iy+1, Nly) )
                Latt%nn_bonds(ii, 2) = Latt%inv_n_list( npbc(ix-1, Nlx), npbc(iy+1, Nly) )
                Latt%nn_bonds(ii, 3) = Latt%inv_n_list( npbc(ix-1, Nlx), npbc(iy-1, Nly) )
                Latt%nn_bonds(ii, 4) = Latt%inv_n_list( npbc(ix+1, Nlx), npbc(iy-1, Nly) )
            enddo
        enddo
! define the nearest neighbors on space-time
        do nt = 1, Ltrot
            do ii = 1, Lq
                iit = Latt%inv_t_list(ii, nt)
                Latt%LT_bonds(iit, 0) = iit
                Latt%LT_bonds(iit, 1) = Latt%inv_t_list(Latt%L_bonds(ii, 1), nt)
                Latt%LT_bonds(iit, 2) = Latt%inv_t_list(Latt%L_bonds(ii, 2), nt)
                Latt%LT_bonds(iit, 3) = Latt%inv_t_list(Latt%L_bonds(ii, 3), nt)
                Latt%LT_bonds(iit, 4) = Latt%inv_t_list(Latt%L_bonds(ii, 4), nt)
                Latt%LT_bonds(iit, 5) = Latt%inv_t_list(ii, npbc(nt+1, Ltrot))
                Latt%LT_bonds(iit, 6) = Latt%inv_t_list(ii, npbc(nt-1, Ltrot))
            enddo
        enddo
	    return
    end subroutine Lattice_make

    subroutine Lattice_clear(this)
        type(SquareLattice), intent(inout) :: this
        if (allocated(this%n_list))        deallocate(this%n_list)
        if (allocated(this%inv_n_list))    deallocate(this%inv_n_list)
        if (allocated(this%bond_list))     deallocate(this%bond_list)
        if (allocated(this%inv_bond_list)) deallocate(this%inv_bond_list)
        if (allocated(this%t_list))        deallocate(this%t_list)
        if (allocated(this%inv_t_list))    deallocate(this%inv_t_list)
        if (allocated(this%L_bonds))       deallocate(this%L_bonds)
        if (allocated(this%nn_bonds))      deallocate(this%nn_bonds)
        if (allocated(this%LT_bonds))      deallocate(this%LT_bonds)
        if (allocated(this%imj))           deallocate(this%imj)
        if (allocated(this%group_1))       deallocate(this%group_1)
        if (allocated(this%group_2))       deallocate(this%group_2)
        if (allocated(this%group_3))       deallocate(this%group_3)
        if (allocated(this%group_4))       deallocate(this%group_4)
        if (allocated(this%bond_bonds))    deallocate(this%bond_bonds)
        if (allocated(this%site_bonds))    deallocate(this%site_bonds)
        if (allocated(this%xk_v))          deallocate(this%xk_v)
        if (allocated(this%aimj_v))        deallocate(this%aimj_v)
        if (allocated(this%k_dot_r))       deallocate(this%k_dot_r)
        if (allocated(ZKRON))              deallocate(ZKRON)
        return
    end subroutine Lattice_clear
    
    subroutine FFT_RtoK_1(gr, gk, Latt)
        complex(kind=8), dimension(Lq), intent(in) :: gr
        complex(kind=8), dimension(Lq), intent(out) :: gk
        class(SquareLattice), intent(in) :: Latt
        integer :: imj, nk
        
        gk = dcmplx(0.d0, 0.d0)
        do imj = 1, Lq
            do nk = 1, Lq
                gk(nk) = gk(nk) + exp( dcmplx(0.d0, Latt%k_dot_r(nk, imj)) ) * gr(imj)
            enddo
        enddo
        gk = gk/dble(Lq)
        return
    end subroutine FFT_RtoK_1
    
    subroutine FFT_RtoK_2(gr, gk, Latt)
        complex(kind=8), dimension(:,:), intent(in) :: gr
        complex(kind=8), dimension(:,:), intent(out) :: gk
        class(SquareLattice), intent(in) :: Latt
        integer :: imj, nk, nf, NN2
        
        NN2 = size(gr, 2)
        if ((size(gr, 1) .ne. Lq) .or. (size(gk, 1) .ne. Lq)) then
            write(6,*), "incorrect matrix size in FFT_RtoK_2"; stop
        endif
        gk = dcmplx(0.d0, 0.d0)
        do nf = 1, NN2
            do imj = 1, Lq
                do nk = 1, Lq
                    gk(nk, nf) = gk(nk, nf) + exp( dcmplx(0.d0, Latt%k_dot_r(nk, imj)) ) * gr(imj, nf)
                enddo
            enddo
        enddo
        gk = gk/dble(Lq)
        return
    end subroutine FFT_RtoK_2
    
    subroutine FFT_RtoK_3(gr, gk, Latt)
        complex(kind=8), dimension(:,:,:), intent(in) :: gr
        complex(kind=8), dimension(:,:,:), intent(out) :: gk
        class(SquareLattice), intent(in) :: Latt
        integer :: imj, nk, nf2, nf3, NN2, NN3
        
        NN2 = size(gr, 2); NN3 = size(gr, 3)
        if ((size(gr, 1) .ne. Lq) .or. (size(gk, 1) .ne. Lq)) then
            write(6,*), "incorrect matrix size in FFT_RtoK_3"; stop
        endif
        gk = dcmplx(0.d0, 0.d0)
        do nf3 = 1, NN3
            do nf2 = 1, NN2
                do imj = 1, Lq
                    do nk = 1, Lq
                        gk(nk, nf2, nf3) = gk(nk, nf2, nf3) + exp( dcmplx(0.d0, Latt%k_dot_r(nk, imj)) ) * gr(imj, nf2, nf3)
                    enddo
                enddo
            enddo
        enddo
        gk = gk/dble(Lq)
        return
    end subroutine FFT_RtoK_3
    
    subroutine FFT_RtoK_4(gr, gk, Latt)
        complex(kind=8), dimension(:,:,:,:), intent(in) :: gr
        complex(kind=8), dimension(:,:,:,:), intent(out) :: gk
        class(SquareLattice), intent(in) :: Latt
        integer :: imj, nk, nf2, nf3, nf4, NN2, NN3, NN4
        
        NN2 = size(gr, 2); NN3 = size(gr, 3); NN4 = size(gr, 4)
        if ((size(gr, 1) .ne. Lq) .or. (size(gk, 1) .ne. Lq)) then
            write(6,*), "incorrect matrix size in FFT_RtoK_4"; stop
        endif
        gk = dcmplx(0.d0, 0.d0)
        do nf4 = 1, NN4
            do nf3 = 1, NN3
                do nf2 = 1, NN2
                    do imj = 1, Lq
                        do nk = 1, Lq
                            gk(nk, nf2, nf3, nf4) = gk(nk, nf2, nf3, nf4) + exp( dcmplx(0.d0, Latt%k_dot_r(nk, imj)) ) * gr(imj, nf2, nf3, nf4)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        gk = gk/dble(Lq)
        return
    end subroutine FFT_RtoK_4
end module MyLattice
