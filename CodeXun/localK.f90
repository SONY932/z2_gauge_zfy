module LocalK_mod
    use Multiply_mod
    implicit none
    
    public
    private :: LocalK_metro, phi_new
    
    type(AccCounter) :: Acc_Kl, Acc_Kt
    real(kind=8), dimension(:,:,:), allocatable :: phi_new
    
contains
    subroutine LocalK_init()
        call Acc_Kl%init()
        call Acc_Kt%init()
        allocate(phi_new(Nspin, Lq, Ltrot))
        return
    end subroutine LocalK_init
    
    subroutine LocalK_clear()
        deallocate(phi_new)
        return
    end subroutine LocalK_clear
    
    subroutine LocalK_reset()
        call Acc_Kl%reset()
        call Acc_Kt%reset()
        phi_new = NsigL_K%phi
        return
    end subroutine LocalK_reset
    
    subroutine LocalK_metro(Gr, iseed, ii, ntau)
        use MyMats
! Arguments:
	    complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gr
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau
!   Local: 
        real(kind=8), external :: ranf
        complex(kind=8) :: Proddet
        complex(kind=8), dimension(Norb, Norb) :: Prod, Prodinv, Gr_local, mat_tmp
        complex(kind=8) :: Vhlp(2, Ndim), Uhlp(Ndim, 2), temp(Ndim, 2), Diff(Ndim, Ndim)
        real(kind=8) :: ratio_fermion, ratio_boson, ratio_re, ratio_re_abs
        real(kind=8) :: random, Xdif, xflip
        integer :: ns, P(Norb), j, no, sign, nl, nr
        real(kind=8), dimension(Nspin) :: vec_new, vec_old

! Local update on a two-component spin vector on space-time (ii, ntau)
        do ns = 1, Nspin
            xflip = ranf(iseed)
            Xdif = dble((xflip - 0.5) * abs(valr0))
            phi_new(ns, ii, ntau) = NsigL_K%phi(ns, ii, ntau) + Xdif
        enddo
        vec_new(:) = phi_new(:, ii, ntau)
        vec_old(:) = NsigL_K%phi(:, ii, ntau)
        do no = 1, Norb
            P(no) = Latt%inv_o_list(ii, no)
        enddo
        if (Latt%b_list(ii, 2) == 1) sign = 1
        if (Latt%b_list(ii, 2) == 2) sign = -1
! Calculate fermionic Metropolis ratio within 2*2 matrix space
        call Op_K%get_delta(vec_old, vec_new, sign) ! update Delta matrix in Op_K
        Prod = dcmplx(0.d0, 0.d0)
        do nr = 1, Norb
            do nl = 1, Norb
                Gr_local(nl, nr) = ZKRON(nl, nr) - Gr(P(nl), P(nr))
            enddo
        enddo
        call mmult(mat_tmp, Op_K%Delta, Gr_local) ! 2*2 matrix multiplication
        do nr = 1, Norb
            do nl = 1, Norb
                Prod(nl, nr) = ZKRON(nl, nr) + mat_tmp(nl, nr)
            enddo
        enddo
        Proddet = Prod(1,1) * Prod(2,2) - Prod(1,2) * Prod(2,1)
        ratio_fermion = real(Proddet * dconjg(Proddet))
! Calculate total Metropolis ratio     
        ratio_boson = NsigL_K%bosonratio(phi_new, ii, ntau, Latt)
        ratio_re = dble(ratio_fermion * ratio_boson)
        ratio_re_abs = abs(ratio_re)
        random = ranf(iseed)
! Upgrade Green's function
        if (ratio_re_abs .gt. random) then
            call Acc_Kl%count(.true.)
            Prodinv(1,1) = Prod(2,2)
            Prodinv(2,2) = Prod(1,1)
            Prodinv(1,2) = - Prod(1,2)
            Prodinv(2,1) = - Prod(2,1)
            Prodinv = Prodinv / Proddet
            Uhlp = dcmplx(0.d0, 0.d0); Vhlp = dcmplx(0.d0, 0.d0)
            temp = dcmplx(0.d0, 0.d0); Diff = dcmplx(0.d0, 0.d0)
! Vhlp(1:2, 1:Ndim) = Del(1:2) * (1 - Grup)(P(1):P(2), 1:Ndim); Uhlp(1:Ndim, 1:2) = Grup(1:Ndim, P(1):P(2))
            do no = 1, Norb
                do j = 1, Ndim
                    Uhlp(j, no) = Gr(j, P(no))
                    Vhlp(no, j) = - Op_K%Delta(no, 1) * Gr(P(1), j) - Op_K%Delta(no, 2) * Gr(P(2), j)
                enddo
                Vhlp(no, P(1)) = Vhlp(no, P(1)) + Op_K%Delta(no, 1)
                Vhlp(no, P(2)) = Vhlp(no, P(2)) + Op_K%Delta(no, 2)
            enddo
            call mmult(temp, Uhlp, Prodinv)
            call mmult(Diff, temp, Vhlp)
            Gr = Gr - Diff ! output Gr in each spin-orbital sector
! Flip: 
            NsigL_K%phi(:, ii, ntau) = phi_new(:, ii, ntau)
        else
            call Acc_Kl%count(.false.)
            phi_new(:, ii, ntau) = NsigL_K%phi(:, ii, ntau)
        endif
        return
    end subroutine LocalK_metro
    
    subroutine LocalK_prop_L(Prop, iseed, nt)
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        integer :: ii
        do ii = Lq, 1, -1
            call LocalK_metro(Prop%Gr, iseed, ii, nt)
        enddo
        call Op_K%mmult_L(Prop%Gr, Latt, NsigL_K%phi, nt, 1)
        call Op_K%mmult_R(Prop%Gr, Latt, NsigL_K%phi, nt, -1)
        call Op_K%mmult_L(Prop%UUL, Latt, NsigL_K%phi, nt, 1)
        return
    end subroutine LocalK_prop_L
    
    subroutine LocalK_prop_R(Prop, iseed, nt)
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        integer :: ii
        call Op_K%mmult_R(Prop%Gr, Latt, NsigL_K%phi, nt, 1)
        call Op_K%mmult_L(Prop%Gr, Latt, NsigL_K%phi, nt, -1)
        do ii = 1, Lq
            call LocalK_metro(Prop%Gr, iseed, ii, nt)
        enddo
        call Op_K%mmult_R(Prop%UUR, Latt, NsigL_K%phi, nt, 1)
        return
    end subroutine LocalK_prop_R
    
    subroutine LocalK_therm(ii, ntau, iseed)
! Arguments: 
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau
! Local: 
        real(kind=8), external :: ranf
        real(kind=8) :: xflip, random, Xdif, ratio_boson
        integer :: ns

        do ns = 1, Nspin
            xflip = ranf(iseed)
            Xdif = dble((xflip - 0.5) * abs(valrt(ns)))
            phi_new(ns, ii, ntau) = NsigL_K%phi(ns, ii, ntau) + Xdif
        enddo
        ratio_boson = NsigL_K%bosonratio(phi_new, ii, ntau, Latt)
        random = ranf(iseed)
        if (ratio_boson .gt. random) then
            call Acc_Kt%count(.true.)
            NsigL_K%phi(:, ii, ntau) = phi_new(:, ii, ntau)
        else
            call Acc_Kt%count(.false.)
            phi_new(:, ii, ntau) = NsigL_K%phi(:, ii, ntau)
        endif
        return
    end subroutine LocalK_therm
end module LocalK_mod