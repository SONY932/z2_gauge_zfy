module NonInteract
    use MyLattice
    implicit none

    public
    private :: def_hamT
    
    type :: OperatorKinetic
        complex(kind=8), dimension(:,:), allocatable :: expT_P, expT_M
    contains
        procedure :: make => opT_make
        procedure :: set => opT_set
        procedure :: mmult_R => opT_mmult_R
        procedure :: mmult_L => opT_mmult_L
        final :: opT_clear
    end type OperatorKinetic
    
contains
    subroutine opT_make(this)
        class(OperatorKinetic), intent(inout) :: this
        allocate(this%expT_P(Ndim, Ndim), this%expT_M(Ndim, Ndim))
        this%expT_P = dcmplx(0.d0, 0.d0); this%expT_M = dcmplx(0.d0, 0.d0)
        return
    end subroutine opT_make
    
    subroutine opT_clear(this)
        type(OperatorKinetic), intent(inout) :: this
        deallocate(this%expT_P, this%expT_M)
        return
    end subroutine opT_clear
    
    subroutine def_hamT(HamT, Latt)
! Arguments: 
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: HamT
        class(SquareLattice), intent(in) :: Latt
! Local: 
        complex(kind=8) :: Z
        integer :: i, ii, i_0, i_n, ix, iy, nf, no
        
        HamT = dcmplx(0.d0, 0.d0)
!  nearest bond hopping
        do no = 1, Norb
            do ii = 1, Lq
                i_0 = Latt%inv_o_list(Latt%L_Bonds(ii, 0), no)
                ix = Latt%n_list(ii, 1)
                iy = Latt%n_list(ii, 2)
                do nf = 1, Nbond
                    i_n = Latt%inv_o_list(Latt%L_Bonds(ii, nf), no)
                    if (nf == 1) then
                        Z = dcmplx( - RT(no, nf), 0.d0) * &
                            &   exp(- dcmplx(0.d0, 1.d0) * 2.d0 * Pi * NB_field * dble(iy)/dble(Lq)) ! adding one flux quanta(NB_field=1) on finite lattice under Landau gauge A=(-By,0,0)
                    elseif (nf == 2) then
                        if (iy .NE. Nly) then
                            Z = dcmplx( - RT(no, nf), 0.d0)
                        else
                            Z = dcmplx( - RT(no, nf), 0.d0) * &
                                &   exp( dcmplx(0.d0, 1.d0) * 2.d0 * Pi * NB_field * dble(ix)/dble(Nlx))
                        endif
                    else 
                        write(6,*) "incorrect nearest neighbor", Nbond, nf; stop
                    endif
                    HamT(i_0, i_n)  = HamT(i_0, i_n) + Z
                    HamT(i_n, i_0)  = HamT(i_n, i_0) + dconjg(Z)
                enddo
            enddo
        enddo
        do no = 1, Norb
            do ii = 1, Lq
                i = Latt%inv_o_list(ii, no)
                HamT(i, i) = HamT(i, i) - dcmplx(mu, 0.d0) ! uniform chemical potential
            enddo
        enddo
        return
    end subroutine def_hamT

    subroutine opT_set(this, Latt)
        use MyMats
! Arguments: 
        class(OperatorKinetic), intent(inout) :: this
        class(SquareLattice), intent(in) :: Latt
! Local: 
!        real(kind=8) :: degen, en_free
        integer :: i, nl, nr
        complex(kind=8), dimension(Ndim, Ndim) :: HamT, Hlp1, Hlp1dag, temp1, temp2
        real(kind=8), dimension(Ndim) :: WC
        complex(kind=8), dimension(Ndim) :: dmat1, dmat2
        
        call def_HamT(HamT, Latt)
        call diag(HamT, Hlp1, WC)

        dmat1 = dcmplx(0.d0, 0.d0)
        dmat2 = dcmplx(0.d0, 0.d0)
        forall( i = 1:Ndim )
            dmat1(i) = dcmplx(exp(- Dtau * WC(i)), 0.d0)
            dmat2(i) = dcmplx(exp( Dtau * WC(i)), 0.d0)
        endforall
        Hlp1dag = dconjg(transpose(Hlp1))
        do nr = 1, Ndim
            do nl = 1, Ndim
                temp1(nl, nr) = Hlp1(nl, nr) * dmat1(nr)
                temp2(nl, nr) = Hlp1(nl, nr) * dmat2(nr)
            enddo
        enddo
        call mmult(this%expT_P, temp1, Hlp1dag) ! output
        call mmult(this%expT_M, temp2, Hlp1dag) ! output
        return
    end subroutine opT_set
    
    subroutine opT_mmult_R(this, Mat, nflag)
        use MyMats
!	In Mat Out exp(-Dtau*T) * Mat for nflag = 1
!	In Mat Out exp( Dtau*T) * Mat for nflag = -1  
        class(OperatorKinetic), intent(in) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        integer, intent(in) :: nflag
        complex(kind=8), dimension(Ndim, Ndim) :: temp
        if (nflag == 1) then
            call mmult(temp, this%expT_P, Mat)
        elseif (nflag == -1) then
            call mmult(temp, this%expT_M, Mat)
        else
            write(6,*) "incorrect nflag in opT_mmult_R"; stop
        endif
        Mat = temp
        return
    end subroutine opT_mmult_R
    
    subroutine opT_mmult_L(this, Mat, nflag)
        use MyMats
!	In Mat Out Mat * exp(-Dtau*T) for nflag = 1
!	In Mat Out Mat * exp( Dtau*T) for nflag = -1
        class(OperatorKinetic), intent(in) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        integer, intent(in) :: nflag
        complex(kind=8), dimension(Ndim, Ndim) :: temp
        if (nflag == 1) then
            call mmult(temp, Mat, this%expT_P)
        elseif (nflag == -1) then
            call mmult(temp, Mat, this%expT_M)
        else
            write(6,*) "incorrect nflag in opT_mmult_L"; stop
        endif
        Mat = temp
        return
    end subroutine opT_mmult_L
end module NonInteract