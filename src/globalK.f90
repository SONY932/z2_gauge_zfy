module GlobalK_mod
    use Multiply_mod
    use DQMC_Model_mod
    use ProcessMatrix
    use MyLattice
    use calc_basic
    implicit none

    real(kind=8), parameter :: RATIO_EPS = 1.d-300

    public
    private :: ratioK_fermion

contains
    real(kind=8) function ratioK_fermion(GrU, GrD, sigma_new, ii, ntau) result(ratio_fermion)
! Arguments:
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: GrU, GrD
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_new
        integer, intent(in) :: ii, ntau
! Local:
        complex(kind=8), dimension(2, 2) :: ProdU, ProdD, ProdinvU, ProdinvD
        complex(kind=8), dimension(2, 2) :: GrU_local, GrD_local, matU_tmp, matD_tmp
        complex(kind=8), dimension(Ndim, 2) :: UhlpU, UhlpD
        complex(kind=8), dimension(2, Ndim) :: VhlpU, VhlpD
        complex(kind=8), dimension(Ndim, 2) :: temp
        complex(kind=8), dimension(Ndim, Ndim) :: Diff
        complex(kind=8) :: ProddetU, ProddetD
        real(kind=8) :: S_old, S_new
        integer :: P(2), nl, nr, j

        S_old = NsigL_K%sigma(ii, ntau)
        S_new = sigma_new(ii, ntau)
        if (S_new == S_old) then
            ratio_fermion = 1.d0
            return
        endif

        P(1) = Latt%bond_list(ii, 1)
        P(2) = Latt%bond_list(ii, 2)
        call Op_K%get_delta(S_old, S_new)

        ProdU = dcmplx(0.d0, 0.d0)
        ProdD = dcmplx(0.d0, 0.d0)
        do nr = 1, 2
            do nl = 1, 2
                GrU_local(nl, nr) = ZKRON(nl, nr) - GrU(P(nl), P(nr))
                GrD_local(nl, nr) = ZKRON(nl, nr) - GrD(P(nl), P(nr))
            enddo
        enddo

        matU_tmp = matmul(Op_K%Delta, GrU_local)
        matD_tmp = matmul(Op_K%Delta, GrD_local)
        do nr = 1, 2
            do nl = 1, 2
                ProdU(nl, nr) = ZKRON(nl, nr) + matU_tmp(nl, nr)
                ProdD(nl, nr) = ZKRON(nl, nr) + matD_tmp(nl, nr)
            enddo
        enddo

        ProddetU = ProdU(1, 1) * ProdU(2, 2) - ProdU(1, 2) * ProdU(2, 1)
        ProddetD = ProdD(1, 1) * ProdD(2, 2) - ProdD(1, 2) * ProdD(2, 1)
        ratio_fermion = max(abs(ProddetU * ProddetD), RATIO_EPS)

! Update spin-up Green's function
        ProdinvU(1, 1) = ProdU(2, 2)
        ProdinvU(2, 2) = ProdU(1, 1)
        ProdinvU(1, 2) = -ProdU(1, 2)
        ProdinvU(2, 1) = -ProdU(2, 1)
        ProdinvU = ProdinvU / ProddetU
        UhlpU = dcmplx(0.d0, 0.d0)
        VhlpU = dcmplx(0.d0, 0.d0)
        temp = dcmplx(0.d0, 0.d0)
        Diff = dcmplx(0.d0, 0.d0)
        do nl = 1, 2
            do j = 1, Ndim
                UhlpU(j, nl) = GrU(j, P(nl))
                VhlpU(nl, j) = - Op_K%Delta(nl, 1) * GrU(P(1), j) - Op_K%Delta(nl, 2) * GrU(P(2), j)
            enddo
            VhlpU(nl, P(1)) = VhlpU(nl, P(1)) + Op_K%Delta(nl, 1)
            VhlpU(nl, P(2)) = VhlpU(nl, P(2)) + Op_K%Delta(nl, 2)
        enddo
        temp = matmul(UhlpU, ProdinvU)
        Diff = matmul(temp, VhlpU)
        GrU = GrU - Diff

! Update spin-down Green's function
        ProdinvD(1, 1) = ProdD(2, 2)
        ProdinvD(2, 2) = ProdD(1, 1)
        ProdinvD(1, 2) = -ProdD(1, 2)
        ProdinvD(2, 1) = -ProdD(2, 1)
        ProdinvD = ProdinvD / ProddetD
        UhlpD = dcmplx(0.d0, 0.d0)
        VhlpD = dcmplx(0.d0, 0.d0)
        temp = dcmplx(0.d0, 0.d0)
        Diff = dcmplx(0.d0, 0.d0)
        do nl = 1, 2
            do j = 1, Ndim
                UhlpD(j, nl) = GrD(j, P(nl))
                VhlpD(nl, j) = - Op_K%Delta(nl, 1) * GrD(P(1), j) - Op_K%Delta(nl, 2) * GrD(P(2), j)
            enddo
            VhlpD(nl, P(1)) = VhlpD(nl, P(1)) + Op_K%Delta(nl, 1)
            VhlpD(nl, P(2)) = VhlpD(nl, P(2)) + Op_K%Delta(nl, 2)
        enddo
        temp = matmul(UhlpD, ProdinvD)
        Diff = matmul(temp, VhlpD)
        GrD = GrD - Diff

        return
    end function ratioK_fermion

    subroutine GlobalK_prop_L(PropU, PropD, log_ratio_fermion, sigma_new, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        real(kind=8), intent(inout) :: log_ratio_fermion
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_new
        integer, intent(in) :: nt
        integer :: ii
        real(kind=8) :: ratio

        do ii = 2*Lq, 1, -1
            ratio = ratioK_fermion(PropU%Gr, PropD%Gr, sigma_new, ii, nt)
            log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
        enddo
        call Op_K%mmult_L(PropU%Gr, Latt, sigma_new, nt, 1)
        call Op_K%mmult_L(PropD%Gr, Latt, sigma_new, nt, 1)
        call Op_K%mmult_R(PropU%Gr, Latt, sigma_new, nt, -1)
        call Op_K%mmult_R(PropD%Gr, Latt, sigma_new, nt, -1)
        call Op_K%mmult_L(PropU%UUL, Latt, sigma_new, nt, 1)
        call Op_K%mmult_L(PropD%UUL, Latt, sigma_new, nt, 1)
        return
    end subroutine GlobalK_prop_L

    subroutine GlobalK_prop_R(PropU, PropD, log_ratio_fermion, sigma_new, nt)
        class(Propagator), intent(inout) :: PropU, PropD
        real(kind=8), intent(inout) :: log_ratio_fermion
        real(kind=8), dimension(2*Lq, Ltrot), intent(in) :: sigma_new
        integer, intent(in) :: nt
        integer :: ii
        real(kind=8) :: ratio

        call Op_K%mmult_R(PropU%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_R(PropD%Gr, Latt, NsigL_K%sigma, nt, 1)
        call Op_K%mmult_L(PropU%Gr, Latt, NsigL_K%sigma, nt, -1)
        call Op_K%mmult_L(PropD%Gr, Latt, NsigL_K%sigma, nt, -1)
        do ii = 1, 2*Lq
            ratio = ratioK_fermion(PropU%Gr, PropD%Gr, sigma_new, ii, nt)
            log_ratio_fermion = log_ratio_fermion + log(max(ratio, RATIO_EPS))
        enddo
        call Op_K%mmult_R(PropU%UUR, Latt, sigma_new, nt, 1)
        call Op_K%mmult_R(PropD%UUR, Latt, sigma_new, nt, 1)
        return
    end subroutine GlobalK_prop_R
end module GlobalK_mod