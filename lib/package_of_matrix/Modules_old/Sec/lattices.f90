      Module Lattices

        Use Matrix
        Type Lattice
           Integer          :: N, Ns
           Integer, pointer :: list(:,:), invlist(:,:), nnlist(:,:), listk(:,:), invlistk(:,:),  &
                &              a1_p(:), a2_p(:),  L1_p(:), L2_p(:)
           Real (Kind=8), pointer :: b1_p(:), b2_p(:), BZ1_p(:), BZ2_p(:)
        end Type Lattice

        Interface Iscalar
           module procedure Iscalar_II, Iscalar_IR, Iscalar_RR
        end Interface
        Interface npbc
           module procedure npbc_I, npbc_R
        end Interface
        Interface Xnorm
           module procedure Xnorm_I, Xnorm_R
        end Interface
        Interface Fourier_K_to_R
           module procedure FT_K_to_R
        end Interface
        Interface Fourier_R_to_K
           module procedure FT_R_to_K
        end Interface

      Contains
        
        subroutine Make_lattice(L1_p, L2_p, a1_p, a2_p, Latt) 

          ! This is for a general tilted square lattice defined by the vector  
          ! a1, a2 
          ! L1_p = (L1_x, L1_y)  and L2_p = (-L1_y, L1_x)  
          
          Implicit Real (Kind=8) (A-G,O-Z)
          Implicit Integer (H-N)
          
          Integer, dimension(:) :: L1_p, L2_p, a1_p, a2_p
          Type (Lattice) :: Latt

          Real (Kind=8), dimension(:), allocatable :: xk_p, b1_p, b2_p, BZ1_p, BZ2_p
          Integer,       dimension(:), allocatable :: i_p, nd_p, i1_p

          ndim = size(L1_p)
          allocate (Latt%L2_p(ndim), Latt%L1_p(ndim), Latt%a1_p(ndim) , Latt%a2_p(ndim), &
               &    Latt%b1_p(ndim), Latt%b2_p(ndim), Latt%BZ1_p(ndim), Latt%BZ2_p(ndim) ) 
          Zero = 1.E-6
          Latt%L1_p = L1_p
          Latt%L2_p = L2_p
          Latt%a1_p = a1_p
          Latt%a2_p = a2_p
          If (Iscalar(a1_p,a2_p).ne.0) then 
             write(6,*) 'Error in Lattice. The lattice vectors are not orthogonal to each other'
             stop
          endif
          Ns   =  4*nint(Xnorm(L1_p) )
          Latt%Ns = Ns
          Allocate ( Latt%List(Ns*Ns,ndim), Latt%Invlist(-Ns:Ns, -Ns:Ns) )
          Allocate ( i_p(ndim), nd_p(ndim), i1_p(ndim) )

          !Write(6,*) 'Ns, Ndim : ', Ns, Ndim
          !Setting up real space lattice
          nc = 0
          do i1 = -Ns,Ns
             do i2 = -Ns,Ns
                i_p  = i1*a1_p + i2*a2_p
                if     ( Iscalar(i_p,L1_p) .gt. 0                                  .and.   &
                     &   Iscalar(i_p,L2_p) .gt. 0                                  .and.   &
                     &  dble(Iscalar(i_p,L1_p))/Xnorm(L1_p).le.(Xnorm(L1_p)+Zero)  .and.   &
                     &  dble(Iscalar(i_p,L2_p))/Xnorm(L2_p).le.(Xnorm(L2_p)+Zero)        ) &
                     &                                                                   then
                   !write(10,"(I4,2x,I4)") i_p(1), i_p(2)
                   nc = nc + 1
                   Latt%list(nc,1) = i1
                   Latt%list(nc,2) = i2
                   Latt%invlist(i1, i2 ) = nc
                 endif
             enddo
          enddo
          Latt%N = nc

          !Setting up k-space lattice
          Allocate ( b1_p(ndim), b2_p(ndim), xk_p(ndim), BZ1_p(ndim), BZ2_p(ndim) )
          Allocate ( Latt%Listk(Ns*Ns,ndim), Latt%Invlistk(-Ns:Ns, -Ns:Ns) )
          pi   = acos(-1.d0)
          b1_p = 2.0*pi*dble( L1_p ) / ( xnorm(L1_p)**2 )
          b2_p = 2.0*pi*dble( L2_p ) / ( xnorm(L2_p)**2 )
          Latt%b1_p  = b1_p
          Latt%b2_p  = b2_p
          BZ1_p      = 2.d0*pi*dble(a1_p)/(xnorm(a1_p)**2)
          BZ2_p      = 2.d0*pi*dble(a2_p)/(xnorm(a2_p)**2)
          Latt%BZ1_p = BZ1_p
          Latt%BZ2_p = BZ2_p
          nc = 0
          do m = -Ns,Ns
             do n = -Ns,Ns
                xk_p = m * b1_p + n*b2_p
                !write(55,*) n,m, xk_p
                if     ( Iscalar(xk_p,BZ1_p) .gt. Zero                               .and.   &
                 &  Iscalar(xk_p,BZ2_p) .gt. Zero                                    .and.   &
                 &  Iscalar(xk_p,BZ1_p)/(Xnorm(BZ1_p)).le.(2.0*pi/xnorm(a1_p)+Zero)  .and.   &
                 &  Iscalar(xk_p,BZ2_p)/(Xnorm(BZ2_p)).le.(2.0*pi/xnorm(a2_p)+Zero)        ) &
                 &                                                                    then
                   !write(11,"(F14.7,2x,F14.7)")  xk_p(1), xk_p(2)
                   nc = nc + 1
                   Latt%listk(nc,1) = m
                   Latt%listk(nc,2) = n
                   Latt%invlistk(m,n) = nc
                endif
             enddo
          enddo
          If (nc.ne.Latt%N) Then 
             write(6,*) 'Error ', nc, Latt%N
             stop
          endif

          !Setup nnlist
          Allocate ( Latt%nnlist(Ns*Ns,8) )
          
          do nr = 1, Latt%N
             do n = 1,8
                if (n.eq.1) nd_p =  a1_p
                if (n.eq.2) nd_p =  a2_p
                if (n.eq.3) nd_p = -a1_p
                if (n.eq.4) nd_p = -a2_p
                if (n.eq.5) nd_p =   a1_p + a2_p
                if (n.eq.6) nd_p =  -a1_p + a2_p
                if (n.eq.7) nd_p =  -a1_p - a2_p
                if (n.eq.8) nd_p =   a1_p - a2_p
                i_p  = Latt%list(nr,1)*Latt%a1_p + Latt%list(nr,2)*a2_p  + nd_p
                call npbc(i1_p, i_p , Latt%L1_p, Latt%L2_p)
                call npbc(i_p , i1_p, Latt%L1_p, Latt%L2_p)
                nnr1 =  nint (dble(Iscalar(i_p,Latt%a1_p))/Xnorm(Latt%a1_p)**2 )
                nnr2 =  nint (dble(Iscalar(i_p,Latt%a2_p))/Xnorm(Latt%a2_p)**2 )
                nnr  = Latt%invlist(nnr1,nnr2)
                Latt%nnlist(nr,n) = nnr
             enddo
          enddo

          deallocate ( i_p, nd_p, i1_p )
          deallocate ( b1_p, b2_p, xk_p, BZ1_p, BZ2_p )

        end subroutine MAKE_LATTICE

!********
        subroutine npbc_I(nr_p, n_p, L1_p, L2_p) 
      
          Implicit Real (Kind=8) (A-G,O-Z)
          Implicit Integer (H-N)

          integer, dimension(:) ::  nr_p, n_p, L1_p, L2_p
          
          Zero = 1.E-8
          nr_p = n_p 
          X = dble(Iscalar(nr_p,L1_p))/(Xnorm(L1_p)**2)
          if (X .gt. 1.0 + Zero ) nr_p = nr_p - L1_p
          if (X .lt. Zero       ) nr_p = nr_p + L1_p   

          X = dble(Iscalar(nr_p,L2_p))/(Xnorm(L2_p)**2)
          if (X .gt. 1.0 + Zero ) nr_p = nr_p - L2_p
          if (X .lt. Zero       ) nr_p = nr_p + L2_p   

        end subroutine npbc_I


        subroutine npbc_R(nr_p, n_p, L1_p, L2_p) 
      
          Implicit Real (Kind=8) (A-G,O-Z)
          Implicit Integer (H-N)

          Real (Kind=8), dimension(:) ::  nr_p, n_p, L1_p, L2_p
          
          Zero = 1.E-8
          nr_p = n_p 
          X = Iscalar(nr_p,L1_p)/(Xnorm(L1_p)**2)
          if (X .gt. 1.0 + Zero ) nr_p = nr_p - L1_p
          if (X .lt. Zero       ) nr_p = nr_p + L1_p   

          X = Iscalar(nr_p,L2_p)/(Xnorm(L2_p)**2)
          if (X .gt. 1.0 + Zero ) nr_p = nr_p - L2_p
          if (X .lt. Zero       ) nr_p = nr_p + L2_p   

        end subroutine npbc_R
        
!********
        integer function Iscalar_II(i_p, j_p)
          Implicit none
          integer, dimension(:) :: i_p, j_p
          integer i
          
          Iscalar_II = 0
          !write(6,*) size(i_p)
          do i = 1,  size(i_p)
            ! write(6,*) i
             Iscalar_II = Iscalar_II + i_p(i)*j_p(i)
          enddo
        end function Iscalar_II

!********
        Real (Kind=8)  function Iscalar_IR(x_p, j_p)
          Implicit none
          Real (Kind=8), dimension(:) ::  x_p
          integer, dimension(:) ::  j_p
          integer i
          
          Iscalar_IR = 0.d0
          !write(6,*) size(i_p)
          do i = 1,  size(x_p)
            ! write(6,*) i
             Iscalar_IR = Iscalar_IR + x_p(i)*dble(j_p(i))
          enddo
        end function Iscalar_IR
!********

        Real (Kind=8)  function Iscalar_RR(x_p, y_p)
          Implicit none
          Real (Kind=8), dimension(:) ::  x_p, y_p
          integer i
          
          Iscalar_RR = 0.d0
          do i = 1,  size(x_p)
             Iscalar_RR = Iscalar_RR + x_p(i)*y_p(i)
          enddo
        end function Iscalar_RR

!********
        Real (Kind=8) function Xnorm_I(i_p)
          Implicit none
          integer, dimension(:) :: i_p
          integer :: i

          Xnorm_I = 0.d0
          do i = 1,  size(i_p)
             Xnorm_I = Xnorm_I + dble(i_p(i)*i_p(i))
          enddo
          Xnorm_I = sqrt(Xnorm_I)
        end function Xnorm_I

!********
        Real (Kind=8) function Xnorm_R(x_p)
          Implicit none
          Real (Kind=8), dimension(:) :: x_p
          integer :: i

          Xnorm_R = 0.d0
          do i = 1,  size(x_p)
             Xnorm_R = Xnorm_R + x_p(i)*x_p(i)
          enddo
          Xnorm_R = sqrt(Xnorm_R)
        end function Xnorm_R

!********
        subroutine Print_latt(Latt)
          
          Implicit Real (Kind=8) (A-G,O-Z)
          Implicit Integer (H-N)
          
          Type (Lattice) :: Latt
          Integer          :: i_p(2),nd_p(2)
          Real    (Kind=8) :: x_p(2)

          Open (Unit=56,file="Real_space_latt", status = "unknown")
          Open (Unit=57,file="K_space_latt", status = "unknown")
          Open (Unit=58,file="nn_latt", status = "unknown")
          do n = 1, Latt%n
             i_p = Latt%list(n,1)*Latt%a1_p + Latt%list(n,2)*Latt%a2_p 
             write(56,"(I6,2x,I6)") i_p(1), i_p(2)
             x_p = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p 
             write(57,"(F14.7,2x,F14.7)") x_p(1), x_p(2)
             write(58,*)
             write(58,"('I :',I6,2x,I6)") i_p(1), i_p(2)
             do nn = 1,8
                if (nn.eq.1) nd_p =  Latt%a1_p
                if (nn.eq.2) nd_p =  Latt%a2_p
                if (nn.eq.3) nd_p = -Latt%a1_p
                if (nn.eq.4) nd_p = -Latt%a2_p
                if (nn.eq.5) nd_p =   Latt%a1_p + Latt%a2_p
                if (nn.eq.6) nd_p =  -Latt%a1_p + Latt%a2_p
                if (nn.eq.7) nd_p =  -Latt%a1_p - Latt%a2_p
                if (nn.eq.8) nd_p =   Latt%a1_p - Latt%a2_p

                nnr = Latt%nnlist(n,nn)
                !Write(6,*) 'nnr : ', nnr
                i_p = Latt%list(nnr,1)*Latt%a1_p + Latt%list(nnr,2)*Latt%a2_p 
                write(58,"('I+(',I4,',',I4,')=',2x,I6,2x,I6)") nd_p(1),nd_p(2),i_p(1), i_p(2)
             enddo
          enddo
          close(56)
          close(57)
          close(58)
          
        end subroutine Print_latt

!******* 
        subroutine FT_K_to_R( Xin_K, Xout_R, Latt)
          
          Implicit Real (Kind=8) (A-G,O-Z)
          Implicit Integer (H-N)
          
          Type (Lattice)                             :: Latt
          Type (Mat_R ), Dimension(:,:)              :: Xin_K, Xout_R 
          Real (Kind=8), Dimension(:,:), allocatable :: X_MAT
          Real (Kind=8)                              :: XK_p(2)
          Integer                                    :: IR_p(2)

          Ltrot  = size(Xin_K,2  ) 
          norb   = size(Xin_K(1,1)%el,1) 
          LQ     = Latt%N

          !Write(6,*) 'Ltrot, norb ', Ltrot, norb
          !Write(6,*) Xin_K(1,1)%el(1,1)
          !Write(6,*) Xin_K(Latt%N,Ltrot)%el(1,1)
          
          allocate ( X_MAT(norb,norb) )

          
          do nt = 1,Ltrot
             do nr = 1,LQ
                IR_p =  Latt%list(nr,1)*Latt%a1_p + Latt%list(nr,2)*Latt%a2_p  
                X_MAT = 0.d0
                do nk = 1,LQ
                   XK_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
                   X_MAT = X_MAT + cos(Iscalar(XK_p,IR_p))*Xin_K(nk,nt)%el
                enddo
                Xout_R(nr,nt)%el = X_MAT/dble(LQ)
             enddo
          enddo

          deallocate(X_Mat)
        end subroutine FT_K_to_R


        subroutine FT_R_to_K( Xin_R, Xout_K, Latt)
          
          Implicit Real (Kind=8) (A-G,O-Z)
          Implicit Integer (H-N)
          
          Type (Lattice)                             :: Latt
          Type (Mat_R ), Dimension(:,:)              :: Xin_R, Xout_K 
          Real (Kind=8), Dimension(:,:), allocatable :: X_MAT
          Real (Kind=8)                              :: XK_p(2)
          Integer                                    :: IR_p(2)

          Ltrot  = size(Xin_R,2  ) 
          norb   = size(Xin_R(1,1)%el,1) 
          LQ     = Latt%N

          !Write(6,*) 'Ltrot, norb ', Ltrot, norb
          !Write(6,*) Xin_R(1,1)%el(1,1)
          !Write(6,*) Xin_R(Latt%N,Ltrot)%el(1,1)
          
          allocate ( X_MAT(norb,norb) )

          
          do nt = 1,Ltrot
             do nk = 1,LQ
                XK_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
                X_MAT = 0.d0
                do nr = 1,LQ
                   IR_p =  Latt%list(nr,1)*Latt%a1_p + Latt%list(nr,2)*Latt%a2_p  
                   X_MAT = X_MAT + cos(Iscalar(XK_p,IR_p))*Xin_R(nr,nt)%el
                enddo
                Xout_K(nk,nt)%el = X_MAT
             enddo
          enddo

          deallocate(X_Mat)
        end subroutine FT_R_to_K

        
      end Module Lattices

      
