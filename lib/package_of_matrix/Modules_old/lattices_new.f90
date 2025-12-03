      Module Lattices_new

        Use Matrix
        Type Lattice
           Integer          :: N, Ns
           Integer, pointer :: list(:,:), invlist(:,:), nnlist(:,:,:), listk(:,:), &
                &              invlistk(:,:),  &
                &              a1_p(:), a2_p(:),  L1_p(:), L2_p(:), imj(:,:)
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
           module procedure FT_K_to_R, FT_K_to_R_Mat, FT_K_to_R_C 

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

          Real (Kind=8), dimension(:), allocatable :: xk_p, b1_p, b2_p, BZ1_p, BZ2_p, b_p
          Integer,       dimension(:), allocatable :: i_p, nd_p, i1_p, a_p, imj_p
          
          Real (Kind=8), allocatable :: Mat(:,:), Mat_inv(:,:)

          ndim = size(L1_p)
          allocate (Latt%L2_p(ndim), Latt%L1_p(ndim), Latt%a1_p(ndim) , Latt%a2_p(ndim), &
               &    Latt%b1_p(ndim), Latt%b2_p(ndim), Latt%BZ1_p(ndim), Latt%BZ2_p(ndim) )
          Zero = 1.E-6
          Latt%L1_p = L1_p
          Latt%L2_p = L2_p
          Latt%a1_p = a1_p
          Latt%a2_p = a2_p

          
          Allocate ( i_p(ndim), nd_p(ndim), i1_p(ndim) , a_p(ndim), imj_p(ndim) )
          ! Count the number of lattice points. 
          Ns   =  4 * nint(Xnorm(L1_p) )
          Ns1  =  4 * nint(Xnorm(L2_p) )
          if (Ns1 .gt. Ns) Ns = Ns1
          nc = 0
          L  = 0
          do i1 = -Ns,Ns
             do i2 = -Ns,Ns
                i_p  = i1*a1_p + i2*a2_p
                L_f = 1
                do i = 1,4
                   if (i.eq.1) a_p = L2_p 
                   if (i.eq.2) a_p = L1_p  
                   if (i.eq.3) a_p = L2_p - L1_p 
                   if (i.eq.4) a_p = L2_p + L1_p
                   if  (  Iscalar( i_p, a_p )    .le.  xnorm(a_p)**2/2.d0 + Zero   .and.   &
                        & Iscalar( i_p, a_p )    .ge. -xnorm(a_p)**2/2.d0 + Zero    ) then
                      L_f = L_f * 1
                   else
                      L_f = 0
                   endif
                enddo
                if (L_f .eq. 1) then   
                   nc = nc + 1
                   if (abs(i1).gt.L)  L = abs(i1)
                   if (abs(i2).gt.L)  L = abs(i2)
                endif
             enddo
          enddo
          LQ = nc
          Latt%Ns = LQ
          Latt%N  = LQ
          L       = 2*L
          Write(6,*) L, LQ
          Allocate ( Latt%List(LQ,ndim), Latt%Invlist(-L:L, -L:L ) )

          !Setting up real space lattice 
          nc = 0
          do i1 = -L,L
             do i2 = -L,L
                i_p  = i1*a1_p + i2*a2_p
                L_f = 1
                do i = 1,4
                   if (i.eq.1) a_p = L2_p 
                   if (i.eq.2) a_p = L1_p  
                   if (i.eq.3) a_p = L2_p - L1_p 
                   if (i.eq.4) a_p = L2_p + L1_p
                   if  (  Iscalar( i_p, a_p )    .le.  xnorm(a_p)**2/2.d0 + Zero   .and.   &
                        & Iscalar( i_p, a_p )    .ge. -xnorm(a_p)**2/2.d0 + Zero    ) then
                      L_f = L_f * 1
                   else
                      L_f = 0
                   endif
                enddo
                if (L_f .eq. 1) then   
                   !write(10,"(I4,2x,I4)") i_p(1), i_p(2)
                   nc = nc + 1
                   Latt%list(nc,1) = i1
                   Latt%list(nc,2) = i2
                   Latt%invlist(i1, i2 ) = nc
                endif
             enddo
          enddo

          !Setting up k-space lattice
          Allocate ( b1_p(ndim), b2_p(ndim), xk_p(ndim), BZ1_p(ndim), BZ2_p(ndim),b_p(ndim) )
          Allocate ( Latt%Listk(LQ,ndim), Latt%Invlistk(-L:L, -L:L) )
          pi   = acos(-1.d0)
          ! Setup the 2X2 matrix to determine  BZ1_p, BZ2_p
          Allocate ( Mat(2 , 2), Mat_inv( 2 , 2 ) ) 
          Mat(1,1) = dble(a1_p(1))
          Mat(1,2) = dble(a1_p(2))
          Mat(2,1) = dble(a2_p(1))
          Mat(2,2) = dble(a2_p(2))
          X = Mat(1,1)*Mat(2,2) - Mat(2,1)*Mat(1,2)
          Mat_inv(1,1) =  Mat(2,2)/X
          Mat_inv(2,2) =  Mat(1,1)/X
          Mat_inv(1,2) = -Mat(1,2)/X
          Mat_inv(2,1) = -Mat(2,1)/X
          BZ1_p(1)      = 2.d0*pi*Mat_inv(1,1)
          BZ1_p(2)      = 2.d0*pi*Mat_inv(2,1)
          BZ2_p(1)      = 2.d0*pi*Mat_inv(1,2)
          BZ2_p(2)      = 2.d0*pi*Mat_inv(2,2)
          Latt%BZ1_p = BZ1_p
          Latt%BZ2_p = BZ2_p
          

          X =  2.d0*pi / ( Iscalar(BZ1_p,L1_p) * Iscalar(BZ2_p,L2_p) -   &
               &           Iscalar(BZ2_p,L1_p) * Iscalar(BZ1_p,L2_p)   )
          X = abs(X)

          b1_p = X*( Iscalar(BZ2_p,L2_p) * BZ1_p - Iscalar(BZ1_p,L2_p) * BZ2_p ) 
          b2_p = X*( Iscalar(BZ1_p,L1_p) * BZ2_p - Iscalar(BZ2_p,L1_p) * BZ1_p )

          Deallocate ( Mat,  Mat_inv ) 


          Latt%b1_p  = b1_p
          Latt%b2_p  = b2_p

          !write(6,*) 'BZ1: ', BZ1_p(1), BZ1_p(2)
          !write(6,*) 'BZ2: ', BZ2_p(1), BZ2_p(2)
          !write(6,*) 'b1 : ', b1_p(1), b1_p(2)
          !write(6,*) 'b2 : ', b2_p(1), b2_p(2)
          nc = 0
          do m = -L,L
             do n = -L,L
                xk_p = dble(m)* b1_p + dble(n) * b2_p
                L_f = 1
                do i = 1,4
                   if (i.eq.1) b_p = BZ2_p 
                   if (i.eq.2) b_p = BZ1_p  
                   if (i.eq.3) b_p = BZ2_p - BZ1_p 
                   if (i.eq.4) b_p = BZ2_p + BZ1_p
                   if  (  Iscalar( xk_p, b_p )    .le.  xnorm(b_p)**2/2.d0 + Zero   .and.   &
                        & Iscalar( xk_p, b_p )    .ge. -xnorm(b_p)**2/2.d0 + Zero    ) then
                      L_f = L_f * 1
                   else
                      L_f = 0
                   endif
                enddo
                if (L_f .eq. 1) then   
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
          Allocate ( Latt%nnlist(LQ,-1:1,-1:1) )
          
          do nr = 1, Latt%N
             do nd1 = -1,1
                do nd2 = -1,1
                   nd_p = nd1*a1_p + nd2*a2_p
                   i_p  = Latt%list(nr,1)*Latt%a1_p + Latt%list(nr,2)*a2_p  + nd_p
                   call npbc(i1_p, i_p , Latt%L1_p, Latt%L2_p)
                   call npbc(i_p , i1_p, Latt%L1_p, Latt%L2_p)
                   nnr1 =  nint ( Iscalar(Latt%BZ1_p,i_p) / (2.d0*pi) )
                   nnr2 =  nint ( Iscalar(Latt%BZ2_p,i_p) / (2.d0*pi) )
                   nnr  = Latt%invlist(nnr1,nnr2)
                   Latt%nnlist(nr,nd1,nd2) = nnr
                enddo
             enddo
          enddo

          !Setup imj 
          If (LQ  .lt. 1000 ) then 
             Allocate ( Latt%imj(LQ,LQ) )
             
             do nr = 1, Latt%N
                i_p = Latt%list(nr,1)*Latt%a1_p + Latt%list(nr,2)*a2_p  + nd_p
                do nr1 = 1,Latt%N
                   i1_p = Latt%list(nr1,1)*Latt%a1_p + Latt%list(nr1,2)*a2_p  + nd_p
                   imj_p = i_p - i1_p
                   call npbc(i1_p  , imj_p , Latt%L1_p, Latt%L2_p)
                   call npbc(imj_p , i1_p, Latt%L1_p, Latt%L2_p)
                   imj_1 =  nint ( Iscalar(Latt%BZ1_p,imj_p) / (2.d0*pi) )
                   imj_2 =  nint ( Iscalar(Latt%BZ2_p,imj_p) / (2.d0*pi) )
                   imj   = Latt%invlist(imj_1,imj_2)
                   Latt%imj(nr,nr1) = imj 
                enddo
             enddo
          endif
          
          deallocate ( i_p, nd_p, i1_p , imj_p)
          deallocate ( b1_p, b2_p, xk_p, BZ1_p, BZ2_p )


        end subroutine MAKE_LATTICE

!********
        subroutine npbc_I(nr_p, n_p, L1_p, L2_p) 
      
          Implicit Real (Kind=8) (A-G,O-Z)
          Implicit Integer (H-N)

          integer, dimension(:) ::  nr_p, n_p, L1_p, L2_p
          
          integer, dimension(:), allocatable :: x_p

          Zero = 1.E-8
          nr_p = n_p 
          ndim = size(nr_p)

          allocate (x_p(ndim))

          do  i = 1,4
             if (i.eq.1) x_p = L2_p
             if (i.eq.2) x_p = L1_p
             if (i.eq.3) x_p = L2_p - L1_p
             if (i.eq.4) x_p = L2_p + L1_p
             
             X = dble(Iscalar(nr_p,x_p))/(Xnorm(x_p)**2)
             if (X .gt.   0.5+Zero  ) nr_p = nr_p - x_p
             if (X .lt.  -0.5+Zero  ) nr_p = nr_p + x_p   
          enddo
          
          deallocate(x_p)

        end subroutine npbc_I


        subroutine npbc_R(nr_p, n_p, L1_p, L2_p) 
      
          Implicit Real (Kind=8) (A-G,O-Z)
          Implicit Integer (H-N)

          Real (Kind=8), dimension(:) ::  nr_p, n_p, L1_p, L2_p
          
          Real (Kind=8), dimension(:), allocatable :: x_p
          ndim = size(nr_p)

          allocate (x_p(ndim))
          Zero = 1.E-8
          nr_p = n_p 

          do  i = 1,4
             if (i.eq.1) x_p = L2_p
             if (i.eq.2) x_p = L1_p
             if (i.eq.3) x_p = L2_p - L1_p
             if (i.eq.4) x_p = L2_p + L1_p
             
             X =  Iscalar(nr_p,x_p)/(Xnorm(x_p)**2)
             if (X .gt.   0.5+Zero  ) nr_p = nr_p - x_p
             if (X .lt.  -0.5+Zero  ) nr_p = nr_p + x_p   
          enddo
          
          deallocate(x_p)

        end subroutine npbc_R

!********
        integer Function Inv_K(XK_P,Latt) 
          
          Implicit None
          Real (Kind=8)  :: XK_P(2)
          Type (Lattice) :: Latt
          
          Real (Kind=8)     :: XK1_P(2), X, Zero
          Integer :: nk 

          Zero = 1.D-10
          nk = 1
          do 
             XK1_P = Latt%listk(nk,1)*latt%b1_p + Latt%listk(nk,2)*latt%b2_p
             if (Xnorm(XK1_P - XK_P)  < Zero ) then
                Inv_K = nk
                exit
             elseif (nk < Latt%N) then 
                nk = nk + 1
             else
                write(6,*) 'Error in Inv_K Lattice_new'
                stop
             endif
          enddo

        end Function Inv_K
        
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

          Open (Unit=55,file="Latt_info", status = "unknown")
          write(55,*) ' Reciprocal vector 1: ', Latt%BZ1_p(1), Latt%BZ1_p(2)
          write(55,*) ' Reciprocal vector 2: ', Latt%BZ2_p(1), Latt%BZ2_p(2)
          write(55,*) ' Latt       vector 1: ', Latt%a1_p(1), Latt%a1_p(2)
          write(55,*) ' Latt       vector 2: ', Latt%a2_p(1), Latt%a2_p(2)
          close(55)
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
             do nd1 = -1,1
                do nd2 = -1,1
                   nd_p =   nd1*Latt%a1_p + nd2*Latt%a2_p
                   nnr = Latt%nnlist(n,nd1,nd2)
                   !Write(6,*) 'nnr : ', nnr
                   i_p = Latt%list(nnr,1)*Latt%a1_p + Latt%list(nnr,2)*Latt%a2_p 
                   write(58,"('I+(',I4,',',I4,')=',2x,I6,2x,I6)") nd_p(1),nd_p(2),i_p(1), i_p(2)
                enddo
             enddo
          enddo
          close(56)
          close(57)
          close(58)
        end subroutine Print_latt
!******* 
        subroutine FT_K_to_R_Mat( Xin_K, Xout_R, Latt)
          
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
        end subroutine FT_K_to_R_Mat



        subroutine FT_K_to_R( Xin_K, Xout_R, Latt)
          
          Implicit Real (Kind=8) (A-G,O-Z)
          Implicit Integer (H-N)
          
          Type (Lattice)                             :: Latt
          Real (Kind=8), Dimension(:,:)              :: Xin_K, Xout_R 
          Real (Kind=8)                              :: XK_p(2)
          Integer                                    :: IR_p(2)

          Ltrot  = size(Xin_K,2  ) 
          LQ     = Latt%N

          !Write(6,*) 'Ltrot, norb ', Ltrot, norb
          !Write(6,*) Xin_K(1,1)%el(1,1)
          !Write(6,*) Xin_K(Latt%N,Ltrot)%el(1,1)
          

          do nt = 1,Ltrot
             do nr = 1,LQ
                IR_p =  Latt%list(nr,1)*Latt%a1_p + Latt%list(nr,2)*Latt%a2_p  
                X_MAT = 0.d0
                do nk = 1,LQ
                   XK_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
                   X_MAT = X_MAT + cos(Iscalar(XK_p,IR_p))*Xin_K(nk,nt)
                enddo
                Xout_R(nr,nt) = X_MAT/dble(LQ)
             enddo
          enddo

        end subroutine FT_K_to_R


        subroutine FT_K_to_R_C( Xin_K, Xout_R, Latt)
          
          Implicit Real (Kind=8) (A-G,O-Z)
          Implicit Integer (H-N)
          
          Type (Lattice)                             :: Latt
          Complex (Kind=8), Dimension(:,:)           :: Xin_K, Xout_R 
          Complex (Kind=8)                           :: Z
          Real    (Kind=8)                           :: XK_p(2)
          Integer                                    :: IR_p(2)

          Ltrot  = size(Xin_K,2  ) 
          LQ     = Latt%N

          !Write(6,*) 'Ltrot, norb ', Ltrot, norb
          !Write(6,*) Xin_K(1,1)%el(1,1)
          !Write(6,*) Xin_K(Latt%N,Ltrot)%el(1,1)
          

          do nt = 1,Ltrot
             do nr = 1,LQ
                IR_p =  Latt%list(nr,1)*Latt%a1_p + Latt%list(nr,2)*Latt%a2_p  
                Z = cmplx(0.d0,0.d0)
                do nk = 1,LQ
                   XK_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
                   Z = Z + cmplx(cos(Iscalar(XK_p,IR_p)),0.d0)*Xin_K(nk,nt)
                enddo
                Xout_R(nr,nt) = Z/cmplx(dble(LQ),0.d0)
             enddo
          enddo

        end subroutine FT_K_to_R_C


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

        
      end Module Lattices_new

      
