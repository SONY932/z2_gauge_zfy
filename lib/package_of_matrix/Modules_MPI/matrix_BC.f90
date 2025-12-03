     MODULE Matrix_BC

       Use Matrix
       Interface Bcast_mat
          module procedure Bcast_mat_2d_R, Bcast_mat_1d_R
       end Interface
       Interface Reduce_sum_mat 
          module procedure Reduce_sum_mat_2d_R, Reduce_sum_mat_1d_R
       end Interface

       Contains

         subroutine  Bcast_mat_2d_R( data_mat, N_BC ) 

           Implicit Real (KIND=8) (A-G,O-Z)
           Implicit Integer (H-N)

           include 'mpif.h'

           type (Mat_R), Dimension(:,:) :: data_mat

           Real (Kind=8), Dimension(:,:), allocatable :: Collect
           
           INTEGER  STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
           
           N_x = Size(data_mat,1) 
           N_y = Size(data_mat,2)
           N_dim = data_mat(1,1)%dim 
           Allocate  ( Collect(N_x,N_y) ) 
           Write(6,*) " In Bcast :: ", N_x, N_y, N_dim
           
           do no1 = 1, N_dim
              do no2 = 1, N_dim
                 If (Irank.eq.N_BC ) Then
                    do nx = 1,N_x
                       do ny = 1,N_y
                          collect(nx,ny) = data_mat(nx,ny)%el(no1,no2)
                       enddo
                    enddo
                 endif
                 Ntmp = N_x*N_y
                 call MPI_BCAST(collect,Ntmp,MPI_REAL8,N_BC,MPI_COMM_WORLD,IERR)
                 if (Irank.ne. N_BC) Then
                    do nx = 1,N_x
                       do ny = 1,N_y
                          data_mat(nx,ny)%el(no1,no2) =   collect(nx,ny) 
                       enddo
                    enddo
                 endif
              enddo
           enddo

           deallocate ( Collect ) 
           
         end subroutine Bcast_mat_2d_R


!++++++++++++

         subroutine  Bcast_mat_1d_R( data_mat, N_BC ) 

           Implicit Real (KIND=8) (A-G,O-Z)
           Implicit Integer (H-N)

           include 'mpif.h'

           type (Mat_R), Dimension(:) :: data_mat

           Real (Kind=8), Dimension(:), allocatable :: Collect
           
           INTEGER  STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
           
           N_x   = Size(data_mat,1) 
           N_dim = data_mat(1)%dim 
           Allocate  ( Collect(N_x) ) 
           Write(6,*) " In Bcast :: ", N_x,  N_dim
           
           do no1 = 1, N_dim
              do no2 = 1, N_dim
                 If (Irank.eq.N_BC ) Then
                    do nx = 1,N_x
                       collect(nx) = data_mat(nx)%el(no1,no2)
                    enddo
                 endif
                 Ntmp = N_x
                 call MPI_BCAST(collect,Ntmp,MPI_REAL8,N_BC,MPI_COMM_WORLD,IERR)
                 if (Irank.ne. N_BC) Then
                    do nx = 1,N_x
                       data_mat(nx)%el(no1,no2) =   collect(nx) 
                    enddo
                 endif
              enddo
           enddo

           deallocate ( Collect ) 
           
         end subroutine Bcast_mat_1d_R

!++++++++++++


         subroutine  Reduce_sum_mat_2d_R( data_mat ) 
           
           ! Averages over all nodes and places result into data_mat. 

           Implicit Real (KIND=8) (A-G,O-Z)
           Implicit Integer (H-N)

           include 'mpif.h'

           type (Mat_R), Dimension(:,:) :: data_mat
           
           Real (Kind=8), Dimension(:,:), allocatable :: Collect, Collect1
           
           INTEGER  STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
           
           N_x = Size(data_mat,1) 
           N_y = Size(data_mat,2)
           N_dim = data_mat(1,1)%dim 
           Allocate  ( Collect(N_x,N_y), Collect1(N_x,N_y) ) 
           Write(6,*) " In Bcast :: ", N_x, N_y, N_dim
           
           do no1 = 1, N_dim
              do no2 = 1, N_dim
                 do nx = 1,N_x
                    do ny = 1,N_y
                       collect(nx,ny) = data_mat(nx,ny)%el(no1,no2)
                    enddo
                 enddo
                 Ntmp = N_x*N_y
                 CALL MPI_REDUCE(Collect,Collect1,Ntmp,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
                 do nx = 1,N_x
                    do ny = 1,N_y
                       data_mat(nx,ny)%el(no1,no2) =   Collect1(nx,ny) / dble(Isize)
                    enddo
                 enddo
              enddo
           enddo

           deallocate ( Collect, Collect1 ) 
           
         end subroutine Reduce_sum_mat_2d_R


!++++++++++++


         subroutine  Reduce_sum_mat_1d_R( data_mat ) 
           
           ! Averages over all nodes and places result into data_mat. 

           Implicit Real (KIND=8) (A-G,O-Z)
           Implicit Integer (H-N)

           include 'mpif.h'

           type (Mat_R), Dimension(:) :: data_mat
           
           Real (Kind=8), Dimension(:), allocatable :: Collect, Collect1
           
           INTEGER  STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
           
           N_x   = Size(data_mat,1) 
           N_dim = data_mat(1)%dim 
           Allocate  ( Collect(N_x), Collect1(N_x) ) 
           Write(6,*) " In Bcast :: ", N_x, N_dim
           
           do no1 = 1, N_dim
              do no2 = 1, N_dim
                 do nx = 1,N_x
                    collect(nx) = data_mat(nx)%el(no1,no2)
                 enddo
                 Ntmp = N_x
                 CALL MPI_REDUCE(Collect,Collect1,Ntmp,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
                 do nx = 1,N_x
                    data_mat(nx)%el(no1,no2) =   Collect1(nx) / dble(Isize)
                 enddo
              enddo
           enddo
           deallocate ( Collect, Collect1 ) 
           
         end subroutine Reduce_sum_mat_1d_R


       end MODULE Matrix_BC

