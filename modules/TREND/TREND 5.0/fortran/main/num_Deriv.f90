!
!Developed by Sven Pohl - 15.03.2019
!This is fortran module for the calculation of numerical derivatives
!The modules uses an external function to calculate the the function derivatives


      module numerical_differentiation

      use module_all_types
      implicit none

      double precision , dimension(11), parameter :: order_1_1 =  (/-1D0/2D0,0D0,1D0/2,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/)
      double precision , dimension(11), parameter :: order_1_2 =  (/ 1D0/12D0 , -2D0/3D0 , 0D0 , 2D0/3D0 ,-1D0/12D0 , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_1_3 =  (/ -1D0/60D0 , 3D0/20D0 , -3D0/4D0 , 0D0 , 3D0/4D0 ,-3D0/20D0 , 1D0/60D0 , 0D0 , 0D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_1_4 = (/ 1D0/280D0 , -4D0/105D0 , 1D0/5D0 , -4D0/5D0 , 0D0 , 4D0/5D0 , -1D0/5D0 , 4D0/105D0 , -1D0/280D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_1_5 = (/ 0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0 ,0D0 /)
      double precision , dimension(11,5), parameter :: order_1 = (/order_1_1,order_1_2,order_1_3,order_1_4,order_1_5/)

      double precision , dimension(11), parameter :: order_2_1 = (/ 1D0 , -2D0 , 1D0 , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_2_2 = (/ -1D0/12D0 , 4D0/3D0 , -5D0/2D0 , 4D0/3D0 , -1D0/12D0 ,0D0 , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_2_3 = (/ 1D0/90D0 , -3D0/20D0 , 3D0/2D0 , -49D0/18D0 , 3D0/2D0 , -3D0/20D0 , 1D0/90D0   , 0D0 , 0D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_2_4 = (/ -1D0/560D0 , 8D0/315D0 , -1D0/5D0 , 8D0/5D0 , -205D0/72D0 , 8D0/5D0 , -1D0/5D0 , 8D0/315D0 , -1D0/560D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_2_5 = (/ 0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, 0D0 , 0D0 /)
      double precision , dimension(11,5), parameter :: order_2 = (/order_2_1,order_2_2,order_2_3,order_2_4,order_2_5/)

      double precision , dimension(11), parameter :: order_3_1 = (/ 0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0 /)
      double precision , dimension(11), parameter :: order_3_2 = (/ -1D0/2D0 , 1D0 , 0D0 , -1D0 , 1D0/2D0 , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_3_3 = (/ 1D0/8D0 , -1D0 , 13D0/8D0 , 0D0 , -13D0/8D0 , 1D0 , -1D0/8D0 , 0D0 , 0D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_3_4 = (/ -7D0/240D0 , 3D0/10D0 , -169D0/120D0 , 61D0/30D0 , 0D0 , -61D0/30D0 , 169D0/120D0 , -3D0/10D0 , 7D0/240D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_3_5 = (/ 0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0 /)
      double precision , dimension(11,5), parameter :: order_3 = (/order_3_1,order_3_2,order_3_3,order_3_4,order_3_5/)

      double precision , dimension(11), parameter :: order_4_1 = (/ 0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0 /)
      double precision , dimension(11), parameter :: order_4_2 = (/ 1D0 , -4D0 , 6D0 , -4D0 , 1D0  , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_4_3 = (/ -1D0/6D0 , 2D0 , -13D0/2D0 , 28D0/3D0 , -13D0/2D0 , 2D0 , -1D0/6D0 , 0D0 , 0D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_4_4 = (/ 7D0/240D0 , -2D0/5D0 , 169D0/60D0 ,   -122D0/15D0 , 91D0/8D0 , -122D0/15D0 , 169D0/60D0 , -2D0/5D0 , 7D0/240D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_4_5 = (/ 0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0 /)
      double precision , dimension(11,5), parameter :: order_4 = (/order_4_1,order_4_2,order_4_3,order_4_4,order_4_5/)


      double precision , dimension(11), parameter :: order_5_1 =(/ 0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/)
      double precision , dimension(11), parameter :: order_5_2 =(/ 0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/)
      double precision , dimension(11), parameter :: order_5_3 =(/ -1D0/2D0 , 2D0 , -5D0/2D0 , 0D0 , 5D0/2D0 , -2D0 , 1D0/2D0   , 0D0 , 0D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_5_4 =(/ 1D0/6D0 , -3D0/2D0 , 13D0/3D0 , -29D0/6D0 , 0D0 , 29D0/6D0 , -13D0/3D0 , 3D0/2D0 , -1D0/6D0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_5_5 =(/ -13D0/288D0 , 19D0/36D0 , -87D0/32D0 , 13D0/2D0 , -323D0/48D0 , 0D0 , 323D0/48D0 , -13D0/2D0 , 87D0/32D0 , -19D0/36D0 , 13D0/288D0 /)
      double precision , dimension(11,5), parameter :: order_5 = (/order_5_1,order_5_2,order_5_3,order_5_4,order_5_5/)

      double precision , dimension(11), parameter :: order_6_1 =(/ 0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/)
      double precision , dimension(11), parameter :: order_6_2 =(/ 0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/)
      double precision , dimension(11), parameter :: order_6_3 = (/ 1D0 , -6D0 , 15D0 , -20D0 , 15D0 , -6D0 , 1D0   , 0d0 , 0d0 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_6_4 = (/ -1D0/4D0 , 3D0 , -13D0 , 29D0 , -75D0/2D0 , 29D0 , -13D0 , 3D0 , -1D0/4 , 0D0 , 0D0 /)
      double precision , dimension(11), parameter :: order_6_5 = (/ 13D0/240D0 , -19D0/24D0 , 87D0/16D0 , -39D0/2D0 , 323D0/8D0 , -1023D0/20D0 , 323D0/8D0 , -39D0/2D0 , 87D0/16D0 , -19D0/24D0 , 13D0/240D0 /)
      double precision , dimension(11,5), parameter :: order_6 = (/order_6_1,order_6_2,order_6_3,order_6_4,order_6_5/)

      !MATRIX WITH ALL COEFFICIENTS
      !(nr,prec,deriv)
      double precision , dimension(11,5,6), parameter ::deriv_fac = (/order_1,order_2,order_3,order_4,order_5,order_6/)

      double precision, parameter:: n_del = 1.D-2
      double precision, parameter:: n_del2 = 1.D-2
      

      double precision, dimension(11), parameter :: stenc_cat_1 = (/-n_del,0d0,n_del,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
      double precision, dimension(11), parameter :: stenc_cat_2 = (/-2d0*n_del,-n_del,0d0,1d0*n_del,2d0*n_del,0d0,0d0,0d0,0d0,0d0,0d0/)
      double precision, dimension(11), parameter :: stenc_cat_3 = (/-3d0*n_del,-2d0*n_del,-n_del,0d0,n_del,2d0*n_del,3d0*n_del,0d0,0d0,0d0,0d0/)
      double precision, dimension(11), parameter :: stenc_cat_4 = (/-4d0*n_del,-3d0*n_del,-2d0*n_del,-n_del,0d0,n_del,2d0*n_del,3d0*n_del,4d0*n_del,0d0,0d0/)
      double precision, dimension(11), parameter :: stenc_cat_5 = (/-5d0*n_del,-4d0*n_del,-3d0*n_del,-2d0*n_del,-n_del,0d0,n_del,2d0*n_del,3d0*n_del,4d0*n_del,5d0*n_del/)

      double precision, dimension(11,5), parameter:: stenc_from_cat = (/ stenc_cat_1,stenc_cat_2,stenc_cat_3,stenc_cat_4,stenc_cat_5 /)

      contains

      double precision function generic_numerical_derivative(gl,f,x,y,dx,dy,nr_stenc,prec_cat)
 

      implicit none

      double precision, external  :: f !function which should be derived
      double precision :: x,y !input variables
      integer:: dx,dy !direction of differentiation
      !modify variables you may need for your external function!
      type(type_gl):: gl

      !------------------------------!
      integer:: loop_d
      integer:: nr_stenc !Number of stencil points
      integer:: prec_cat
      !( 3 Stencil -> Cat. 1)
      !( 5 Stencil -> Cat. 2)
      !( 7 Stencil -> Cat. 3)
      !( 9 Stencil -> Cat. 4)
      !(11 Stencil -> Cat. 5)
      ! Available are: 1st derivative(Cat. 1 2 3 4    )
      !                2st derivative(Cat. 1 2 3 4    )
      !                3th derivative(Cat.   2 3 4    )
      !                4th derivative(Cat.   2 3 4    )
      !                5th derivative(Cat.     3 4 5  )
      !                6th derivative(Cat.     3 4 5  )
      ! Dx number of derivatives in x direction
      ! Dy number of derivatives in y direction

      generic_numerical_derivative = 0.d0

      if(dx .ne. 0 .and. dy .eq. 0) then !derivative with respect to x variable
          do loop_d=1,nr_stenc
              generic_numerical_derivative = deriv_fac(loop_d,prec_cat,dx)*(f(gl, x+stenc_from_cat(loop_d,prec_cat),y)) + generic_numerical_derivative
          end do
          generic_numerical_derivative=generic_numerical_derivative/(n_del)**dx
          return
      elseif(dy .ne. 0 .and. dx .eq. 0) then !derivative with respect to y variable
          do loop_d=1,nr_stenc
              generic_numerical_derivative = deriv_fac(loop_d,prec_cat,dy)*(f(gl,x,y+stenc_from_cat(loop_d,prec_cat))) + generic_numerical_derivative
          end do
          generic_numerical_derivative=generic_numerical_derivative/(n_del)**dy
          return
      elseif(dy .eq. 1 .and. dx .eq. 1) then !mixed derivative in x and y direction
          generic_numerical_derivative= ( f(gl,x-n_del,y) + f(gl,x+n_del,y) + f(gl,x,y-n_del)+ f(gl,x,y+n_del) - 4d0*f(gl,x,y))/(n_del**2)
      end if

      end function

      end module