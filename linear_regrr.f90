program ntolog
      implicit none
      character(len=80):: filename,outfile
      integer :: i,j,ierror,n,dat,end_n,n_old
      integer,parameter:: runit=10
      real(kind=8),allocatable,dimension(:) :: x,y,lnX,lnY
      real(kind=8) :: sum_x,sum_y,sum_x2,sum_xy,sum_y2
      real(kind=8) :: sum_lnX,sum_lnY,lnX2,lnXY,lnY2,sum_lnX2,sum_lnY2,sum_lnXY
      real(kind=8) :: slope,intercept,regr
      real(kind=8) :: alpha,diff,r2_d,r2_alph 
      write(*,*) "Input file name"
      read(*,*) filename
      write(*,*) "Output file name"
      read(*,*) outfile
      write(*,*) "Please provide the extent of data to be fitted for calculations in ns."
      read(*,*) dat
!-------------------------------------------------------------------------------
      call read_lines(filename,n)
      write(*,*) "This ",filename," contains ",n," lines."
      write(*,*) "Doing fitting on 1 ns data in segments of 100 ps."
!-------------------------------------------------------------------------------
   !   sum_x = 0.0 ; sum_y = 0.0 ; sum_x2 = 0.0 ; sum_xy = 0.0
   !   sum_y2 = 0.0 ; sum_lnY2 = 0.0
   !   sum_lnX = 0.0 ; sum_lnY = 0.0 ; sum_lnX2 = 0.0 ; sum_lnXY = 0.0
!-------------------------------------------------------------------------------
      dat = dat * 1000 +40   ! dat *1000 ps data
      open(unit=12,file=outfile,status='replace',action='write',iostat=ierror)
      open(unit=runit,file=filename,status='old',action='read',iostat=ierror)
      allocate(x(dat),y(dat))
      allocate(lnX(dat),lnY(dat))

rdat: do j = 1,dat
      read(runit,*,iostat=ierror) x(j) , y(j)
      lnY(j) = log(y(j)) ; lnX(j) = log(x(j))
      write(13,100,iostat=ierror) lnX(j) , lnY(j)
      enddo rdat
      n = 1
      end_n = 40
      do j = 1,3! total 10 segments of 100ps each
!-------------------------------------------------------------------------------
      sum_x = 0.0 ; sum_y = 0.0 ; sum_x2 = 0.0 ; sum_xy = 0.0
      sum_y2 = 0.0 ; sum_lnY2 = 0.0
      sum_lnX = 0.0 ; sum_lnY = 0.0 ; sum_lnX2 = 0.0 ; sum_lnXY = 0.0
!-------------------------------------------------------------------------------
      do i = n,end_n ! updation after fitting of each 100ps data
     ! read(runit,*,iostat=ierror) x(j) , y(j)
      sum_x = sum_x + x(i) ; sum_y = sum_y + y(i)
     ! write(*,*) j,i,x(i),y(i),sum_x,sum_y
      sum_x2 = sum_x2 + (x(i)*x(i))
      sum_y2 = sum_y2 + (y(i)*y(i))
      sum_xy = sum_xy + (x(i)*y(i))
     ! lnY(i) = log(y(i)) ; lnX(i) = log(x(i))
      if (x(i).ne.0.0 .and. y(i).ne.0.0) then
      sum_lnX = sum_lnX + lnX(i) ; sum_lnY = sum_lnY + lnY(i)
     ! write(*,*) "log========================================"
     ! write(*,*) j,i,lnX(i),lnY(i),sum_lnX,sum_lnY
     ! write(*,*) "==========================================="
      lnX2 = lnX(i)*lnX(i) ; lnXY = lnX(i)*lnY(i) ; lnY2 = lnY(i)*lnY(i)
      sum_lnX2 = sum_lnX2 + lnX2
      sum_lnY2 = sum_lnY2 + lnY2
      sum_lnXY = sum_lnXY + lnXY
     ! write(13,100,iostat=ierror) lnX(i) , lnY(i)
      endif
      enddo
     ! n = end_n + 1
     ! end_n = end_n +100
     ! if (end_n .gt. 1000) exit   ! do while loop
     ! enddo
!-------------------------------------------------------------------------------
      call fit_and_test(sum_x,sum_y,sum_x2,sum_xy,sum_y2,slope,intercept,regr)
      r2_d = regr*regr
      diff = slope
      call fit_and_test(sum_lnX,sum_lnY,sum_lnX2,sum_lnXY,sum_lnY2,slope,intercept,regr)
      r2_alph = regr*regr
      alpha = slope
      write(12,101) j,diff,r2_d,alpha,r2_alph
      !n_old = n
      n = end_n + 1
      end_n = end_n +500 
      if (end_n .gt. 1040) exit   ! do while loop
      enddo
!-------------------------------------------------------------------------------
      close(10)
      close(12)
100  format(3x,f10.4,6x,f10.6)
101  format(1x,i3,4(3x,E10.4))
! endif openif
 deallocate(x,y,lnX,lnY)
end program  ntolog
!===============================================================================
      subroutine fit_and_test(sum_1,sum_2,sum_1sq,sum_12,sum_2sq,m,b,r)
!-------------------------------------------------------------------------------
      real(kind=8) :: n_data,sum_1,sum_2,sum_1sq,sum_12,sum_2sq,m,b,r
      real(kind=8) :: product_12,product_11,product_22,product_n_2sq
      real(kind=8) :: product_n_sum12,product_n_1sq,sub_n12_12,sub_n1sq_11
!-------------------------------------------------------------------------------
      n_data = 100
      product_12 = sum_1*sum_2 ; product_11 = sum_1*sum_1
      product_n_2sq = n_data*sum_2sq ; product_22 = sum_2*sum_2
      product_n_sum12 = n_data*sum_12 ; product_n_1sq = n_data*sum_1sq
      sub_n12_12 = product_n_sum12 - product_12
      sub_n1sq_11 = product_n_1sq - product_11
!-------------------------------------------------------------------------------
      m = (sub_n12_12) / (sub_n1sq_11) ! slope
      b = ((sum_2)*product_11) - ((sum_1)*(sum_12)) / sub_n1sq_11 ! y-intercept
      r = sub_n12_12 / sqrt((sub_n1sq_11)*(product_n_2sq-product_22))
!-------------------------------------------------------------------------------
              end subroutine fit_and_test
!===============================================================================
      subroutine read_lines(fname,nlines)
!-------------------------------------------------------------------------------
      integer:: nlines,ierro
      character(len=80) :: fname
!===============================================================================
      nlines = 0
!      runit = 10
      open(unit=13,file=fname,status='old',action='read',iostat=ierro)
      openif: if (ierro==0) then
      do
        read(13,*,iostat=ierro)
        if (ierro.ne.0) exit
        nlines = nlines + 1
      enddo
      endif openif
      close(13)
      write(*,*) nlines
              end subroutine read_lines
!end program ntolog
