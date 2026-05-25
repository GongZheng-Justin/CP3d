!  gfortran -ffree-line-length-none -cpp -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fbacktrace -Wall -lm -g -O3 -mcmodel=large -funroll-loops -floop-optimize -fcheck=all refine_stl.f90 -o refine


module m_divide_trangle

contains

subroutine divide_triangle(n_triangle, area, area_target, p1_i, p2_i, p3_i, p1_o, p2_o, p3_o)
   implicit none
   integer, intent(out) :: n_triangle
   real(8), intent(inout) :: area
   real(8), intent(in) :: area_target
   real(8), intent(in) :: p1_i(3), p2_i(3), p3_i(3)
   real(8), pointer, intent(out) :: p1_o(:,:), p2_o(:,:), p3_o(:,:)

   ! local
   integer :: level, i_level, k, k1, k2, i_long, kt
   real(8) :: p1_t(3), p2_t(3), p3_t(3), p_middle(3), l_12, l_23, l_31
   
   n_triangle = 1
   if (associated(p1_o)) deallocate(p1_o)
   if (associated(p2_o)) deallocate(p2_o)
   if (associated(p3_o)) deallocate(p3_o)
   if (area <= area_target) then
      allocate(p1_o(3, 1))
      allocate(p2_o(3, 1))
      allocate(p3_o(3, 1))
      p1_o(:, 1) = p1_i
      p2_o(:, 1) = p2_i
      p3_o(:, 1) = p3_i
      return
   endif
   
   level = 0
   do
      if (area <= area_target) exit
      area = area / 2.0_8
      n_triangle = n_triangle * 2
      level = level + 1
   enddo
   allocate(p1_o(3, n_triangle))
   allocate(p2_o(3, n_triangle))
   allocate(p3_o(3, n_triangle))
   p1_o(:, 1) = p1_i
   p2_o(:, 1) = p2_i
   p3_o(:, 1) = p3_i
   
   kt = 1
   do i_level = 1, level
      do k = 1, kt
         k1 = k;  k2 = k + kt
         p1_t = p1_o(:, k);  p2_t = p2_o(:, k);  p3_t = p3_o(:, k)
         
         l_12 = (p1_t(1) - p2_t(1)) ** 2 + (p1_t(2) - p2_t(2)) ** 2 + (p1_t(3) - p2_t(3)) ** 2
         l_23 = (p2_t(1) - p3_t(1)) ** 2 + (p2_t(2) - p3_t(2)) ** 2 + (p2_t(3) - p3_t(3)) ** 2
         l_31 = (p3_t(1) - p1_t(1)) ** 2 + (p3_t(2) - p1_t(2)) ** 2 + (p3_t(3) - p1_t(3)) ** 2

         i_long = 3
         if (l_12 > l_23) then
            if (l_12 > l_31) i_long = 1
         else
            if (l_23 > l_31) i_long = 2
         endif
         if (i_long == 1) then
            if (p1_t(1) > p2_t(1)) then
               p_middle = (p1_t + p2_t) / 2.0_8
            else
               p_middle = (p2_t + p1_t) / 2.0_8
            endif
            p1_o(:, k1) = p1_t;  p2_o(:, k1) = p_middle;  p3_o(:, k1) = p3_t
            p1_o(:, k2) = p3_t;  p2_o(:, k2) = p_middle;  p3_o(:, k2) = p2_t          
         elseif (i_long == 2) then
            if (p2_t(1) > p3_t(1)) then
               p_middle = (p2_t + p3_t) / 2.0_8
            else
               p_middle = (p3_t + p2_t) / 2.0_8
            endif
            p1_o(:, k1) = p2_t;  p2_o(:, k1) = p_middle;  p3_o(:, k1) = p1_t
            p1_o(:, k2) = p1_t;  p2_o(:, k2) = p_middle;  p3_o(:, k2) = p3_t
         else
            if (p3_t(1) > p1_t(1)) then
               p_middle = (p3_t + p1_t) / 2.0_8
            else
               p_middle = (p1_t + p3_t) / 2.0_8
            endif
            p1_o(:, k1) = p3_t;  p2_o(:, k1) = p_middle;  p3_o(:, k1) = p2_t
            p1_o(:, k2) = p2_t;  p2_o(:, k2) = p_middle;  p3_o(:, k2) = p1_t
         endif
      enddo
      kt = kt * 2
   enddo
endsubroutine divide_triangle

endmodule m_divide_trangle

program main
   use m_divide_trangle
   implicit none

   character(len = *), parameter :: case_name = 'xijiao'
   real(8) :: area_target = 0.8_8 * (0.1_8 / 200.0_8 )**2
   real(8), parameter :: x_off = 0.05_8, y_off = 0.10_8, z_off = 0.05_8
   !character(len = *), parameter :: case_name = 'block'
   !real(8) :: area_target = 0.8_8 * ((12.0_8 / 600) * (1.0_8 / 100) * (6.0_8 / 600) )**(2.0_8 / 3.0_8)
   !real(8), parameter :: x_off = 0.00_8, y_off = 0.00_8, z_off = 0.00_8

   type :: stl_ele
     sequence
     real(4) :: normal(3)
     real(4) :: p1(3)
     real(4) :: p2(3)
     real(4) :: p3(3)
     integer(2) :: marker = int(0, 2)
   endtype stl_ele
   
   character(len = 1) :: dummy1, dummy2
   character(len = 80) :: string
   type(stl_ele) :: element

   integer :: ui, uo, ierr, it, i, n_triangle, k, n_tot
   real(8) :: normal(3), p1(3), p2(3), p3(3), area
   real(8), pointer :: p1_o(:,:), p2_o(:,:), p3_o(:,:)
   
   p1_o => null()
   p2_o => null()
   p3_o => null()
   
   open (newunit = ui, file = case_name // '.stl', form = 'formatted', status = 'old', &
         action = 'read', iostat = ierr)
   open (newunit = uo, file = case_name // '_out.stl', form = 'unformatted', status = 'replace', &
         action = 'write', iostat = ierr, access = 'stream')
   read (unit = ui, fmt = '(A)', iostat = ierr) string
   write(unit = uo, iostat = ierr) string
   write(unit = uo, iostat = ierr) int(0, 4)
         
   i = 1
   n_tot = 0
   do
      i = i + 1
      it = mod(i, 7)
      select case (it)
      case (2)
         read (unit = ui, fmt = *, iostat = ierr) dummy1, dummy2, normal
      case (4)
         read (unit = ui, fmt = *, iostat = ierr) dummy1, p1
      case (5)
         read (unit = ui, fmt = *, iostat = ierr) dummy1, p2
      case (6)
         read (unit = ui, fmt = *, iostat = ierr) dummy1, p3
         p1(1) = p1(1) + x_off;  p1(2) = p1(2) + y_off;  p1(3) = p1(3) + z_off
         p2(1) = p2(1) + x_off;  p2(2) = p2(2) + y_off;  p2(3) = p2(3) + z_off
         p3(1) = p3(1) + x_off;  p3(2) = p3(2) + y_off;  p3(3) = p3(3) + z_off    
                  
         call clc_area(p1, p2, p3, area)
         call divide_triangle(n_triangle, area, area_target, p1, p2, p3, p1_o, p2_o, p3_o)
         n_tot = n_tot + n_triangle

         do k = 1, n_triangle
            p1 = p1_o(:, k);  p2 = p2_o(:, k);  p3 = p3_o(:, k)            
            element%normal = real(normal, 4)
            element%p1 = real(p1, 4)
            element%p2 = real(p2, 4)
            element%p3 = real(p3, 4)
            write(unit = uo, iostat = ierr) element
            
         enddo
      case default

         read (unit = ui, fmt = *, iostat = ierr)
      endselect
      if (ierr /= 0) exit
   enddo

   write(unit = uo, pos = 81) n_tot
   close(unit = ui, iostat = ierr)
   close(unit = uo, iostat = ierr)

   !print*, 'OK', n_line, n_triangle
endprogram main


subroutine clc_area(tp1, tp2, tp3, area)
   implicit none
   real(8), intent(in) :: tp1(3), tp2(3), tp3(3)
   real(8), intent(out) :: area
     
   ! locals
   real(8) :: pt1(3), pt2(3), pt3(3)

   pt1 = tp2 - tp1;  pt2 = tp3 - tp1
   pt3(1) = pt1(2) * pt2(3) - pt1(3) * pt2(2)
   pt3(2) = pt1(3) * pt2(1) - pt1(1) * pt2(3)
   pt3(3) = pt1(1) * pt2(2) - pt1(2) * pt2(1)
   area = 0.5_8 * sqrt(pt3(1) * pt3(1) + pt3(2) * pt3(2) + pt3(3) * pt3(3))  
endsubroutine clc_area
