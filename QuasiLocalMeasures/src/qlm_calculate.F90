#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_calculate (CCTK_ARGUMENTS)
  use cctk
  use qlm_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer   :: my_proc
  integer   :: hn
  integer, parameter :: qlm_root_proc = 0
  
  character :: msg*1000
  character :: slabel*2, ilabel*8
  character(len=200) :: odir
  integer   :: nchars
 
  logical   :: did_allocate
  logical   :: do_surface
  logical   :: owns_surface
  
  did_allocate = .false.
  
  my_proc   = CCTK_MyProc (cctkGH)

  if (calculate_every <= 0) return

  if (mod(cctk_iteration, calculate_every) /= 0) then
     if (my_proc == qlm_root_proc .and. veryverbose/=0) then
        write (msg, '("Skipping quasi-local quantities at iteration ",i8)') cctk_iteration
        call CCTK_INFO (msg)
     end if
     return
  end if

  if (num_surfaces > 0) then
     call allocate_variables (int(maxntheta), int(maxnphi))
     did_allocate = .true.
  end if

  do hn = 1, num_surfaces

     owns_surface = my_proc == qlm_root_proc
     do_surface = .true.

     ! start calculations already? 
     if (cctk_time < begin_qlm_calculations_after(hn)) do_surface = .false.
     
     if (owns_surface .and. (verbose/=0 .or. veryverbose/=0)) then
        if (do_surface) then
           write (msg, '("Calculating quasi-local quantities for surface ",i4)') hn - 1
        else
           write (msg, '("Skipping quasi-local quantities for surface ",i4)') hn - 1
        end if
        call CCTK_INFO (msg)
     end if
     
     if (owns_surface .and. do_surface) then
        if (surface_index(hn) == -1 .and. CCTK_EQUALS(surface_name(hn), "")) then
           qlm_calc_error(hn) = 1
           qlm_have_valid_data(hn) = 0
           qlm_have_killing_vector(hn) = 0
           do_surface = .false.
        end if
     end if
     
     if (owns_surface .and. do_surface) then
        call qlm_import_surface (CCTK_PASS_FTOF, hn)
        if (qlm_calc_error(hn) /= 0) do_surface = .false.
     endif
     
     if (owns_surface .and. do_surface) then
        call qlm_set_coordinates (CCTK_PASS_FTOF, hn)
     end if
     
     call qlm_interpolate (CCTK_PASS_FTOF, hn, owns_surface .and. do_surface)
     
     if (owns_surface .and. do_surface) then
        if (qlm_calc_error(hn) /= 0) goto 9999
        
        call qlm_calc_tetrad (CCTK_PASS_FTOF, hn)
        call qlm_calc_newman_penrose (CCTK_PASS_FTOF, hn)
        call qlm_calc_weyl_scalars (CCTK_PASS_FTOF, hn)
        call qlm_calc_twometric (CCTK_PASS_FTOF, hn)
        if (CCTK_EQUALS(killing_vector_method, "axial")) then
           call qlm_killing_axial (CCTK_PASS_FTOF, hn)
        else if (CCTK_EQUALS(killing_vector_method, "eigenvector")) then
           call qlm_killing_transport (CCTK_PASS_FTOF, hn)
           if (qlm_calc_error(hn) /= 0) goto 9999
           if (qlm_have_killing_vector(hn) /= 0) then
              call qlm_killing_normalise (CCTK_PASS_FTOF, hn)
           end if
        else if (CCTK_EQUALS(killing_vector_method, "gradient")) then
           call qlm_killing_gradient (CCTK_PASS_FTOF, hn)
           if (qlm_have_killing_vector(hn) /= 0) then
              call qlm_killing_normalise (CCTK_PASS_FTOF, hn)
           end if
        else
           call CCTK_WARN (0, "internal error")
        end if
        if (qlm_calc_error(hn) /= 0) goto 9999
        if (qlm_have_killing_vector(hn) /= 0) then
           call qlm_killing_test (CCTK_PASS_FTOF, hn)
           call qlm_calc_coordinates (CCTK_PASS_FTOF, hn)
        end if
        call qlm_calc_3determinant (CCTK_PASS_FTOF, hn)
        call qlm_analyse (CCTK_PASS_FTOF, hn)
        if (qlm_have_killing_vector(hn) /= 0) then
           call qlm_multipoles (CCTK_PASS_FTOF, hn)
           call qlm_multipoles_normalise (CCTK_PASS_FTOF, hn)
        end if

        if (output_vtk_every /= 0) then
           if (mod (cctk_iteration, output_vtk_every) == 0) then
              write (slabel, '(I2.2)') hn
              write (ilabel, '(I8.8)') cctk_iteration
              call CCTK_ParameterValString (nchars, "out_dir", "IOUtil", odir) 
              call qlm_output_vtk (CCTK_PASS_FTOF, hn, &
                   odir(1:nchars) // '/surface' // slabel // '_' // ilabel // &
                   '.vtk')
           end if
        end if
        
9999    continue
        
        if (qlm_timederiv_order(hn) < 2) then
           call CCTK_WARN (2, "There were not enough past time levels available for accurate calculations")
        end if
     end if
     
  end do
 
  if (did_allocate) then
     ! release 2D arrays
     call deallocate_variables
  end if
 
  call qlm_broadcast (cctkGH)

  if (veryverbose/=0) then
     call CCTK_INFO ("Done.")
  end if
  
end subroutine qlm_calculate
