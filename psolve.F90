! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine calculates new particle concentrations.
!!
!!   The basic form from which the solution is derived is
!!   ( new_value - old_value ) / dtime = source_term - loss_rate*new_value
!!
!!  Modified  Sep-1997  (McKie)
!!    New particle concentrations due to coagulation processes
!!    were moved to the csolve routine.  Csolve is called to
!!    update particle concentrations due to coagulation.
!!    This new psolve now updates particle concentrations due
!!    to the faster calcs of the non-coag microphysical processes.
!!
!! @author Eric Jensen, Bill McKie
!! @version Oct-1995, Sep-1997
subroutine psolve(carma, cstate, iz, ibin, ielem, rc, dtime_in)


  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_planet_mod
  use carma_condensate_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(in)                  :: iz      !! z index
  integer, intent(in)                  :: ibin    !! bin index
  integer, intent(in)                  :: ielem   !! element index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure
  real(kind=f), intent(in)             :: dtime_in !! thread-local dtime for this substep

  ! Local declarations
  integer                        :: igroup    ! group index
  integer                        :: iepart
  integer                        :: igto
  integer                        :: iz_no_sed
  real(kind=f)                   :: ppd       ! particle prodocution rate
  real(kind=f)                   :: pc_nonuc  ! particles - no nucleation
  real(kind=f)                   :: pls       ! particle loss rate
  real(kind=f)                   :: sed_rate
  real(kind=f)                   :: rnuclgtot
  real(kind=f)                   :: dsed


  ! Define current group & particle number concentration element indices
  igroup = igelem(ielem)
  iepart = ienconc(igroup)

  if(do_grow) then

    ! Compute total production rate
    ppd = rnucpe(ibin,ielem,iz) + rhompe(ibin,ielem,iz)  + growpe(ibin,ielem,iz)  + evappe(ibin,ielem,iz) + phochemprod(iz,ibin,ielem)
    !ppd = rnucpe(ibin,ielem,iz)/ xyzmet(iz) + rhompe(ibin,ielem,iz)/ xyzmet(iz)  + growpe(ibin,ielem,iz)/ xyzmet(iz)  + evappe(ibin,ielem,iz)/ xyzmet(iz)  + phochemprod(iz,ibin,ielem)/ xyzmet(iz)
    !write(*,*) 'psolve ppd', iz, ibin, ielem, rnucpe(ibin,ielem,iz), rhompe(ibin,ielem,iz), growpe(ibin,ielem,iz), evappe(ibin,ielem,iz)
    !write(*,*) iz, ibin, ielem, rnucpe(ibin,ielem,iz)
    !write(*,*), iz,ibin,ielem,phochemprod(iz,ibin,ielem)
    ! Sum up nucleation loss rates
    rnuclgtot = rnuclgsum(ibin,igroup,iz)

!    if (rnucpe(ibin,ielem,iz) > 1.e-50_f) then
!    	write(*,*) ibin, ielem, rnucpe(ibin,ielem,iz), 'rnucpe'
!    end if
!
!    if (evappe(ibin,ielem,iz) > 1.e-50_f) then
!    	write(*,*) ibin, ielem, evappe(ibin,ielem,iz), 'evappe'
!    end if
!
!    if (growpe(ibin,ielem,iz) > 1.e-50_f) then
!    	write(*,*) ibin, ielem, growpe(ibin,ielem,iz), 'growpe'
!    end if

    ! Compute total loss rate
    pls = rnuclgtot + growlg(ibin,igroup,iz) + evaplg(ibin,igroup,iz)
    !pls = rnuclg(ibin,igroup,iz,iz)/ xyzmet(iz) + growlg(ibin,igroup,iz)/ xyzmet(iz) + evaplg(ibin,igroup,iz)/ xyzmet(iz)
        !write(*,*) 'psolve pls', iz, ibin, ielem,rnuclgtot,  growlg(ibin,igroup,iz),  evaplg(ibin,igroup,iz), pc(iz, ibin, ielem)

    ! Update net particle number concentration during current timestep
    ! due to production and loss rates.
    pc(iz,ibin,ielem) = (pc(iz,ibin,ielem) + dtime_in * ppd) / (ONE + pls * dtime_in)
        !write(*,*) 'psolve pc', iz, ibin, ielem,ppd,pls,pc(iz, ibin, ielem)

   ! if (ibin .eq. 1) then
!	write(*,*) 'psolve', evaplg(1,ielem,iz), iz, ielem
 !   endif

    ! Figure out how many particles were produced from nucleation. This is just
    ! for statistics and is done as a total for the step, not per substep.

    !pc_nonuc = (pc(iz,ibin,ielem) + dtime_in * (ppd - &
     !   rnucpe(ibin,ielem,iz) - rhompe(ibin,ielem,iz))) / (ONE + (pls) * dtime_in)
    !pc_nucl(iz,ibin,ielem) = pc_nucl(iz,ibin,ielem) + (pc(iz,ibin,ielem) - pc_nonuc)

    pc_nucl(iz,ibin,ielem) = (pc_nucl(iz,ibin,ielem) + dtime_in * (rnucpe(ibin,ielem,iz) + rhompe(ibin,ielem,iz))) / (ONE + rnuclgtot * dtime_in)

    pc_psolve(iz,ibin,ielem) = pc(iz,ibin,ielem)                                        !PETER
  end if

  ! Prevent particle concentrations from dropping below SMALL_PC
  call smallconc(carma, cstate, iz, ibin, ielem, rc)

  rnucpeup(ibin,ielem,iz) = rnucpe(ibin,ielem,iz)	!PETER

  !  Return to caller with new particle number concentrations.
  return
end
