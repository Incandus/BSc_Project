!!****if* source/Simulation/SimulationMain/WindTunnel/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the wind tunnel with a step problem
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_pAmbient    Initial ambient pressure
!!  sim_rhoAmbient  Initial ambient density
!!  sim_windVel     Inflow velocity (parallel to x-axis)
!!  gamma           the Gamma EOS thing
!!  smallp          minimum for pressure
!!  smallx          minimum for abundance
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"

  call RuntimeParameters_get('eos_singlespeciesa', eos_singlespeciesa)
  call RuntimeParameters_get('eos_singlespeciesz', eos_singlespeciesz)
  call RuntimeParameters_get('sim_pResev', sim_pResev)
  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_TAmbient', sim_TAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('sim_dcrit', sim_dcrit)
  call RuntimeParameters_get('sim_dexit', sim_dexit)
  call RuntimeParameters_get('sim_lopt', sim_lopt)
  call RuntimeParameters_get('sim_angle', sim_angle)
  call RuntimeParameters_get('sim_obslen', sim_obslen)

  call Logfile_stamp("initializing for gasjet + step", 'run_init')
  write (*,*) "flash:  initializing for gas jet + step"

end subroutine Simulation_init
