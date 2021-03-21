!!****if* source/Simulation/SimulationMain/WindTunnel/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Mach 3 wind tunnel
!!  problem.
!!
!!  References:  Emery, A. E., 1968, JCP, 2, 306
!!               Woodward, P. and Colella, P., 1984, JCP, 54, 115
!!
!!
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY: sim_TAmbient, sim_rhoAmbient, sim_gamma, &
     &  sim_smallP, sim_smallX, sim_pAmbient, sim_pResev, eos_singlespeciesa, &
     & eos_singlespeciesz, sim_dcrit, sim_dexit, sim_lopt, sim_angle, sim_obslen
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkPtr, Grid_releaseBlkPtr


  implicit none

#include "constants.h"
#include "Flash.h"
  real     ::  r, phi
  integer,intent(IN) :: blockID
  real,pointer :: solnData(:,:,:,:)
  integer :: istat
  real :: rho_zone, velx_zone, vely_zone, velz_zone, temp_zone, &
       ener_zone, ekin_zone, eint_zone, pres_zone
  real, parameter :: Pi = 3.1415927

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis  
  integer  ::  i, j, k, n, jLo, jHi
  integer  ::  ii, jj, kk
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac
  logical :: gcell = .true.
  real :: dexit, dcrit, lopt, angle, obslen, ypos, xpos, line, anglerad



!===============================================================================







! In this problem the initial conditions are spatially uniform.
  
  
  !rho_zone = sim_rhoAmbient
  temp_zone = sim_TAmbient
  pres_zone = sim_pAmbient
  rho_zone = (pres_zone/temp_zone)/(8.3107e7/eos_singlespeciesa)  !the factor here is kb/(A*amu) in cgs units
  
  vely_zone = 0.0
  velx_zone = 0.0
  velz_zone = 0.0


  ! Compute the gas energy and set the gamma-values needed for
  ! the equation of state.
  ekin_zone = 0.5 * (velx_zone**2 + vely_zone**2 + velz_zone**2)

  eint_zone = pres_zone / (sim_gamma-1.)
  eint_zone = eint_zone / rho_zone
  ener_zone = eint_zone + ekin_zone
  ener_zone = max(ener_zone, sim_smallP)


  call Grid_getBlkPtr(blockID, solnData, CENTER)
#if NSPECIES > 0
  solnData(SPECIES_BEGIN,:,:,:) =  1.0-(NSPECIES-1)*sim_smallX
  solnData(SPECIES_BEGIN+1:SPECIES_END,:,:,:) =     sim_smallX
#endif

  ! store the variables in the block's unk data
  solnData(DENS_VAR,:,:,:) = rho_zone
  solnData(TEMP_VAR,:,:,:) = temp_zone
!  solnData(PRES_VAR,:,:,:) = pres_zone
  solnData(ENER_VAR,:,:,:) = ener_zone
#ifdef EINT_VAR
  solnData(EINT_VAR,:,:,:) = eint_zone
#endif
  solnData(GAMC_VAR,:,:,:) = sim_gamma
  solnData(GAME_VAR,:,:,:) = sim_gamma


  solnData(VELX_VAR,:,:,:) = velx_zone
  solnData(VELY_VAR,:,:,:) = vely_zone
  solnData(VELZ_VAR,:,:,:) = velz_zone

  call Grid_releaseBlkPtr(blockID, solnData, CENTER)
 
 
   call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat); xCoord = 0.0
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat); yCoord = 0.0
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat); zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
  !
  !     For each cell
  !  
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_getBlkPtr(blockId,solnData)
#endif
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     ! Find a real difference between z's if problem is >= 3D
     if (NDIM > 2) then
        if (k .eq. 1) then
           dzz = zCoord(2) - zCoord(1) 
        else
           dzz = zCoord(k) - zCoord(k-1) 
        endif
     ! Otherwise this problem is <= 2D, so dzz is meaningless
     else
       dzz = 0.0
     endif
     zz = zCoord(k)
     
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        ! Find a real difference between y's if problem is >= 2D
        if (NDIM > 1) then
           if (j .eq. 1) then
              dyy = yCoord(2) - yCoord(1) 
           else
              dyy = yCoord(j) - yCoord(j-1) 
           endif
        ! Otherwise this problem is <= 1D, so dyy is meaningless
        else
          dyy = 0.0
        endif
        yy = yCoord(j)
        

        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
           xx = xCoord(i)
           if (i .eq. 1) then
              dxx = xCoord(2) - xCoord(1) 
           else
              dxx = xCoord(i) - xCoord(i-1) 
           endif
           
           

           
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k

!!!!!!! THIS BLOCK INITIALISES THE INITIAL CONDITIONS
#ifdef FL_NON_PERMANENT_GUARDCELLS
			
			#if NSPECIES > 0
			  solnData(SPECIES_BEGIN,:,:,:) =  1.0-(NSPECIES-1)*sim_smallX
			  solnData(SPECIES_BEGIN+1:SPECIES_END,:,:,:) =     sim_smallX
			#endif

		  ! store the variables in the block's unk data
		  solnData(DENS_VAR,:,:,:) = rho_zone
		  solnData(TEMP_VAR,:,:,:) = temp_zone
		!  solnData(PRES_VAR,:,:,:) = pres_zone
		  solnData(ENER_VAR,:,:,:) = ener_zone
#ifdef EINT_VAR
		  solnData(EINT_VAR,:,:,:) = eint_zone
#endif
		  solnData(GAMC_VAR,:,:,:) = sim_gamma
		  solnData(GAME_VAR,:,:,:) = sim_gamma


		  solnData(VELX_VAR,:,:,:) = velx_zone
		  solnData(VELY_VAR,:,:,:) = vely_zone
		  solnData(VELZ_VAR,:,:,:) = velz_zone

#else
           if (NSPECIES > 0) then
              ! putting in the value of the default species
              call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                   axis, 1.0e0-(NSPECIES-1)*sim_smallX)


              !if there is only one species, this loop will not execute
              do n = SPECIES_BEGIN+1, SPECIES_END

                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                      axis, sim_smallX)


              enddo
           end if

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho_zone)
           !call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, ener_zone)    
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velx_zone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vely_zone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velz_zone)


!!!!!! THIS BLOCK INPUTS BOUNDARIES (ie sets BDRY_VAR to be 1.0)   !!! IMPORTANT: FOR BOUNDARIES TO WORK YOU HAVE TO USE THE UNSPLIT HYDRO SOLVER
#ifdef BDRY_VAR           
         
            
            
           !dcrit = 0.1 !throat diameter
           !dexit = 1.0 !exit diameter
           !lopt = 1.5 !length 
           dcrit = sim_dcrit
           dexit = sim_dexit
           lopt = sim_lopt
	   angle = sim_angle
	   obslen = sim_obslen
	    
	   anglerad = angle*Pi/180
           ypos = lopt + obslen*COS(anglerad)
           xpos = dexit/2 - obslen*SIN(anglerad)


           r = dcrit/2 + yCoord(j)*(dexit-dcrit)/(2*lopt)
           !line = lopt + ((xCoord(i) - (dexit/2)) * (ypos - lopt)/(xpos-dexit*0.5))
	   line = ((yCoord(j) - lopt) * ((xpos-dexit*0.5) / (ypos - lopt)) + (dexit * 0.5))


           if ( yCoord(j) <= lopt .and. xCoord(i) >= r ) then
              call Grid_putPointData(blockId, CENTER, BDRY_VAR, EXTERIOR, axis, 1.0)
              call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho_zone*1000)
           else if  (yCoord(j) <= ypos .and. yCoord(j) >= lopt .and. xCoord(i) >= line) then
              call Grid_putPointData(blockId, CENTER, BDRY_VAR, EXTERIOR, axis, 1.0)
              call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho_zone*1000)
           else
              call Grid_putPointData(blockId, CENTER, BDRY_VAR, EXTERIOR, axis, -1.0)
           end if

#endif

#endif
        enddo
     enddo
  enddo
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_releaseBlkPtr(blockID, solnData)
 
 
 
 
#endif
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

 
 
 
 

  return
end subroutine Simulation_initBlock



