!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single incompressible flow solver and corresponding time tracker
   type(incomp),      public :: fs
   type(timetracker), public :: time
   type(sgsmodel),    public :: sgs
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: SR
   
   !> Fluid viscosity
   real(WP) :: visc
   
contains
   
   
   !> Function that localizes the top (y+) of the domain
   function top_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function top_locator
   
   
   !> Function that localizes the bottom (y-) of the domain
   function bottom_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function bottom_locator
   
   
   !> Function that localizes the bottom (y-) of the domain
   !> This version limits the zone to a small patch
   function bottom_locator_small(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin.and.sqrt(pg%xm(i)**2+pg%zm(k)**2).le.pg%dx(i)) isIn=.true.
   end function bottom_locator_small
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR(1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Create an incompressible flow solver with bconds
      create_solver: block
         use ils_class,    only: rbgs,amg,pcg_amg,pcg_parasail,gmres,gmres_pilut,smg,pfmg,pcg
         use incomp_class, only: dirichlet,convective,neumann,clipped_neumann
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Define boundary conditions
         call fs%add_bcond(name='outflow',type=dirichlet      ,locator=bottom_locator_small,face='y',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='inflow' ,type=clipped_neumann,locator=   top_locator,face='y',dir=+1,canCorrect=.true. )
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_amg,implicit_ils=gmres)
      end block create_solver
      
      
      ! Create an LES model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs
         
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=fs%cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use incomp_class, only: bcond
         type(bcond), pointer :: my_bc
         integer :: n,i,j,k
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Apply Dirichlet at the bottom
         call fs%get_bcond('outflow',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            fs%V(i,j,k)=-1.0_WP
         end do
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
         call fs%correct_mfr()
         ! Calculate interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block initialize_velocity
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='test')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('div',fs%div)
         call ens_out%add_scalar('sgs_visc',sgs%visc)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Reset here fluid properties
         fs%visc=visc
         
         ! Turbulence modeling
         call fs%get_strainrate(Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)
         resU=fs%rho
         call sgs%get_visc(dt=time%dtold,rho=resU,Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)
         where (sgs%visc.lt.-fs%visc)
            sgs%visc=-fs%visc
         end where
         !fs%visc=fs%visc+sgs%visc
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%correct_mfr()
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho
            fs%V=fs%V-time%dt*resV/fs%rho
            fs%W=fs%W-time%dt*resW/fs%rho
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation
