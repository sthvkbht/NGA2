!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP,SP
   use geometry,          only: cfg,D
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use df_class,          only: dfibm
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use pardata_class,     only: pardata
   implicit none
   private
   
   !> Get an an incompressible solver, pressure solver, and corresponding time tracker
   type(incomp) :: fs
   type(dfibm)  :: df
   type(fft2d)  :: ps
   type(timetracker) :: time
   
   !> Implicit solver
   logical :: use_implicit
   type(ddadi) :: vs
   
   !> SGS model
   logical :: use_sgs
   type(sgsmodel) :: sgs
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(event) :: ens_evt
   type(ensight) :: ens_out
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,ibmfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP) :: visc,omega

   
contains
   
  !> Function that localizes the left (x-) of the domain
  function left_of_domain(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    logical :: isIn
    isIn=.false.
    if (i.eq.pg%imin) isIn=.true.
  end function left_of_domain

  !> Function that localizes the left (x+) of the domain
   function right_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
    end function right_of_domain

   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max time',time%tmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
         use incomp_class, only: dirichlet,slip
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Define boundary conditions
         call fs%add_bcond(name='bottom',type=dirichlet,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='top'   ,type=slip     ,locator=right_of_domain,face='x',dir=+1,canCorrect=.false.)
         ! Configure pressure solver
         ps=fft2d(cfg=cfg,name='Pressure',nst=7)
         ! Check if implicit velocity solver is used
         call param_read('Use implicit solver',use_implicit)
         if (use_implicit) then
            ! Configure implicit solver
            vs=ddadi(cfg=cfg,name='Velocity',nst=7)
            ! Finish flow solver setup
            call fs%setup(pressure_solver=ps,implicit_solver=vs)
         else
            ! Finish flow solver setup
            call fs%setup(pressure_solver=ps)
         end if
      end block create_flow_solver
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)) 
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
       end block allocate_work_arrays


       ! Initialize our direct forcing solver
       initialize_df: block
         use mathtools, only: twoPi,arctan
         real(SP), allocatable :: normal(:,:), v1(:,:), v2(:,:), v3(:,:)
         real(SP), allocatable :: centroid(:,:), area(:)
         real(SP) :: a(3),b(3),c(3)
         integer :: i,j,k,np,nx,nt,iunit,ierr
         integer(kind=2) :: attr
         real(WP) :: Dp,Lp,dx,x,theta,r
         real(WP) :: scale_factor,shift_x
         character(len=80) :: stlfile,header
         call param_read('STL file',stlfile)
         df=dfibm(cfg=cfg,name='IBM')
         df%can_move=.true.
         call param_read('Angular speed',omega)
         ! Initialize marker particles
         if (df%cfg%amRoot) then
            ! Read the STL file
            open(newunit=iunit,file=trim(stlfile),access="stream",form="unformatted",action="read",status="old",iostat=ierr)
            read(iunit) header
            read(iunit) np
            allocate(normal(3,np),v1(3,np),v2(3,np),v3(3,np))
            allocate(centroid(3,np),area(np))
            do i=1,np
               read(iunit) normal(:,i)
               read(iunit) v1(:,i)
               read(iunit) v2(:,i)
               read(iunit) v3(:,i)
               read(iunit) attr
               a=v1(:,i); b=v2(:,i); c=v3(:,i)
               centroid(:,i) = (a+b+c)/3.0_SP
               area(i) = 0.5_SP * sqrt( &
                    ((b(2)-a(2))*(c(3)-a(3)) - (b(3)-a(3))*(c(2)-a(2)))**2 + &
                    ((b(3)-a(3))*(c(1)-a(1)) - (b(1)-a(1))*(c(3)-a(3)))**2 + &
                    ((b(1)-a(1))*(c(2)-a(2)) - (b(2)-a(2))*(c(1)-a(1)))**2 )
            end do
            a(1)=0.5_SP*(MINVAL(centroid(1,:))+MAXVAL(centroid(1,:)))
            a(2)=0.5_SP*(MINVAL(centroid(2,:))+MAXVAL(centroid(2,:)))
            a(3)=0.5_SP*(MINVAL(centroid(3,:))+MAXVAL(centroid(3,:)))
            print*, "COM is ",a(1),a(2),a(3)
            print*, "Max and min area is ", maxval(area), minval(area)
            do i=1,np
               centroid(:,i) = centroid(:,i) - a
            end do
            close(iunit)

            ! Initialize marker particles
            call df%resize(np)
            scale_factor=1.0_WP
            shift_x=0.50_WP*df%cfg%x(df%cfg%imax+1)
            do i=1,np
               df%p(i)%pos(3) = real(centroid(1,i), WP) * scale_factor
               df%p(i)%pos(2) = real(centroid(2,i), WP) * scale_factor
               df%p(i)%pos(1) = real(centroid(3,i), WP) * scale_factor + shift_x
               ! Set various parameters for the marker
               df%p(i)%id  =1
               df%p(i)%vel =0.0_WP
               ! Assign element area
               df%p(i)%dA=real(area(i),WP)*scale_factor*scale_factor
               ! Assign outward normal vector
               df%p(i)%norm(3) = real(normal(1,i),WP)
               df%p(i)%norm(2) = real(normal(2,i),WP)
               df%p(i)%norm(1) = real(normal(3,i),WP)
               ! Locate the particle on the mesh
               df%p(i)%ind=df%cfg%get_ijk_global(df%p(i)%pos,[df%cfg%imin,df%cfg%jmin,df%cfg%kmin])
               ! Activate the particle
               df%p(i)%flag=0
            end do
            deallocate(v1,v2,v3,normal,centroid,area)
         endif

         ! Communicate
         call df%sync()

         ! Get initial volume fraction
         call df%update_VF()

         ! All processes initialize IBM objects
         call df%setup_obj()

         if (df%cfg%amRoot) then
            print*,"===== Direct Forcing Setup Description ====="
            print*,'Number of marker particles', df%np
            print*,'Number of IBM objects', df%nobj
         end if
       end block initialize_df


       ! Create partmesh object for marker particle output
       create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=1,nvec=1,name='ibm')
         pmesh%varname(1)='area'
         pmesh%vecname(1)='velocity'
         call df%update_partmesh(pmesh)
         do i=1,df%np_
            pmesh%var(1,i)=df%p(i)%dA
            pmesh%vec(:,1,i)=df%p(i)%vel
         end do
       end block create_pmesh


       ! Initialize our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         ! Initial fields
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         call fs%get_div()
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         call param_read('Use SGS model',use_sgs)
         if (use_sgs) sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='mixer')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('markers',pmesh)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('Gib',cfg%Gib)
         call ens_out%add_scalar('VF',df%VF)
         call ens_out%add_scalar('pressure',fs%P)
         if (use_sgs) call ens_out%add_scalar('visc_sgs',sgs%visc)
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
         ! Create IBM monitor
         ibmfile=monitor(amroot=df%cfg%amRoot,name='ibm')
         call ibmfile%add_column(time%n,'Timestep number')
         call ibmfile%add_column(time%t,'Time')
         call ibmfile%add_column(df%VFmin,'VF min')
         call ibmfile%add_column(df%VFmax,'VF max')
         call ibmfile%add_column(df%Fx,'Fx')
         call ibmfile%add_column(df%Fy,'Fy')
         call ibmfile%add_column(df%Fz,'Fz')
         call ibmfile%add_column(df%np,'Marker number')
         call ibmfile%add_column(df%Umin,'Marker Umin')
         call ibmfile%add_column(df%Umax,'Marker Umax')
         call ibmfile%add_column(df%Vmin,'Marker Vmin')
         call ibmfile%add_column(df%Vmax,'Marker Vmax')
         call ibmfile%add_column(df%Wmin,'Marker Wmin')
         call ibmfile%add_column(df%Wmax,'Marker Wmax')
         call ibmfile%write()
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

         ! Update IBM
         call df%rotate(dt=time%dtmid,omega=omega,axis='x')
         
         ! Turbulence modeling
         if (use_sgs) then
            sgs_modeling: block
               use sgsmodel_class, only: vreman
               resU=fs%rho
               call fs%get_gradu(gradU)
               call sgs%get_visc(type=vreman,dt=time%dtold,rho=resU,gradu=gradU)
               where (cfg%Gib.gt.0.0_WP) sgs%visc=0.0_WP
               fs%visc=visc+sgs%visc
            end block sgs_modeling
         end if
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW

            ! Finish update
            if (use_implicit) then
               ! Form implicit residuals
               call fs%solve_implicit(time%dt,resU,resV,resW)
               ! Apply these residuals
               fs%U=2.0_WP*fs%U-fs%Uold+resU
               fs%V=2.0_WP*fs%V-fs%Vold+resV
               fs%W=2.0_WP*fs%W-fs%Wold+resW
            else
               ! Apply these residuals
               fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho
               fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho
               fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho
            end if

            ! Apply IB forcing to enforce BC at the pipe walls
            ibforcing: block
              integer :: i,j,k
              do k=fs%cfg%kmin_,fs%cfg%kmax_
                 do j=fs%cfg%jmin_,fs%cfg%jmax_
                    do i=fs%cfg%imin_,fs%cfg%imax_
                       fs%U(i,j,k)=fs%U(i,j,k)*sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                       fs%V(i,j,k)=fs%V(i,j,k)*sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
                       fs%W(i,j,k)=fs%W(i,j,k)*sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
                    end do
                 end do
              end do
              call fs%cfg%sync(fs%U)
              call fs%cfg%sync(fs%V)
              call fs%cfg%sync(fs%W)
            end block ibforcing

            ! Apply direct forcing on propellor
            ibm_correction: block
              integer :: i,j,k
              resU=fs%rho
              call df%get_source(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU)
              do k=fs%cfg%kmin_,fs%cfg%kmax_
                 do j=fs%cfg%jmin_,fs%cfg%jmax_
                    do i=fs%cfg%imin_,fs%cfg%imax_
                       fs%U(i,j,k)=fs%U(i,j,k)+sum(fs%itpr_x(:,i,j,k)*df%srcU(i-1:i,j,k))
                       fs%V(i,j,k)=fs%V(i,j,k)+sum(fs%itpr_y(:,i,j,k)*df%srcV(i,j-1:j,k))
                       fs%W(i,j,k)=fs%W(i,j,k)+sum(fs%itpr_z(:,i,j,k)*df%srcW(i,j,k-1:k))
                    end do
                 end do
              end do
              call fs%cfg%sync(fs%U)
              call fs%cfg%sync(fs%V)
              call fs%cfg%sync(fs%W)
            end block ibm_correction

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
         if (ens_evt%occurs()) then
            update_pmesh: block
              integer :: i
              call df%update_partmesh(pmesh)
              do i=1,df%np_
                 pmesh%var(1,i)=df%p(i)%dA
                 pmesh%vec(:,1,i)=df%p(i)%vel
              end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call df%get_max()
         call mfile%write()
         call cflfile%write()
         call ibmfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,gradU)
      
   end subroutine simulation_final
   
end module simulation
