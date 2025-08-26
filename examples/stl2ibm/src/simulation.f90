!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use incomp_class,      only: incomp
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Ensight postprocessing
   type(ensight)  :: ens_out
   
   public :: simulation_init,simulation_run,simulation_final
   
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      implicit none
      
      
      ! Create Ensight output
      create_ensight: block
         
         ! Create ensight output for ply surface
         ens_out=ensight(cfg=cfg,name='stl2ibm')
         
         ! Add variables to output
         call ens_out%add_scalar('Gib',cfg%Gib)
         call ens_out%write_data(0.0_WP)

         
      end block create_ensight
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
   end subroutine simulation_final
   
end module simulation
