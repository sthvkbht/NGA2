!> AMR I/O class: Checkpoint and restart using AMReX VisMF
!> Registration-based: solvers/users register fields, write all at once
!> Read piecemeal: read header first, then individual fields
module amrio_class
   use precision,      only: WP
   use string,         only: str_medium
   use iso_c_binding,  only: c_ptr,c_char,c_int,c_null_char
   use amrgrid_class,  only: amrgrid
   use amrdata_class,  only: amrdata
   use amrex_interface, only: amrmfab_vismf_write,amrmfab_vismf_read, &
      amrcheckpoint_prebuild_dirs,amrcheckpoint_mfab_prefix
   implicit none
   private

   public :: amrio

   !> Linked list node for registered data
   type :: data_node
      type(amrdata), pointer :: ptr => null()     !< Pointer to amrdata
      character(len=str_medium) :: name           !< Name for this data
      type(data_node), pointer :: next => null()  !< Next in list
   end type data_node

   !> Linked list node for scalar metadata
   type :: scalar_node
      character(len=str_medium) :: name           !< Name for this scalar
      real(WP) :: value                           !< Scalar value
      type(scalar_node), pointer :: next => null() !< Next in list
   end type scalar_node

   !> AMR I/O handler for checkpoints
   type :: amrio
      class(amrgrid), pointer :: amr => null()    !< Pointer to AMR grid
      type(data_node), pointer :: first => null() !< First registered data
      integer :: ndata = 0                        !< Number of registered data
      type(scalar_node), pointer :: first_scalar => null() !< First registered scalar
      integer :: nscalar = 0                      !< Number of registered scalars
   contains
      procedure :: initialize                     !< Initialize with grid and I/O aggregation
      procedure :: add_data                       !< Register a data field
      procedure :: add_scalar                     !< Register a scalar value
      procedure :: write                          !< Write all registered to checkpoint
      procedure :: read_header                    !< Read checkpoint header (time, step, fields)
      procedure :: read_data                      !< Read a specific field from checkpoint
      procedure :: get_scalar                     !< Get a scalar value from checkpoint
      procedure :: finalize                       !< Clean up registered data list
   end type amrio

contains

   !> Initialize I/O handler with AMR grid
   subroutine initialize(this, amr, nfiles)
      use amrex_interface, only: amrvismf_set_noutfiles
      implicit none
      class(amrio), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      integer, intent(in) :: nfiles  !< Number of output files (1=aggregated, nranks=file-per-rank)
      this%amr => amr
      this%ndata = 0
      nullify(this%first)
      call amrvismf_set_noutfiles(nfiles)
   end subroutine initialize


   !> Register an amrdata field for checkpointing
   subroutine add_data(this, data, name)
      implicit none
      class(amrio), intent(inout) :: this
      type(amrdata), target, intent(in) :: data
      character(len=*), intent(in) :: name
      type(data_node), pointer :: new_node, current

      ! Create new node
      allocate(new_node)
      new_node%ptr => data
      new_node%name = trim(name)
      nullify(new_node%next)

      ! Add to end of list
      if (.not.associated(this%first)) then
         this%first => new_node
      else
         current => this%first
         do while (associated(current%next))
            current => current%next
         end do
         current%next => new_node
      end if
      this%ndata = this%ndata + 1
   end subroutine add_data


   !> Register a scalar value for checkpointing (dt, ambient pressure, etc.)
   subroutine add_scalar(this, name, value)
      implicit none
      class(amrio), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), intent(in) :: value
      type(scalar_node), pointer :: new_node, current

      ! Create new node
      allocate(new_node)
      new_node%name = trim(name)
      new_node%value = value
      nullify(new_node%next)

      ! Add to end of list
      if (.not.associated(this%first_scalar)) then
         this%first_scalar => new_node
      else
         current => this%first_scalar
         do while (associated(current%next))
            current => current%next
         end do
         current%next => new_node
      end if
      this%nscalar = this%nscalar + 1
   end subroutine add_scalar


   !> Write all registered data to checkpoint directory
   subroutine write(this, dirname, time, step)
      use messager, only: log
      use parallel, only: MPI_REAL_WP
      use mpi_f08,  only: MPI_BCAST,MPI_INTEGER
      implicit none
      class(amrio), intent(inout) :: this
      character(len=*), intent(in) :: dirname      !< Checkpoint directory name
      real(WP), intent(in) :: time                 !< Simulation time
      integer, intent(in) :: step                  !< Step number

      integer :: lev, ierr, iunit
      character(len=256) :: mfab_path
      character(len=str_medium) :: header_file
      type(data_node), pointer :: current
      type(scalar_node), pointer :: scurrent

      ! Create directory hierarchy
      call amrcheckpoint_prebuild_dirs(trim(dirname)//c_null_char, &
         'Level_'//c_null_char, this%amr%nlevels)

      ! Write Header file (root only)
      if (this%amr%amRoot) then
         header_file = trim(dirname)//'/Header'
         open(newunit=iunit, file=trim(header_file), status='replace', action='write')
         write(iunit,'(A)') 'NGA2 AMR Checkpoint v3'  ! v3 includes scalar metadata
         write(iunit,'(I0)') this%amr%nlevels - 1  ! finest_level (0-indexed)
         write(iunit,'(I0)') step
         write(iunit,'(ES23.16)') time
         ! Write number of registered fields and their names
         write(iunit,'(I0)') this%ndata
         current => this%first
         do while (associated(current))
            write(iunit,'(A)') trim(current%name)
            current => current%next
         end do
         ! Write number of registered scalars and their name/value pairs
         write(iunit,'(I0)') this%nscalar
         scurrent => this%first_scalar
         do while (associated(scurrent))
            write(iunit,'(A,1X,ES23.16)') trim(scurrent%name), scurrent%value
            scurrent => scurrent%next
         end do
         ! Write refinement ratio per level
         do lev = 0, this%amr%maxlvl
            write(iunit,'(I0,A,I0)') lev, ' ', this%amr%rref(min(lev, this%amr%maxlvl-1))
         end do
         close(iunit)
         call log('Wrote checkpoint Header: '//trim(header_file))
      end if

      ! Write MultiFab data for each registered field
      current => this%first
      do while (associated(current))
         do lev = 0, this%amr%maxlvl
            call amrcheckpoint_mfab_prefix(mfab_path, len(mfab_path), lev, &
               trim(dirname)//c_null_char, 'Level_'//c_null_char, &
               trim(current%name)//c_null_char)
            call amrmfab_vismf_write(current%ptr%mf(lev)%p, trim(mfab_path)//c_null_char)
         end do
         current => current%next
      end do

      if (this%amr%amRoot) call log('Wrote checkpoint: '//trim(dirname)//' ('// &
         trim(adjustl(itoa(this%ndata)))//' fields)')
   contains
      function itoa(i) result(str)
         integer, intent(in) :: i
         character(len=12) :: str
         write(str,'(I0)') i
      end function itoa
   end subroutine write


   !> Read checkpoint header to get time, step, and available field names
   !> Also reads scalar metadata into internal list for get_scalar access
   subroutine read_header(this, dirname, time, step, nfields, fieldnames)
      use messager, only: log,die
      use parallel, only: MPI_REAL_WP
      use mpi_f08,  only: MPI_BCAST,MPI_INTEGER,MPI_CHARACTER
      implicit none
      class(amrio), intent(inout) :: this
      character(len=*), intent(in) :: dirname                  !< Checkpoint directory
      real(WP), intent(out) :: time                            !< Simulation time
      integer, intent(out) :: step                             !< Step number
      integer, intent(out), optional :: nfields                !< Number of fields
      character(len=str_medium), allocatable, intent(out), optional :: fieldnames(:)

      integer :: ierr, iunit, finest_level, nf, ns, i
      character(len=256) :: line
      character(len=str_medium) :: header_file, sname
      character(len=str_medium), allocatable :: names(:)
      real(WP) :: svalue
      real(WP), allocatable :: svalues(:)
      character(len=str_medium), allocatable :: snames(:)
      type(scalar_node), pointer :: new_node, current

      header_file = trim(dirname)//'/Header'

      ! Clear any existing scalar list
      current => this%first_scalar
      do while (associated(current))
         new_node => current%next
         deallocate(current)
         current => new_node
      end do
      nullify(this%first_scalar)
      this%nscalar = 0

      ! Root reads header
      ns = 0
      if (this%amr%amRoot) then
         open(newunit=iunit, file=trim(header_file), status='old', action='read', iostat=ierr)
         if (ierr .ne. 0) call die('Cannot open checkpoint Header: '//trim(header_file))
         read(iunit,'(A)') line  ! Title (v2 or v3)
         read(iunit,*) finest_level
         read(iunit,*) step
         read(iunit,*) time
         read(iunit,*) nf
         allocate(names(nf))
         do i = 1, nf
            read(iunit,'(A)') names(i)
         end do
         ! Check if v3 format with scalars
         if (index(line, 'v3') .gt. 0) then
            read(iunit,*) ns
            allocate(snames(ns), svalues(ns))
            do i = 1, ns
               read(iunit,*) snames(i), svalues(i)
            end do
         end if
         close(iunit)
         call log('Read checkpoint Header: '//trim(header_file))
      end if

      ! Broadcast to all ranks
      call MPI_BCAST(step, 1, MPI_INTEGER, 0, this%amr%comm, ierr)
      call MPI_BCAST(time, 1, MPI_REAL_WP, 0, this%amr%comm, ierr)
      call MPI_BCAST(nf, 1, MPI_INTEGER, 0, this%amr%comm, ierr)
      call MPI_BCAST(ns, 1, MPI_INTEGER, 0, this%amr%comm, ierr)

      if (.not.this%amr%amRoot) allocate(names(nf))
      do i = 1, nf
         call MPI_BCAST(names(i), str_medium, MPI_CHARACTER, 0, this%amr%comm, ierr)
      end do

      ! Broadcast and store scalars
      if (ns .gt. 0) then
         if (.not.this%amr%amRoot) allocate(snames(ns), svalues(ns))
         do i = 1, ns
            call MPI_BCAST(snames(i), str_medium, MPI_CHARACTER, 0, this%amr%comm, ierr)
            call MPI_BCAST(svalues(i), 1, MPI_REAL_WP, 0, this%amr%comm, ierr)
         end do
         ! Store in internal list
         do i = 1, ns
            call this%add_scalar(trim(snames(i)), svalues(i))
         end do
         deallocate(snames, svalues)
      end if

      ! Return optional outputs
      if (present(nfields)) nfields = nf
      if (present(fieldnames)) then
         allocate(fieldnames(nf))
         fieldnames = names
      end if
      deallocate(names)
   end subroutine read_header


   !> Read a specific field from checkpoint directory
   subroutine read_data(this, dirname, data, dataname)
      use messager, only: log
      implicit none
      class(amrio), intent(inout) :: this
      character(len=*), intent(in) :: dirname      !< Checkpoint directory
      type(amrdata), intent(inout) :: data         !< Data to read into
      character(len=*), intent(in) :: dataname     !< Name of field to read

      integer :: lev
      character(len=256) :: mfab_path

      ! Read MultiFab data for each level
      ! NOTE: The amrdata must already be allocated with correct BoxArray/DistroMap
      do lev = 0, this%amr%maxlvl
         call amrcheckpoint_mfab_prefix(mfab_path, len(mfab_path), lev, &
            trim(dirname)//c_null_char, 'Level_'//c_null_char, &
            trim(dataname)//c_null_char)
         call amrmfab_vismf_read(data%mf(lev)%p, trim(mfab_path)//c_null_char)
      end do

      if (this%amr%amRoot) call log('Read checkpoint field: '//trim(dataname))
   end subroutine read_data


   !> Get a scalar value by name (after read_header has been called)
   subroutine get_scalar(this, name, value, found)
      use messager, only: warn
      implicit none
      class(amrio), intent(in) :: this
      character(len=*), intent(in) :: name
      real(WP), intent(out) :: value
      logical, intent(out), optional :: found
      type(scalar_node), pointer :: current

      current => this%first_scalar
      do while (associated(current))
         if (trim(current%name) .eq. trim(name)) then
            value = current%value
            if (present(found)) found = .true.
            return
         end if
         current => current%next
      end do

      ! Not found
      value = 0.0_WP
      if (present(found)) then
         found = .false.
      else
         call warn('Scalar not found in checkpoint: '//trim(name))
      end if
   end subroutine get_scalar


   !> Finalize: clean up registered data and scalar lists
   subroutine finalize(this)
      implicit none
      class(amrio), intent(inout) :: this
      type(data_node), pointer :: current, next
      type(scalar_node), pointer :: scurrent, snext

      ! Clean up data list
      current => this%first
      do while (associated(current))
         next => current%next
         deallocate(current)
         current => next
      end do
      nullify(this%first)
      this%ndata = 0

      ! Clean up scalar list
      scurrent => this%first_scalar
      do while (associated(scurrent))
         snext => scurrent%next
         deallocate(scurrent)
         scurrent => snext
      end do
      nullify(this%first_scalar)
      this%nscalar = 0

      nullify(this%amr)
   end subroutine finalize

end module amrio_class
