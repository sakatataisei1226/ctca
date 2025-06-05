! program myesses_example
!     use m_myesses, only: t_ESSimulator, new_ESSimulator
!     implicit none

!     type(t_ESSimulator) :: essimulator

!     integer :: ierr
!     integer :: nnode, myid

!     ! Initialize MPI
!     call MPI_Init(ierr)
!     call MPI_Comm_size(MCW, nnode, ierr)
!     call MPI_Comm_rank(MCW, myid, ierr)

!     ! Initialize HDF5
!     call h5fortran_init

!     essimulator = new_ESSimulator(nnode, myid)

!     call essimulator%initialize()

!     call essimulator%run(snapshot=snap)

!     call essimulator%finalize()

!     ! Finalize HDF5
!     call h5fortran_finalize

!     ! Finalize MPI
!     call hdffinalize()
!     call MPI_Finalize(mpierr)
! end program

! module m_myesses
!     use m_h5fortran, only: t_h5file, h5fortran_init, h5fortran_finalize
! #define MCW MPI_COMM_WORLD
!     implicit none

!     type t_ESSimulator
!         integer :: istep
!         double precision :: t
!         integer :: nnode, myid
!     contains
!         procedure :: initialize => essimulator_initialize
!         procedure :: run => essimulator_run
!         procedure :: update => essimulator_update
!         procedure :: show_status => essimulator_show_status
!         procedure :: finalize => essimulator_finalize
!     end type

!     private
!     public t_ESSimulator
!     public new_ESSimulator
! contains

!     function new_ESSimulator(nnode, myid) result(obj)
!         integer, intent(in) :: nnode
!         integer, intent(in) :: myid
!         type(t_ESSimulator) :: obj

!         obj%nnode = nnode
!         obj%myid = myid
!     end function

!     subroutine essimulator_run(self)
!         type(t_ESSimulator), intent(inout) :: self
!         integer :: ierr

!         ! Run Main Loop
!         do self%istep = 1, nstep
!             self%t = self%t + dt

!             ! Show Simulation Status
!             call self%show_status()

!             ! Update Simulation Step (TODO)
!             call self%update()
!         end do
!     end subroutine

!     subroutine essimulator_initialize(self)
!         type(t_ESSimulator), intent(inout) :: self

!         ! Input Parameters (TODO)

!         ! Initialize Simulation (TODO)
!     end subroutine

!     subroutine essimulator_update(self)
!         type(t_ESSimulator), intent(inout) :: self
!     end subroutine

!     subroutine show_simulation_status(self)
!         type(t_ESSimulator), intent(inout) :: self

!         if (self%myid == 0) then
!             print *, '**** step ---------', self%istep
!         end if
!     end subroutine

!     subroutine essimulator_finalize(self)
!         type(t_ESSimulator), intent(inout) :: self

!     end subroutine

! end module
