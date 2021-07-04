! 2014-01-27: Added calculation of 2nd order derivative of a before the
!             first step in the constrained case
! 2014-01-31: Removed above calculation (did not give better results)
!             by commenting out
! 2014-02-09: Added debug subroutines print_matrix and print_vector
! 2014-02-28: Added gena_stats in order to calculate the average and
!             maximal newton steps
! 2014-03-04: Added opts%const_mass_matrix in order to explicitly use
!             a constant mass matrix
! 2014-03-04: Added opts%diag_mass_matrix in order to explicitly use
!             a diagonal mass matrix
! 2014-03-04: Added opts%banded_iteration_matrix in order to exploit
!             a possible banded structure of the iteration matrix St.
!             (only used for unconstrained systems)
! 2014-03-04: Modificated the Newton method to only evaluate the
!             iteration matrix St once and use it until convergence
!             TODO: * condition for recalculating St
!
! 2014-03-05: Added switch to recalculate the iteration matrix St in
!             every Newton step
! 2014-03-07: Fixed bug where gena_const_M or gena_const_diag_M would be
!             allocated more than once
!             Added subroutine gena_print_stats
!             Added stats%time, integration time measurement and
!             adjusted gena_print_stats
! 2014-03-15: Added subroutine gena_cleanup
! 2014-04-01: Added a standard call of the output function directly
!             after its initialization call. Thus, the initialization
!             is decoupled from outputting the initial conditions.
! 2015-02-05: Added support for 2014-03-04 in the constrained case
! 2015-03-12: Added pertubing initial values in order to overcome the
!             order reduction in the first steps for unconstrained case
! 2015-03-12: Changed the initialisation of vd1=0.0 to vd1=vd in
!             gena_solveTimeStep
! 2015-03-18: Added pertubing initial values in order to overcome the
!             order reduction in the first steps for constrained case
!             TODO: Activate pertubing should be standard
!
! 2015-04-14: Added subroutine print_vector_int
!             Added support for banded iteration matrix in the
!             constrained case using a vector jour for ordering the
!             iteration matrix. The vector jour must be user-provided.
! 2015-04-21: Implemented curvature term gena_Z
! 2015-04-23: Added Lagrange multipliers to the gena_problem structure
!             Renamed this%ng to this%sizel
!             Changed the interface to use the this%sizeX variables
!             instead of size(this%X)
!             In the time step subroutines, changed some occurences of
!             size(X) to sX (which was introduced by an associate block)
!             Also removed the internal variables ng and sv and
!             replaced them with said associate variables
!             Added calculation of this%vd to gena_pertube
!             Added calculation of this%vd to gena_pertubeConstrained
!                 TODO: Debug pertube
! 2015-04-28: Moved the first output subroutine call after the
!             pertubation of initial values
! 2015-05-03: Reverted  2015-03-12, because for small step sizes the
!             solution to testprob5 would converge to a different
!             solution. WE NEED TO KNOW WHY! TODO
!             Same applies to the constrained case (testprob4, 4a)
! 2015-05-25: Corrected gena_num_Kt
!             renamed gena_pertube and gena_pertubeConstrained to
!             gena_calcInitial and gena_calcInitialConstrained and
!             these subroutines are called unconditionally in
!             gena_integrate
!             moved the if (this%opts%pertube == 1) into the latter
! 2015-05-26: Implemented exploiting band structure in gena_num_Ct,
!             gena_num_Kt and gena_num_Kt_lambda
! 2015-05-28: Fixed exploiting band structure in the funcitons above
!             Added counting of calls of gena_g and gena_B
!             Added experimental support for the stabilized index-2
!             system. Hashtag: STAB2
! 2015-06-17: Changed the calculation of the initial values in the
!             stabilized index-2 system from the choice
!                 \dot v_{n+1}^{(0)}=0
!             to
!                 a_{n+1}^{(0)} = a_n.
!             The update for a_{n+1} after convergence of the newton
!             iteration was changed accordingly.
! 2015-06-23: Applied the changed from 2015-06-17 to the ODE case and
!             the normal index-3 case.
!             Changed the criterion to stop the newton iteration from
!             testing, if the residuals are small to if the increments
!             are small.
!             The first is used in BrulsCardonaArnold12 and the latter
!             in CISM15.pdf ! TODO
!              TODO: Testen
! 2015-07-08: Changed gena_num_Ct, gena_num_Kt and gena_num_Kt_lambda.
!             The step size for the finite differences was changed from
!                 h = sqrt(this%opts%atol)
!             to
!                 h = max( abs( x(i) )*1.0e-8_8,  1.0e-12_8)
!             with x=v or x=q resp.
! 2015-07-21: Changed calculation of the pertubed initial values
! 2015-07-22: Added macro STNUM in order to calculate the iteration
!             matrix in the stabilized index-2 formulation via finite
!             differences. It is only meant for debugging. Note that
!             gena_M must be correct.
! 2015-07-22: Fixed error in the calculation of the initial value for
!             Dq in stabilized index-2 formulation
!             (added - B^T(q)*eta)
! 2015-09-24: Fixed error in the convergence test for the stab. index-2
!             formulation. (dividing by (sv+2*sl) instead of (sv+sl))
! 2015-09-28: Added support for reordering iteration matrix in the
!             stabilized index-2 case
! 2015-09-28: Removed hashtag STAB2. Stabilized index-2 is no longer
!             experimental. The parameter this%opts%stab2 takes the
!             place of the macro STAB2.
! 2015-09-28: Removed hashtags NOCT and NOKT. The parameters
!             this%opts%no_Ct and this%opts%no_Kt take their places.
!
! TODO: Anpassen der testprobs auf die Änderungen von 2013-03-04
! TODO: Anpassen der testprobs auf die Änderungen von 2015-04-21


module gena
   implicit none

   ! definition of integrator options
   type  :: gena_options
      ! system is constrained?
      integer :: constrained = 0
      ! use stabilized index two formulation
      integer :: stab2 = 0
      ! mass matrix is constant?
      integer :: const_mass_matrix = 0
      ! mass matrix is a diagonal matrix? If == 1 then gena_diag_M
      ! is used in order to calculate the mass matrix rather than
      ! gena_M
      integer :: diag_mass_matrix = 0
      ! iteration matrix is banded? (in contrained case: an ordered iteration matrix)
      integer :: banded_iteration_matrix = 0
      ! if iteration matrix is banded: number of subdiagonals
      integer :: nr_subdiag = 0
      ! if iteration matrix is banded: number of superdiagonals
      integer :: nr_superdiag = 0
      ! vector, that is used to order the iteration matrix in order to obtain a banded structure in the constrained case
      integer, dimension(:), allocatable  :: jour
      ! recalculate the iteration matrix in every Newton step
      integer :: recalc_iteration_matrix = 0
      ! pertube initial values?
      integer :: pertube = 0 ! DEBUG
      ! parameter s for pertubing initial values
      real(8) :: pertube_s = 1.0_8 ! TODO: Better value?
      ! use numerical approximation for Ct and Kt
      integer :: use_num_Ct = 1
      integer :: use_num_Kt = 1
      ! omit Ct and Kt resp. in the iteration matrix
      integer :: no_Ct = 0
      integer :: no_Kt = 0
      ! variables for error control of newton method (absolute and relativ tolerance)
      real(8)  :: atol = 1.0e-10_8
      real(8)  :: rtol = 1.0e-8_8
!      real(8) :: tolr   = 1.0e-5_8
!      real(8) :: tolphi = 1.0e-8_8 ! irrelevant for this%opts%constrained == 0
!#ifdef STAB2
!      real(8) :: tolB   = 1.0e-5_8
!#endif
      integer :: imax   = 5
      ! integration span
      real(8) :: t0 = 0.0_8
      real(8) :: te = 1.0_8
      integer :: nsteps = 100
   end type gena_options

   ! definition of integrator statistics
   type  :: gena_statistics
      ! current number of newton steps TODO: private
      integer  :: newt_steps_curr = 0
      ! number of newton steps
      integer  :: newt_steps_sum = 0
      ! maximum number of newton steps
      integer  :: newt_steps_max = 0
      ! average number of newton steps
      real(8)  :: newt_steps_avg = 0.0_8
      ! number of calls
      integer  :: ngcalls = 0
      integer  :: nBcalls = 0
      ! integration time
      real(8)  :: time = 0.0_8
   end type gena_statistics

   ! definition of abstract problem type
   type, abstract :: gena_problem
      ! variables that define the state of the integrator
      real(8)                             :: t = 0.0_8
      integer                             :: sizeq = 0
      real(8), dimension(:), allocatable  :: q
      integer                             :: sizev = 0
      real(8), dimension(:), allocatable  :: v
      real(8), dimension(:), allocatable  :: vd
      real(8), dimension(:), allocatable  :: a  ! TODO: private
      integer                             :: sizel = 0 ! number of constraints, irrelevant for this%opts%constrained == 0
      real(8), dimension(:), allocatable  :: l ! Lagrange multipliers, only needed in the constrained case
      real(8), dimension(:), allocatable  :: eta ! Auxiliar variables eta, only needed in the stabilized index-2 case
      ! integrator options
      type(gena_options)      :: opts  ! TODO: rename to gena_opts
      ! solver statistics
      type(gena_statistics)   :: gena_stats
      ! constant mass matrix (if opts%const_mass_matrix == 1)
      real(8), dimension(:,:), allocatable   :: gena_const_M
      ! constant diagonal mass matrix (if opts%diag_mass_matrix == 1, also)
      real(8), dimension(:), allocatable     :: gena_const_diag_M
      ! internal variables, momentan für numerische dissipation = 9/10 TODO: private
      real(8) :: gena_alpha_m =  8.0_8/ 19.0_8
      real(8) :: gena_alpha_f =  9.0_8/ 19.0_8
      real(8) :: gena_beta    =100.0_8/361.0_8
      real(8) :: gena_gamma   = 21.0_8/ 38.0_8
      !! numerische dissipation = 0.5
      !real(8) :: gena_alpha_m =  0.0_8/19.0_8
      !real(8) :: gena_alpha_f =  1.0_8/ 3.0_8
      !real(8) :: gena_beta    =  4.0_8/ 9.0_8
      !real(8) :: gena_gamma   =  5.0_8/ 6.0_8
   ! definition of deferred procedures
   contains
         ! function $M$ in $M(q) \dot v = -g(q,v,t)$
      procedure(gena_M),              deferred :: gena_M
         ! function $M$ in $M(q) \dot v = -g(q,v,t)$, diagonal case
      procedure(gena_diag_M),         deferred :: gena_diag_M
         ! function $g$ in $M(q) \dot v = -g(q,v,t)$
      procedure(gena_g),              deferred :: gena_g
         ! operation in the lie space $q_n * \exp(h\cdot\widetilde{\Delta q_n})$
      procedure(gena_qlpexphDqtilde), deferred :: gena_qlpexphDqtilde
         ! operation in the lie algebra $inversetilde([\tilde{v},\tilde{w}])$, may be a dummy, if system is unconstrained
      procedure(gena_itlbtvtw),         deferred :: gena_itlbtvtw
#ifdef extra_integrate
         ! special operation: exp(4/3*log(x1*x2^(-1)))*x2
      procedure(gena_exp43logx1ix2x2), deferred :: gena_exp43logx1ix2x2
#endif
         ! tilde operator in the lie group (may be a dummy function, if gena_num_Kt is not used)  ! TODO: This is not used anymore
      procedure(gena_tilde),          deferred :: gena_tilde ! TODO
         ! tangent damping matrix $C_t$
      procedure(gena_Ct),             deferred :: gena_Ct
         ! tangent stiffness matrix $K_t$
      procedure(gena_Kt),             deferred :: gena_Kt
         ! tangent stiffness matrix $K_t$ in the constrained case
         ! (may depend on the Lagrange multipliers)
      procedure(gena_Kt_lambda),      deferred :: gena_Kt_lambda
         ! tangent operator of the exponential map
      procedure(gena_Tg),             deferred :: gena_Tg
         ! norm    ! TODO: this is actually not used anymore
      procedure(gena_norm),           deferred :: gena_norm   ! TODO
         ! subroutine for output, is called after every integration step
      procedure(gena_outputFunction), deferred :: gena_outputFunction
         ! subroutine to initialize the problem
      procedure(gena_init),           deferred :: gena_init
         ! function $\Phi(q)$
      procedure(gena_phi),            deferred :: gena_phi
         ! function $B(q)$, where $B(q)v = d/dq \Phi(q) * (DL_q(e) v)$
      procedure(gena_b),              deferred :: gena_b
         ! function $F(q,v) = d/(dq) (B(q)v) * (DL_q(e) v)$
      procedure(gena_Z),              deferred :: gena_Z
         ! function Z, where $Zw = Z(q(Dq)) (v(Dq + BT(q)eta), T w)$
      procedure(gena_matZ),           deferred :: gena_matZ
         ! subroutine to calculate correct initial values for vd and a
      procedure                  :: gena_calcInitial => gena_calcInitial ! TODO: private
         ! subroutine to calculate correct initial values for vd and a, also apply pertubing, if needed
      procedure                  :: gena_calcInitialConstrained => gena_calcInitialConstrained
         ! subroutine to integrate one time step
      procedure                  :: gena_solveTimeStep => gena_solveTimeStep ! TODO: private
         ! subroutine to integrate one time step in a constrained system
      procedure                  :: gena_solveConstrainedTimeStep => gena_solveConstrainedTimeStep
         ! subroutine to integrate one time step in a constrained system
      procedure                  :: gena_solveConstrainedTimeStep_stab2 => gena_solveConstrainedTimeStep_stab2
         ! subroutine for integration
      procedure                  :: gena_integrate => gena_integrate
#ifdef extra_integrate
      procedure                  :: gena_extra_integrate => gena_extra_integrate
#endif
         ! subroutine for numerical approximation of Ct
      procedure                  :: gena_num_Ct => gena_num_Ct
         ! subroutine for numerical approximation of Kt
      procedure                  :: gena_num_Kt => gena_num_Kt
         ! subroutine for numerical approximation of Kt in the constrained case
      procedure                  :: gena_num_Kt_lambda => gena_num_Kt_lambda
         ! In order to use the numerical approximations of Ct and Kt
         ! respectively, the procedures gena_Ct and gena_Kt as well as
         ! their constrained pendants gen_CT_lambda and gena_Kt_lambde
         ! have to be referenced to gena_num_Ct and gena_num_Kt resp.
      procedure                  :: gena_print_stats => gena_print_stats
         ! Clean up the gena_problem object (only internal variables
         ! are reset, gena_alpha_m, gena_alpha_f, gena_beta and
         ! gena_gamma are NOT reset)
      procedure                  :: gena_cleanup => gena_cleanup
#ifdef STNUM
      procedure                  :: St_Phi => St_Phi
      procedure                  :: St_num => St_num
#endif
   end type gena_problem

   ! abstract interface of the problem
   abstract interface
      pure function gena_M(this, q) result(rslt)
         import                                          :: gena_problem
         ! input
         class(gena_problem),   intent(in)               :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev,this%sizev)       :: rslt
      end function gena_M

      pure function gena_diag_M(this, q) result(rslt)
         import                                          :: gena_problem
         ! input
         class(gena_problem),   intent(in)               :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev)                  :: rslt
      end function gena_diag_M

      pure function gena_g(this, q, v, t) result(rslt)
         import                              :: gena_problem
         ! input
         class(gena_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8), dimension(:), intent(in)   :: v
         real(8),               intent(in)   :: t
         ! result
         real(8), dimension(this%sizev)      :: rslt
      end function gena_g

      pure function gena_qlpexphDqtilde(this, q, h, dq) result(rslt)
         import                              :: gena_problem
         ! input
         class(gena_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8),               intent(in)   :: h
         real(8), dimension(:), intent(in)   :: dq
         ! result
         real(8), dimension(this%sizeq)      :: rslt
      end function gena_qlpexphDqtilde

      pure function gena_itlbtvtw(this, v, w) result(rslt)
         import                              :: gena_problem
         ! input
         class(gena_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         real(8), dimension(:), intent(in)   :: w
         ! result
         real(8), dimension(this%sizev)      :: rslt
      end function gena_itlbtvtw

#ifdef extra_integrate
      pure function gena_exp43logx1ix2x2(this, x1, x2) result(rslt)
         import      :: gena_problem
         ! input
         class(gena_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: x1
         real(8), dimension(:), intent(in)   :: x2
         ! result
         real(8), dimension(this%sizeq)      :: rslt
      end function gena_exp43logx1ix2x2
#endif

      pure function gena_tilde(this, v) result(rslt)
         import                              :: gena_problem
         ! input
         class(gena_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         ! result
         real(8), dimension(this%sizeq)      :: rslt
      end function gena_tilde

      pure function gena_Ct(this, q, v, t) result(rslt)
         import                                    :: gena_problem
         ! input
         class(gena_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function gena_Ct

      pure function gena_Kt(this, q, v, vd, t) result(rslt)
         import                                    :: gena_problem
         ! input
         class(gena_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8), dimension(:), intent(in)         :: vd
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function gena_Kt

      pure function gena_Kt_lambda(this, q, v, vd, l, t) result(rslt)
         import                                    :: gena_problem
         ! input
         class(gena_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8), dimension(:), intent(in)         :: vd
         real(8), dimension(:), intent(in)         :: l
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function gena_Kt_lambda

      pure function gena_Tg(this, h, dq) result(rslt)
         import                                    :: gena_problem
         ! input
         class(gena_problem),   intent(in)         :: this
         real(8),               intent(in)         :: h
         real(8), dimension(:), intent(in)         :: dq
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function gena_Tg

      pure function gena_norm(this, v) result(rslt)
         import                  :: gena_problem
         ! input
         class(gena_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         ! result
         real(8)                             :: rslt
      end function gena_norm

      subroutine gena_outputFunction(this,info)
         import                           :: gena_problem
         ! input
         class(gena_problem), intent(in)  :: this
         integer,             intent(in)  :: info
      end subroutine gena_outputFunction

      subroutine gena_init(this)
         import                              :: gena_problem
         ! input/output
         class(gena_problem), intent(inout)  :: this
      end subroutine gena_init

      pure function gena_phi(this,q) result(rslt)
         import                                       :: gena_problem
         ! input
         class(gena_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         ! result
         real(8), dimension(this%sizel)               :: rslt
      end function gena_phi

      pure function gena_B(this,q) result(rslt)
         import                                       :: gena_problem
         ! input
         class(gena_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         ! result
         real(8), dimension(this%sizel,this%sizev)    :: rslt
      end function gena_B

      pure function gena_Z(this,q,v) result(rslt)
         import                                       :: gena_problem
         ! input
         class(gena_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         real(8), dimension(:),            intent(in) :: v
         ! result
         real(8), dimension(this%sizel)               :: rslt
      end function gena_Z

      pure function gena_matZ(this,q,v,T) result(rslt)
         import                                       :: gena_problem
         ! input
         class(gena_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         real(8), dimension(:),            intent(in) :: v
         real(8), dimension(:,:),          intent(in) :: T
         ! result
         real(8), dimension(this%sizel, this%sizev)   :: rslt
      end function gena_matZ

   end interface

   contains

#ifdef STNUM

   pure function St_num(this,dq,h,l,eta,t) result(rslt)
      implicit none
      ! input
      class(gena_problem),    intent(in)  :: this
      real(8), dimension(:),  intent(in)  :: dq
      real(8),                intent(in)  :: h
      real(8), dimension(:),  intent(in)  :: l
      real(8), dimension(:),  intent(in)  :: eta
      real(8),                intent(in)  :: t
      ! result
      real(8), dimension(this%sizev+2*this%sizel,this%sizev+2*this%sizel)  :: rslt
      ! internal
      real(8), dimension(this%sizev)      :: wdq
      real(8), dimension(this%sizel)      :: wl
      integer                             :: i
      real(8), dimension(this%sizev+2*this%sizel)  :: St_Phi0
      !
      St_Phi0 = this%St_Phi(dq, h, h*l, eta, t)
      !
      do i = 1,this%sizev
         wdq    = 0.0_8
         wdq(i) = max( abs(dq(i))*1.0e-8_8,  1.0e-12_8)
         rslt(:,i) = (this%St_Phi(dq+wdq,h,h*l, eta, t) - St_Phi0)/ wdq(i)
      end do
      !
      do i = 1,this%sizel
         wl    = 0.0_8
         wl(i) = max( abs(h*l(i))*1.0e-8_8,  1.0e-12_8)
         rslt(:,this%sizev+i) = (this%St_Phi(dq,h,h*l + wl, eta, t) - St_Phi0)/ wl(i)
      end do
      !
      do i = 1,this%sizel
         wl    = 0.0_8
         wl(i) = max( abs(eta(i))*1.0e-8_8,  1.0e-12_8)
         rslt(:,this%sizev+this%sizel+i) = (this%St_Phi(dq,h,h*l, eta + wl, t) - St_Phi0)/ wl(i)
      end do
   end function St_num

   pure function St_Phi(this,dq,h,hl,eta,t) result(rslt)
      implicit none
      ! input
      class(gena_problem),    intent(in)  :: this
      real(8), dimension(:),  intent(in)  :: dq
      real(8),                intent(in)  :: h
      real(8), dimension(:),  intent(in)  :: hl
      real(8), dimension(:),  intent(in)  :: eta
      real(8),                intent(in)  :: t
      ! result
      real(8), dimension(this%sizev+2*this%sizel)  :: rslt
      ! internal
      real(8), dimension(this%sizeq)            :: q
      real(8), dimension(this%sizev)            :: v
      real(8), dimension(this%sizev)            :: vd
      real(8), dimension(this%sizev)            :: BTeta
      real(8), dimension(this%sizel,this%sizev) :: B
      !
      associate (                         &
            alm => this%gena_alpha_m,     &
            alf => this%gena_alpha_f,     &
            bet => this%gena_beta,        &
            gam => this%gena_gamma,       &
            sv  => this%sizev,            &
            sl  => this%sizel             )
         BTeta = matmul(transpose(this%gena_B(this%q)), eta)

         q  = this%gena_qlpexphDqtilde(this%q, h, dq)
         v  = gam/bet*(dq + BTeta) + (1-gam/bet)*this%v + h*(1-gam/(2*bet))*this%a
         vd = (1-alm)/(bet*(1-alf)) * ( ((dq + BTeta) - this%v)/h - 0.5_8*this%a ) + (this%a - alf*this%vd)/(1-alf)

         B = this%gena_B(q)

         rslt(      1:      sv) = h*(matmul(this%gena_M(q),vd) + this%gena_g(q,v,t)) + matmul(transpose(B),hl)
         rslt(   sv+1:   sv+sl) = 1/h*this%gena_Phi(q)
         rslt(sv+sl+1:sv+sl+sl) = matmul(B, v)
      end associate
   end function St_Phi
#endif

   ! subroutine for pertubing initial values
   subroutine gena_calcInitial(this, h)
      implicit none
      class(gena_problem), intent(inout)  :: this  ! problem object
      real(8),             intent(in   )  :: h     ! step size
      ! internal variables
      real(8), dimension(this%sizeq, 2)            :: qpmsh  ! $q_{\pm sh}$
      real(8), dimension(this%sizev, 2)            :: vpmsh  ! $v_{\pm sh}$
      real(8), dimension(this%sizev, this%sizev,2) :: M      ! mass matrix
      real(8), dimension(this%sizev, 2)            :: vdpmsh ! $\dot v_{\pm sh}$
      integer, dimension(this%sizev)               :: ipiv   ! pivot vector for dgesv
      integer                                      :: info   ! info flag for dgesv
      !

      associate ( t  => this%t,            &
                  q  => this%q,            &
                  v  => this%v,            &
                  vd => this%vd,           &
                  a  => this%a,            &
                  sv => this%sizev,        &
                  s  => this%opts%pertube_s)
         ! calculate $\dot v(t_0)$
         if (this%opts%diag_mass_matrix == 1) then
            if (this%opts%const_mass_matrix == 1) then
               vd = -this%gena_g(q, v, t)/this%gena_const_diag_M
            else
               vd = -this%gena_g(q, v, t)/this%gena_diag_M(q)
            end if
            ! count calls
            this%gena_stats%ngcalls = this%gena_stats%ngcalls + 1
         else
            if (this%opts%const_mass_matrix == 1) then
               M(:,:,1) = this%gena_const_M
            else
               M(:,:,1) = this%gena_M(q)
            end if
            ! we need to solve a linear equation, first calulate the rhs
            vd = -this%gena_g(q, v, t)
            ! count calls
            this%gena_stats%ngcalls = this%gena_stats%ngcalls + 1
            ! then solve the system
            call dgesv(         &! solve the System A*X=B and save the result in B
                       sv,      &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                       1,       &! number of right hand sides (=size(B,2))
                       M(:,:,1),&! matrix A
                       sv,      &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                       ipiv,    &! integer pivot vector; it is not needed
                       vd,      &! matrix B
                       sv,      &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                       info)     ! integer information flag
            ! Now vd actually contains $\dot v(t_0)$
            if (info .ne. 0)  print*, "dgesv sagt info=", info ! TODO
         end if

         ! calculate $q_{\pm sh}$
         qpmsh(:,1) = this%gena_qlpexphDqtilde(q, h, &
                        -s * h * v   +   0.5_8 * s**2 * h**2 * vd)
         qpmsh(:,2) = this%gena_qlpexphDqtilde(q, h, &
                         s * h * v   +   0.5_8 * s**2 * h**2 * vd)

         ! calculate $v_{\pm sh}$
         vpmsh(:,1) = v   -   s * h * vd
         vpmsh(:,2) = v   +   s * h * vd

         ! calculate the rhs for solving the linear system M*vdpmsh = -g
         vdpmsh(:,1) = -this%gena_g(qpmsh(:,1), vpmsh(:,1), t)
         vdpmsh(:,2) = -this%gena_g(qpmsh(:,2), vpmsh(:,2), t)
         ! count calls
         this%gena_stats%ngcalls = this%gena_stats%ngcalls + 2
         ! after call of dgesv this variable will be overwritten and will
         ! actually contain $\dot v_{\pm sh}$

         if (this%opts%diag_mass_matrix == 1) then
            if (this%opts%const_mass_matrix == 1) then
               vdpmsh(:,1) = vdpmsh(:,1) / this%gena_const_diag_M
               vdpmsh(:,2) = vdpmsh(:,2) / this%gena_const_diag_M
            else
               vdpmsh(:,1) = vdpmsh(:,1) / this%gena_diag_M(qpmsh(:,1))
               vdpmsh(:,2) = vdpmsh(:,2) / this%gena_diag_M(qpmsh(:,2))
            end if
         else
            if (this%opts%const_mass_matrix == 1) then
               M(:,:,1) = this%gena_const_M
               call dgesv(         &! solve the System A*X=B and save the result in B
                          sv,      &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                          2,       &! number of right hand sides (=size(B,2))
                          M(:,:,1),&! matrix A
                          sv,      &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                          ipiv,    &! integer pivot vector; it is not needed
                          vdpmsh,  &! matrix B
                          sv,      &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                          info)     ! integer information flag
               ! Now vdpmsh actually contains $\dot v_{\pm sh}$
               if (info .ne. 0)  print*, "dgesv sagt info=", info ! TODO
            else
               M(:,:,1) = this%gena_M(qpmsh(:,1))
               M(:,:,2) = this%gena_M(qpmsh(:,2))
               call dgesv(            &! solve the System A*X=B and save the result in B
                          sv,         &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                          1,          &! number of right hand sides (=size(B,2))
                          M(:,:,1),   &! matrix A
                          sv,         &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                          ipiv,       &! integer pivot vector; it is not needed
                          vdpmsh(:,1),&! matrix B
                          sv,         &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                          info)        ! integer information flag
               if (info .ne. 0)  print*,  "dgesv sagt info=", info ! TODO
               call dgesv(            &! solve the System A*X=B and save the result in B
                          sv,         &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                          1,          &! number of right hand sides (=size(B,2))
                          M(:,:,2),   &! matrix A
                          sv,         &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                          ipiv,       &! integer pivot vector; it is not needed
                          vdpmsh(:,2),&! matrix B
                          sv,         &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                          info)        ! integer information flag
               if (info .ne. 0)  print*,  "dgesv sagt info=", info ! TODO
            end if
         end if

         ! Calculate a better approximation for a
         a = vd + (this%gena_alpha_m - this%gena_alpha_m) *      &
                (vdpmsh(:,2) - vdpmsh(:,1))/(2*s)

      end associate
   end subroutine gena_calcInitial

   subroutine gena_calcInitialConstrained(this, h) !TODO TODO
      implicit none
      class(gena_problem), intent(inout)  :: this  ! problem object
      real(8),             intent(in   )  :: h     ! step size
      ! internal variables
      integer                                                              :: i
      real(8), dimension(this%sizev+this%sizel)                            :: vdl
      real(8), dimension(this%sizeq, 2)                                    :: qpmsh   ! $q_{\pm sh}$
      real(8), dimension(this%sizev, 2)                                    :: vpmsh   ! $v_{\pm sh}$
      real(8), dimension(this%sizev+this%sizel, this%sizev+this%sizel,2)   :: MBB     ! Matrix
      real(8), dimension(this%sizev+this%sizel, 2)                         :: vdlpmsh ! $(\dot v_{\pm sh} \\ \lambda_{\pm sh})$
      real(8), dimension(this%sizel, this%sizev, 2)                        :: Bqpmsh  ! $B(q_{\pm sh})$
      integer, dimension(this%sizev+this%sizel)                            :: ipiv    ! pivot vector for dgesv
      integer, dimension(this%sizev+this%sizel)                            :: ipiv0   ! pivot vector for dgesv for LU factorization of MBB0, this is actually needed
      integer                                                              :: info    ! info flag for dgesv
      real(8), dimension(this%sizev)                                       :: hvdd    ! $\ddot v(t_0)$
      real(8), dimension(this%sizev)                                       :: lq0ph   ! $l_0^q/h$
      real(8), dimension(this%sizev+this%sizel, this%sizev+this%sizel)     :: MBB0    ! Matrix
      real(8), dimension(this%sizel, this%sizev)                           :: B0      ! $B(q_0)$
      real(8), dimension(this%sizev+this%sizel)                            :: DvDl    ! $(\Delta_0^v \\ \Delta_0^\lambda)$
      !

      associate ( t  => this%t,            &
                  q  => this%q,            &
                  v  => this%v,            &
                  vd => this%vd,           &
                  a  => this%a,            &
                  l  => this%l,            &
                  sv => this%sizev,        &
                  sl => this%sizel,        &
                  s  => this%opts%pertube_s)
         ! calulate $B(q_0)$
         B0 = this%gena_b(q)
         ! count calls
         this%gena_stats%nBcalls = this%gena_stats%nBcalls + 1

         ! calculate $MBB0$
         if (this%opts%diag_mass_matrix == 1) then
            ! Set the mass matrix part to zero beforehand
            MBB0(1:sv, 1:sv) = 0.0_8
            if (this%opts%const_mass_matrix == 1) then
               MBB0(1:sv, 1) = this%gena_const_diag_M
            else
               MBB0(1:sv, 1) = this%gena_diag_M(q)
            end if
            do concurrent (i=2:sv)
               MBB0(i,i) = MBB0(i,1)
               MBB0(i,1) = 0.0_8
            end do
         else
            if (this%opts%const_mass_matrix == 1) then
               MBB0(1:sv, 1:sv) = this%gena_const_M
            else
               MBB0(1:sv, 1:sv) = this%gena_M(q)
            end if
         end if
         MBB0(sv+1:sv+sl, 1:sv) =           B0
         MBB0(1:sv, sv+1:sv+sl) = transpose(B0)
         MBB0(sv+1:sv+sl, sv+1:sv+sl) = 0.0_8

         ! we need to solve a linear equation, first calulate the rhs
         vdl(1:sv)       = -this%gena_g(q, v, t)
         vdl(sv+1:sv+sl) = -this%gena_Z(q, v)
         ! count calls
         this%gena_stats%ngcalls = this%gena_stats%ngcalls + 1
         ! then solve the system
         call dgesv(          &! solve the System A*X=B and save the result in B
                     sv+sl,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                     1,       &! number of right hand sides (=size(B,2))
                     MBB0,    &! matrix A
                     sv+sl,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                     ipiv0,   &! integer pivot vector; it is not needed
                     vdl,     &! matrix B
                     sv+sl,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                     info)     ! integer information flag
         ! Now vdl actually contains $\dot v(t_0)$ and $\lambda(t_0)
         if (info .ne. 0)  print*, "dgesv sagt info=", info ! TODO

         ! apply the calculated values
         vd = vdl(   1:   sv)
         l  = vdl(sv+1:sv+sl)

         ! calulate $q_{\pm sh}$
         qpmsh(:,1) = this%gena_qlpexphDqtilde(q, h, &
                        -s * v   +   0.5_8 * s**2 * h * vd)  ! DEBUG
         qpmsh(:,2) = this%gena_qlpexphDqtilde(q, h, &
                         s * v   +   0.5_8 * s**2 * h * vd)  ! DEBUG

         ! calculate $v_{\pm sh}$
         vpmsh(:,1) = v   -   s * h * vd
         vpmsh(:,2) = v   +   s * h * vd

         ! calculate $B(q_{\pm sh})$
         Bqpmsh(:,:,1) = this%gena_b(qpmsh(:,1))
         Bqpmsh(:,:,2) = this%gena_b(qpmsh(:,2))
         ! count calls
         this%gena_stats%nBcalls = this%gena_stats%nBcalls + 2

         ! calculate the rhs for solving the linear system MBB*vdlpmsh = [-g \\ -Z]
         vdlpmsh(1:sv, 1) = -this%gena_g(qpmsh(:,1), vpmsh(:,1), this%t - s*h)
         vdlpmsh(sv+1:sv+sl, 1) = -this%gena_Z(qpmsh(:,1), vpmsh(:,1))
         vdlpmsh(1:sv, 2) = -this%gena_g(qpmsh(:,2), vpmsh(:,2), this%t + s*h)
         vdlpmsh(sv+1:sv+sl, 2) = -this%gena_Z(qpmsh(:,2), vpmsh(:,2))
         ! after call of dgesv this variable will be overwritten and will
         ! actually contain $\dot v_{\pm sh}$
         ! count calls
         this%gena_stats%ngcalls = this%gena_stats%ngcalls + 2

         ! use constant mass matrix, if possible
         if (this%opts%diag_mass_matrix == 1) then
            ! set MBB to zero in beforehand
            MBB(1:sv, 1:sv, :) = 0.0_8
            if (this%opts%const_mass_matrix == 1) then
               MBB(1:sv, 1, 1) = this%gena_const_diag_M
               MBB(1:sv, 1, 2) = this%gena_const_diag_M
            else
               MBB(1:sv, 1, 1) = this%gena_diag_M(qpmsh(:,1))
               MBB(1:sv, 1, 2) = this%gena_diag_M(qpmsh(:,2))
            end if
            ! shift the diagonal onto the actual diagonal
            do concurrent (i=2:sv)
               MBB(i,i, 1) = MBB(i,1, 1)
               MBB(i,1, 1) = 0.0_8
               MBB(i,i, 2) = MBB(i,1, 2)
               MBB(i,1, 2) = 0.0_8
            end do
         else
            if (this%opts%const_mass_matrix == 1) then
               MBB(1:sv, 1:sv, 1) = this%gena_const_M
               MBB(1:sv, 1:sv, 2) = this%gena_const_M
            else
               MBB(1:sv, 1:sv, 1) = this%gena_M(qpmsh(:,1))
               MBB(1:sv, 1:sv, 2) = this%gena_M(qpmsh(:,2))
            end if
         end if

         ! insert B
         MBB(1:sv, sv+1:sv+sl, 1) = transpose(Bqpmsh(:,:,1))
         MBB(sv+1:sv+sl, 1:sv, 1) =           Bqpmsh(:,:,1)
         MBB(1:sv, sv+1:sv+sl, 2) = transpose(Bqpmsh(:,:,2))
         MBB(sv+1:sv+sl, 1:sv, 2) =           Bqpmsh(:,:,2)

         ! insert zeros
         MBB(sv+1:sv+sl, sv+1:sv+sl, :) = 0.0_8
         !! DEBUG
         !call print_matrix(MBB(:,:,1),'MBBmsh')
         !call print_matrix(MBB(:,:,2),'MBBpsh')
         !call print_vector(vdlpmsh(:,1),'rhsmsh')
         !call print_vector(vdlpmsh(:,2),'rhspsh')
         !! GUBED
         ! solve for vdlpmsh
         call dgesv(             &! solve the System A*X=B and save the result in B
                    sv+sl,       &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                    1,           &! number of right hand sides (=size(B,2))
                    MBB(:,:,1),  &! matrix A
                    sv+sl,       &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                    ipiv,        &! integer pivot vector; it is not needed
                    vdlpmsh(:,1),&! matrix B
                    sv+sl,       &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                    info)         ! integer information flag
         if (info .ne. 0)  print*,  "dgesv sagt info=", info, __FILE__ ! TODO
         call dgesv(             &! solve the System A*X=B and save the result in B
                    sv+sl,       &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                    1,           &! number of right hand sides (=size(B,2))
                    MBB(:,:,2),  &! matrix A
                    sv+sl,       &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                    ipiv,        &! integer pivot vector; it is not needed
                    vdlpmsh(:,2),&! matrix B
                    sv+sl,       &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                    info)         ! integer information flag
         if (info .ne. 0)  print*,  "dgesv sagt info=", info ! TODO

         !! DEBUG
         !call print_vector(vdlpmsh(:,1),'vdlmsh')
         !call print_vector(vdlpmsh(:,2),'vdlpsh')
         !! GUBED

         ! calculate $\ddot v(t_0)$
         hvdd = (vdlpmsh(1:sv,2) - vdlpmsh(1:sv,1))/(2*s)

         !! DEBUG
         !call print_vector(hvdd/h, 'vdd')
         !! GUBED

         ! Calculate a better approximation for a
         a = vd + (this%gena_alpha_m - this%gena_alpha_f) * hvdd

         ! pertube initial values for v, if wanted and applicable
         if (this%opts%pertube == 1 .and. this%opts%stab2 .ne. 1) then

            ! calculate local truncation error divided by h
            !lq0 = (h**2)/6 * (                                                                  &
            !           (1 - 6*this%gena_beta - 3*(this%gena_alpha_m - this%gena_alpha_f))*hvdd  &
            !         + h*this%gena_itlbtvtw(v, vd)/2     )
            lq0ph = (1 - 6*this%gena_beta - 3*(this%gena_alpha_m - this%gena_alpha_f))/6 * h * hvdd &
                      +  h**2 * this%gena_itlbtvtw(v, vd) / 12 ! DEBUG

            ! calculate the rhs for solving the linear system MBB0*DvDl = [0 \\ B0*lq0/h]
            DvDl(1:sv) = 0.0_8
            DvDl(sv+1:sv+sl) = matmul(B0, lq0ph) ! DEBUG
            ! after call of dgesv this variable will be overwritten and will
            ! actually contain $(\Delta_0^v \\ \Delta_0^\lambda)$

            ! solve for DvDl (MBB0 has been LU-factorized by dgesv before)
            call dgetrs(                &! solve the system A*X=B and save the result in B
                        'No transpose', &! A is not a transposed matrix
                        sv+sl,          &! number of linear equations
                        1,              &! number of right hand sides
                        MBB0,           &! matrix A
                        sv+sl,          &! leading dimension of A, in this case is equal to the number of linear equations
                        ipiv0,          &! integer pivot vector; it is not really needed
                        DvDl,           &! matrix B
                        sv+sl,          &! leading dimension of B, in this case is equal to the number of linear equations
                        info)            ! integer information flag
            !call dgesv(                      &! solve the System A*X=B and save the result in B
            !           sv+sl,&! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
            !           1,                    &! number of right hand sides (=size(B,2))
            !           MBB0,                 &! matrix A
            !           sv+sl,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
            !           ipiv,                 &! integer pivot vector; it is not needed
            !           DvDl,                 &! matrix B
            !           sv+sl,&! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
            !           info)                  ! integer information flag
            if (info .ne. 0)  print*,  "dgesv sagt info=", info ! TODO

            ! calculate correctly pertubed (initial) velocity
            v = v + DvDl(   1:   sv)
            !v = [4.61536480720402_8, &
            !             0.0_8,        &
            !             0.0_8,        &
            !             0.0_8,        &
            !           150.0_8,        &
            ! -4.61635233894283_8]
            l = l + DvDl(sv+1:sv+sl)

            !! DEBUG
            !call print_vector(lq0ph*h, 'lq0')
            !call print_vector(DvDl(1:sv), 'Deltav')
            !call print_vector(v, 'v_per')
            !call print_vector(a, 'a_per')
            !! GUBED
         end if

      end associate
   end subroutine gena_calcInitialConstrained

   ! subroutine for integrating one time step
   subroutine gena_solveConstrainedTimeStep(this,t1)
      implicit none
      class(gena_problem), intent(inout)  :: this  ! problem object
      real(8)            , intent(in   )  :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                                      :: i     ! for iteration
      integer, dimension(this%sizev+this%sizel)    :: ipiv  ! needed for dgesv
      integer                                      :: info  ! needed for dgesv

      ! intenal real variables
      real(8)                                                           :: h       ! step size
      real(8), dimension(this%sizev)                                    :: vd1     ! $\dot v_{n+1}$
      !real(8), dimension(this%sizev)                                    :: a1      ! $a_{n+1}$ not needed 2015-06-23
      real(8), dimension(this%sizev)                                    :: v1      ! $v_{n+1}$
      real(8), dimension(this%sizev)                                    :: Dq      ! $\Delta q_n$
      real(8), dimension(this%sizeq)                                    :: q1      ! $q_{n+1}}
      real(8), dimension(this%sizev+this%sizel)                         :: res     ! $res(\dots)$
      real(8), dimension(this%sizev+this%sizel,this%sizev+this%sizel)   :: St      ! $S_t$
      real(8), dimension(:,:), allocatable                              :: Stb     ! permuted $S_t$ in banded format
      !real(8)                                                           :: beta_d  ! $\beta'$ not needed 2015-06-23
      !real(8)                                                           :: gamma_d ! $\gamma'$ not needed 2015-06-23
      real(8), dimension(this%sizev+this%sizel)                         :: Dxl     ! $(\Delta x, \Delta \lambda)^T$

      ! internal logical
      logical                                                           :: converged

      ! associate construct for better readability
      associate (alf => this%gena_alpha_f,    &
                 alm => this%gena_alpha_m,    &
                 bet => this%gena_beta,       &
                 gam => this%gena_gamma,      &
                 sv  => this%sizev,           &
                 sl  => this%sizel,           &
                 v   => this%v,               &
                 vd  => this%vd,              &
                 a   => this%a,               &
                 q   => this%q,               &
                 t   => this%t,               &
                 l   => this%l)
         ! calculation of step size $h$
         h = t1 - t

#ifdef ZEROINIT
         Dq = 0.0_8
         v1 = (1 - gam/bet)*v + h*(1 - gam/(2*bet))*a
         vd1 = (1-alm)/(bet*(1-alf))*( (Dq - v)/h - 0.5*a)  +  (a - alf*vd)/(1-alf)
         l = 0.0_8
#else
         ! initialization of the new values $\dot v_{n+1}$, $a_{n+1}$ and $v_{n+1}$ and the value $\Delta q_n$ and $\lambda_{n+1}$
         !vd1 = 0.0_8 ! vd ! 0.0_8 DEBUG
         !a1  = (alf*vd - alm*a)/(1.0_8 - alm)
         !v1  = v + h*(1.0_8 - gam)*a + gam*h*a1
         !Dq  = v + (0.5_8 - bet)*h*a + bet*h*a1
         !!l   = l ! DEBUG
         !l   = 0.0_8
         Dq  = v + 0.5*h*a
         v1  = v +     h*a
         vd1 = (a - alf*vd)/(1-alf)
#endif


         ! loop of the newton method
         converged = .false.
         do i=1,this%opts%imax
            !! DEBUG
            !print *, t, i
            !print *, Dq
            !! GUBED
            !print *, 'Newton: ', i, this%opts%imax DEBUG
            ! calculate the new value $q_{n+1}$
            q1  = this%gena_qlpexphDqtilde(q,h,Dq)

            ! caluclate the residue $res$
            res(1:sv) = this%gena_g(q1,v1,t1) + matmul(transpose(this%gena_B(q1)),l)
            ! count calls
            this%gena_stats%ngcalls = this%gena_stats%ngcalls + 1
            this%gena_stats%nBcalls = this%gena_stats%nBcalls + 1
            !
            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  res(1:sv) = res(1:sv) + this%gena_const_diag_M*vd1
               else
                  res(1:sv) = res(1:sv) + matmul(this%gena_const_M,vd1)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  res(1:sv) = res(1:sv) + this%gena_diag_M(q1)*vd1
               else
                  res(1:sv) = res(1:sv) + matmul(this%gena_M(q1),vd1)
               end if
            end if
            ! scale
            res(1:sv) = h*res(1:sv)
            res(sv+1:sv+sl) = this%gena_phi(q1)/h

            !! check if norm of the residue is sufficiently small to stop the iteration
            !if ((this%gena_norm(res(1:sv))       < this%opts%tolr  )  .and. &
            !    (this%gena_norm(res(sv+1:sv+sl)) < this%opts%tolphi)) then
            !   exit
            !end if

            ! solve the linear System $S_t \cdot \Delta xl = -res$
            Dxl = -res ! the result will be saved in the right hand side's spot

            if (i == 1 .or. this%opts%recalc_iteration_matrix == 1) then
               ! calculate iteration matrix
               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     St(1:sv,1:sv) = 0.0_8
                     forall (i=1:sv)
                        St(    i   ,     i   ) = this%gena_const_diag_M(i) * (1-alm)/(bet*(1-alf))
                     end forall
                  else
                     St(   1:sv   ,   1:sv   ) = this%gena_const_M * (1-alm)/(bet*(1-alf))
                  end if
               else
                  if (this%opts%diag_mass_matrix == 1) then
                     St(1:sv,1:sv) = 0.0_8
                     St(   1:sv   ,   1   ) = this%gena_diag_M(q1) * (1-alm)/(bet*(1-alf))
                     forall (i=2:sv)
                        St(i,i) = St(i,1)
                        St(i,1) = 0.0_8
                     end forall
                  else
                     St(   1:sv   ,   1:sv   ) = this%gena_M(q1) * (1-alm)/(bet*(1-alf))
                  end if
               end if
               St(   1:sv   ,sv+1:sv+sl) = transpose(this%gena_B(q1)) ! TODO: dont calculate twice
               St(sv+1:sv+sl,   1:sv   ) = matmul(this%gena_B(q1), this%gena_Tg(h,Dq))
               St(sv+1:sv+sl,sv+1:sv+sl) = 0.0_8
               ! count calls
               this%gena_stats%nBcalls = this%gena_stats%nBcalls + 2

               if (this%opts%no_Ct == 0) then
                  if (this%opts%use_num_Ct == 1) then
                     St(1:sv,1:sv) = St(1:sv,1:sv) + this%gena_num_Ct(q1,v1,t1) * h*gam/bet
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + this%sizev + 1
                     end if
                  else
                     St(1:sv,1:sv) = St(1:sv,1:sv) + this%gena_Ct(q1,v1,t1) * h*gam/bet
                  end if
               end if

               if (this%opts%no_Kt == 0) then
                  if (this%opts%use_num_Kt == 1) then
                     St(1:sv,1:sv) = St(1:sv,1:sv) + matmul(this%gena_num_Kt_lambda(q1,v1,vd1,l,t1),this%gena_Tg(h,Dq)) * h*h ! TODO: dgemm oder sowas
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                        this%gena_stats%nBcalls = this%gena_stats%nBcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + this%sizev + 1
                        this%gena_stats%nBcalls = this%gena_stats%nBcalls + this%sizev + 1
                     end if
                  else
                     St(1:sv,1:sv) = St(1:sv,1:sv) + matmul(this%gena_Kt_lambda(q1,v1,vd1,l,t1),this%gena_Tg(h,Dq)) * h*h ! TODO: dgemm oder sowas
                  end if
               end if

               !! DEBUG !
               !print *, "t=", t
               !print *, "cond=", mycond(St)
               !print *, "maxstep=", this%gena_stats%newt_steps_max
               !!print *,
               !call print_vector(q,'qn')
               !call print_vector(v,'vn')
               !call print_vector(vd,'vdn')
               !call print_vector(a,'an')
               !call print_vector(l,'ln')
               !call print_vector(res,'resn')
               !call print_matrix(St,'Stn')
               !errorstop "mopp"
               !if (t > 2.5) then
               !   call print_matrix(St,'St')
               !   !call print_vector_int(this%opts%jour,'jour')
               !   stop
               !end if
               !! GUBED !

               !!DEBUG
               ! call print_matrix(St, 'Stb')
               !!GUBED

                        !!print *, St ! DEBUG
               ! Actually solve the linear System $S_t \cdot \Delta xl = -res$
               if (this%opts%banded_iteration_matrix == 1) then
#ifdef DEBUG_PRINT_ITERATION_MATRIX_AT
                  if (t >= DEBUG_PRINT_ITERATION_MATRIX_AT) then
                     call print_matrix(St,'St')
                     call print_vector_int(this%opts%jour,'jour')
                     stop "Successfully printed iteration matrix"
                  end if
#endif
                  ! TODO
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag)
                     if (.not. allocated(Stb)) then
                        allocate(Stb(lo*2+hi+1, sv+sl))
                     end if
                     !Stb = 0.0_8 ! DEBUG
                     forall (i=1:sv+sl)
                        Stb((max(1, i-hi)-i+lo+1+hi):(min(sv+sl,i+lo)-i+lo+1+hi),i) &
                            = St(this%opts%jour(max(1, i-hi):min(sv+sl, i + lo)), this%opts%jour(i))
                     end forall
                     Dxl = Dxl(this%opts%jour)
                     call dgbsv(           &! solve the system A*X=B and save the result in B
                                 sv+sl,    &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,       &!
                                 hi,       &!
                                 1,        &! number of right hand sides (=size(B,2))
                                 Stb,      &! matrix A
                                 2*lo+hi+1,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,     &! integer pivot vector; it is not needed
                                 Dxl,      &! matrix B
                                 sv+sl,    &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)      ! integer information flag
                     Dxl(this%opts%jour) = Dxl
                     !call print_vector(Dxl,'Dxl_b') ! DEBUG
                     !stop ! DEBUG
                  end associate
               else
                  call dgesv(       &! solve the System A*X=B and save the result in B
                              sv+sl,&! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,    &! number of right hand sides (=size(B,2))
                              St,   &! matrix A
                              sv+sl,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv, &! integer pivot vector; it is not needed
                              Dxl,  &! matrix B
                              sv+sl,&! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)  ! integer information flag
                  !call print_vector(Dxl,'Dxl') ! DEBUG
                  !stop ! DEBUG
                  ! TODO: Errors abfragen
               end if

               if (info.ne.0) print *, "dgesv sagt info=", info

               ! update $\dot v_{n+1}$, $v_{n+1}$ and $\Delta q_n$
            else ! (St was not (re)calculated)
               ! Actually solve the linear System $S_t \cdot \Delta xl = -res$
               if (this%opts%banded_iteration_matrix == 1) then
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag, &
                             uStb => ubound(Stb,1))
                     Dxl = Dxl(this%opts%jour) !TODO
                     call dgbtrs(               &
                                 'No transpose',&! solve the System A*X=B and save the result in B
                                 sv+sl,         &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,            &!
                                 hi,            &!
                                 1,             &! number of right hand sides (=size(B,2))
                                 Stb,           &! matrix A
                                 2*lo+hi+1,     &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,          &! integer pivot vector; it is not needed
                                 Dxl,           &! matrix B
                                 sv+sl,         &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)           ! integer information flag
                     Dxl(this%opts%jour) = Dxl
                  end associate
               else
                  call dgetrs(                &
                              'No transpose', & ! solve the System A*X=B and save the result in B
                              sv + sl,        & ! number of linear equations (=size(A,1))  ! Vorsicht: double precision muss real(8) sein
                              1,              & ! number of right hand sides (=size(B,2))
                              St,             & ! matrix A
                              sv + sl,        & ! leading dimesnion of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv,           & ! integer pivot vector; holds the pivoting that was calculated in the dgesv call
                              Dxl,            & ! matrix B
                              sv + sl,        & ! leading dimension of B, in thie case is equal to the number of linear equations (=size(A,1))
                              info)             ! integer information flag
               end if
               ! DEBUG
               if (info.ne.0) print*, 'dgetrs sagt info=', info
               ! GUBED
            end if

            !! DEBUG
            !call print_vector(Dxl,'Dxln')
            !stop
            !! GUBED

            ! Check for convergence
            if ( ( sum((Dxl(   1:   sv) / (this%opts%atol + this%opts%rtol*abs(Dq )))**2 )  &
                  +sum((Dxl(sv+1:sv+sl) / (this%opts%atol + this%opts%rtol*abs(h*l)))**2 ) )&
                   / (sv+sl)   <=  1.0_8 ) then
               converged = .true.
            end if

            ! update $\dot v_{n+1}$, $v_{n+1}$, $\Delta q_n$ and $\lambda_{n+1}$
#define REL *1
            Dq  = Dq  + Dxl(1:sv)                              REL
            v1  = v1  + gam/bet * Dxl(1:sv)                    REL
            vd1 = vd1 + (1-alm)/(bet*(1-alf)*h) * Dxl(1:sv)    REL
            l   = l   + Dxl(sv+1:sv+sl)/h                      REL
#undef REL

            ! update solver stats
            this%gena_stats%newt_steps_curr = i

            ! exit loop, if converged
            if ( converged ) exit

         end do ! end newton method

         ! DEBUG
         if (this%gena_stats%newt_steps_curr >= this%opts%imax) then
            print *, ' ' ! DEBUG
            print *, "Newton did not converge at t=", t
            !! DEBUG !
            call this%gena_print_stats
            print *, "Iteration matrix:"
            call print_matrix(St,'St')
            errorstop
            !! GUBED !
         end if


         ! The integration step is done; now the new values have to be
         ! applied to the problem object, and the new value for a can
         ! be calculated
         a  = ((1-alf)*vd1 + alf*vd - alm*a) / (1-alm)
         t  = t1
         q  = q1
         v  = v1
         vd = vd1

      end associate

   end subroutine gena_solveConstrainedTimeStep

   ! subroutine for integrating one time step
   subroutine gena_solveConstrainedTimeStep_stab2(this,t1)
      implicit none
      class(gena_problem), intent(inout)  :: this  ! problem object
      real(8)            , intent(in   )  :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                                      :: i     ! for iteration
      integer, dimension(this%sizev+2*this%sizel)  :: ipiv  ! needed for dgesv
      integer                                      :: info  ! needed for dgesv

      ! intenal real variables
      real(8)                                                              :: h       ! step size
      real(8), dimension(this%sizev)                                       :: vd1     ! $\dot v_{n+1}$
      real(8), dimension(this%sizev)                                       :: v1      ! $v_{n+1}$
      real(8), dimension(this%sizev)                                       :: Dq      ! $\Delta q_n$
      real(8), dimension(this%sizeq)                                       :: q1      ! $q_{n+1}}
      real(8), dimension(this%sizev+2*this%sizel)                          :: res     ! $res(\dots)$
      real(8), dimension(this%sizev+2*this%sizel,this%sizev+2*this%sizel)  :: St      ! $S_t$
      real(8), dimension(this%sizev,this%sizev)                            :: Mstar   ! $M^\ast$
      real(8), dimension(this%sizel,this%sizev)                            :: Bqn     ! $B(q_n)$
      real(8), dimension(this%sizel,this%sizev)                            :: B       ! $B$
      real(8), dimension(:,:), allocatable                                 :: Stb     ! permuted $S_t$ in banded format
      real(8), dimension(this%sizev+2*this%sizel)                          :: Dxleta  ! $(\Delta x, \Delta \lambda, \Delta \eta)^T$

      ! internal logical
      logical                                                              :: converged
#ifdef DEBUG2
      integer :: sigma, sigmas
#endif

      ! associate construct for better readability
      associate (alf => this%gena_alpha_f,    &
                 alm => this%gena_alpha_m,    &
                 bet => this%gena_beta,       &
                 gam => this%gena_gamma,      &
                 sv  => this%sizev,           &
                 sl  => this%sizel,           &
                 v   => this%v,               &
                 vd  => this%vd,              &
                 a   => this%a,               &
                 q   => this%q,               &
                 t   => this%t,               &
                 l   => this%l,               &
                 eta => this%eta              )
         ! calculation of step size $h$
         h = t1 - t

         ! calculate B(q_n)
         Bqn = this%gena_B(q)
         ! count calls
         this%gena_stats%nBcalls = this%gena_stats%nBcalls + 1

!#define MY_T 3.28_8
#define MY_T 3.2_8

#ifdef ZEROINIT
         Dq = 0.0_8
         v1 = (1 - gam/bet)*v + h*(1 - gam/(2*bet))*a
         vd1 = (1-alm)/(bet*(1-alf))*( (Dq - v)/h - 0.5*a)  +  (a - alf*vd)/(1-alf)
         l = 0.0_8
         eta = 0.0_8
#else
#ifdef DEBUG1
         ! DEBUG
         if (t < MY_T) then
            Dq = 0.0_8
            v1 = (1 - gam/bet)*v + h*(1 - gam/(2*bet))*a
            vd1 = (1-alm)/(bet*(1-alf))*( (Dq - v)/h - 0.5*a)  +  (a - alf*vd)/(1-alf)
            l = 0.0_8
            eta = 0.0_8
         else
            print *, 't = ', t
         end if
         ! GUBED
#endif

         ! initialization of the new values $\dot v_{n+1}$, $a_{n+1}$ and $v_{n+1}$ and the value $\Delta q_n$
         ! choice \dot v_{n+1}^{(0)} = 0, probably bad
         !vd1 =  0.0_8 ! vd ! 0.0_8 DEBUG DEBUG DEBUG
         !a1  = (vd - alm*a)/(1.0_8 - alm) !(alf*vd - alm*a)/(1.0_8 - alm) DEBUG
         !v1  = v + h*(1.0_8 - gam)*a + gam*h*a1
         !Dq  = v + (0.5_8 - bet)*h*a + bet*h*a1 !+matmul(Bqn,eta), but eta=0
         ! choice a_{n+1}^{(0)} = a_n, probably better

         Dq  = v + 0.5*h*a - matmul(transpose(Bqn),eta)
         v1 = v + h*a
         vd1 = (a - alf*vd)/(1-alf)

         !Dq = v + 0.5*h*a
         !v1  = gam/bet*(Dq + matmul(transpose(Bqn),eta)) + (1 - gam/bet)*v + h*(1 - gam/(2*bet))*a
         !vd1 = (1-alm)/(bet*(1-alf))*( (Dq + matmul(transpose(Bqn),eta) - v)/h - 0.5*a)  +  (a - alf*vd)/(1-alf)


         !Dq = v + h*a
         !v1  = gam/bet*(Dq + matmul(transpose(Bqn),eta)) + (1 - gam/bet)*v + h*(1 - gam/(2*bet))*a
         !vd1 = (1-alm)/(bet*(1-alf))*( (Dq + matmul(transpose(Bqn),eta) - v)/h - 0.5*a)  +  (a - alf*vd)/(1-alf)

         !v1  = gam/bet*Dq + (1 - gam/bet)*v + h*(1 - gam/(2*bet))*a  ! is the same (in exact arithmetics) as


         !!!!!!v1 = v + h*a + gam/bet*matmul(transpose(Bqn),eta)
         !vd1 = (1-alm)/(bet*(1-alf))*( (Dq - v)/h - 0.5*a)  +  (a - alf*vd)/(1-alf) ! is the same (in exact arithmetics) as

         !!!!!!!!vd1 = (a - alf*vd)/(1-alf) + (1-alm)/((1-alf)*bet*h)*matmul(transpose(Bqn),eta)
         !vd1 = (1-alm)/(bet*(1-alf))*( (Dq + matmul(transpose(Bqn),eta) - v)/h - 0.5*a)  +  (a - alf*vd)/(1-alf)

         ! initialization of and $\lambda_{n+1}$, $\eta_{n+1}
         !l   = l    ! this is not necessary to do, obviously
         !eta = eta  ! dito
#endif

         ! loop of the newton method
         converged = .false.
         do i=1,this%opts%imax
#ifdef DEBUG1
            ! DEBUG
            if (t>=MY_T) then
               write(*,'(A)',advance='NO') 'Dq'
#ifdef ZEROINIT
               write(*,'(A)',advance='NO') 'zero'
#else
               write(*,'(A)',advance='NO') 'ours'
#endif
               write(*,'(A,I5,A)',advance='NO') '(', i, ',:) ...'
               call print_vector(Dq,'')

               write(*,'(A)',advance='NO') 'hl'
#ifdef ZEROINIT
               write(*,'(A)',advance='NO') 'zero'
#else
               write(*,'(A)',advance='NO') 'ours'
#endif
               write(*,'(A,I5,A)',advance='NO') '(', i, ',:) ...'
               call print_vector(h*l,'')

               write(*,'(A)',advance='NO') 'eta'
#ifdef ZEROINIT
               write(*,'(A)',advance='NO') 'zero'
#else
               write(*,'(A)',advance='NO') 'ours'
#endif
               write(*,'(A,I5,A)',advance='NO') '(', i, ',:) ...'
               call print_vector(eta,'')
            end if
            ! GUBED
#endif
            !! DEBUG
            !print *, '------------------------------------------------------'
            !print *, 't', t, 'i', i
            !print *, 'eta', eta
            !print *, 'Bqn', Bqn
            !print *, 'Dq', Dq
            !print *, 'v1', v1
            !print *, 'vd1', vd1
            !print *, ''
            !! GUBED
            !print *, 'Newton: ', i, this%opts%imax DEBUG

#ifdef DEBUG2
#define MY_I 5
            if (t>= MY_T .and. i >= MY_I) then
               sigmas = 500
            else
               sigmas = 0
            end if

            do sigma=1,sigmas+1
               if (sigma == 2) then
                  print *, 'clear J p sigma Fx'
                  print *, 'sigma(1)=0'
                  call print_matrix(St,'J')
                  call print_vector(Dxleta,'p')
               end if

               if (sigma > 1) then
                  write(*,'(A,I5,A)',advance='no') 'Fx(:,', sigma-1, ')...'
                  call print_vector(-res,'')

                  !if (sigma > 1e-4 * sigmas) STOP
                  if (sigma == sigmas + 1) STOP

                  print *, 'sigma(', sigma, ')=', (sigma-1.0_8)/(sigmas-1), ';'

                  Dq  = Dq  + 1/(sigmas-1.0_8) * Dxleta(1:sv)
                  l   = l   + 1/(sigmas-1.0_8) * Dxleta(sv   +1:sv+sl  )/h
                  eta = eta + 1/(sigmas-1.0_8) * Dxleta(sv+sl+1:sv+sl+sl)
               end if

               v1  = gam/bet*(Dq + matmul(transpose(Bqn),eta)) + (1 - gam/bet)*v + h*(1 - gam/(2*bet))*a
               vd1 = (1-alm)/(bet*(1-alf))*( (Dq + matmul(transpose(Bqn),eta) - v)/h - 0.5*a)  +  (a - alf*vd)/(1-alf)
#endif

            ! calculate the new value $q_{n+1}$
            q1  = this%gena_qlpexphDqtilde(q,h,Dq)

            ! calculate new B(q1)
            B = this%gena_B(q1)
            ! count calls
            this%gena_stats%nBcalls = this%gena_stats%nBcalls + 1

            ! caluclate the residue $res_h$
            res(1:sv) = this%gena_g(q1,v1,t1) + matmul(transpose(B),l)
            ! count calls
            this%gena_stats%ngcalls = this%gena_stats%ngcalls + 1
            !
            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  res(1:sv) = res(1:sv) + this%gena_const_diag_M*vd1
               else
                  res(1:sv) = res(1:sv) + matmul(this%gena_const_M,vd1)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  res(1:sv) = res(1:sv) + this%gena_diag_M(q1)*vd1
               else
                  res(1:sv) = res(1:sv) + matmul(this%gena_M(q1),vd1)
               end if
            end if
            res(1:sv) = res(1:sv) * h
            res(sv+1:sv+sl) = this%gena_phi(q1) / h
            res(sv+sl+1:sv+sl+sl) = matmul(B, v1)


            !! check if norm of the residue is sufficiently small to stop the iteration
            !if ((this%gena_norm(res(1:sv))             < this%opts%tolr  )  .and. &
            !    (this%gena_norm(res(sv+1:sv+sl))       < this%opts%tolphi)  .and. &  ! TODO
            !    (this%gena_norm(res(sv+sl+1:sv+sl+sl)) < this%opts%tolB  )    ) then
            !   exit
            !end if

#ifdef DEBUG2
               if (sigma .ne. 1) cycle
#endif

            ! solve the linear System $S_t \cdot \Delta xl = -res$
            Dxleta = -res ! the result will be saved in the right hand side's spot


            if (i == 1 .or. this%opts%recalc_iteration_matrix == 1) then
#ifdef STNUM
               St = this%St_num(Dq,h,l,eta,t1)
#else
               ! calculate iteration matrix
               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     Mstar = 0.0_8
                     forall (i=1:sv)
                        Mstar(i, i) = this%gena_const_diag_M(i)
                     end forall
                  else
                     Mstar = this%gena_const_M
                  end if
               else
                  if (this%opts%diag_mass_matrix == 1) then
                     Mstar = 0.0_8
                     Mstar(1:sv, 1) = this%gena_diag_M(q1)
                     forall (i=2:sv)
                        Mstar(i,i) = Mstar(i,1)
                        Mstar(i,1) = 0.0_8
                     end forall
                  else
                     Mstar = this%gena_M(q1)
                  end if
               end if
               Mstar = Mstar*(1-alm)/(bet*(1-alf))

               if (this%opts%no_Ct == 0) then
                  ! Add damping matrix to Mstar
                  if (this%opts%use_num_Ct == 1) then
                     Mstar = Mstar + h*gam/bet * this%gena_num_Ct(q1,v1,t1)
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + this%sizev + 1
                     end if
                  else
                     Mstar = Mstar + h*gam/bet * this%gena_Ct(q1,v1,t1)
                  end if
               end if

               St(      1:sv,             1:sv      ) = Mstar
               St(      1:sv      ,    sv+1:sv+sl   ) = transpose(B)
               St(sv   +1:sv+sl   ,       1:sv      ) = matmul(B, this%gena_Tg(h,Dq))
               St(sv   +1:sv+sl   ,    sv+1:sv+sl   ) = 0.0_8
#ifdef NOZ
               St(sv+sl+1:sv+sl+sl,       1:sv      ) = gam/bet*B
#else
               St(sv+sl+1:sv+sl+sl,       1:sv      ) = gam/bet*B + h*this%gena_matZ(q1,v1,this%gena_Tg(h,Dq))
#endif
               St(sv+sl+1:sv+sl+sl, sv   +1:sv+sl   ) = 0.0_8
               St(sv+sl+1:sv+sl+sl, sv+sl+1:sv+sl+sl) = gam/bet*matmul(B,transpose(Bqn))
               St(sv   +1:sv+sl   , sv+sl+1:sv+sl+sl) = 0.0_8
               St(      1:sv      , sv+sl+1:sv+sl+sl) = matmul(Mstar,transpose(Bqn))

               if (this%opts%no_Kt == 0) then
                  if (this%opts%use_num_Kt == 1) then
                     St(1:sv,1:sv) = St(1:sv,1:sv) + h*h*matmul(this%gena_num_Kt_lambda(q1,v1,vd1,l,t1),this%gena_Tg(h,Dq)) ! TODO: dgemm oder sowas
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                        this%gena_stats%nBcalls = this%gena_stats%nBcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + this%sizev + 1
                        this%gena_stats%nBcalls = this%gena_stats%nBcalls + this%sizev + 1
                     end if
                  else
                     St(1:sv,1:sv) = St(1:sv,1:sv) + h*h*matmul(this%gena_Kt_lambda(q1,v1,vd1,l,t1),this%gena_Tg(h,Dq)) ! TODO: dgemm oder sowas
                  end if
               end if
! endif STNUM
#endif

               !!DEBUG
               !if (t1 > 0.1_8) then
               !if ( i >= 2) then
               !print *, i
               !call print_matrix(St,'St')
               !   pause
               !   call print_vector(Dxleta,'Dxleta')
               !end if
               !if (t1 > 0.11_8)  stop
               !!GUBED

               !! DEBUG !
               !print *, "t=", t, "cond=", mycond(St), "nmax=", this%gena_stats%newt_steps_max
               !print *,
               !if (t > 4.5_8) then
               !   call print_matrix(St,'St')
               !   stop
               !end if
               !call print_vector_int(this%opts%jour,'jour')
               !stop
               !! GUBED !

                        !!print *, St ! DEBUG
               ! Actually solve the linear System $S_t \cdot \Delta xleta = -res$
               if (this%opts%banded_iteration_matrix == 1) then
#ifdef DEBUG_PRINT_ITERATION_MATRIX_AT
                  if (t >= DEBUG_PRINT_ITERATION_MATRIX_AT) then
                     call print_matrix(St,'St')
                     call print_vector_int(this%opts%jour,'jour')
                     stop "Successfully printed iteration matrix"
                  end if
#endif
                  ! TODO
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag)
                     if (.not. allocated(Stb)) then
                        allocate(Stb(lo*2+hi+1, sv+2*sl))
                     end if
                     !Stb = 0.0_8 ! DEBUG
                     forall (i=1:sv+2*sl)
                        Stb((max(1, i-hi)-i+lo+1+hi):(min(sv+sl,i+lo)-i+lo+1+hi),i) &
                            = St(this%opts%jour(max(1, i-hi):min(sv+2*sl, i + lo)), this%opts%jour(i))
                     end forall
                     Dxleta = Dxleta(this%opts%jour)
                     call dgbsv(           &! solve the system A*X=B and save the result in B
                                 sv+2*sl,  &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,       &!
                                 hi,       &!
                                 1,        &! number of right hand sides (=size(B,2))
                                 Stb,      &! matrix A
                                 2*lo+hi+1,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,     &! integer pivot vector; it is not needed
                                 Dxleta,   &! matrix B
                                 sv+2*sl,  &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)      ! integer information flag
                     Dxleta(this%opts%jour) = Dxleta
                     !call print_vector(Dxl,'Dxl_b') ! DEBUG
                     !stop ! DEBUG
                  end associate
               else
                  !!! DEBUG
                  !if (t > 5.3_8) then
                  !   call print_matrix(St,'St')
                  !   call print_vector_int(this%opts%jour,'jour')
                  !   stop
                  !end if
                  !!! GUBED
                  call dgesv(          &! solve the System A*X=B and save the result in B
                              sv+sl+sl,&! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,       &! number of right hand sides (=size(B,2))
                              St,      &! matrix A
                              sv+sl+sl,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv,    &! integer pivot vector; it is not needed
                              Dxleta,  &! matrix B
                              sv+sl+sl,&! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)     ! integer information flag
                  !call print_vector(Dxl,'Dxl') ! DEBUG
                  !stop ! DEBUG
                  ! TODO: Errors abfragen
               end if

               if (info.ne.0) print *, "dgesv sagt info=", info

               ! update $\dot v_{n+1}$, $v_{n+1}$ and $\Delta q_n$
            else ! (St was not (re)calculated)
               ! Actually solve the linear System $S_t \cdot \Delta xleta = -res$
               if (this%opts%banded_iteration_matrix == 1) then
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag, &
                             uStb => ubound(Stb,1))
                     Dxleta = Dxleta(this%opts%jour) !TODO
                     call dgbtrs(               &
                                 'No transpose',&! solve the System A*X=B and save the result in B
                                 sv+2*sl,       &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,            &!
                                 hi,            &!
                                 1,             &! number of right hand sides (=size(B,2))
                                 Stb,           &! matrix A
                                 2*lo+hi+1,     &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,          &! integer pivot vector; it is not needed
                                 Dxleta,        &! matrix B
                                 sv+2*sl,       &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)           ! integer information flag
                     Dxleta(this%opts%jour) = Dxleta
                  end associate
               else
                  call dgetrs(                &
                              'No transpose', & ! solve the System A*X=B and save the result in B
                              sv + sl + sl,   & ! number of linear equations (=size(A,1))  ! Vorsicht: double precision muss real(8) sein
                              1,              & ! number of right hand sides (=size(B,2))
                              St,             & ! matrix A
                              sv + sl + sl,   & ! leading dimesnion of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv,           & ! integer pivot vector; holds the pivoting that was calculated in the dgesv call
                              Dxleta,         & ! matrix B
                              sv + sl + sl,   & ! leading dimension of B, in thie case is equal to the number of linear equations (=size(A,1))
                              info)             ! integer information flag
               end if
               ! DEBUG
               if (info.ne.0) print*, 'dgetrs sagt info=', info
               ! GUBED
            end if

#ifdef DEBUG2
            end do
#endif

#ifdef DEBUG1
            ! DEBUG
            if (t>=MY_T) then
               write(*,'(A)',advance='NO') 'Dxleta'
#ifdef ZEROINIT
               write(*,'(A)',advance='NO') 'zero'
#else
               write(*,'(A)',advance='NO') 'ours'
#endif
               write(*,'(A,I5,A)',advance='NO') '(', i, ',:) ...'
               call print_vector(Dxleta,'')
            end if
            ! GUBED
#endif

            ! Check for convergence
            if ( ( sum((Dxleta(      1:      sv) / (this%opts%atol + this%opts%rtol*abs(Dq )))**2 )  &
                  +sum((Dxleta(sv+   1:   sv+sl) / (this%opts%atol + this%opts%rtol*abs(h*l)))**2 )  &
                  +sum((Dxleta(sv+sl+1:sv+sl+sl) / (this%opts%atol + this%opts%rtol*abs(eta)))**2 ) )&
                   / (sv+2*sl)   <=  1.0_8 ) then
               converged = .true.
            end if

            !!DEBUG
            !print *, 'norm:', norm2(Dxleta)
            !!GUBED

            ! update $\dot v_{n+1}$, $v_{n+1}$, $\Delta q_n$ and $\lambda_{n+1}$, $\eta_{n+1}$
! DEBUGGETY
#define REL *1
            Dq  = Dq  + Dxleta(1:sv)                                                                             REL
            v1  = v1  + gam/bet*(Dxleta(1:sv) + matmul(transpose(Bqn),Dxleta(sv+sl+1:sv+sl+sl)))                 REL
            vd1 = vd1 + (1-alm)/(bet*(1-alf)*h)*(Dxleta(1:sv) + matmul(transpose(Bqn),Dxleta(sv+sl+1:sv+sl+sl))) REL
            l   = l   + Dxleta(sv   +1:sv+sl  )/h                                                                REL
            eta = eta + Dxleta(sv+sl+1:sv+sl+sl)                                                                 REL
#undef REL
            !v1  = gam/bet*(Dq + matmul(transpose(Bqn),eta)) + (1 - gam/bet)*v + h*(1 - gam/(2*bet))*a
            !vd1 = (1-alm)/(bet*(1-alf))*( (Dq + matmul(transpose(Bqn),eta) - v)/h - 0.5*a)  +  (a - alf*vd)/(1-alf)

            !v1  = v1  + gam/bet*( - matmul(transpose(Bqn),eta))
            !vd1 = vd1 + (1-alm)/(bet*(1-alf)*h)*(- matmul(transpose(Bqn),eta))

            !Dq  = Dq  + Dxleta(1:sv)
            !l   = l   + Dxleta(sv   +1:sv+sl  )/h
            !eta = eta + Dxleta(sv+sl+1:sv+sl+sl)

            !v1  = v1  + gam/bet*(Dxleta(1:sv) + matmul(transpose(B),eta))
            !vd1 = vd1 + (1-alm)/(bet*(1-alf)*h)*(Dxleta(1:sv) + matmul(transpose(B),eta))

            !Bqn = B


            ! update solver stats
            this%gena_stats%newt_steps_curr = i

            ! exit loop, if converged
            if ( converged ) exit

         end do ! end newton method

#ifdef DEBUG1
         ! DEBUG
         if (t>=MY_T) then
            i = this%gena_stats%newt_steps_curr+1
            write(*,'(A)',advance='NO') 'Dq'
#ifdef ZEROINIT
            write(*,'(A)',advance='NO') 'zero'
#else
            write(*,'(A)',advance='NO') 'ours'
#endif
            write(*,'(A,I5,A)',advance='NO') '(', i, ',:) ...'
            call print_vector(Dq,'')

            write(*,'(A)',advance='NO') 'hl'
#ifdef ZEROINIT
            write(*,'(A)',advance='NO') 'zero'
#else
            write(*,'(A)',advance='NO') 'ours'
#endif
            write(*,'(A,I5,A)',advance='NO') '(', i, ',:) ...'
            call print_vector(h*l,'')

            write(*,'(A)',advance='NO') 'eta'
#ifdef ZEROINIT
            write(*,'(A)',advance='NO') 'zero'
#else
            write(*,'(A)',advance='NO') 'ours'
#endif
            write(*,'(A,I5,A)',advance='NO') '(', i, ',:) ...'
            call print_vector(eta,'')

            write(*,'(A)',advance='NO') 'res'
#ifdef ZEROINIT
            write(*,'(A)',advance='NO') 'zero'
#else
            write(*,'(A)',advance='NO') 'ours'
#endif
            write(*,'(A)',advance='NO') '...'
            call print_vector(res,'')

         end if
         if (t >= MY_T) STOP
         ! GUBED
#endif

         ! DEBUG
         if (this%gena_stats%newt_steps_curr >= this%opts%imax) then
            print *, ' ' ! DEBUG
            print *, "Newton did not converge at t=", t
            !! DEBUG !
            call this%gena_print_stats
            print *, "Iteration matrix:"
            !call print_matrix(St,'St')
            errorstop
            !! GUBED !
         end if

         ! The integration step is done; now the new values have to be
         ! applied to the problem object, and the new value for a can
         ! be calculated
         a  = ((1-alf)*vd1 + alf*vd - alm*a)/(1-alm)
         t  = t1
         q  = q1
         v  = v1
         vd = vd1

      end associate

   end subroutine gena_solveConstrainedTimeStep_stab2

   ! subroutine for integrating one time step
   subroutine gena_solveTimeStep(this,t1)
      implicit none
      class(gena_problem), intent(inout)        :: this  ! problem object
      real(8)            , intent(in   )        :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                                   :: i     ! for iteration
      integer, dimension(this%sizev)            :: ipiv  ! needed for dgesv
      integer                                   :: info  ! needed for dgesv

      ! intenal real variables
      real(8)                                   :: h       ! step size
      real(8), dimension(this%sizev)            :: vd1     ! $\dot v_{n+1}$
      !real(8), dimension(this%sizev)            :: a1      ! $a_{n+1}$ not needed 2015-06-23
      real(8), dimension(this%sizev)            :: v1      ! $v_{n+1}$
      real(8), dimension(this%sizev)            :: Dq      ! $\Delta q_n$
      real(8), dimension(this%sizeq)            :: q1      ! $q_{n+1}}
      real(8), dimension(this%sizev)            :: res     ! $res$
      real(8), dimension(this%sizev,this%sizev) :: St      ! $S_t$
      real(8), dimension(:,:), allocatable      :: Stb     ! $S_t$ in banded format
      !real(8)                                   :: beta_d  ! $\beta'$  not needed 2015-06-23
      !real(8)                                   :: gamma_d ! $\gamma'$ not needed 2015-06-23
      real(8), dimension(this%sizev)            :: Dx      ! $\Delta x$
#ifdef ITERATIVE_RELAXATION
      real(8)                                   :: relgam
      real(8)                                   :: relalpha
      real(8)                                   :: resrel(this%sizev)
      real(8)                                   :: Dqrel(this%sizev)
      real(8)                                   :: qrel1(this%sizeq)
      real(8)                                   :: vrel1(this%sizev)
      real(8)                                   :: vdrel1(this%sizev)
      real(8)                                   :: nDx, nDxold
#endif

      ! internal logical
      logical                                   :: converged


#ifdef ITERATIVE_RELAXATION
      relgam = 1.0_8
#endif

      ! associate construct for better readability
      associate (alf => this%gena_alpha_f,    &
                 alm => this%gena_alpha_m,    &
                 bet => this%gena_beta,       &
                 gam => this%gena_gamma,      &
                 sv  => this%sizev,           &
                 v   => this%v,               &
                 vd  => this%vd,              &
                 a   => this%a,               &
                 q   => this%q,               &
                 t   => this%t)
         ! calculation of step size $h$
         h = t1 - t

#ifdef ZEROINIT
         Dq = 0.0_8
         v1 = (1 - gam/bet)*v + h*(1 - gam/(2*bet))*a
         vd1 = (1-alm)/(bet*(1-alf))*( (Dq - v)/h - 0.5*a)  +  (a - alf*vd)/(1-alf)
#else
         ! initialization of the new values $\dot v_{n+1}$, $a_{n+1}$ and $v_{n+1}$ and the value $\Delta q_n$
         !vd1 = 0.0_8 ! vd !0.0_8 ! DEBUG TODO: hier vd? Blöde Idee, aber warum???
         !a1  = (alf*vd - alm*a)/(1.0_8 - alm)
         !v1  = v + h*(1.0_8 - gam)*a + gam*h*a1
         !Dq  = v + (0.5_8 - bet)*h*a + bet*h*a1
         Dq  = v + 0.5*h*a
         v1  = v + h*a
         vd1 = (a - alf*vd)/(1-alf)
#endif

         ! loop of the newton method
         converged = .false.
         do i=1,this%opts%imax
            ! calculate the new value $q_{n+1}$
            q1  = this%gena_qlpexphDqtilde(q,h,Dq)
            ! caluclate the residue $res$
            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  res = this%gena_const_diag_M*vd1 + this%gena_g(q1,v1,t1)
               else
                  res = matmul(this%gena_const_M,vd1) + this%gena_g(q1,v1,t1)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  res = this%gena_diag_M(q1)*vd1 + this%gena_g(q1,v1,t1)
               else
                  res = matmul(this%gena_M(q1),vd1) + this%gena_g(q1,v1,t1)
               end if
            end if
            ! scale
            res = h*res
            ! count calls
            this%gena_stats%ngcalls = this%gena_stats%ngcalls + 1

            !! DEBUG
            !print *, this%gena_norm(res)

            !! check if norm of the residue is sufficiently small to stop the iteration
            !if (this%gena_norm(res) < this%opts%tolr) then
            !   exit
            !end if

            Dx = -res ! the result will be saved in the right hand side's spot

            if (i == 1  .or. this%opts%recalc_iteration_matrix == 1) then
               ! calculate iteration matrix
               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     St = 0.0_8
                     forall (i=1:sv)
                        St(i,i) = this%gena_const_diag_M(i) * (1-alm)/(bet*(1-alf))
                     end forall
                  else
                     St = this%gena_const_M * (1-alm)/(bet*(1-alf))
                  end if
               else
                  if (this%opts%diag_mass_matrix == 1) then
                     St = 0.0_8
                     St(:,1) = this%gena_diag_M(q1) * (1-alm)/(bet*(1-alf))
                     forall (i=2:sv)
                        St(i,i) = St(i,1)
                        St(i,1) = 0.0_8
                     end forall
                  else
                     St = this%gena_M(q1) * (1-alm)/(bet*(1-alf))
                  end if
               end if

               if (this%opts%no_Ct == 0) then
                  if (this%opts%use_num_Ct == 1) then
                     St = St + this%gena_num_Ct(q1,v1,t1) * h*gam/bet
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + this%sizev + 1
                     end if
                  else
                     St = St + this%gena_Ct(q1,v1,t1) * h*gam/bet
                  end if
               end if

               if (this%opts%no_Kt == 0) then
                  if (this%opts%use_num_Kt == 1) then
                     St = St + matmul(this%gena_num_Kt(q1,v1,vd1,t1),this%gena_Tg(h,Dq)) * h*h! TODO: dgemm oder sowas: Idee: Man kann per compilerflags einschalten, dass für matmul direkt eine blas-routine benutzt wird
                     ! count calls
#ifdef NOBANDEDDER
                     if (.false.) then
#else
                     if (this%opts%banded_iteration_matrix == 1) then
#endif
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + &
                           this%opts%nr_subdiag + this%opts%nr_superdiag + 1
                     else
                        this%gena_stats%ngcalls = this%gena_stats%ngcalls + this%sizev + 1
                     end if
                  else
                     St = St + matmul(this%gena_Kt(q1,v1,vd1,t1),this%gena_Tg(h,Dq)) * h*h ! TODO: dgemm oder sowas
                  end if
               end if

               !! DEBUG !
               !print *, 't=', t
               !call print_vector(q,'q')
               !call print_vector(v,'v')
               !call print_vector(vd,'vd')
               !call print_vector(a,'a')
               !!call print_matrix(St,'Stn')
               !call print_vector(res, 'res')
               !!call print_vector(this%gena_const_diag_M,'Mn')
               !!print *, "t=", t
               !!if (t > 2.0_8) then
               !!   call print_matrix(St,'St')
               !!   print *, "cond=", mycond(St)
               !!   exit
               !!endif
               !!print *,
               !!stop
               !! GUBED !

               ! solve the linear System $S_t \cdot \Delta x = -res$
               if (this%opts%banded_iteration_matrix == 1) then
#ifdef DEBUG_PRINT_ITERATION_MATRIX_AT
                  if (t >= DEBUG_PRINT_ITERATION_MATRIX_AT) then
                     call print_matrix(St,'St')
                     call print_vector_int(this%opts%jour,'jour')
                     stop "Successfully printed iteration matrix"
                  end if
#endif
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag, &
                             uStb => ubound(Stb,1))
                     if (.not. allocated(Stb)) then
                        allocate(Stb(lo*2+hi+1,this%sizev))
                     end if
                     forall (i=1:this%sizev)
                        Stb((max(1, i-hi)-i+lo+1+hi):(min(this%sizev,i+lo)-i+lo+1+hi),i) &
                            = St(max(1, i-hi):min(this%sizev, i + lo),i)
                     end forall
                     ! TODO: Conditional jump or move depends on unititialised value(s) in the next line?
                     call dgbsv(       &! solve the system A*X=B and save the result in B
                                 sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,   &!
                                 hi,   &!
                                 1,    &! number of right hand sides (=size(B,2))
                                 Stb,  &! matrix A
                                 2*lo+hi+1, &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv, &! integer pivot vector; it is not needed
                                 Dx,   &! matrix B
                                 sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)  ! integer information flag
                  end associate
               else
                  call dgesv(       &! solve the System A*X=B and save the result in B
                              sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,    &! number of right hand sides (=size(B,2))
                              St,   &! matrix A
                              sv,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv, &! integer pivot vector; it is not needed
                              Dx,   &! matrix B
                              sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)  ! integer information flag
               end if
               ! TODO: Errors abfragen
               !!! DEBUG
               !   print *, 'St = ...'
               !   print *, '[', St(1,1), ',', St(1,2), ';'
               !   print *, ' ', St(2,1), ',', St(2,2), '];'
               !   print *, 'mres = ...'
               !   print *, '[', -res(1), ',', -res(2), '];'
               !   print *, 'Dx = ...'
               !   print *, '[', Dx(1), ',', Dx(2), '];'
               !   print *, ''
               if (info.ne.0) print *, 'dgesv sagt info=', info
               !! GEBUD
            else
               ! solve the linear System $S_t \cdot \Delta x = -res$
               if (this%opts%banded_iteration_matrix == 1) then
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag  )
                     call dgbtrs(      &
                                 'No transpose', &! solve the System A*X=B and save the result in B
                                 sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,   &!
                                 hi,   &!
                                 1,    &! number of right hand sides (=size(B,2))
                                 Stb,  &! matrix A
                                 2*lo+hi+1, &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv, &! integer pivot vector; it is not needed
                                 Dx,   &! matrix B
                                 sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)  ! integer information flag
                  end associate
               else
                  call dgetrs(      &
                              'No transpose', &! solve the System A*X=B and save the result in B
                              sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,    &! number of right hand sides (=size(B,2))
                              St,   &! matrix A
                              sv,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv, &!
                              Dx,   &! matrix B
                              sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)  ! integer information flag
               end if
               ! DEBUG
               if (info.ne.0) print *, 'dgesv sagt info=', info
               ! GEBUD
            end if

#ifdef ITERATIVE_RELAXATION

            do while (.true.)
               Dqrel  = Dq  + relgam * Dx
               vrel1  = v1  + relgam * gam/bet * Dx
               vdrel1 = vd1 + relgam * (1-alm)/(bet*(1-alf)*h)*Dx
               qrel1  = this%gena_qlpexphDqtilde(q,h,Dqrel)

               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     resrel = this%gena_const_diag_M*vdrel1 + this%gena_g(qrel1,vrel1,t1)
                  else
                     resrel = matmul(this%gena_const_M,vdrel1) + this%gena_g(qrel1,vrel1,t1)
                  end if
               else
                  if (this%opts%diag_mass_matrix == 1) then
                     resrel = this%gena_diag_M(qrel1)*vdrel1 + this%gena_g(qrel1,vrel1,t1)
                  else
                     resrel = matmul(this%gena_M(qrel1),vdrel1) + this%gena_g(qrel1,vrel1,t1)
                  end if
               end if
               ! scale
               resrel = h*resrel

               if (norm2(resrel) < (1.0_8 - 1.0e-4_8*relgam)*norm2(res)) exit
               if (relgam < (2.0_8**(-30))) then
                  print *, "Could not find relaxation"
                  !pause
                  Dqrel  = Dq  + Dx
                  vrel1  = v1  + gam/bet * Dx
                  vdrel1 = vd1 + (1-alm)/(bet*(1-alf)*h)*Dx
                  qrel1  = this%gena_qlpexphDqtilde(q,h,Dq)
                  exit
               end if

               relgam = relgam/2
            end do

            ! DEBUG
            print *, ''
            print *, 't      = ', t
            print *, 'i      = ', i
            print *, 'relgam = ', nint(log(relgam)/log(2.0_8))
            print *, 'Dx     = ', norm2(Dx)
            print *, 'resrel = ', norm2(resrel)
            relalpha = norm2(Dx) / nDxold
            print *, 'relalp = ', relalpha
            print *, 'err    = ', sum((Dx / (this%opts%atol + this%opts%rtol*abs(Dq)))**2 ) / (sv)
            print *, 'scaler = ', relalpha/(1-relalpha)*sum((Dx / (this%opts%atol + this%opts%rtol*abs(Dq)))**2 ) / (sv)
            ! GUBED

            Dq = Dqrel
            v1 = vrel1
            vd1 = vdrel1

            nDx = norm2(Dx)

            ! Check for convergence
            if (i == 1) then
               ! Check for convergence
               if (sum((Dx / (this%opts%atol + this%opts%rtol*abs(Dq)))**2 )  &
                      / (sv)   <=  1.0_8 ) then
                  converged = .true.
                  print *, '################################################'
               end if
            else
               relalpha = nDx / nDxold
               if (relalpha < 1.0_8) then
                  if ( relalpha/(1-relalpha) *&
                          sum((Dx / (this%opts%atol + this%opts%rtol*abs(Dq)))**2 )  &
                         / (sv)   <=  1.0_8 ) then
                     converged = .true.
                     print *, '################################################'
                  end if
               end if
            end if

            relgam = 1.0_8 !min(1.0_8, 4*relgam)

            nDxold = nDx

#endif

            !! DEBUG
            !call print_vector(Dx,'Dxn')
            !stop
            !! GUBED



#ifndef ITERATIVE_RELAXATION
            ! Check for convergence
            if ( sum((Dx / (this%opts%atol + this%opts%rtol*abs(Dq)))**2 )  &
                  / (sv)   <=  1.0_8 ) then
               converged = .true.
            end if

! DEBUG!!
#define REL *1
            ! update $\dot v_{n+1}$, $v_{n+1}$ and $\Delta q_n$
            Dq  = Dq  + Dx                          REL
            v1  = v1  + gam/bet * Dx                REL
            vd1 = vd1 + (1-alm)/(bet*(1-alf)*h)*Dx  REL
#undef REL
! GUBED!!
#endif

            ! update solver stats
            this%gena_stats%newt_steps_curr = i

            ! exit loop, if converged
            if ( converged ) exit

         end do ! end newton method

         !! DEBUG
         if (this%gena_stats%newt_steps_curr >= this%opts%imax) then
            print *, ' ' ! DEBUG
            print *, "Newton did not converge at t=", t
            !!! DEBUG !
            !call this%gena_print_stats
            !print *, "Iteration matrix:"
            !call print_matrix(St,'St')
            errorstop
            !!! GUBED !
         end if


         ! The integration step is done; now the new values have to be
         ! applied to the problem object, and the new value for a can
         ! be calculated
         a  = ((1-alf)*vd1 + alf*vd - alm*a) / (1-alm)
         t  = t1
         q  = q1
         v  = v1
         vd = vd1


      end associate

   end subroutine gena_solveTimeStep

   !subroutine gena_trap_solveTimeStep(this, h)
   !
   !end subroutine gena_trap_solveTimeStep

   ! subroutine to integrate
   subroutine gena_integrate(this)
      implicit none
      class(gena_problem), intent(inout)  :: this

      integer  :: n  ! needed for iteration

      real(8)                             :: h  ! step size $h$
      real(8)                             :: t1 ! next time $t_{n+1}$
      !real(8)                             :: alf_orig ! needed in order to calculate a 2nd order
      !real(8)                             :: alm_orig ! approximation for a at the beginning
      !real(8)                             :: bet_orig ! of the algorithm (in the constrained
      !real(8)                             :: gam_orig ! case).
      !real(8), dimension(:), allocatable  :: q_orig
      !real(8), dimension(:), allocatable  :: v_orig
      !real(8), dimension(:), allocatable  :: vd_orig
      real(4)                             :: times(2), time ! for dtime


      ! problem initialization
      call this%gena_init()

      ! TODO: More error checking
      if (this%opts%constrained == 1             .and. &
          this%opts%banded_iteration_matrix == 1 .and. &
          .not. allocated(this%opts%jour)               ) then
         print *, '[ERROR] opts%jour not allocated, aborting'
         errorstop
      end if
      ! TODO: t0 und t vergleichen (müssen gleich sein)

      ! initialize output function
      call this%gena_outputFunction(0)

      ! Calculate step size $h$
      h = (this%opts%te - this%opts%t0)/this%opts%nsteps

      ! Set stats of solver to zero
      this%gena_stats%newt_steps_curr = 0
      this%gena_stats%newt_steps_sum = 0
      this%gena_stats%newt_steps_max = 0
      this%gena_stats%newt_steps_avg = 0
      this%gena_stats%ngcalls = 0
      this%gena_stats%nBcalls = 0

      ! if mass matrix is constant, calculate it
      if (this%opts%const_mass_matrix == 1) then
         if (this%opts%diag_mass_matrix == 1) then
            if (.not. allocated(this%gena_const_diag_M)) then
               allocate(this%gena_const_diag_M(this%sizev))
            end if
            this%gena_const_diag_M = this%gena_diag_M(this%q)
         else
            if (.not. allocated(this%gena_const_M)) then
               allocate(this%gena_const_M(this%sizev,this%sizev))
            end if
            this%gena_const_M = this%gena_M(this%q)
         end if
      end if

      ! start stopwatch
      time = dtime(times)

      if (this%opts%constrained == 0) then

         ! calculate initial values
         call this%gena_calcInitial(h)

         ! output for the first time
         call this%gena_outputFunction(1)

         ! integration loop
         do n=1,this%opts%nsteps
            ! Calculate the next time $t_{n+1}$
            t1 = this%opts%t0 + n*h
            ! solve time step
            call this%gena_solveTimeStep(t1)
            ! update solver stats
            this%gena_stats%newt_steps_sum = this%gena_stats%newt_steps_sum &
               + this%gena_stats%newt_steps_curr
            this%gena_stats%newt_steps_avg = real(this%gena_stats%newt_steps_sum,8)/n
            !this%gena_stats%newt_steps_avg =           &
            !    (this%gena_stats%newt_steps_avg*(n-1) + &
            !    this%gena_stats%newt_steps_curr)/n
            this%gena_stats%newt_steps_max = max( &
               this%gena_stats%newt_steps_max,    &
               this%gena_stats%newt_steps_curr  )
            ! output normally
            call this%gena_outputFunction(1)
         end do

      else
!! commented out on ! 2014-01-31
!!         ! calculate a 2nd order approximation of a
!!         alf_orig = this%gena_alpha_f
!!         alm_orig = this%gena_alpha_m
!!         bet_orig = this%gena_beta
!!         gam_orig = this%gena_gamma
!!         q_orig   = this%q
!!         v_orig   = this%v
!!         vd_orig  = this%vd
!!         this%gena_alpha_f = (alf_orig+alm_orig)/2.0_8 ! TODO: Besserer Wert?
!!         this%gena_alpha_m = this%gena_alpha_f
!!         this%gena_gamma   = 1.0_8/2.0_8
!!         this%gena_beta    = 1.0_8/2.0_8 ! TODO: Besserer Wert?
!!         this%a = vd_orig
!!         call this%gena_solveConstrainedTimeStep( &
!!                  this%opts%t0+(alm_orig-alf_orig)*h)
!!         this%gena_alpha_f = alf_orig
!!         this%gena_alpha_m = alm_orig
!!         this%gena_beta    = bet_orig
!!         this%gena_gamma   = gam_orig
!!         this%q            = q_orig
!!         this%v            = v_orig
!!         this%a            = this%vd
!!         this%vd           = vd_orig
!!         this%t            = this%opts%t0

         ! calculate correct initial values (and pertube them, if wanted)
         call this%gena_calcInitialConstrained(h)

         ! output for the first time
         call this%gena_outputFunction(1)

         ! integration loop
         do n=1,this%opts%nsteps
            ! Calculate the next time $t_{n+1}$
            t1 = this%opts%t0 + n*h

            !! DEBUG
            !if (this%t > 2.67_8) then
            !   call print_vector(this%q,'qn')
            !   call print_matrix(this%gena_B(this%q),'Be')
            !   call print_matrix(gena_num_B(this,this%q),'Bd')
            !   stop
            !end if
            !! GUBED

            !! DEBUG
            !if (n > 500) then
            !   exit
            !end if
            !! GUBED

            if (this%opts%stab2 == 1) then
               ! solve time step with the stabilized index-2 system
               call this%gena_solveConstrainedTimeStep_stab2(t1)
            else
               ! solve time step
               call this%gena_solveConstrainedTimeStep(t1)
            end if

            ! update solver stats
            this%gena_stats%newt_steps_avg =           &
               (this%gena_stats%newt_steps_avg*(n-1) + &
                this%gena_stats%newt_steps_curr)/n
            this%gena_stats%newt_steps_max = max( &
               this%gena_stats%newt_steps_max,    &
               this%gena_stats%newt_steps_curr  )
            ! output normally
            call this%gena_outputFunction(1)
         end do
      end if

      ! stop stopwatch
      this%gena_stats%time = dtime(times) - time

      ! output to terminate
      call this%gena_outputFunction(99)

   end subroutine gena_integrate

   pure function gena_num_Ct(this, q, v, t) result(rslt)
      ! input
      class(gena_problem),   intent(in)         :: this
      real(8), dimension(:), intent(in)         :: q
      real(8), dimension(:), intent(in)         :: v
      real(8),               intent(in)         :: t
      ! result
      real(8), dimension(this%sizev,this%sizev) :: rslt
      ! internal
      integer                                   :: i
      integer                                   :: j
      integer                                   :: diags
      real(8), dimension(this%sizev)            :: w
      real(8), dimension(this%sizev)            :: g0
      !
      ! Ct is the derivative of r wrt v. But v only appears in g, so
      ! it is sufficient to calculate the derivative of g.
      !
      ! Calculate g ! TODO: man könnte g0 durchschleifen oder so
      g0 = this%gena_g(q,v,t)
#ifdef NOBANDEDDER
      if (.false.) then
#else
      ! if the iteration matrix is banded, Ct must be banded, as well
      if (this%opts%banded_iteration_matrix == 1) then
#endif
         ! the number of total off-diagonals
         diags = this%opts%nr_subdiag+this%opts%nr_superdiag
         ! most of the result will be zero
         rslt = 0.0_8
         ! loop over the first diags columns of Ct
         do i=1,min(diags+1,this%sizev)
            ! Construct vector w
            w = 0.0_8
            do j=i, this%sizev, 1+diags
               w(j) = max( abs(v(j))*1.0e-8_8,  1.0e-12_8)
            end do
            ! calculate finite differences
            rslt(:,i) = this%gena_g(q,v+w,t) - g0
            ! move the parts of the finite difference, that belog to
            ! an other column of Ct, also divide by h=w(k)
            do j=i, this%sizev, 1+diags
               rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), j) =      &
                     rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i)  &
                     / w(j)
               if (i .ne. j) then
                  rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i) = 0.0_8
               end if
            end do
         end do
      else
         ! set w to zero
         w = 0.0_8
         ! loop over the columns of Ct
         do i=1,this%sizev
            if (.not. i == 1) w(i-1) = 0.0_8
            w(i) = max( abs(v(i))*1.0e-8_8,  1.0e-12_8)
            rslt(:,i) = (this%gena_g(q,v+w,t) - g0)/w(i)
         end do
      end if
   end function gena_num_Ct

   pure function gena_num_Kt(this, q, v, vd, t) result(rslt)
      ! input
      class(gena_problem),   intent(in)         :: this
      real(8), dimension(:), intent(in)         :: q
      real(8), dimension(:), intent(in)         :: v
      real(8), dimension(:), intent(in)         :: vd
      real(8),               intent(in)         :: t
      ! result
      real(8), dimension(this%sizev,this%sizev) :: rslt
      ! internal
      integer                                   :: i
      integer                                   :: j
      integer                                   :: diags
      real(8), dimension(this%sizev)            :: w
      real(8), dimension(this%sizev)            :: g0
      !
      ! Kt is the derivative of r wrt q. WE ASSUME THAT THE MASS MATRIX
      ! M DOES NOT DEPEND ON Q, OTHERWISE THIS APPROXIMATION WILL BE
      ! FALSE! TODO! TODO!
      ! Thus, it is sufficient to calculate the derivative of g.
      !
      ! Calculate g ! TODO: man könnte g0 durchschleifen oder so
      g0 = this%gena_g(q,v,t)
#ifdef NOBANDEDDER
      if (.false.) then
#else
      ! if the iteration matrix is banded, Kt must be banded, as well
      if (this%opts%banded_iteration_matrix == 1) then
#endif
         ! most of Kt is zero
         rslt = 0.0_8
         ! the number of total off-diagonals
         diags = this%opts%nr_subdiag+this%opts%nr_superdiag
         ! loop over the first diags columns of Kt
         do i=1,min(diags+1,this%sizev)
            ! Construct vector w
            w = 0.0_8
            do j=i, this%sizev, 1+diags
               w(j) = max( abs(q(j))*1.0e-8_8,  1.0e-12_8)
            end do
            ! calculate finite differences
            rslt(:,i) = this%gena_g(this%gena_qlpexphDqtilde(q,1.0_8, w),v,t) - g0
            ! move the parts of the finite difference, that belong to
            ! an other column of Kt, also divide by h=w(k)
            do j=i, this%sizev, 1+diags
               rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), j) =     &
                     rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i) &
                       / w(j)
               if (i .ne. j) then
                  rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i) = 0.0_8
               end if
            end do
         end do
      else
         ! set w to zero
         w = 0.0_8
         ! loop over the columns of Ct
         do i=1,this%sizev
            if (.not. i == 1) w(i-1) = 0.0_8
            w(i) = max( abs(q(i))*1.0e-8_8,  1.0e-12_8)
            rslt(:,i) = (this%gena_g(this%gena_qlpexphDqtilde(q,1.0_8, w),v,t) - g0) / w(i)
            !w    = 0.0_8
            !w(i) = h
            !rslt(:,i) = (this%gena_g(q+this%gena_tilde(w),v,t) - g0)/h
         end do
      end if
   end function gena_num_Kt

! DEBUG
   pure function gena_num_B(this, q) result(rslt)
      ! input
      class(gena_problem),   intent(in)         :: this
      real(8), dimension(:), intent(in)         :: q
      ! result
      real(8), dimension(this%sizel,this%sizev) :: rslt
      ! internal
      integer                                   :: i
      real(8), dimension(this%sizev)            :: w
      real(8), dimension(this%sizev)            :: Phi0
      !
      !
      ! Calculate Phi
      Phi0 = this%gena_Phi(q)

      ! set w to zero
      w = 0.0_8
      ! loop over the columns of Ct
      do i=1,this%sizev
         if (.not. i == 1) w(i-1) = 0.0_8
         w(i) = max( abs(q(i))*1.0e-8_8,  1.0e-12_8)
         rslt(:,i) = (this%gena_Phi(this%gena_qlpexphDqtilde(q,1.0_8, w)) - Phi0) / w(i)
         !w    = 0.0_8
         !w(i) = h
         !rslt(:,i) = (this%gena_g(q+this%gena_tilde(w),v,t) - g0)/h
      end do
   end function gena_num_B
! GUBED

   pure function gena_num_Kt_lambda(this, q, v, vd, l, t) result(rslt)
      ! input
      class(gena_problem),            intent(in)   :: this
      real(8), dimension(:),          intent(in)   :: q
      real(8), dimension(:),          intent(in)   :: v
      real(8), dimension(:),          intent(in)   :: vd
      real(8), dimension(:),          intent(in)   :: l
      real(8),                        intent(in)   :: t
      ! result
      real(8), dimension(this%sizev,this%sizev)    :: rslt
      ! internal
      integer                                      :: i
      integer                                      :: j
      integer                                      :: diags
      real(8), dimension(this%sizev)               :: w
      real(8), dimension(this%sizev)               :: g0
      !
      ! Kt is the derivative of r wrt q. WE ASSUME THAT THE MASS MATRIX
      ! M DOES NOT DEPEND ON Q, OTHERWISE THIS APPROXIMATION WILL BE
      ! FALSE! TODO! TODO!
      ! Thus, it is sufficient to calculate the derivative of g.
      !
      ! Calculate g ! TODO: man könnte g0 durchschleifen oder so
      g0 = this%gena_g(q,v,t) + matmul(transpose(this%gena_B(q)),l)
#ifdef NOBANDEDDER
      if (.false.) then
#else
      ! if the iteration matrix is banded, Kt must be banded, as well
      if (this%opts%banded_iteration_matrix == 1) then
#endif
         ! most of Kt is zero
         rslt = 0.0_8
         ! the number of total off-diagonals
         diags = this%opts%nr_subdiag+this%opts%nr_superdiag
         ! loop over the first diags columns of Kt
         do i=1,min(diags+1,this%sizev)
            ! Construct vector w
            w = 0.0_8
            do j=i, this%sizev, 1+diags
               w(j) = max( abs(q(j))*1.0e-8_8,  1.0e-12_8)
            end do
            ! calculate finite differences
            rslt(:,i) = this%gena_g(this%gena_qlpexphDqtilde(q,1.0_8, w),v,t) &
               + matmul(transpose(this%gena_B(this%gena_qlpexphDqtilde(q,1.0_8,w))),l)- g0
            ! move the parts of the finite difference, that belog to
            ! an other column of Kt, also divide by h=w(k)
            do j=i, this%sizev, 1+diags
               rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), j) =     &
                     rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i) &
                     / w(j)
               if (i .ne. j) then
                  rslt(max(1, j-this%opts%nr_superdiag):min(j+this%opts%nr_subdiag, this%sizev), i) = 0.0_8
               end if
            end do
         end do
      else
         ! set w to zero
         w = 0.0_8
         ! loop over the columns of Ct
         do i=1,this%sizev
            if (.not. i == 1)  w(i-1) = 0.0_8
            w(i) = max( abs(q(i))*1.0e-8_8,  1.0e-12_8)
            !rslt(:,i) = (this%gena_g(q+this%gena_tilde(w),v,t) + matmul(transpose(this%gena_B(q+this%gena_tilde(w))),l) - g0)/h
            !w(i) = 1.0_8 ! TODO
            rslt(:,i) = (this%gena_g(this%gena_qlpexphDqtilde(q,1.0_8,w),v,t)  &
             + matmul(transpose(this%gena_B(this%gena_qlpexphDqtilde(q,1.0_8,w))),l) - g0)/w(i)
         end do
      end if
   end function gena_num_Kt_lambda

   subroutine gena_print_stats(this)
      ! input
      class(gena_problem), intent(in)  :: this
      !
      print *, 'time:          ', this%gena_stats%time
      print *, '#calls of g:   ', this%gena_stats%ngcalls
      print *, '#calls of B:   ', this%gena_stats%nBcalls
      print *, 'newt_steps_max:', this%gena_stats%newt_steps_max
      print *, 'newt_steps_avg:', this%gena_stats%newt_steps_avg
   end subroutine

   subroutine print_matrix(A,Aname)
      implicit none
      ! input
      real(8), dimension(:,:) :: A
      character(len=*)       :: Aname
      ! internal
      integer                 :: i
      !
      print *, '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print *, trim(Aname), ' = ['
      do i=1,ubound(A,1) - 1
         print *, A(i,:), ';'
      end do
      print *, A(ubound(A,1),:), '];'
   end subroutine print_matrix

   subroutine print_vector_int(A,Aname)
      implicit none
      ! input
      integer, dimension(:)   :: A
      character(len=*)        :: Aname
      !
      print *, '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print *, trim(Aname), ' = ['
      print *, A, '];'
   end subroutine print_vector_int

   subroutine print_vector(A,Aname)
      implicit none
      ! input
      real(8), dimension(:)   :: A
      character(len=*)        :: Aname
      !
      print *, '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print *, trim(Aname), ' = ['
      print *, A, '];'
   end subroutine print_vector

   subroutine gena_cleanup(this)
      implicit none
      ! input/output
      class(gena_problem), intent(inout)  :: this
      !
      this%opts%constrained = 0
      this%opts%const_mass_matrix = 0
      this%opts%diag_mass_matrix = 0
      this%opts%banded_iteration_matrix = 0
      this%opts%nr_subdiag = 0
      this%opts%nr_superdiag = 0
      this%opts%recalc_iteration_matrix = 0
      this%opts%use_num_Ct = 1
      this%opts%use_num_Kt = 1
      this%opts%atol = 1.0e-10_8
      this%opts%rtol = 1.0e-8_8
      this%opts%imax   = 5
      this%opts%t0 = 0.0_8
      this%opts%te = 1.0_8
      this%opts%nsteps = 100
      !
      this%gena_stats%newt_steps_curr = 0
      this%gena_stats%newt_steps_max = 0
      this%gena_stats%newt_steps_avg = 0.0_8
      this%gena_stats%ngcalls = 0
      this%gena_stats%nBcalls = 0
      this%gena_stats%time = 0.0_8
      !
      this%t = 0.0_8
      this%sizeq = 0
      this%sizev = 0
      this%sizel = 0
      !
      if (allocated(this%q))  deallocate(this%q)
      if (allocated(this%v))  deallocate(this%v)
      if (allocated(this%vd)) deallocate(this%vd)
      if (allocated(this%a))  deallocate(this%a)
      if (allocated(this%l))  deallocate(this%l)
      if (allocated(this%eta)) deallocate(this%eta)
      if (allocated(this%gena_const_M)) deallocate(this%gena_const_M)
      if (allocated(this%gena_const_diag_M)) deallocate(this%gena_const_diag_M)
      if (allocated(this%opts%jour)) deallocate(this%opts%jour)
   end subroutine gena_cleanup

#ifdef extra_integrate
   subroutine gena_extra_integrate(this)
      implicit none
      class(gena_problem), intent(inout)  :: this

      integer  :: n  ! needed for iteration

      real(8)                             :: h, hold  ! step size $h$
      real(8)                             :: t0
      real(8)                             :: tol ! Error tolerance
      real(8), allocatable, dimension(:)  :: thatq
      real(8), allocatable, dimension(:)  :: thatv
      real(8), allocatable, dimension(:)  :: thata
      real(8), allocatable, dimension(:)  :: savedq
      real(8), allocatable, dimension(:)  :: savedv
      real(8), allocatable, dimension(:)  :: saveda
      real(8), allocatable, dimension(:)  :: olda
      real(8), allocatable, dimension(:)  :: errq
      real(8), allocatable, dimension(:)  :: errv
      real(8)                             :: errind


      print *, "bin da"

      ! problem initialization
      call this%gena_init()

      ! Allocate Variable
      allocate(thatq (this%sizeq))
      allocate(thatv (this%sizev))
      allocate(thata (this%sizev))
      allocate(savedq(this%sizeq))
      allocate(savedv(this%sizev))
      allocate(saveda(this%sizev))
      allocate(olda  (this%sizev))
      allocate(errq  (this%sizeq))
      allocate(errv  (this%sizev))

      ! TODO: t0 und t vergleichen (müssen gleich sein)

      ! output for the first time
      call this%gena_outputFunction(0)

      ! output the initial conditions
      call this%gena_outputFunction(1)

      ! TODO: Support for Stats

      ! if mass matrix is constant, calculate it
      if (this%opts%const_mass_matrix == 1) then
         if (this%opts%diag_mass_matrix == 1) then
            if (.not. allocated(this%gena_const_diag_M)) then
               allocate(this%gena_const_diag_M(this%sizev))
            end if
            this%gena_const_diag_M = this%gena_diag_M(this%q)
         else
            if (.not. allocated(this%gena_const_M)) then
               allocate(this%gena_const_M(this%sizev,this%sizev))
            end if
            this%gena_const_M = this%gena_M(this%q)
         end if
      end if

      ! Use 10**(-this%opts%nsteps) as tolerance TODO
      tol = 10.0_8**(-this%opts%nsteps)

      ! Set initial time
      t0 = this%t

      ! Calculate initial step size TODO
      h = 10.0_8**(-this%opts%nsteps/2.0_8)

      olda = this%a

      if (this%opts%constrained == 0) then

         n = 0

         ! Start integration loop
         do while (.true.)

            n = n+1

            ! Save current configuration for possible recalculation
            savedq = this%q
            savedv = this%v
            saveda = this%a

            ! calculate appropriate a
            this%a = this%a + (this%gena_alpha_f - this%gena_alpha_m)*0.5_8 &
                                *(this%a - olda)

            ! Integrate with two half-sized steps
            call this%gena_solveTimeStep(t0 + h/2.0_8)
            call this%gena_solveTimeStep(t0 + h)

            ! Save results in that
            thatq = this%q
            thatv = this%v

            ! Set this back
            this%t = t0
            this%q = savedq
            this%v = savedv
            this%a = saveda

            ! Integrate with one step
            call this%gena_solveTimeStep(t0 + h)

            ! Richardson-Extrapolation for the velocity
            this%v = (4.0_8*thatv - this%v)/3.0_8
            errv = this%v - thatv

            ! Richardson-Extrapolation for the position
            this%q = this%gena_exp43logx1ix2x2(thatq,this%q)
            errq = this%q - thatq

            ! Calculate error indicator
            errind = max(maxval(abs(errq)/(3*(tol + norm2(thatq)*tol))), &
                         maxval(abs(errv)/(3*(tol + norm2(thatv)*tol))))
            print *, 't0     = ', t0
            print *, 'h      = ', h
            print *, 'errind = ', errind

            ! Assign new step size
            if (errind > 0.0_8) then
               hold = h
               h = min(max((1/errind)**(1.0_8/3.0_8)*0.9_8, 0.2_8), 1.5_8)*h ! safety factor 0.9

               ! Recalculate a
               this%a = this%a + (this%gena_alpha_f -this%gena_alpha_m) &
                                   *(h/hold - 1.0_8)*(this%a - olda)
            end if

            ! no endless loops, please
            if (n >= 100000 ) then
               exit
            end if

            ! h may not be too small
            if (h < 1e-16) then
               print *, 'h too small'
               exit
            end if


            if (errind >= 1.0_8) then
               print *, 'not ok'
               print *, ''
               this%t = t0
               this%q = savedq
               this%v = savedv
               this%q = saveda

               cycle
            end if

            print *, 'goo'
            print *, ''

            ! Renew t0
            t0 = this%t

            ! output normally
            call this%gena_outputFunction(1)

            ! end of integration?
            if (this%t > this%opts%te) then
               exit
            end if

         end do

      else
         ! TODO
      end if

      ! output to terminate
      call this%gena_outputFunction(99)

   end subroutine gena_extra_integrate
#endif

#if 1
   function mycond(A) result(rslt)
      implicit none
      ! input
      real(8), intent(in)  :: A(:,:)
      ! result
      real(8)              :: rslt
      ! internal
      real(8)              :: AA(size(A,1),size(A,2))
      real(8)              :: S(min(size(A,1), size(A,2)))
      real(8), allocatable :: work(:)
      integer              :: lwork
      real(8)              :: U(size(A,1),size(A,1))
      real(8)              :: VT(size(A,2),size(A,2))
      integer              :: iwork(8*min(size(A,1), size(A,2)))
      integer              :: info
      !
      AA = A
      lwork = 3*min(size(A,1),size(A,2)) + max(max(size(A,1),size(A,2)),7*min(size(A,1),size(A,2)))
      allocate(work(lwork))
      call DGESDD('A',       &!JOBZ, calculate singular vectors TODO: Will fail if they are not calculated
                  size(A,1), &!M,
                  size(A,2), &!N,
                  AA,        &!A,
                  size(A,1), &!LDA,
                  S,         &!S,
                  U,         &!U,
                  size(A,1), &!LDU,
                  VT,        &!VT,
                  size(A,2), &!LDVT,
                  work,      &!WORK,
                  -1,        &!LWORK,
                  iwork,     &!IWORK
                  info)      !INFO )
      if (info /= 0) then
         errorstop "In cond: dgesdd did not succeed"
      end if
      AA = A
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      call DGESDD('A',       &!JOBZ, calculate singular vectors TODO: Will fail if they are not calculated
                  size(A,1), &!M,
                  size(A,2), &!N,
                  AA,        &!A,
                  size(A,1), &!LDA,
                  S,         &!S,
                  U,         &!U,
                  size(A,1), &!LDU,
                  VT,        &!VT,
                  size(A,2), &!LDVT,
                  work,      &!WORK,
                  lwork,     &!LWORK,
                  iwork,     &!IWORK
                  info)      !INFO )
      if (info /= 0) then
         errorstop "In cond: dgesdd did not succeed"
      end if
      if (abs(S(min(size(A,1),size(A,2)))) < 1.0e-16_8) then
         errorstop "In cond: Matrix A is singular"
      end if
      rslt = S(1)/S(min(size(A,1),size(A,2)))
   end function mycond
#endif

end module gena
