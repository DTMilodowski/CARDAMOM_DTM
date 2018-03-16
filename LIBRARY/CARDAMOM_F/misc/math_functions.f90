
module math_functions

  !!!!!!!!!!!
  ! Module contains functions needed to mathematical calculations in CARDAMOM
  !!!!!!!!!!!

  implicit none

  ! assume default private
  private

  ! make explicit bits we want others to see
  public :: randn, std, idum, random_normal, random_uniform, rnstrt

  !!!!!!!!!!!
  ! Subroutines rand(), narray() and rnstrt() are from:
  !!!!!!!!!!!

  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2000-09-10  Time: 16:37:48
  ! Latest revision - 16 January 2003

  ! FORTRAN 77 version of "ran_array"
  ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
  !       including the MODifICATIONS made in the 9th printing (2002)
  ! ********* see the book for explanations and caveats! *********
  ! Author: Steve Kifowit
  ! http://ourworld.compuserve.com/homepages/steve_kifowit
  ! with modifications by Alan Miller to rnarry and rnstrt based upon
  ! Knuth's code.

  ! For Donald Knuth's Fortran 77 versions, go to:
  ! http://www-cs-faculty.stanford.edu/~knuth/programs
  ! Look for frng.f and frngdb.f

  ! NOTE: that minimum number of values to be returned is 100
  !!!!!!!!!!!

  !!!!!!!!!!!
  ! Function randn()
  !!!!!!!!!!!

  ! Code is a modified version of that found in the uniform distribution generator from Numerical receipes
  ! Modified to give 0-1 unform on input of value 0 and normal distribution (mean = 0 sd = 1) on input of 1
  ! Modified by TLS

  !!!!!!!!!!!

  ! randn() related seed value
  double precision :: idum

  ! rand(), narray(), rnstrt() related
  integer, parameter  :: kk=100, ll=37, mm=2**30, tt=70, kkk=kk+kk-1
  integer, save       :: ranx(kk)

  contains

  !
  !--------------------------------------------------------------------
  !
  double precision function std(a,n)

    ! Function to determine the standard deviation
    ! inputs are the vector of values and number of values included

    implicit none

    ! declare inputs
    integer, intent(in) :: n ! number of values in vector
    double precision, dimension(n), intent(in) :: a

    ! declare local variables
    double precision mean, sq_diff_sum, diff, variance
    integer i

    ! if no length has been returned then provide value which ensures crash (i.e.
    ! infinity)
    if (n == 0) then
        std=0d0
        write(*,*) "no sample size has been provided to std function"
        return
    endif

    ! first calculate the mean
    mean=sum(a)/dble(n)

    ! ensure zero values
    diff = 0d0 ; sq_diff_sum = 0d0

    ! calculate cumulative square difference
    do i = 1, n
       diff = a(i)-mean
       sq_diff_sum = sq_diff_sum + diff**2d0
    end do

    ! calculate the variance
    variance = sq_diff_sum/dble(n-1)

    ! return the standard deviation
    std = sqrt(variance)

    return

  end function std
  !
  !------------------------------------------------------------------
  !
  double precision function randn(option)

    ! from Numerical Receipes p271 Press et al., 1986 2nd Edition Chapter 7,
    ! Random Numbers function returns real random number between 0-1 based on an initial start
    ! point (ran1).
    ! The start point (default = -1) is reinitialised every time the model runs
    ! providing the same distribution each run
    ! To ensure random numbers each time use the sum of the current system time

    ! modified based on blooms C code to alter range of random numbers

    implicit none
    integer IA,IM,IQ,IR,NTAB,NDIV,option
    double precision AM,EPS,RNMX,const,r1,r2,pi
    parameter(IA=16807,IM=2147483647,AM=1d0/dble(IM),IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2d-30,RNMX=1d0-EPS)
    integer j,k,iv(NTAB),iy
    SAVE iv,iy
    DATA iv /NTAB*0/
    DATA iy /0/

    const=1d0
    pi=3.141592653589793d0

    if (option == 0) then
      if (idum < 0d0 .or. iy == 0) then
          idum = max(-idum,const)
          do j = (NTAB+8), 1, -1
             k = nint(idum/dble(IQ))
             idum=dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
             if (idum < 0d0) idum = idum + dble(IM)
             if (j < NTAB) iv(j) = nint(idum)
          enddo
          iy = iv(1)
      endif
      k = nint(idum/dble(IQ))
      idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
      if (idum < 0d0) idum = idum+dble(IM)
      j = 1+iy/NDIV
      iy = iv(j)
      iv(j) = nint(idum)

      ! output now
      randn=min(AM*dble(iy), RNMX)
      return

    else

      if (idum < 0d0 .or. iy == 0) then
          idum = max(-idum,const)
          do j = (NTAB+8),1,-1
             k = nint(idum/dble(IQ))
             idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
             if (idum < 0d0) idum = idum+dble(IM)
             if (j < NTAB) iv(j) = nint(idum)
          enddo
          iy = iv(1)
      endif
      k=nint(idum)/IQ
      idum=dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
      if (idum < 0d0) idum = idum+dble(IM)
      j = 1+iy/NDIV
      iy = iv(j)
      iv(j) = nint(idum)
      r1 = max(min(AM*dble(iy), RNMX),1d-30)

      if (idum < 0d0 .or. iy == 0) then
          idum = max(-idum,const)
          do j = (NTAB+8),1,-1
             k = nint(idum)/IQ
             idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
             if (idum < 0d0) idum = idum+dble(IM)
             if (j < NTAB) iv(j) = nint(idum)
          enddo
          iy = iv(1)
      endif
      k = nint(idum)/IQ
      idum = dble(IA)*(idum-dble(k*IQ))-dble(IR*k)
      if (idum < 0d0) idum = idum + dble(IM)
      j = 1+iy/NDIV
      iy = iv(j)
      iv(j) = nint(idum)
      r2 = max(min(AM*dble(iy), RNMX),1d-30)

      ! output now
      randn = sqrt(-2d0*log(r1)) * cos(2d0*pi*r2)
      return

    endif

  end function
  !
  !--------------------------------------------------------------------
  !
  subroutine random_normal(uniform, random_length, uniform_random_vector, fn_val)

    ! Generate a random normal deviate using the polar method.
    ! Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for
    ! generating normal variables',
    ! Siam Rev., vol.6, 260-264, 1964.

    implicit none

    ! arguments
    integer, intent(in) :: random_length
    integer, intent(inout) :: uniform
    double precision, intent(out) :: fn_val
    double precision, dimension(random_length), intent(inout) :: uniform_random_vector

    ! Local variables
    double precision :: u, sumsq
    double precision, save :: v, sln
    logical, save   :: second = .false.
    double precision, parameter :: one = 1d0, vsmall = tiny( one )

    if (second) then
        ! If second, use the second random number generated on last call
        second = .false.
        fn_val = v*sln

    else

        ! First call; generate a pair of random normals
        second = .true.
        do
!           CALL RANDOM_NUMBER( u )
!           CALL RANDOM_NUMBER( v )
           u = uniform_random_vector(uniform) ; uniform = uniform + 1
           v = uniform_random_vector(uniform) ; uniform = uniform + 1
           u = scale( u, 1 ) - one
           v = scale( v, 1 ) - one
           sumsq = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
           if (uniform >= random_length) then
               call random_uniform(uniform_random_vector,random_length)
               uniform = 1
           endif
           if (sumsq < one) exit
        end do
        sln = sqrt(-scale( log(sumsq), 1 ) / sumsq)
        fn_val = u*sln
    end if

    return

  end subroutine random_normal
  !
  !--------------------------------------------------------------------
  !
  subroutine random_uniform(u, n)

    ! Generate an array of n double precision values between 0 and 1.

    integer, intent(in)  :: n ! number of random values wanted
    double precision, intent(out) :: u(n) ! output vector

    ! Local array
    integer, allocatable, dimension(:)  :: aa

    ! allocate memory
    allocate(aa(n))

    call rnarry(aa, n)
    u(1:n) = scale( dble(aa), -30)

    ! tidy
    deallocate(aa)

    return

  end subroutine random_uniform
  !
  !--------------------------------------------------------------------
  !
  subroutine rnarry(aa, n)

    ! Generate an array of n integers between 0 and 2^30-1.

    integer, intent(in)   :: n
    integer, intent(out)  :: aa(n)

    ! Local variables
    integer  :: j

    aa(1:kk) = ranx(1:kk)
    do j = kk + 1, n
       aa(j) = aa(j-kk) - aa(j-ll)
       if (aa(j) < 0) aa(j) = aa(j) + mm
    end do
    do j = 1, ll
       ranx(j) = aa(n+j-kk) - aa(n+j-ll)
       if (ranx(j) < 0) ranx(j) = ranx(j) + mm
    end do
    do j = ll+1, kk
       ranx(j) = aa(n+j-kk) - ranx(j-ll)
       if (ranx(j) < 0) ranx(j) = ranx(j) + mm
    end do

    return

  end subroutine rnarry
  !
  !--------------------------------------------------------------------
  !
  subroutine rnstrt(seed)

    ! Initialize integer array ranx using the input seed.

    integer, intent(in)  :: seed

    ! Local variables
    integer  :: x(kkk), j, ss, sseed, t

    if (seed < 0) then
        sseed = mm - 1 - mod(-1-seed, mm)
    else
        sseed = mod(seed, mm)
    end if
    ss = sseed - mod(sseed,2) + 2
    do j = 1, kk
       x(j) = ss
       ss = ishft(ss, 1)
       if (ss >= mm) ss = ss - mm + 2
    end do
    x(kk+1:kkk) = 0
    x(2) = x(2)+1
    ss = sseed
    t = tt - 1
10  do j = kk, 2, -1
       x(j+j-1) = x(j)
    end do
    do j = kkk, kk + 1, -1
       x(j-(kk-ll)) = x(j-(kk-ll)) - x(j)
       if (x(j-(kk-ll)) < 0) x(j-(kk-ll)) = x(j-(kk-ll)) + mm
       x(j-kk) = x(j-kk) - x(j)
       if (x(j-kk) < 0) x(j-kk) = x(j-kk) + mm
    end do
    if (mod(ss,2) == 1) then
        do j = kk, 1, -1
           x(j+1) = x(j)
        end do
        x(1) = x(kk+1)
        x(ll+1) = x(ll+1) - x(kk+1)
        if (x(ll+1) < 0) x(ll+1) = x(ll+1) + mm
    end if
    if (ss /= 0) THEN
        ss = ishft(ss, -1)
    else
        t = t - 1
    end if
    if (t > 0) GO TO 10

    do j = 1, ll
       ranx(j+kk-ll) = x(j)
    end do
    do j = ll+1, kk
       ranx(j-ll) = x(j)
    end do

    do j = 1, 10
       call rnarry(x,kkk)
    end do

    return
  end subroutine rnstrt
  !
  !--------------------------------------------------------------------
  !
end module math_functions
