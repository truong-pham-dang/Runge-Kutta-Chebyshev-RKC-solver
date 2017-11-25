
c                             Example A
c
c  This example shows how to use RKC.  It solves a system of ODEs that
c  arise from semi-discretization of the reaction-diffusion equation
c
c            U  = U   + (1 - U)*U**2    for   t >= 0,   0 <= x <= 10
c             t    xx
c
c  Dirichlet boundary conditions specify U(0,t) and U(10,t) for all t >= 0
c  and the initial values U(x,0) are specified. These values are taken from
c  an analytical solution that is evaluated in sol(x,t) so that the numerical
c  solution can be compared to a known solution.
c
c  A semi-discretization of the PDE is obtained by choosing a set of points
c  {x_i} in [0, 10] and approximating U(x_i,t) by a function y_i(t).  Here
c  neqn+2 equally spaced points x_i are used for neqn = 99.  When the second
c  partial derivative of U with respect to x is approximated by central
c  differences, a system of neqn ODEs is obtained for the y_i(t).  The
c  initial values y_i(0) are given by U(x_i,0) = sol(x_i,0).
c
c  A common way to present the computed results is to plot approximations
c  to U(x,tout) on [0, 10] for a selection of times tout.  This example
c  shows how to compute approximations to the y_i(t) at these specific times.
c  They are written to an output file for plotting; the most convenient way
c  to do this will depend on the system and the plotting package used.
c
c  Because an analytical solution U(x,t) is available, the maximum error
c  of the approximation to U(x,tend) is computed and displayed.  Here
c  tend = 15.  It should be appreciated that this error has two parts,
c  one the error made by RKC in the time integration and the other from
c  the spacial discretization.  Some statistics about the integration are
c  also displayed.
c
      integer           neqn,nout
      parameter        (neqn=99,nout=4)
      integer           info(4),idid
      double precision  t,tend,rtol,atol
      double precision  y(neqn),work(8+5*neqn)
      integer           i,next
      double precision  dx,delta,sol,tout(nout),yout(neqn),
     &                  truey,error
      integer           nfe,nsteps,naccpt,nrejct,nfesig,maxm
      common /rkcdid/   nfe,nsteps,naccpt,nrejct,nfesig,maxm
      external          f
      character*5       cnext
c------------------------------------------------------------
c  Specify the interval of integration in time.  The initial
c  values of the solution at mesh points are provided by the
c  analytical solution sol(x,t).
c------------------------------------------------------------
      t = 0d0
      tend = 15d0
      dx = 10d0/(neqn+1)
      do i = 1, neqn
        y(i) = sol(i*dx,t)
      enddo

c----------------------------------------------------------------------
c  Initialize the output:  Define the times at which solutions are to
c  be reported for plotting and output the number of these times.  A
c  solution is computed on a mesh of equally spaced points in [0, 10].
c  Output the number of points in the mesh and then the mesh itself.
c  Because tout(1) = t, output the initial values for the neqn solution
c  components along with the values given at 0 and 10.
c----------------------------------------------------------------------
      delta = (tend - t)/(nout-1)
      do i = 1,nout
        tout(i) = t + (i-1)*delta
      enddo
      open(10,file='exaout')
      write(10,'(i10)') nout,neqn+2
      write(10,*) 'Vector of space:'
      do i = 0,neqn+1
        write(10,'(e10.4)') i*dx
      enddo
      next = 1

      write(cnext,'(i5.5)') next
      open(unit = next, file = 'solution_'//trim(cnext)//'.dat',
     & status = 'replace')
      write(next,*) 'x, u, u_e'
      do i = 0,neqn+1
        if (i == 0) then
            write(next,*) i*dx, sol(0d0,tout(next)), sol(0d0,tout(next))
        else if (i == neqn + 1) then
            write(next,*) i*dx, sol(10d0,tout(next)), sol(10d0,tout(next))
        else
            write(next,*) i*dx, y(i), sol(i*dx,tout(next))
        endif
      enddo
      close(unit = next)

      write(10,*) 'Vector of initial solution at time ',tout(next)
      write(10,'(e10.4)') sol(0d0,tout(next))
      write(10,'(e10.4)') y
      write(10,'(e10.4)') sol(10d0,tout(next))
      next = next + 1
c--------------------------------------------------
c  To compute results at specific times, the code
c  must return after each step.  Common choices for
c  info(*) have value 0.
c  info(1) = 0 -- return after each step.
c  info(2) = 0 -- RKC computes the spectral radius.
c  info(3) = 0 -- the Jacobian may not be constant.
c  info(4) = 0 -- ATOL is a scalar.
c--------------------------------------------------
      info(1) = 0
      info(2) = 0
      info(3) = 0
      info(4) = 0
c-------------------------
c  Specify the tolerances.
c-------------------------
      rtol = 1d-4
      atol = rtol
c-----------------------------
c  Initialize the integration.
c-----------------------------
      idid = 0
c---------------------
c  Take a single step:
c---------------------
40    continue
      call rkc(neqn,f,y,t,tend,rtol,atol,info,work,idid)
c-------------------------------------------------------------
c  Was the step successful?  If not, quit with an explanation.
c-------------------------------------------------------------
      if(idid .gt. 2) then
        write(*,*) ' Failed at t = ',t,' with idid = ',idid
        stop
      endif
c--------------------------------------------------------------
c  To get output at specific points, step towards TEND with RKC
c  until the integration passes the next output point. Compute
c  a result at the point using RKCINT.  There might be several
c  output points in the span of a single step by RKC.
c--------------------------------------------------------------
50    continue
      if(t .ge. tout(next)) then
        call rkcint(work,tout(next),yout)
        write(10,*) 'Vector of solution at time ',tout(next)
        write(10,'(e10.4)') sol(0d0,tout(next))
        write(10,'(e10.4)') yout
        write(10,'(e10.4)') sol(10d0,tout(next))

        write(cnext,'(i5.5)') next
        open(unit = next, file = 'solution_'//trim(cnext)//'.dat',
     & status = 'replace')
        write(next,*) 'x, u, u_e'
        do i = 0,neqn+1
            if (i == 0) then
                write(next,*) i*dx, sol(0d0,tout(next)), sol(0d0,tout(next))
            else if (i == neqn + 1) then
                write(next,*) i*dx, sol(10d0,tout(next)), sol(10d0,tout(next))
            else
                write(next,*) i*dx, yout(i), sol(i*dx,tout(next))
            endif
        enddo
        close(unit = next)

        next = next + 1
        if(next .le. nout) goto 50
      endif
c--------------------------------------
c  Monitor the cost of the integration.
c--------------------------------------
      if(nsteps .ge. 5000) then
        write(*,*) ' Quit because of too much work.'
      endif
c-------------------------------------
c  If not done yet, take another step.
c-------------------------------------
      if(idid .eq. 2) goto 40
c------------------------------------------------------
c  Done.  Compute the error and report some statistics.
c------------------------------------------------------
      error = 0d0
      do i = 1, neqn
        truey = sol(i*dx,t)
        error = max(error,abs(y(i) - truey))
      enddo
      write(*,'(/a,d8.1,a,f6.1,a,d8.2)') ' With rtol = atol =',rtol,
     &  ', the maximum error at tend =',tend,' was',error
      write(*,'(a,i5,a)') ' The integration cost',nfe,
     &  ' function evaluations.'
      write(*,'(a,i4,a,i3,a)') ' There were',nsteps,' steps (',
     &  nrejct,' rejected).'
      write(*,'(a,i4)') ' The maximum number of stages used was',maxm
      end

      double precision function sol(x,t)
c------------------------------------------------------------
c  An analytical solution to the reaction-diffusion equation.
c------------------------------------------------------------
      double precision x,t
      double precision v,z
      v = sqrt(0.5d0)
      z = x - v*t
      sol = 1d0/(1d0 + exp(v*z))
      return
      end

      subroutine f(neqn,t,y,dy)
c---------------------------------------------------------------
c  Semi-discretization of reaction-diffusion equation by central
c  differences.  The analytical solution sol(x,t) is used for
c  Dirichlet boundary conditions at x = 0 and x = 10.
c---------------------------------------------------------------
      integer           neqn
      double precision  t,y(neqn),dy(neqn)
      integer           i
      double precision  dx,dxsq,sol
c
      dx = 10d0/(neqn+1)
      dxsq = dx**2
      dy(1) = (sol(0d0,t)- 2d0*y(1) + y(2))/dxsq +
     &                                   (1d0 - y(1))*y(1)**2
      do 10 i = 2,neqn-1
         dy(i) = (y(i-1) - 2d0*y(i) + y(i+1))/dxsq +
     &                                   (1d0 - y(i))*y(i)**2
10    continue
      dy(neqn) = (y(neqn-1)- 2d0*y(neqn) + sol(10d0,t))/dxsq +
     &                                   (1d0 - y(neqn))*y(neqn)**2
      return
      end

      double precision function spcrad(neqn,t,y)
c--------------------------
c  This is a dummy routine.
c--------------------------
      integer          neqn
      double precision t,y(neqn)
      spcrad = 0d0
      return
      end




