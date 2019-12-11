C     fcn is the name of the user-supplied subroutine which
C     calculates the functions. fcn must be declared
C     in an external statement in the user calling
C     program, and should be written as follows.
C
C     subroutine fcn(n, x, fvec, iflag)
C       integer n, iflag
C       double precision x(n), fvec(n)
C       ----------
C       calculate the functions at x and
C       return this vector in fvec.
C       ---------
C       return
C     end
C
C     the value of iflag should not be changed by fcn unless
C     the user wants to terminate execution of hybrd1.
C     in this case set iflag to a negative integer.
C
C     n is a positive integer input variable set to the number
C       of functions and variables.
C
C     x is an array of length n. on input x must contain
C       an initial estimate of the solution vector. on output x
C       contains the final estimate of the solution vector.
C
C     fvec is an output array of length n which contains
C       the functions evaluated at the output x.
C
C     tol is a nonnegative input variable. termination occurs
C       when the algorithm estimates that the relative error
C       between x and the solution is at most tol.
C
C     returns an integer output variable. if the user has
C       terminated execution, info is set to the (negative)
C       value of iflag. see description of fcn. otherwise,
C       info is set as follows.
C       info = 0   improper input parameters.
C       info = 1   algorithm estimates that the relative error
C                  between x and the solution is at most tol.
C       info = 2   number of calls to fcn has reached or exceeded
C                  200*(n+1).
C       info = 3   tol is too small. no further improvement in
C                  the approximate solution x is possible.
C       info = 4   iteration is not making good progress.

      subroutine root(n,x,tol,info)
      integer n,info
      double precision tol
      double precision x(n)
      double precision fvec(n),wa((n*(3*n+13))/2)
Cf2py intent(hide) :: n
Cf2py intent(in, out, copy), depend(n) :: x
Cf2py intent(in) :: tol
Cf2py intent(out) :: info

      call hybrd1(n,x,fvec,tol,info,wa,(n*(3*n+13))/2)
      end subroutine

      subroutine fcn(n,x,fvec,iflag)
      integer n, iflag
      double precision x(n),fvec(n)
      external fun
Cf2py intent(hide) :: n
Cf2py intent(in, copy), depend(n) :: x
Cf2py intent(out, copy), depend(n) :: fvec
Cf2py intent(in, out) :: iflag
Cf2py intent(callback, hide) :: fun

      call fun(n, x, iflag, fvec)
      end subroutine

      subroutine hybrd1(n,x,fvec,tol,info,wa,lwa)
      integer n,info,lwa
      double precision tol
      double precision x(n),fvec(n),wa(lwa)

      integer index,j,lr,maxfev,ml,mode,mu,nfev,nprint
      double precision epsfcn,factor,one,xtol,zero
      data factor,one,zero /1.0d2,1.0d0,0.0d0/
      info = 0
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. tol .lt. zero .or. lwa .lt. (n*(3*n + 13))/2)
     *   go to 20
c
c     call hybrd.
c
      maxfev = 200*(n + 1)
      xtol = tol
      ml = n - 1
      mu = n - 1
      epsfcn = zero
      mode = 2
      do 10 j = 1, n
         wa(j) = one
   10    continue
      nprint = 0
      lr = (n*(n + 1))/2
      index = 6*n + lr
      call hybrd(n,x,fvec,xtol,maxfev,ml,mu,epsfcn,wa(1),mode,
     *           factor,nprint,info,nfev,wa(index+1),n,wa(6*n+1),lr,
     *           wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info .eq. 5) info = 4
   20 continue
      return
c
c     last card of subroutine hybrd1.
c
      end

      subroutine hybrd(n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,
     *                 mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,
     *                 qtf,wa1,wa2,wa3,wa4)
      integer n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr
      double precision xtol,epsfcn,factor
      double precision x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr),
     *                 qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)

      integer i,iflag,iter,j,jm1,l,msum,ncfail,ncsuc,nslow1,nslow2
      integer iwa(1)
      logical jeval,sing
      double precision actred,delta,epsmch,fnorm,fnorm1,one,pnorm,
     *                 prered,p1,p5,p001,p0001,ratio,sum,temp,xnorm,
     *                 zero
      double precision dpmpar,enorm
      data one,p1,p5,p001,p0001,zero
     *     /1.0d0,1.0d-1,5.0d-1,1.0d-3,1.0d-4,0.0d0/
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      info = 0
      iflag = 0
      nfev = 0
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. xtol .lt. zero .or. maxfev .le. 0
     *    .or. ml .lt. 0 .or. mu .lt. 0 .or. factor .le. zero
     *    .or. ldfjac .lt. n .or. lr .lt. (n*(n + 1))/2) go to 300
      if (mode .ne. 2) go to 20
      do 10 j = 1, n
         if (diag(j) .le. zero) go to 300
   10    continue
   20 continue
c
c     evaluate the function at the starting point
c     and calculate its norm.
c
      iflag = 1
      call fcn(n,x,fvec,iflag)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = enorm(n,fvec)
c
c     determine the number of calls to fcn needed to compute
c     the jacobian matrix.
c
      msum = min0(ml+mu+1,n)
c
c     initialize iteration counter and monitors.
c
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
c
c     beginning of the outer loop.
c
   30 continue
         jeval = .true.
c
c        calculate the jacobian matrix.
c
         iflag = 2
         call fdjac1(n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,
     *               wa2)
         nfev = nfev + msum
         if (iflag .lt. 0) go to 300
c
c        compute the qr factorization of the jacobian.
c
         call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
c
c        on the first iteration and if mode is 1, scale according
c        to the norms of the columns of the initial jacobian.
c
         if (iter .ne. 1) go to 70
         if (mode .eq. 2) go to 50
         do 40 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) .eq. zero) diag(j) = one
   40       continue
   50    continue
c
c        on the first iteration, calculate the norm of the scaled x
c        and initialize the step bound delta.
c
         do 60 j = 1, n
            wa3(j) = diag(j)*x(j)
   60       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta .eq. zero) delta = factor
   70    continue
c
c        form (q transpose)*fvec and store in qtf.
c
         do 80 i = 1, n
            qtf(i) = fvec(i)
   80       continue
         do 120 j = 1, n
            if (fjac(j,j) .eq. zero) go to 110
            sum = zero
            do 90 i = j, n
               sum = sum + fjac(i,j)*qtf(i)
   90          continue
            temp = -sum/fjac(j,j)
            do 100 i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
  100          continue
  110       continue
  120       continue
c
c        copy the triangular factor of the qr factorization into r.
c
         sing = .false.
         do 150 j = 1, n
            l = j
            jm1 = j - 1
            if (jm1 .lt. 1) go to 140
            do 130 i = 1, jm1
               r(l) = fjac(i,j)
               l = l + n - i
  130          continue
  140       continue
            r(l) = wa1(j)
            if (wa1(j) .eq. zero) sing = .true.
  150       continue
c
c        accumulate the orthogonal factor in fjac.
c
         call qform(n,n,fjac,ldfjac,wa1)
c
c        rescale if necessary.
c
         if (mode .eq. 2) go to 170
         do 160 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  160       continue
  170    continue
c
c        beginning of the inner loop.
c
  180    continue
c
c           if requested, call fcn to enable printing of iterates.
c
            if (nprint .le. 0) go to 190
            iflag = 0
            if (mod(iter-1,nprint) .eq. 0) call fcn(n,x,fvec,iflag)
            if (iflag .lt. 0) go to 300
  190       continue
c
c           determine the direction p.
c
            call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
c
c           store the direction p and x + p. calculate the norm of p.
c
            do 200 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  200          continue
            pnorm = enorm(n,wa3)
c
c           on the first iteration, adjust the initial step bound.
c
            if (iter .eq. 1) delta = dmin1(delta,pnorm)
c
c           evaluate the function at x + p and calculate its norm.
c
            iflag = 1
            call fcn(n,wa2,wa4,iflag)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = enorm(n,wa4)
c
c           compute the scaled actual reduction.
c
            actred = -one
            if (fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
c
c           compute the scaled predicted reduction.
c
            l = 1
            do 220 i = 1, n
               sum = zero
               do 210 j = i, n
                  sum = sum + r(l)*wa1(j)
                  l = l + 1
  210             continue
               wa3(i) = qtf(i) + sum
  220          continue
            temp = enorm(n,wa3)
            prered = zero
            if (temp .lt. fnorm) prered = one - (temp/fnorm)**2
c
c           compute the ratio of the actual to the predicted
c           reduction.
c
            ratio = zero
            if (prered .gt. zero) ratio = actred/prered
c
c           update the step bound.
c
            if (ratio .ge. p1) go to 230
               ncsuc = 0
               ncfail = ncfail + 1
               delta = p5*delta
               go to 240
  230       continue
               ncfail = 0
               ncsuc = ncsuc + 1
               if (ratio .ge. p5 .or. ncsuc .gt. 1)
     *            delta = dmax1(delta,pnorm/p5)
               if (dabs(ratio-one) .le. p1) delta = pnorm/p5
  240       continue
c
c           test for successful iteration.
c
            if (ratio .lt. p0001) go to 260
c
c           successful iteration. update x, fvec, and their norms.
c
            do 250 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
               fvec(j) = wa4(j)
  250          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  260       continue
c
c           determine the progress of the iteration.
c
            nslow1 = nslow1 + 1
            if (actred .ge. p001) nslow1 = 0
            if (jeval) nslow2 = nslow2 + 1
            if (actred .ge. p1) nslow2 = 0
c
c           test for convergence.
c
            if (delta .le. xtol*xnorm .or. fnorm .eq. zero) info = 1
            if (info .ne. 0) go to 300
c
c           tests for termination and stringent tolerances.
c
            if (nfev .ge. maxfev) info = 2
            if (p1*dmax1(p1*delta,pnorm) .le. epsmch*xnorm) info = 3
            if (nslow2 .eq. 5) info = 4
            if (nslow1 .eq. 10) info = 5
            if (info .ne. 0) go to 300
c
c           criterion for recalculating jacobian approximation
c           by forward differences.
c
            if (ncfail .eq. 2) go to 290
c
c           calculate the rank one modification to the jacobian
c           and update qtf if necessary.
c
            do 280 j = 1, n
               sum = zero
               do 270 i = 1, n
                  sum = sum + fjac(i,j)*wa4(i)
  270             continue
               wa2(j) = (sum - wa3(j))/pnorm
               wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
               if (ratio .ge. p0001) qtf(j) = sum
  280          continue
c
c           compute the qr factorization of the updated jacobian.
c
            call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
            call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
            call r1mpyq(1,n,qtf,1,wa2,wa3)
c
c           end of the inner loop.
c
            jeval = .false.
            go to 180
  290    continue
c
c        end of the outer loop.
c
         go to 30
  300 continue
c
c     termination, either normal or user imposed.
c
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0) call fcn(n,x,fvec,iflag)
      return
c
c     last card of subroutine hybrd.
c
      end

      subroutine fdjac1(n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,
     *                  wa1,wa2)
      integer n,ldfjac,iflag,ml,mu
      double precision epsfcn
      double precision x(n),fvec(n),fjac(ldfjac,n),wa1(n),wa2(n)

      integer i,j,k,msum
      double precision eps,epsmch,h,temp,zero
      double precision dpmpar
      data zero /0.0d0/
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      eps = dsqrt(dmax1(epsfcn,epsmch))
      msum = ml + mu + 1
      if (msum .lt. n) go to 40
c
c        computation of dense approximate jacobian.
c
         do 20 j = 1, n
            temp = x(j)
            h = eps*dabs(temp)
            if (h .eq. zero) h = eps
            x(j) = temp + h
            call fcn(n,x,wa1,iflag)
            if (iflag .lt. 0) go to 30
            x(j) = temp
            do 10 i = 1, n
               fjac(i,j) = (wa1(i) - fvec(i))/h
   10          continue
   20       continue
   30    continue
         go to 110
   40 continue
c
c        computation of banded approximate jacobian.
c
         do 90 k = 1, msum
            do 60 j = k, n, msum
               wa2(j) = x(j)
               h = eps*dabs(wa2(j))
               if (h .eq. zero) h = eps
               x(j) = wa2(j) + h
   60          continue
            call fcn(n,x,wa1,iflag)
            if (iflag .lt. 0) go to 100
            do 80 j = k, n, msum
               x(j) = wa2(j)
               h = eps*dabs(wa2(j))
               if (h .eq. zero) h = eps
               do 70 i = 1, n
                  fjac(i,j) = zero
                  if (i .ge. j - mu .and. i .le. j + ml)
     *               fjac(i,j) = (wa1(i) - fvec(i))/h
   70             continue
   80          continue
   90       continue
  100    continue
  110 continue
      return
c
c     last card of subroutine fdjac1.
c
      end