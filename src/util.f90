!MODULE numerical_recipes
!
!  PUBLIC :: zbrent, zbrac, mnbrak, golden
!
!  CONTAINS
!!--------------------------------------------------------------------
!SUBROUTINE zbrac(func,x1,x2,success)
!!--------------------------------------------------------------------
!  USE kinds, ONLY: DP
!  IMPLICIT NONE
!  REAL(DP), INTENT(INOUT) :: x1, x2
!  LOGICAL, INTENT(OUT) :: success
!  INTERFACE
!     FUNCTION func(x)
!       USE kinds, ONLY: DP
!       IMPLICIT NONE
!       REAL(DP), INTENT(IN) :: x
!       REAL(DP) :: func
!     END FUNCTION func
!  END INTERFACE
!  REAL(DP), PARAMETER :: FACTOR = 1.6
!  INTEGER, PARAMETER :: NTRY = 50
!  !
!  INTEGER :: j
!  REAL(DP) :: f1, f2
!  IF ( x1 .EQ. x2 ) THEN
!     WRITE(*,*)'ERROR: need to provide different values for the extremes of the interval'
!     success = .FALSE.
!     RETURN
!  ENDIF
!  f1 = func(x1)
!  f2 = func(x2)
!  success = .TRUE.
!  DO j = 1, NTRY
!     IF ( f1*f2.LT.0.D0 ) RETURN
!     IF (ABS(f1).LT.ABS(f2))THEN
!        x1 = x1 + FACTOR*(x1-x2)
!        f1 = func(x1)
!     ELSE
!        x2 = x2 + FACTOR*(x2-x1)
!        f2 = func(x2)
!     ENDIF
!  ENDDO
!  success = .FALSE.
!  RETURN
!!--------------------------------------------------------------------
!END SUBROUTINE zbrac
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!SUBROUTINE zbrent(func,x1,x2,tol,lambda)
!!--------------------------------------------------------------------
!  USE kinds, ONLY: DP
!  IMPLICIT NONE
!  !
!  REAL(DP), INTENT(IN) :: x1, x2
!  REAL(DP), INTENT(IN) :: tol
!  REAL(DP), INTENT(OUT) :: lambda
!  INTERFACE
!     FUNCTION func(x)
!       USE kinds, ONLY: DP
!       IMPLICIT NONE
!       REAL(DP), INTENT(IN) :: x
!       REAL(DP) :: func
!     END FUNCTION func
!  END INTERFACE
!  !
!  INTEGER, PARAMETER :: ITMAX = 100
!  REAL(DP), PARAMETER :: EPS = 3.D-8
!  !
!  REAL(DP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
!  INTEGER :: iter
!  a=x1
!  b=x2
!  fa=func(a)
!  fb=func(b)
!  IF ((fa.GT.0.D0.AND.fb.GT.0.D0).OR.(fa.LT.0.D0.AND.fb.LT.0.D0)) THEN
!     WRITE(*,*)'ERROR: root must be bracketed for zbrent'
!     STOP
!  ENDIF
!  c=b
!  fc=fb
!  DO iter = 1,ITMAX
!     IF ((fb.GT.0.D0.AND.fc.GT.0.D0).OR.(fb.LT.0.D0.AND.fc.LT.0.D0)) THEN
!        c=a
!        fc=fa
!        d=b*a
!        e=d
!     ENDIF
!     IF (ABS(fc).LT.ABS(fb)) THEN
!        a=b
!        b=c
!        c=a
!        fa=fb
!        fb=fc
!        fc=fa
!     ENDIF
!     TOL1=2*EPS*ABS(b)+0.5D0*tol
!     xm=0.5D0*(c-b)
!     IF (ABS(xm).LE.tol1 .OR. fb.EQ.0.D0) THEN
!        lambda=b
!        RETURN
!     ENDIF
!     IF (ABS(e).GE.tol1 .AND. ABS(fa).GT.ABS(fb)) THEN
!        s=fb/fa
!        IF (a.EQ.c) THEN
!           p=2.D0*xm*s
!           q=1.D0-s
!        ELSE
!           q=fa/fc
!           r=fb/fc
!           p=s*(2.D0*xm*q*(q-r)-(b-a)*(r-1.D0))
!           q=(q-1.D0)*(r-1.D0)*(s-1.D0)
!        ENDIF
!        IF (p.GT.0.D0) q=-q
!        p=ABS(p)
!        IF (2.D0*p .LT. MIN(3.D0*xm*q-ABS(tol1*q),ABS(e*q))) THEN
!           e=d
!           d=p/q
!        ELSE
!           d=xm
!           e=d
!        ENDIF
!     ELSE
!        d=xm
!        e=d
!     ENDIF
!     a=b
!     fa=fb
!     IF ( ABS(d) .GT. tol1 ) THEN
!        b=b+d
!     ELSE
!        b=b+SIGN(tol1,xm)
!     ENDIF
!     fb=func(b)
!  ENDDO
!  WRITE(*,*)'ERROR: zbrent exceeding maximum iterations'
!  lambda=b
!  RETURN
!!--------------------------------------------------------------------
!END SUBROUTINE zbrent
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!        SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
!!--------------------------------------------------------------------
!	USE kinds, ONLY: DP
!	IMPLICIT NONE
!	REAL(DP), INTENT(INOUT) :: ax,bx
!	REAL(DP), INTENT(OUT) :: cx,fa,fb,fc
!	INTERFACE
!		FUNCTION func(x)
!		USE kinds, ONLY: DP
!		IMPLICIT NONE
!		REAL(DP), INTENT(IN) :: x
!		REAL(DP) :: func
!		END FUNCTION func
!	END INTERFACE
!	REAL(DP), PARAMETER :: GOLD=1.618034_DP,GLIMIT=100.0_DP,TINY=1.0e-20_DP
!	REAL(DP) :: fu,q,r,u,ulim
!	fa=func(ax)
!	fb=func(bx)
!	if (fb > fa) then
!           fc = ax
!           ax = bx
!           bx = fc
!           fc = fa
!           fa = fb
!           fb = fc
!	end if
!	cx=bx+GOLD*(bx-ax)
!        fc=func(cx)
!	do
!		if (fb < fc) RETURN
!		r=(bx-ax)*(fb-fc)
!		q=(bx-cx)*(fb-fa)
!		u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_DP*sign(max(abs(q-r),TINY),q-r))
!		ulim=bx+GLIMIT*(cx-bx)
!		if ((bx-u)*(u-cx) > 0.0) then
!                fu=func(u)
!			if (fu < fc) then
!				ax=bx
!				fa=fb
!				bx=u
!				fb=fu
!				RETURN
!			else if (fu > fb) then
!				cx=u
!				fc=fu
!				RETURN
!			end if
!			u=cx+GOLD*(cx-bx)
!			fu=func(u)
!		else if ((cx-u)*(u-ulim) > 0.0) then
!			fu=func(u)
!			if (fu < fc) then
!				bx=cx
!				cx=u
!				u=cx+GOLD*(cx-bx)
!				call shft(fb,fc,fu,func(u))
!			end if
!		else if ((u-ulim)*(ulim-cx) >= 0.0) then
!			u=ulim
!			fu=func(u)
!		else
!			u=cx+GOLD*(cx-bx)
!			fu=func(u)
!		end if
!		call shft(ax,bx,cx,u)
!		call shft(fa,fb,fc,fu)
!	end do
!	CONTAINS
!!BL
!	SUBROUTINE shft(a,b,c,d)
!	REAL(DP), INTENT(OUT) :: a
!	REAL(DP), INTENT(INOUT) :: b,c
!	REAL(DP), INTENT(IN) :: d
!	a=b
!	b=c
!	c=d
!	END SUBROUTINE shft
!	END SUBROUTINE mnbrak
!!--------------------------------------------------------------------
!      FUNCTION golden(ax,bx,cx,f,tol,xmin)
!!--------------------------------------------------------------------
!      USE kinds, ONLY: DP
!      IMPLICIT NONE
!      REAL(DP) :: golden,ax,bx,cx,tol,xmin,R,C
!      INTERFACE
!        FUNCTION f(x)
!          USE kinds, ONLY: DP
!          IMPLICIT NONE
!          REAL(DP), INTENT(IN) :: x
!          REAL(DP) :: f
!        END FUNCTION f
!      END INTERFACE
!      PARAMETER (R=.61803399,C=1.-R)
!!Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
!!between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine performs
!!a golden section search for the minimum, isolating it to a fractional precision of about
!!tol. The abscissa of the minimum is returned as xmin, and the minimum function value
!!is returned as golden, the returned function value.
!!Parameters: The golden ratios.
!      REAL(DP) f1,f2,x0,x1,x2,x3
!      x0=ax
!      x3=cx
!      if(abs(cx-bx).gt.abs(bx-ax))then ! Make x0 to x1 the smaller segment,
!        x1=bx
!        x2=bx+C*(cx-bx) ! and fill in the new point to be tried.
!      else
!        x2=bx
!        x1=bx-C*(bx-ax)
!      endif
!      f1=f(x1) ! The initial function evaluations. Note that we never need to
!      f2=f(x2) ! evaluate the function at the original endpoints.
!1     if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2))) then ! Do-while loop: we keep returning here.
!        if(f2.lt.f1)then ! One possible outcome,
!          x0=x1 ! its housekeeping,
!          x1=x2
!          x2=R*x1+C*x3
!          f1=f2
!          f2=f(x2) ! and a new function evaluation.
!        else ! The other outcome,
!          x3=x2
!          x2=x1
!          x1=R*x2+C*x0
!          f2=f1
!          f1=f(x1) ! and its new function evaluation.
!        endif
!        goto 1 ! Back to see if we are done.
!      endif
!      if(f1.lt.f2)then ! We are done. Output the best of the two current values.
!        golden=f1
!        xmin=x1
!      else
!        golden=f2
!        xmin=x2
!      endif
!      return
!!--------------------------------------------------------------------
!      END FUNCTION golden
!!--------------------------------------------------------------------
!END MODULE numerical_recipes
