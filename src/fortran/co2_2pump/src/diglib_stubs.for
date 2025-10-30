C diglib_stubs.f — minimal no-op DIGLIB replacements

      SUBROUTINE BEEP()
      RETURN
      END

      SUBROUTINE BGNPLT()
      RETURN
      END

      SUBROUTINE CLEAR_TEXT()
      RETURN
      END

      SUBROUTINE CURVE(X, Y, N)
      INTEGER N
      REAL X(*), Y(*)
      RETURN
      END

      SUBROUTINE DEVSEL(IDEV)
      INTEGER IDEV
      RETURN
      END

      SUBROUTINE ENDPLT()
      RETURN
      END

      INTEGER FUNCTION GET_CHAR()
      GET_CHAR = -1
      RETURN
      END

      LOGICAL FUNCTION GOODCS()
      GOODCS = .TRUE.
      RETURN
      END

      SUBROUTINE GSCOLR(ICOL)
      INTEGER ICOL
      RETURN
      END

      SUBROUTINE GSDRAW(X, Y)
      REAL X, Y
      RETURN
      END

      SUBROUTINE GSLTYP(ITYPE)
      INTEGER ITYPE
      RETURN
      END

      SUBROUTINE GSMOVE(X, Y)
      REAL X, Y
      RETURN
      END

      SUBROUTINE GSPSTR(X, Y, TEXT)
      REAL X, Y
      CHARACTER*(*) TEXT
      RETURN
      END

      SUBROUTINE INIT_MENU()
      RETURN
      END

      SUBROUTINE INIT_TEXT()
      RETURN
      END

      SUBROUTINE MAPIT(X, Y)
      REAL X, Y
      RETURN
      END

      SUBROUTINE MAPSET(X0, X1, Y0, Y1)
      REAL X0, X1, Y0, Y1
      RETURN
      END

      SUBROUTINE MINMAX(A, N, AMIN, AMAX)
      INTEGER N
      REAL A(*), AMIN, AMAX
      IF (N .GT. 0) THEN
         AMIN = A(1)
         AMAX = A(1)
      ENDIF
      RETURN
      END

      SUBROUTINE NAMDEV(NAME)
      CHARACTER*(*) NAME
      RETURN
      END

      SUBROUTINE RLSDEV()
      RETURN
      END

      SUBROUTINE SCALE2(SX, SY)
      REAL SX, SY
      RETURN
      END

      SUBROUTINE SELDEV(IDEV)
      INTEGER IDEV
      RETURN
      END
