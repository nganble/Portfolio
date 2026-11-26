#ifndef TEST_CASES_H
#define TEST_CASES_H

// =====================================================================
// Test cases: coefficients rho, m (COEFF) and exact solutions (SOLUTION)
//   Choose at compile time, e.g.  -DCOEFF=VAR -DSOLUTION=POLY
//   COEFF    : CONS (isotropic, rho = 1, m = I) | VAR (rho, m variable; m symmetric)
//   SOLUTION : POLY (manufactured) | TRIG | PULSE | PULSETRIG (plane waves, f = 0)
//
// Note: this file is included INSIDE main(), after the domain [xa,xb] x [ya,yb]
//       is defined, because omega is evaluated at (xa,ya). It defines
//       RHO, M11, M12, M21, M22 (+ derivatives), UTRUE, UTRUET, FORCE,
//       and the names coeffName, msName used in the output.
// =====================================================================

	// Define the option for coefficients
	#define CONS 1
	#define VAR   2
	#ifndef COEFF
		#define COEFF CONS
		// #define COEFF VAR
	#endif
	
	#if COEFF == CONS
		const char *coeffName = "const";
		#define RHO(x, y)   (1.0)
		#define M11(x, y)   (1.0)
		#define M11x(x, y)  (0.0)
		#define M12(x, y)   (0.0)
		#define M12x(x, y)  (0.0)
		#define M21(x, y)   (0.0)
		#define M21y(x, y)  (0.0)
		#define M22(x, y)   (1.0)
		#define M22y(x, y)  (0.0)

	#elif COEFF == VAR
		const char *coeffName = "var";
		// symmetric (m12 = m21) and positive definite on (-1,1)^2
		const Real a11 = 0.3;
		const Real a12 = 0.2;
		const Real a22 = 0.6;
		#define RHO(x, y)  (0.1 * (sin((x) + (y)) + 1.1))
		#define M11(x, y)  (a11 * ((x) + (y)) + 2.0)
		#define M11x(x, y) (a11)
		#define M12(x, y)  (0.5 + a12 * ((x) + (y)))
		#define M12x(x, y) (a12)
		#define M21(x, y)  (M12(x, y))
		#define M21y(x, y) (a12)
		#define M22(x, y)  (a22 * (2.0 * (x) + (y)) + 3.0)
		#define M22y(x, y) (a22)

	#endif

	// Define ms solution
	#define TRIG  1
	#define POLY  2
	#define PULSE 3
	#define PULSETRIG 4
	#ifndef SOLUTION
		#define SOLUTION POLY
		// #define SOLUTION TRIG
		// #define SOLUTION PULSE
		// #define SOLUTION PULSETRIG
	#endif

	#if SOLUTION == TRIG
		const char *msName = "trig";
		const Real kx = 1;
		const Real ky = 1;
		const Real omega = sqrt(M11(xa,ya)*kx*kx + (M12(xa,ya) + M21(xa,ya))*kx*ky +
														M22(xa,ya)*ky*ky)/sqrt(RHO(xa,ya));
		#define UTRUE(x, y, t)   (3.0 * cos(kx * (x) + ky * (y) - omega * (t)))
		#define UTRUET(x, y, t)  (3.0 * omega * sin(kx * (x) + ky * (y) - omega * (t)))
		#define FORCE(x, y, t)   (0.0)
	#elif SOLUTION == POLY
		const char *msName = "poly";
		// Coefficients for manufactured solution
		static const Real a0 = .5, a1 = .7, a2 = .9;  // Time coefficients
		static const Real b0 = .5, b1 = .8, b2 = .6;  // Space coefficients in x
		static const Real c0 = .3, c1 = .25, c2 = .2; // Space coefficients in y
		#define UTRUE(x,y,t)   ( (b0 + (x)*( b1 + (x)*b2 ))*(c0 + (y)*( c1 + (y)*c2 ))*( a0 + (t)*( a1 + (t)*a2 )) )
		#define UTRUET(x,y,t)  ( ( b0 + (x)*( b1 +(x)*b2 ))*(c0 + (y)*( c1 + (y)*c2 ))*( a1 + 2.*(t)*a2) )
		#define UTRUEX(x,y,t)  ( ( b1 + 2.*(x)*b2)*(c0 + (y)*( c1 + (y)*c2 ))*( a0 + (t)*( a1 + (t)*a2 )) )
		#define UTRUEY(x,y,t)  ( (b0 + (x)*( b1 + (x)*b2 ))*(c1 + 2.*(y)*c2)*( a0 + (t)*( a1 + (t)*a2 )) )
		#define UTRUETT(x,y,t) ( (b0 + (x)*( b1 +(x)*b2 ))*(c0 + (y)*( c1 + (y)*c2 ))*(2.*a2 + 0.*(t)) )
		#define UTRUEXY(x,y,t) ( (b1 + 2.*(x)*b2)*(c1 + 2.*(y)*c2)*( a0 + (t)*( a1 + (t)*a2 )) )
		#define UTRUEXX(x,y,t) ( ( 2.*b2 )*(c0 + (y)*( c1 + (y)*c2 ))*( a0 + (t)*( a1 + (t)*a2 )) )
		#define UTRUEYY(x,y,t) ( (b0 + (x)*( b1 + (x)*b2 ))*( 2.*c2 )*( a0 + (t)*( a1 + (t)*a2 )) )
		#define FORCE(x,y,t)   (RHO(x, y) * UTRUETT(x, y, t) -\
													 (M11(x, y) * UTRUEXX(x, y, t) + M22(x, y) * UTRUEYY(x, y, t)) - \
												 	 (M12(x, y) + M21(x, y)) * UTRUEXY(x, y, t) - \
												 	 (M11x(x, y) + M21y(x, y)) * UTRUEX(x, y, t) - \
												 	 (M12x(x, y) + M22y(x, y)) * UTRUEY(x, y, t))
	#elif SOLUTION == PULSE
		const char *msName = "pulse";
		// Gaussian pulse: u = 4 exp(-beta s^2), s = kx x + ky y - omega t
		const Real kx = 1;
		const Real ky = 1;
		const Real omega = sqrt(M11(xa,ya)*kx*kx + (M12(xa,ya) + M21(xa,ya))*kx*ky +
								M22(xa,ya)*ky*ky)/sqrt(RHO(xa,ya));
		#define PHASE(x, y, t)  (kx * (x) + ky * (y) - omega * (t))
		const Real beta = 3.0;
		#define UTRUE(x, y, t)   (4.0 * exp(-beta * PHASE(x, y, t) * PHASE(x, y, t)))
		#define UTRUET(x, y, t)  (2.0 * beta * omega * PHASE(x, y, t) * UTRUE(x, y, t))
		#define FORCE(x, y, t)   (0.0)
	#elif SOLUTION == PULSETRIG
		const char *msName = "trigPulse";
		// Modulated pulse: u = 2 exp(-beta s^2) cos(kn s), s = kx x + ky y - omega t
		const Real kx = 1;
		const Real ky = 1;
		const Real omega = sqrt(M11(xa,ya)*kx*kx + (M12(xa,ya) + M21(xa,ya))*kx*ky +
								M22(xa,ya)*ky*ky)/sqrt(RHO(xa,ya));
		#define PHASE(x, y, t)  (kx * (x) + ky * (y) - omega * (t))
		const Real beta = 2.5;
		const Real kn = 2.0;
		#define UTRUE(x, y, t)   (2.0 * exp(-beta * PHASE(x, y, t) * PHASE(x, y, t)) * cos(kn * PHASE(x, y, t)))
		#define UTRUET(x, y, t)  (2.0 * beta * omega * PHASE(x, y, t) * UTRUE(x, y, t) + \
								 kn * omega * 2.0 * exp(-beta * PHASE(x, y, t) * PHASE(x, y, t)) * sin(kn * PHASE(x, y, t)))
		#define FORCE(x, y, t)   (0.0)
	#endif

#endif
