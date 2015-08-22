// Package seawater contains utility functions  to compute oceanographic parameters.
package seawater

import (
	"errors"
	"math"
)

// Density of Standard Mean Ocean Water (Pure Water) using EOS 1980.
// Input:
//   temperature [°C (ITS-90)]
// Return:
// dens(t) : density  [kg m :sup:`3`]
// Examples :
// Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22.
// t = sw_T90conv([0, 0, 30, 30, 0, 0, 30, 30])
// sw.smow(t)
// array([ 999.842594  ,  999.842594  ,  995.65113374,  995.65113374,
//         999.842594  ,  999.842594  ,  995.65113374,  995.65113374])
// References :
// Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
// computation of fundamental properties of seawater. UNESCO Tech. Pap. in
// Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
// http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf
// Millero, F.J. and  Poisson, A. International one-atmosphere equation
// of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
// doi:10.1016/0198-0149(81)90122-9
func sw_smow(T float64) float64 {
	const a0, a1, a2, a3, a4, a5 = 999.842594, 6.793952e-2, -9.095290e-3, 1.001685e-4,
		-1.120083e-6, 6.536332e-9

	dens := a0 + (a1+(a2+(a3+(a4+a5*T)*T)*T)*T)*T
	return dens
}

// Convert ITS-90 temperature to IPTS-68
// T68  = T90 * 1.00024
// Parameters:
//   array of float64 temperature [°C (ITS-90)]
// Return:
//   array of float64 temperature [°C (IPTS-68)]
// Notes:
// The International Practical Temperature Scale of 1968 (IPTS-68) need to be
// correct to the ITS-90. This linear transformation is accurate within
// 0.5 °C for conversion between IPTS-68 and ITS-90 over the
// oceanographic temperature range.
// Examples:
//   T68conv(19.995201151723585)
//   20.0
// References:
// Saunders, P. M., 1991: The International Temperature Scale of 1990,
// ITS-90. WOCE Newsletter, No. 10, WOCE International Project Office,
// Southampton, United Kingdom, 10.
func sw_T68conv(T90 float64) float64 {
	return T90 * 1.00024
}

//func sw_T68conv(T90 ...float64) []float64 {
//	out := make([]float64, len(T90))
//	for i, v := range T90 {
//		out[i] = v * 1.00024
//	}
//	return out
//}

// Convert IPTS-68 or IPTS-48 to temperature to ITS-90.
// T48 apply to all data collected prior to 31/12/1967.
// T68 apply to all data collected between 01/10/1968 and 31/12/1989.
// T90 = T68 / 1.00024
// T90 = T48 - (4.4e-6) * T48 * (100-T48) ) / 1.00024
// Parameters:
//   temperature [°C (IPTS-68) or (IPTS-48)]
//   t_type : string, optional 'T68' (default) or 'T48'
// Return:
// T90 : array_like temperature [°C (ITS-90)]
// Notes :
// The International Practical Temperature Scale of 1968 (IPTS-68) need to be
// correct to the ITS-90. This linear transformation is accurate within
// 0.5 °C for conversion between IPTS-68 and ITS-90 over the
// oceanographic temperature range.
// Examples:
// sw_T90conv(20.004799999999999)
// 20.0
// sw_T90conv(20., 'T48')
// 19.988162840918179
// References:
// Saunders, P. M., 1991: The International Temperature Scale of 1990,
// ITS-90. WOCE Newsletter, No. 10, WOCE International Project Office,
// Southampton, United Kingdom, 10.
// International Temperature Scales of 1948, 1968 and 1990, an ICES
// note, available from http://www.ices.dk/ocean/procedures/its.htm
func sw_T90conv(t float64, t_type string) float64 {
	var T90 float64
	if t_type == "T68" {
		T90 = t / 1.00024
	} else if t_type == "T48" {
		T90 = (t - 4.4e-6*t*(100-t)) / 1.00024
	} else {
		errors.New(`Unrecognized temperature type.  Try "T68" or "T48"`)
	}
	return T90
}

// Density of Sea Water at atmospheric pressure using
// UNESCO 1983 (EOS 1980) polynomial.
// Parameters:
// S = salinity   [psu (PSS-78)]
// T = temperature [degree C (ITS-90)]
// Return:
// density  [kg/m^3] of salt water with properties S,T
func sw_dens0(S, T float64) float64 {
	const b0, b1, b2, b3, b4 = 8.24493e-1, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9
	const c0, c1, c2 = -5.72466e-3, +1.0227e-4, -1.6546e-6
	const d0 = 4.8314e-4

	dens := sw_smow(T) + (b0+(b1+(b2+(b3+b4*T)*T)*T)*T)*S + (c0+(c1+c2*T)*T)*S*math.Sqrt(S) + d0*S*S
	return dens
}

// Secant Bulk Modulus (K) of Sea Water using Equation of state 1980.
// UNESCO polynomial implementation.
// Parameters:
// S = salinity    [psu      (PSS-78)]
// T = temperature [degree C (ITS-90)]
// P = pressure    [db]
// Return:
// Secant Bulk Modulus  [bars]
func sw_seck(S, T, P float64) float64 {

	P = P / 10.0 // convert from db to atmospheric pressure units

	// Pure water terms of the secant bulk modulus at atmos pressure.
	// UNESCO eqn 19 p 18

	const h3, h2, h1, h0 = -5.77905E-7, +1.16092E-4, +1.43713E-3, +3.239908 //[-0.1194975]

	AW := h0 + (h1+(h2+h3*T)*T)*T

	const k2, k1, k0 = 5.2787E-8, -6.12293E-6, +8.50935E-5 //[+3.47718E-5]

	BW := k0 + (k1+k2*T)*T

	const e4, e3, e2, e1, e0 = -5.155288E-5, +1.360477E-2, -2.327105, +148.4206, 19652.21
	//[-1930.06]

	KW := e0 + (e1+(e2+(e3+e4*T)*T)*T)*T // eqn 19

	//--------------------------------------------------------------------
	// SEA WATER TERMS OF SECANT BULK MODULUS AT ATMOS PRESSURE.
	//--------------------------------------------------------------------
	const j0 = 1.91075E-4

	const i2, i1, i0 = -1.6078E-6, -1.0981E-5, 2.2838E-3

	SR := math.Sqrt(S)

	A := AW + (i0+(i1+i2*T)*T+j0*SR)*S

	const m2, m1, m0 = 9.1697E-10, +2.0816E-8, -9.9348E-7

	B := BW + (m0+(m1+m2*T)*T)*S // eqn 18

	const f3, f2, f1, f0 = -6.1670E-5, +1.09987E-2, -0.603459, +54.6746

	const g2, g1, g0 = -5.3009E-4, +1.6483E-2, +7.944E-2

	K0 := KW + (f0+(f1+(f2+f3*T)*T)*T+
		(g0+(g1+g2*T)*T)*SR)*S // eqn 16

	K := K0 + (A+B*P)*P // eqn 15
	return K
}

// Sound Velocity in sea water using UNESCO 1983 polynomial.
// Parameters:
// S = salinity    [psu      (PSS-78)]
// T = temperature [degree C (ITS-90)]
// P = pressure    [db]
// Return:
// sound velocity  [m/s]
func sw_svel(S, T, P float64) float64 {

	P = P / 10 // convert db to bars as used in UNESCO routines

	//------------
	// eqn 34 p.46
	//------------
	const c00, c01, c02, c03, c04, c05 = 1402.388, 5.03711, -5.80852e-2, 3.3420e-4, -1.47800e-6, 3.1464e-9

	const c10, c11, c12, c13, c14 = 0.153563, 6.8982e-4, -8.1788e-6, 1.3621e-7, -6.1185e-10

	const c20, c21, c22, c23, c24 = 3.1260e-5, -1.7107e-6, 2.5974e-8, -2.5335e-10, 1.0405e-12

	const c30, c31, c32 = -9.7729e-9, 3.8504e-10, -2.3643e-12

	Cw := c00 + c01*T + c02*T*T + c03*T*T*T + c04*T*T*T*T + c05*T*T*T*T*T +
		(c10+c11*T+c12*T*T+c13*T*T*T+c14*T*T*T*T)*P +
		(c20+c21*T+c22*T*T+c23*T*T*T+c24*T*T*T*T)*
			P*P + (c30+c31*T+c32*T*T)*P*P*P

	//------------
	// eqn 35. p.47
	//-------------
	const a00, a01, a02, a03, a04 = 1.389, -1.262e-2, 7.164e-5, 2.006e-6, -3.21e-8

	const a10, a11, a12, a13, a14 = 9.4742e-5, -1.2580e-5, -6.4885e-8, 1.0507e-8, -2.0122e-10

	const a20, a21, a22, a23 = -3.9064e-7, 9.1041e-9, -1.6002e-10, 7.988e-12

	const a30, a31, a32 = 1.100e-10, 6.649e-12, -3.389e-13

	A := a00 + a01*T + a02*T*T + a03*T*T*T + a04*T*T*T*T + (a10+a11*T+a12*T*T+a13*T*T*T+a14*T*T*T*T)*
		P + (a20+a21*T+a22*T*T+a23*T*T*T)*P*
		P + (a30+a31*T+a32*T*T)*P*P*P

	//------------
	// eqn 36 p.47
	//------------
	const b00, b01, b10, b11 = -1.922e-2, -4.42e-5, 7.3637e-5, 1.7945e-7

	B := b00 + b01*T + (b10+b11*T)*P

	//------------
	// eqn 37 p.47
	//------------
	const d00, d10 = 1.727e-3, -7.9836e-6

	D := d00 + d10*P

	//------------
	// eqn 33 p.46
	//------------
	svel := Cw + A*S + B*S*math.Sqrt(S) + D*S*S
	return svel
}

// Calculates depth in metres from pressure in dbars.
// Parameters:
// P = pressure    [db]
// LAT = Latitude in decimal degress north [-90..+90]
// Return:
// depth [metres]
func sw_dpth(P, LAT float64) float64 {

	const c1, c2, c3, c4, gam_dash = +9.72659, -2.2512E-5, +2.279E-10, -1.82E-15, 2.184e-6
	DEG2RAD := (math.Pi / 180)

	LAT = math.Abs(LAT)
	X := math.Sin(LAT * DEG2RAD) // convert to radians
	X = X * X
	bot_line := 9.780318*(1.0+(5.2788E-3+2.36E-5*X)*X) + gam_dash*0.5*P
	top_line := (((c4*P+c3)*P+c2)*P + c1) * P
	DEPTHM := top_line / bot_line
	return DEPTHM
}

// Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
// Parameters:
// S = salinity    [psu      (PSS-78)]
// T = temperature [degree C (ITS-90)]
// P = pressure    [db]
// Return:
// adiabatic temperature gradient [degree_C/db]
func sw_adtg(S, T, P float64) float64 {

	const a0, a1, a2, a3 = 3.5803E-5, +8.5258E-6, -6.836E-8, 6.6228E-10

	const b0, b1 = +1.8932E-6, -4.2393E-8

	const c0, c1, c2, c3 = +1.8741E-8, -6.7795E-10, +8.733E-12, -5.4481E-14

	const d0, d1 = -1.1351E-10, 2.7759E-12

	const e0, e1, e2 = -4.6206E-13, +1.8676E-14, -2.1687E-16

	ADTG := a0 + (a1+(a2+a3*T)*T)*T + (b0+b1*T)*(S-35) +
		((c0+(c1+(c2+c3*T)*T)*T)+(d0+d1*T)*
			(S-35))*P + (e0+(e1+e2*T)*T)*P*P
	return ADTG
}

// Calculates potential temperature as per UNESCO 1983 report.
// Parameters:
// S = salinity    [psu      (PSS-78)]
// T = temperature [degree C (ITS-90)]
// P = pressure    [db]
// Return:
// reference pressure  [db]
func sw_ptmp(S, T, P, PR float64) float64 {

	// theta1
	del_P := PR - P
	del_th := del_P * sw_adtg(S, T, P)
	th := T + 0.5*del_th
	q := del_th

	// theta2
	del_th = del_P * sw_adtg(S, th, P+0.5*del_P)
	th = th + (1-1/math.Sqrt(2))*(del_th-q)
	q = (2-math.Sqrt(2))*del_th + (-2+3/math.Sqrt(2))*q

	// theta3
	del_th = del_P * sw_adtg(S, th, P+0.5*del_P)
	th = th + (1+1/math.Sqrt(2))*(del_th-q)
	q = (2+math.Sqrt(2))*del_th + (-2-3/math.Sqrt(2))*q

	// theta4
	del_th = del_P * sw_adtg(S, th, P+del_P)
	PT := th + (del_th-2*q)/6
	return PT
}

// Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.
// Parameters:
// S = salinity    [psu      (PSS-78)]
// T = temperature [degree C (ITS-90)]
// P = pressure    [db]
// Return:
// density  [kg/m^3]
func sw_dens(S, T, P float64) float64 {

	densP0 := sw_dens0(S, T)
	K := sw_seck(S, T, P)
	P = P / 10.0 // convert from db to atm pressure unit
	dens := densP0 / (1 - P/K)
	return dens
}

// Calculates depth in metres from pressure in dbars.
// Parameters:
// S = salinity    [psu      (PSS-78)]
// T = temperature [degree C (ITS-90)]
// P = pressure    [db]
// Return:
// svel = sound velocity  [m/s]
func sw_sigmat(S, T, P float64) float64 {

	dens := sw_dens(S, T, 0)
	return dens - 1000.0
}

// Calculates depth in metres from pressure in dbars.
// Parameters:
// S = salinity    [psu      (PSS-78)]
// T = temperature [degree C (ITS-90)]
// P = pressure    [db]
// Return:
// svel = sound velocity  [m/s]
func sw_sigmateta(S, T, P float64) float64 {

	dens := sw_dens(S, sw_ptmp(S, T, P, 0), 0)
	return dens - 1000.0
}

// Specific Volume Anomaly calculated as
// svan = 1/sw_dens(s,t,p) - 1/sw_dens(35,0,p)
// Parameters:
// S = salinity    [psu      (PSS-78)]
// T = temperature [degree C (ITS-90)]
// P = pressure    [db]
// Return:
// Specific Volume Anomaly  [m^3 kg^-1]
func sw_svan(S, T, P float64) float64 {

	svan := (1 / sw_dens(S, T, P)) - (1 / sw_dens(35, 0, P))
	return svan
}

// Calculates salinity
// Parameters:
// C = conductivity []
// T = temperature [degree C (ITS-90)]
// P = pressure    [db]
// Return:
// salinity
func sw_sal(C, T, P float64) float64 {

	// salinity constants
	const A1, A2, A3, B1, B2, B3, B4, C0, C1, C2, C3, C4 = 2.070E-5, -6.370E-10, 3.989E-15, 3.426E-2,
		4.464E-4, 4.215E-1, -3.107E-3, 6.766097E-1,
		2.00564E-2, 1.104259E-4, -6.9698E-7, 1.0031E-9

		// constants for salinity calculation
	a := [6]float64{0.0080, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081}
	b := [6]float64{0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144}

	if C <= 0 {
		return 0
	}
	C = (C * 10) / 42.914 // convert Siemens/meters to mmhos/cm
	r := 1 + (P*(A1+P*(A2+P*A3)))/(1+B1*T+B2*T*T+B3*C+B4*C*T)
	r = C / (r * (C0 + T*(C1+T*(C2+T*(C3+T*C4)))))

	var sum1, sum2 float64 = 0, 0
	var t float64
	for i, v := range a {
		t = math.Pow(r, (float64(i) / 2.0))
		sum1 += v * t
		sum2 += b[i] * t
	}
	T -= 15
	return (sum1 + sum2*T/(1+0.0162*T))
}
