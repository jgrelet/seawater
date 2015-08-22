// seawater_test
package seawater

import (
	_ "fmt"
	"math"
	"testing"
)

const S1, T1, P1 = 35.0, 30.0, 0.0
const D1, Sigma_t1, Sigma_teta1, Svel1, Potential_Temp1, Specific_anomaly1, Depth1 = 1021.729, 21.729, 21.729, 1545.595, 30.000, 6.071e-06, 0.000
const Lat = 4.0

const S2, T2, P2 = 34.67, 2.48, 10035.0
const D2, Sigma_t2, Sigma_teta2, Svel2, Potential_Temp2, Specific_anomaly2, Depth2 = 1070.136, 27.668, 27.764, 1633.179, 1.242, 8.352e-07, 9758.558

const S3, T3, P3, C3 = 35.1554, 26.9900 * 1.00024, 27.000, 5.538891
const S4, T4, P4, C4 = 35.7918, 18.1986 * 1.00024, 71.000, 4.705818

var T90 = []float64{9.997601, 14.996401, 19.995201, 29.992802}
var T68 = []float64{10.0, 15.0, 20.0, 30.0}

func TestSw_T68conv(t *testing.T) {

	for i, v := range T90 {
		v = sw_T68conv(v)
		v = toFixed(v, 0)
		if v != T68[i] {
			t.Errorf("Expected %f, got %f", T68[i], v)
		}
	}
}

func TestSw_T90conv(t *testing.T) {

	for i, v := range T68 {
		v = sw_T90conv(v, "T68")
		v = toFixed(v, 6)
		if v != T90[i] {
			t.Errorf("Expected %f, got %f", T90[i], v)
		}
	}
	v := sw_T90conv(20.0, "T48")
	v = toFixed(v, 14)
	if v != 19.988162840918179 {
		t.Errorf("Expected %f, got %f", 19.988162840918179, v)
	}
}

//func TestSw_smow(t *testing.T) {
//}

func TestSw_dens(t *testing.T) {

	v := sw_dens(S1, T1, P1)
	v = toFixed(v, 3)
	if v != D1 {
		t.Errorf("Expected %f, got %f", D1, v)
	}
	v = sw_dens(S2, T2, P2)
	v = toFixed(v, 3)
	if v != D2 {
		t.Errorf("Expected %f, got %f", D2, v)
	}
}

func TestSw_sal(t *testing.T) {

	v := sw_sal(C3, T3, P3)
	v = toFixed(v, 4)
	if v != S3 {
		t.Errorf("Expected %f, got %f", S3, v)
	}
	v = sw_sal(C4, T4, P4)
	v = toFixed(v, 4)
	if v != S4 {
		t.Errorf("Expected %f, got %f", S4, v)
	}
}

func TestSw_sigmat(t *testing.T) {

	v := sw_sigmat(S1, T1, P1)
	v = toFixed(v, 3)
	if v != Sigma_t1 { // 21.729
		t.Errorf("Expected %f, got %f", Sigma_t1, v)
	}
	v = sw_sigmat(S2, T2, P2)
	v = toFixed(v, 3)
	if v != Sigma_t2 { //  27.668
		t.Errorf("Expected %f, got %f", Sigma_t2, v)
	}
}

func TestSw_sigmateta(t *testing.T) {

	v := sw_sigmateta(S1, T1, P1)
	v = toFixed(v, 3)
	if v != Sigma_teta1 {
		t.Errorf("Expected %6.3f, got %6.3f", Sigma_teta1, v)
	}
	v = sw_sigmateta(S2, T2, P2)
	v = toFixed(v, 3)
	if v != Sigma_teta2 {
		t.Errorf("Expected %6.3f, got %6.3f", Sigma_teta2, v)
	}
}

func TestSw_svel(t *testing.T) {

	v := sw_svel(S1, T1, P1)
	v = toFixed(v, 3)
	if v != Svel1 {
		t.Errorf("Expected %6.3f, got %6.3f", Svel1, v)
	}
	v = sw_svel(S2, T2, P2)
	v = toFixed(v, 3)
	if v != Svel2 {
		t.Errorf("Expected %6.3f, got %6.3f", Svel2, v)
	}
}

func TestSw_ptmp(t *testing.T) {

	v := sw_ptmp(S1, T1, P1, 0)
	v = toFixed(v, 3)
	if v != Potential_Temp1 {
		t.Errorf("Expected %6.3f, got %6.3f", Potential_Temp1, v)
	}
	v = sw_ptmp(S2, T2, P2, 0)
	v = toFixed(v, 3)
	if v != Potential_Temp2 {
		t.Errorf("Expected %6.3f, got %6.3f", Potential_Temp2, v)
	}
}

func TestSw_svan(t *testing.T) {

	v := sw_svan(S1, T1, P1)
	v = toFixed(v, 9)
	if v != Specific_anomaly1 {
		t.Errorf("Expected %6.6g, got %6.6g", Specific_anomaly1, v)
	}
	v = sw_svan(S2, T2, P2)
	v = toFixed(v, 10)
	if v != Specific_anomaly2 {
		t.Errorf("Expected %6.6g, got %6.6g", Specific_anomaly2, v)
	}
}

func TestSw_dpth(t *testing.T) {

	v := sw_dpth(P1, Lat)
	v = toFixed(v, 3)
	if v != Depth1 {
		t.Errorf("Expected %6.3f, got %6.3f", Depth1, v)
	}
	v = sw_dpth(P2, Lat)
	v = toFixed(v, 3)
	if v != Depth2 {
		t.Errorf("Expected %6.3f, got %6.3f", Depth2, v)
	}
}

// I'm just starting in Go and found it surprising that it has neither a
// "toFixed" function (as in JavaScript), which would accomplish what you want,
// nor even a "round" function.
// I picked up a one-liner round function from elsewhere, and also made
// toFixed() which depends on round():
// from http://stackoverflow.com/
// How can we truncate float64 type to a particular precision in golang?
// Usage:
// fmt.Println(toFixed(1.2345678, 0))  // 1.0
// fmt.Println(toFixed(1.2345678, 1))  // 1.2
// fmt.Println(toFixed(1.2345678, 2))  // 1.23
// fmt.Println(toFixed(1.2345678, 3))  // 1.235 (rounded up)

func round(num float64) int {
	return int(num + math.Copysign(0.5, num))
}

func toFixed(num float64, precision int) float64 {
	output := math.Pow(10, float64(precision))
	return float64(round(num*output)) / output
}
