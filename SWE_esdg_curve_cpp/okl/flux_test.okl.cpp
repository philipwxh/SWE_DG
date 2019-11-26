#define USE_DOUBLE 1

#if USE_DOUBLE
#define dfloat double
#define dfloat4 double4
#define LOGDF log
#define POWDF pow
#define EXPDF exp
#define HALF .5
#define ONE 1.0
#define TWO 2.0
#define ZERO 0.0
#define TOL 1.e-3
#define C1 .3333333333333
#define C2 .2
#define C3 0.142857142857143
#define C4 .1111111111111111
#define TOL_ZERO 1.e-14

#else
#define dfloat float
#define dfloat4 float4
#define LOGDF logf
#define POWDF powf
#define EXPDF expf
#define HALF .5f
#define ONE 1.f
#define TWO 2.f
#define ZERO 0.f;
#define TOL 1.fe-3
#define C1 .333333333f
#define C2 .2f
#define C3 0.142857142857143f
#define C4 .1111111111111111f
#define TOL_ZERO 1.fe-14
#endif

// helpful functions
// take the average of two numbers
#define avg( a, b ) HALF * ( a + b )

void SWE2d_Vflux( const dfloat hL, const dfloat hR, const dfloat uL, 
				  const dfloat uR, const dfloat vL, const dfloat vR,
                        dfloat &FxV1, 	 dfloat &FyV1,    dfloat &FxV2,
                        dfloat &FyV2, 	 dfloat &FxV3,    dfloat &FyV3 ){
	/* inputs:
	hL: 	SWE water height on the left
	hR: 	SWE water height on the right
	uL: 	SWE water velocity in x direction on the left
	uR:  	SWE water velocity in x direction on the right
	vL:  	SWE water velocity in y direction on the left
	vR:  	SWE water velocity in y direction on the right
	FxV123: volume flux for x direction for component 1, 2 and 3
	FyV123: volume flux for y direction for component 1, 2 and 3
	other:
	g:  SWE gravitational constant
	*/
	//     h_avg:     {{h}}
	const dfloat h_avg = avg( hL, hR );
	//     u_avg:     {{u}}
	const dfloat u_avg = avg( uL, uR );
	//     v_avg:     {{v}}
	const dfloat v_avg = avg( vL, vR );
	//     half_g_h_sq_avg:  .5 *  g  *    {{h^2}}
	const dfloat half_g_h_sq_avg = p_g_half * avg( hL * hL, hR * hR );

	const dfloat g_h_avg_sq = p_g * h_avg * h_avg;
 // FxV1: {{hu}}
	FxV1 = avg( hL * uL, hR * uR );
 // FyV1: {{hv}}
	FyV1 = avg( hL * vL, hR * vR );
 // FxV2: {{hu}}*{{u}}  + g *{{h}}*{{h}} - .5g{{h^2}}
	FxV2 = FxV1 * u_avg + g_h_avg_sq - half_g_h_sq_avg;
 // FyV2: {{hv}}* {{u}}
	FyV2 = FyV1 * u_avg;
 // FyV2: {{hu}}* {{v}}
	FxV3 = FxV1 * v_avg;
 // FyV3: {{hv}}*{{v}}  + g*{{h}}*{{h}} - .5g{{h^2}}
	FyV3 = FyV1 * v_avg + g_h_avg_sq - half_g_h_sq_avg;
}