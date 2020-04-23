/* Original program from fdm7-3-hybrid.c on 13 November 2010 */
/* Written by Hideo Aochi */
/* Refer and cite Aochi and Madariaga (BSSA, 2003) */
/* Notice at https://github.com/aochihi/FDM-AM2003 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "nr.h"
#include "nrutil.h"
#include "myfdm5-3.h"
#include "mpi.h"

#define PRM "test.prm" 

#define delta 10
#define reflect 0.001
#define f0 1 

main( int argc, char **argv )
{
int 	i, j, k, l, l1, n, is, it, ir, iw, ly;
double 	***vx0,  ***vy0,  ***vz0, 
 	***txx0, ***tyy0, ***tzz0, ***txy0, ***tyz0, ***txz0,
	***sxx,  ***syy,  ***szz,  ***sxy,  ***syz,  ***sxz,
	***fx,   ***fy,   ***fz,
        ***rho,  ***mu,   ***lam,  ***vp,   ***vs,  ***amp;

double  **wx,    **wy,    **wz,    **vmax,  **vazm,  **vinc,
        **wx2,   **wy2,   **wz2, 
	**vel,	 **seisx, **seisy, **seisz;

double  *xhypo,   *yhypo,   *zhypo,   *strike,  *dip,     *rake,    *slip, 
	*laydep,  *vp0,     *vs0,     *rho0,	*q0, 
	*xweight, *yweight, *zweight, 
	*sxbuf,   *sybuf,   *szbuf,   *rxbuf,   *rybuf,   *rzbuf,
        *sxxbuf,  *syybuf,  *szzbuf,  *sxybuf,  *syzbuf,  *sxzbuf,
        *rxxbuf,  *ryybuf,  *rzzbuf,  *rxybuf,  *ryzbuf,  *rxzbuf,
	*vxx,  *vxy,  *vxz,  *vyx,  *vyy,  *vyz,  *vzx,  *vzy,  *vzz,
	*txxx, *txxy, *txxz, *tyyx, *tyyy, *tyyz, *tzzx, *tzzy, *tzzz,
	*txyx, *txyy, *tyzy, *tyzz, *txzx, *txzz,
	*xobs, *yobs, *zobs;

double  ds,  dt,  pi,  xdum, ydum, zdum, 
	vpdum, vsdum, rhodum, 
	bx,  by,  bz,  muxy, muxz, muyz, b1, b2, weight, qtemp, 
 	pxx, pyy, pzz, pxy,  pyz,  pxz,  dsbiem, dtbiem, time,
	gamma0, g01, g02, g10, g11, g12, g20, g21, g22,
	xhypo0, yhypo0, zhypo0, mo, mw, 
	dump0,  dumpx,  dumpy,  dumpz,  dumpx2,  dumpy2,  dumpz2,
	vabs,   vdum,   wabs; 

char 	flname[50], number[5],    flname2[80],  flname3[80], 
	flname4[80], flname5[80], 
	outdir[50], srcfile1[80], srcfile2[80], srcfile3[80], 
	buf[256], buf2[256], 
	char1[30] = "surface", 	
        char2[30] = "cross", 
	char3[30] = "obs", 
	char4[30] = "maximum",
	char5[30] = "total" ;

long	npml, npmlv, npmlt;

int	***idpmlv, ***idpmlt;

int	*ixhypo, *iyhypo, *izhypo, *insrc, *ista;

int	my_rank, np, resultlength, mpmx, imp, icpu, nmax1, nmax2, 
	ISRC, IDUR, IOBS, NLAYER, ixhypo0, iyhypo0, izhypo0, 
	isend, irecv, 
	NDIM, XMIN, XMAX, YMIN, YMAX, ZMIN, TMAX ;

char	pname[MPI_MAX_PROCESSOR_NAME];

MPI_Status status;
MPI_Request req; 
FILE * 	fp1;
FILE *	fp2;
FILE *	fp3;
FILE *	fp4;
FILE *	fp5;
FILE *  fp_in0;
FILE *  fp_in1;
FILE *  fp_in2;
FILE *	fp_in3;

				/* Initiatlization of MPI */
  MPI_Init ( &argc, &argv );
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Get_processor_name(pname, &resultlength); 
/*  fprintf(stdout,"Process %d of %d on %s\n", my_rank, np, pname); */

/*  my_rank = 0;
  np = 1; 
  fprintf(stdout,"Process %d of %d \n", my_rank, np);  */

				/* Model prameters in meter, second */
  fp_in0 = fopen( PRM, "r");
  if ( fp_in0 == NULL ){
    perror ("failed at fopen 0");
    exit(1);
  }
  fscanf ( fp_in0, "%d", &NDIM);
  fscanf ( fp_in0, "%d %d", &XMIN, &XMAX);
  fscanf ( fp_in0, "%d %d", &YMIN, &YMAX);
  fscanf ( fp_in0, "%d", &ZMIN);
  fscanf ( fp_in0, "%d", &TMAX);

  fscanf ( fp_in0, "%s", outdir);
  fscanf ( fp_in0, "%s", srcfile1);
  fscanf ( fp_in0, "%s", srcfile2);
  fscanf ( fp_in0, "%s", srcfile3);
  fscanf ( fp_in0, "%lf %lf", &ds, &dt);
  fscanf ( fp_in0, "%d", &NLAYER);
  laydep = dvector(0, NLAYER-1);  
  vp0 = dvector(0, NLAYER-1);  
  vs0 = dvector(0, NLAYER-1);  
  rho0 = dvector(0, NLAYER-1);  
  q0  = dvector(0, NLAYER-1);  
  for( ly = 0; ly < NLAYER; ly++){
    fscanf ( fp_in0, "%lf %lf %lf %lf %lf", 
	&laydep[ly], &vp0[ly], &vs0[ly], &rho0[ly], &q0[ly]);
  }
  fclose( fp_in0 );

  pi = acos(-1.0);
  dump0 = 3.0*log10(1.0/reflect)/(2.0*ds*delta);

  if ( my_rank == 0 ){
    printf("\nDimension of FDM order ... %i\n", NDIM );
    printf("Parameter File ... %s\n", PRM);
    printf("Source Model based on ... %s\n", srcfile1 );
    printf("Rupture History from ... %s\n", srcfile2 );
    printf("Station Position at ... %s\n", srcfile3 );
    printf("Output directory ... %s\n", outdir );
    printf("spatial grid ds = %f[m]\n", ds);
    printf("time step dt = %f[s]\n", dt);

    printf("\nModel Region (%i:%i, %i:%i, %i:%i)\n",
        XMIN, XMAX, YMIN, YMAX, ZMIN, 0);
    printf("        (%7.2f:%7.2f, %7.2f:%7.2f, %7.2f:%7.2f) [km]\n",
        XMIN*ds/1000., XMAX*ds/1000., YMIN*ds/1000.,
        YMAX*ds/1000., ZMIN*ds/1000., 0.0);

    printf("PML absorbing boundary, dumping %f, width %d, ratio %f\n",
	dump0, delta, reflect);
    printf("structure model\n");
    printf("depth   vp      vs      rho\n");
    for ( ly = 0; ly < NLAYER; ly++ )
      printf("%7.2f %7.2f %7.2f %7.2f %7.5e\n", 
	laydep[ly], vp0[ly], vs0[ly], rho0[ly], q0[ly] );
  }

				/* source position */
  strcpy(flname, srcfile1);
  fp_in1 = fopen( flname, "r");
  if ( fp_in1 == NULL ){
    perror ("failed at fopen 1");
    exit(1);
  }
  fscanf ( fp_in1, "%lf %lf %lf", &xhypo0, &yhypo0, &zhypo0 );
  fscanf ( fp_in1, "%d", &ISRC );

  ixhypo = ivector(0, ISRC-1);
  iyhypo = ivector(0, ISRC-1);
  izhypo = ivector(0, ISRC-1);
  insrc  = ivector(0, ISRC-1);
  xhypo   = dvector(0, ISRC-1);
  yhypo   = dvector(0, ISRC-1);
  zhypo   = dvector(0, ISRC-1);
  strike  = dvector(0, ISRC-1);
  dip     = dvector(0, ISRC-1);
  rake    = dvector(0, ISRC-1);
  slip	  = dvector(0, ISRC-1);
  xweight = dvector(0, ISRC-1);
  yweight = dvector(0, ISRC-1);
  zweight = dvector(0, ISRC-1);

  for( is = 0;  is < ISRC; is++ ){
    fscanf(fp_in1, "%d %lf %lf %lf", &n, &xhypo[is], &yhypo[is], &zhypo[is] );
  }
  fclose( fp_in1 );

  if ( my_rank == 0 ) {
    printf("NUMBER OF SOURCE %d\n", ISRC);
  }
  ixhypo0 = (int)(xhypo0/ds);
  iyhypo0 = (int)(yhypo0/ds);
  izhypo0 = (int)(zhypo0/ds);
				/* source time function */
  strcpy(flname2, srcfile2);
  fp_in2 = fopen( flname2, "r");
  if ( fp_in2 == NULL ){
    perror ("failed at fopen 2");
    exit(1);
  }
  fgets(buf, 255, fp_in2);
  if( my_rank == 0 ) puts(buf);
  fgets(buf, 255, fp_in2);
  if( my_rank == 0 ) puts(buf);

  fscanf ( fp_in2, "%d", &IDUR );
  fscanf ( fp_in2, "%lf %lf", &dsbiem, &dtbiem );
  vel = dmatrix(0, ISRC-1, 0, IDUR-1);

  if ( my_rank == 0 ) {
    printf("Source duration %f sec\n", dtbiem*(IDUR-1));
    printf("fault segment %f m, %f s\n", dsbiem, dtbiem);
  }
  mo = 0.0;
  for( is = 0; is < ISRC; is++ ){
    fscanf ( fp_in2, "%d", &n);
    fscanf ( fp_in2, "%lf %lf %lf", &strike[is], &dip[is], &rake[is]);

    if(my_rank == 0) printf("%d %f %f %f %f %f %f\n", 
	is, xhypo[is], yhypo[is], zhypo[is], strike[is], dip[is], rake[is]);

    strike[is] = strike[is]/180.*pi;
    dip[is] = dip[is]/180.*pi;
    rake[is] = rake[is]/180.*pi;
    slip[is] = 0.;

    for ( it = 0; it < IDUR; it++ ) {
      fscanf ( fp_in2, "%lf", &vel[is][it]);
      slip[is] += vel[is][it]*dtbiem;
      vel[is][it] = dsbiem * dsbiem * vel[is][it]; 
    }
    mo += slip[is];
  }
  fclose( fp_in2 );

  for ( is = 0; is < ISRC; is++ ){

    ixhypo[is] = (int)(xhypo[is]/ds);
    iyhypo[is] = (int)(yhypo[is]/ds);
    izhypo[is] = (int)(zhypo[is]/ds);
    if( xhypo[is] < 0.0 ) ixhypo[is] = ixhypo[is]-1;
    if( yhypo[is] < 0.0 ) iyhypo[is] = iyhypo[is]-1;
    if( zhypo[is] < 0.0 ) izhypo[is] = izhypo[is]-1;

    xweight[is] = xhypo[is]/ds - ixhypo[is];
    yweight[is] = yhypo[is]/ds - iyhypo[is];
    zweight[is] = zhypo[is]/ds - izhypo[is];

/*    printf("%f %f %f %d %f %f %f\n", 
	ixhypo[is]*ds, iyhypo[is]*ds, izhypo[is]*ds, is,
	xweight[is], yweight[is], zweight[is]); */

    insrc[is] = 1;
    if ( ixhypo[is] > XMAX || ixhypo[is] < XMIN || iyhypo[is] > YMAX ||
         iyhypo[is] < YMIN || izhypo[is] > 0 || izhypo[is] < ZMIN ){
      if ( my_rank == 0 && slip[is] != 0.0 ) printf( 
  "Warning: Source %d (%d, %d, %d) (%f8.2, %f8.2, %f8.2)km is not included \n",
        	is+1, ixhypo[is], iyhypo[is], izhypo[is],
        	xhypo[is]/1000., yhypo[is]/1000., zhypo[is]/1000.);
	insrc[is] = 0;
      /*perror ("hypocenter outside of model region");
      exit(1); */
    }
  }

					/* station position */
  fp_in3 = fopen( srcfile3, "r");
  if ( fp_in3 == NULL ){
    perror ("failed at fopen 3");
    exit(1);
  }
  fscanf ( fp_in3, "%d", &IOBS );
  xobs = dvector(0, IOBS-1);
  yobs = dvector(0, IOBS-1);
  zobs = dvector(0, IOBS-1);
  ista = ivector(0, IOBS-1);

  for ( ir = 0; ir < IOBS; ir++ ){
//  fscanf ( fp_in3, "%s %lf %lf %lf", buf, &xobs[ir], &yobs[ir], &zobs[ir] );
    fscanf ( fp_in3, "%s %lf %lf %lf %s", buf, &xobs[ir], &yobs[ir], &zobs[ir], buf2 );

    ista[ir] = 1;
    if ( xobs[ir] <= XMIN*ds || xobs[ir] >= XMAX*ds ||
        yobs[ir] <= YMIN*ds || yobs[ir] >= YMAX*ds ){
      ista[ir] = 0;
      if ( my_rank == 0 )
        printf("Station %d (%f %f) not included\n", ir+1, xobs[ir], yobs[ir]);
    }
    if ( zobs[ir] <= ZMIN*ds || zobs[ir] > 0){
      zobs[ir] = 0.0;
      if ( my_rank == 0 ) printf("receiver depth set to zero %d\n", ir+1);
    }
  }
  fclose( fp_in3 );
  seisx = dmatrix(0, IOBS-1, 0, TMAX-1);
  seisy = dmatrix(0, IOBS-1, 0, TMAX-1);
  seisz = dmatrix(0, IOBS-1, 0, TMAX-1);

  if ( my_rank == 0 ){
    mo = dsbiem * dsbiem * mo;
    mw = (log10(mo) - 9.1)/1.5;
    /* double mw = 4.5; mo = pow(10., 1.5 * mw + 9.1); */
    printf("  Mw = %f; Mo = %e [N m] \n", mw,  mo);

  } 

						/* allocation for MPI*/
  mpmx = (XMAX - XMIN + 2*delta + 1 )/np + 1;
  nmax1 = ( YMAX - YMIN + 2 * delta + 1 ) * ( -ZMIN + delta + 1 ) * 2;
  nmax2 = ( YMAX - YMIN + 2 * delta + 1 ) * ( -ZMIN + delta + 2 ) * 2;

  vx0 = d3tensor(-1, mpmx+2, YMIN-delta-1, YMAX+delta, ZMIN-delta-1, 0); 
  vy0 = d3tensor(-1, mpmx+2, YMIN-delta-1, YMAX+delta, ZMIN-delta-1, 0);
  vz0 = d3tensor(-1, mpmx+2, YMIN-delta-1, YMAX+delta, ZMIN-delta-1, 0);
  idpmlv = i3tensor(-1, mpmx+2, YMIN-delta-1, YMAX+delta, ZMIN-delta-1, 0);
#pragma omp parallel for default (shared) private (imp, j, k)
  for ( imp = -1; imp <= mpmx +2; imp++ ){
    for ( j = YMIN-delta-1; j <= YMAX+delta; j++ ){
      for ( k = ZMIN-delta-1; k <= 0; k++ ){
	vx0[imp][j][k] = 0.0;
	vy0[imp][j][k] = 0.0;
	vz0[imp][j][k] = 0.0;
	idpmlv[imp][j][k] = -1;
      }
    }
  }

  npmlv = 0;
  for ( imp = 1; imp <= mpmx; imp++ ){
    i = imp2i( my_rank, mpmx, imp, XMIN);
    for ( j = YMIN-delta; j <= YMAX+delta; j++ ){ 
      for ( k = ZMIN-delta; k <= 0; k++ ){
        if( i < XMIN || i >= XMAX || j < YMIN || j >= YMAX || k < ZMIN ){
  	  npmlv += 1;
	  idpmlv[imp][j][k] = npmlv;
	}
      }
    }
  }
  printf( "PML BD for velocity %d %d\n", my_rank, npmlv);
  vxx = dvector(1, npmlv);
  vxy = dvector(1, npmlv);
  vxz = dvector(1, npmlv);
  vyx = dvector(1, npmlv);
  vyy = dvector(1, npmlv);
  vyz = dvector(1, npmlv);
  vzx = dvector(1, npmlv);
  vzy = dvector(1, npmlv);
  vzz = dvector(1, npmlv);
#pragma omp parallel for default (shared) private (npml)
  for ( npml = 1; npml <= npmlv; npml++ ){
        vxx[npml] = 0;
        vxy[npml] = 0;
        vxz[npml] = 0;
        vyx[npml] = 0;
        vyy[npml] = 0;
        vyz[npml] = 0;
        vzx[npml] = 0;
        vzy[npml] = 0;
        vzz[npml] = 0;
  }

  fx = d3tensor(1, mpmx, YMIN-1, YMAX+1, ZMIN-1, 1);
  fy = d3tensor(1, mpmx, YMIN-1, YMAX+1, ZMIN-1, 1);
  fz = d3tensor(1, mpmx, YMIN-1, YMAX+1, ZMIN-1, 1);
  sxx = d3tensor(1, mpmx, YMIN-1, YMAX+1, ZMIN-1, 1);
  syy = d3tensor(1, mpmx, YMIN-1, YMAX+1, ZMIN-1, 1);
  szz = d3tensor(1, mpmx, YMIN-1, YMAX+1, ZMIN-1, 1);
  sxy = d3tensor(1, mpmx, YMIN-1, YMAX+1, ZMIN-1, 1);
  sxz = d3tensor(1, mpmx, YMIN-1, YMAX+1, ZMIN-1, 1);
  syz = d3tensor(1, mpmx, YMIN-1, YMAX+1, ZMIN-1, 1);
#pragma omp parallel for default (shared) private (imp, j, k)
  for ( imp = 1; imp <= mpmx; imp++ ){
    for ( j = YMIN-1; j <= YMAX+1; j++ ){
      for ( k = ZMIN-1; k <= 1; k++ ){
	fx[imp][j][k] = 0.0;
	fy[imp][j][k] = 0.0;
	fz[imp][j][k] = 0.0;
	sxx[imp][j][k] = 0.0;
	syy[imp][j][k] = 0.0;
	szz[imp][j][k] = 0.0;
	sxy[imp][j][k] = 0.0;
	syz[imp][j][k] = 0.0;
	sxz[imp][j][k] = 0.0;
      }
    }
  }

  txx0 = d3tensor(-1, mpmx+2, YMIN-delta, YMAX+delta, ZMIN-delta, 1);
  tyy0 = d3tensor(-1, mpmx+2, YMIN-delta, YMAX+delta, ZMIN-delta, 1);
  tzz0 = d3tensor(-1, mpmx+2, YMIN-delta, YMAX+delta, ZMIN-delta, 1);
  txy0 = d3tensor(-1, mpmx+2, YMIN-delta, YMAX+delta, ZMIN-delta, 1);
  tyz0 = d3tensor(-1, mpmx+2, YMIN-delta, YMAX+delta, ZMIN-delta, 1);
  txz0 = d3tensor(-1, mpmx+2, YMIN-delta, YMAX+delta, ZMIN-delta, 1);
  idpmlt = i3tensor(-1, mpmx+2, YMIN-delta, YMAX+delta, ZMIN-delta, 1);
#pragma omp parallel for default (shared) private (imp, j, k)
  for ( imp = -1; imp <= mpmx+2; imp++ ){
    for ( j = YMIN-delta; j <= YMAX+delta; j++ ){
      for ( k = ZMIN-delta; k <= 1; k++ ){
        txx0[imp][j][k] = 0.0;
        tyy0[imp][j][k] = 0.0;
        tzz0[imp][j][k] = 0.0;
        txy0[imp][j][k] = 0.0;
        tyz0[imp][j][k] = 0.0;
        txz0[imp][j][k] = 0.0;
	idpmlt[imp][j][k] = -1;
      }
    }
  }

  npmlt = 0;
  for(imp = 1; imp <= mpmx; imp++){
    i = imp2i( my_rank, mpmx, imp, XMIN);
    for(j = YMIN-delta; j <= YMAX+delta; j++){
      for(k = ZMIN-delta; k <= 0; k++){
        if ( i < XMIN || i >= XMAX || j < YMIN || j >= YMAX || k < ZMIN ){
	  npmlt += 1;
	  idpmlt[imp][j][k] = npmlt;
	}
      }
    }
  }
  printf( "PML BD for stress %d %d\n", my_rank, npmlt);
  txxx = dvector(1, npmlt);
  tyyx = dvector(1, npmlt);
  tzzx = dvector(1, npmlt);
  txyx = dvector(1, npmlt);
  txzx = dvector(1, npmlt);
  txxy = dvector(1, npmlt);
  tyyy = dvector(1, npmlt);
  tzzy = dvector(1, npmlt);
  txyy = dvector(1, npmlt);
  tyzy = dvector(1, npmlt);
  txxz = dvector(1, npmlt);
  tyyz = dvector(1, npmlt);
  tzzz = dvector(1, npmlt);
  tyzz = dvector(1, npmlt);
  txzz = dvector(1, npmlt);
#pragma omp parallel for default (shared) private (npml)
  for ( npml = 1; npml <= npmlt; npml++ ){
        txxx[npml] = 0.0;
        tyyx[npml] = 0.0;
        tzzx[npml] = 0.0;
        txyx[npml] = 0.0;
        txzx[npml] = 0.0;
        txxy[npml] = 0.0;
        tyyy[npml] = 0.0;
        tzzy[npml] = 0.0;
        txyy[npml] = 0.0;
	tyzy[npml] = 0.0;
        txxz[npml] = 0.0;
        tyyz[npml] = 0.0;
        tzzz[npml] = 0.0;
        tyzz[npml] = 0.0;
        txzz[npml] = 0.0;
  } 

  wx = dmatrix(1, mpmx, YMIN, YMAX);
  wy = dmatrix(1, mpmx, YMIN, YMAX);
  wz = dmatrix(1, mpmx, YMIN, YMAX);
  vmax = dmatrix(1, mpmx, YMIN, YMAX);
  vazm = dmatrix(1, mpmx, YMIN, YMAX);
  vinc = dmatrix(1, mpmx, YMIN, YMAX);
#pragma omp parallel for default (shared) private (imp, j)
  for (imp = 1; imp <= mpmx; imp++ ){
    for ( j = YMIN; j <= YMAX; j++ ){
	wx[imp][j] = 0.0;
	wy[imp][j] = 0.0;
	wz[imp][j] = 0.0;
	vmax[imp][j] = 0.0;
	vazm[imp][j] = 0.0;
	vinc[imp][j] = 0.0;
    }
  }
  wx2 = dmatrix(1, mpmx, ZMIN, 0);
  wy2 = dmatrix(1, mpmx, ZMIN, 0);
  wz2 = dmatrix(1, mpmx, ZMIN, 0);
#pragma omp parallel for default (shared) private (imp, k)
  for (imp = 1; imp <= mpmx; imp++ ){
    for ( k = ZMIN; k <= 0; k++ ){
        wx2[imp][k] = 0.0;
        wy2[imp][k] = 0.0;
        wz2[imp][k] = 0.0;
    }
  }

  sxbuf = dvector(0, nmax1-1);
  sybuf = dvector(0, nmax1-1);
  szbuf = dvector(0, nmax1-1);
  rxbuf = dvector(0, nmax1-1);
  rybuf = dvector(0, nmax1-1);
  rzbuf = dvector(0, nmax1-1);
  sxxbuf = dvector(0, nmax2-1);
  syybuf = dvector(0, nmax2-1);
  szzbuf = dvector(0, nmax2-1);
  sxybuf = dvector(0, nmax2-1);
  syzbuf = dvector(0, nmax2-1);
  sxzbuf = dvector(0, nmax2-1);
  rxxbuf = dvector(0, nmax2-1);
  ryybuf = dvector(0, nmax2-1);
  rzzbuf = dvector(0, nmax2-1);
  rxybuf = dvector(0, nmax2-1);         
  ryzbuf = dvector(0, nmax2-1);
  rxzbuf = dvector(0, nmax2-1);

  rho = d3tensor(0, mpmx+1, YMIN-delta-1, YMAX+delta+1, ZMIN-delta-1, 1);
  mu  = d3tensor(0, mpmx+1, YMIN-delta-1, YMAX+delta+1, ZMIN-delta-1, 1);
  lam = d3tensor(0, mpmx+1, YMIN-delta-1, YMAX+delta+1, ZMIN-delta-1, 1);
  vp  = d3tensor(0, mpmx+1, YMIN-delta-1, YMAX+delta+1, ZMIN-delta-1, 1);
  vs  = d3tensor(0, mpmx+1, YMIN-delta-1, YMAX+delta+1, ZMIN-delta-1, 1);
  amp = d3tensor(0, mpmx+1, YMIN-delta-1, YMAX+delta+1, ZMIN-delta-1, 1);
  if ( my_rank == 0 ){
    printf ("depth   rho     vp      vs      mu      lam\n");
  }
  for ( k = ZMIN-delta-1; k <= 1; k++ ){
    zdum = k * ds/1000.;
    for ( imp = 0; imp <= mpmx+1; imp++ ){
      i = imp2i( my_rank, mpmx, imp, XMIN);
      xdum = i * ds/1000.;  
      for ( j = YMIN-delta-1; j <= YMAX+delta+1; j++ ){
        ydum = j * ds/1000.;

        rho[imp][j][k] = rho0[NLAYER-1];               /* kg/m^3 */
        vp[imp][j][k] = vp0[NLAYER-1];                 /* m/s^3 */
        vs[imp][j][k] = vs0[NLAYER-1];
	qtemp = q0[NLAYER-1];
        for (ly = 0; ly < NLAYER-1; ly++){
	  if ( zdum <= laydep[ly] && zdum > laydep[ly+1] ){
	    rho[imp][j][k] = rho0[ly];
	    vp[imp][j][k] = vp0[ly];
	    vs[imp][j][k] = vs0[ly];
	    qtemp = q0[ly];
	  }
	}

	if ( zdum > laydep[0] ){ 		/* shallow part */
	  vs[imp][j][k] = vs0[0];
	  vp[imp][j][k] = vp0[0];
	  rho[imp][j][k] = rho0[0];
	  
	  vsdum = hardrock( -zdum );		/* for Alps, Pyrenees */
	  vpdum = vsdum*sqrt(3.);
	  rhodum = 1741.*pow(vpdum/1000., 0.25);

	  if ( vsdum <= vs[imp][j][k] ) vs[imp][j][k] = vsdum;
	  if ( vpdum <= vp[imp][j][k] ) vp[imp][j][k] = vpdum;
	  if ( rhodum <= rho[imp][j][k] ) rho[imp][j][k] = rhodum;

	}

        mu[imp][j][k]  = vs[imp][j][k]*vs[imp][j][k]*rho[imp][j][k];
        lam[imp][j][k] = vp[imp][j][k]*vp[imp][j][k]*rho[imp][j][k]
                        - 2.0*mu[imp][j][k];

        if ( qtemp <= 0 ){      
          amp[imp][j][k] = 1.0; 
        } else {
          amp[imp][j][k] = exp(-pi*f0*dt/qtemp);
        }
      }
    }
    if ( my_rank == 0 ){
      printf ("%7.2f, %7.2f, %7.2f, %7.2f, %7.2e, %7.2e %7.4e \n", 
	zdum, rho[1][0][k], vp[1][0][k], vs[1][0][k], mu[1][0][k], 
	lam[1][0][k], amp[1][0][k]);
    }
  }

  if ( my_rank == 0 ) {
    printf ("\ninitialization\n"); 
    printf ("grids for each CPU = %i\n", mpmx);
    printf("nmax %i %i\n", nmax1, nmax2);
  }

  MPI_Barrier(MPI_COMM_WORLD); 

				/* iteration */
  for ( l = 1; l <= TMAX; l++ ){	 	
    time = dt * l;
						/* calculation */
				/* t = (l + 1/2) pour vitesse */
    for(imp = 1; imp <= mpmx; imp++){    	
      i = imp2i( my_rank, mpmx, imp, XMIN );
#pragma omp parallel for default (shared) private (j, k, bx, by, bz, npml, dumpx, dumpx2, dumpy, dumpy2, dumpz, dumpz2, xdum, ydum, zdum) 
      for(j = YMIN - delta; j <= YMAX + delta; j++){	
        for(k = ZMIN - delta; k <= 0; k++){

          bx = 0.5*(1.0/rho[imp][j][k] + 1.0/rho[imp+1][j][k]);
          by = 0.5*(1.0/rho[imp][j][k] + 1.0/rho[imp][j+1][k]);
          bz = 0.5*(1.0/rho[imp][j][k] + 1.0/rho[imp][j][k+1]);
						/* fixed boundary */
	  if( i >= XMAX + delta || j == YMAX+delta ){
	    vx0[imp][j][k] = 0.0;
	    vy0[imp][j][k] = 0.0;
	    vz0[imp][j][k] = 0.0;
	  } else if( idpmlv[imp][j][k] >= 1 ){  		/* PML boundary */
	    npml =  idpmlv[imp][j][k];
            dumpfactor(i, XMIN, XMAX, dump0*vp[imp][j][k], &dumpx, &dumpx2 );
            dumpfactor(j, YMIN, YMAX, dump0*vp[imp][j][k], &dumpy, &dumpy2 );
            dumpfactor(k, ZMIN, 10, dump0*vp[imp][j][k], &dumpz, &dumpz2 );

	
	    xdum = vxx[npml];
	    ydum = vxy[npml];
	    zdum = vxz[npml];
	    vxx[npml] = PMLdump (bx, dt, ds, dumpx2,
			  xdum, txx0[imp][j][k], txx0[imp+1][j][k] );
            vxy[npml] = PMLdump (bx, dt, ds, dumpy, 
			  ydum, txy0[imp][j-1][k], txy0[imp][j][k] );
            vxz[npml] = PMLdump (bx, dt, ds, dumpz, 
			  zdum, txz0[imp][j][k-1], txz0[imp][j][k] );
	    vx0[imp][j][k] = vxx[npml] + vxy[npml] + vxz[npml];

	    xdum = vyx[npml];	
	    ydum = vyy[npml];	
	    zdum = vyz[npml];	
	    vyx[npml] = PMLdump (by, dt, ds, dumpx,
			  xdum, txy0[imp-1][j][k], txy0[imp][j][k] );
	    vyy[npml] = PMLdump (by, dt, ds, dumpy2,
			  ydum, tyy0[imp][j][k], tyy0[imp][j+1][k] );
	    vyz[npml] = PMLdump (by, dt, ds, dumpz,
			  zdum, tyz0[imp][j][k-1], tyz0[imp][j][k] );
	    vy0[imp][j][k] = vyx[npml] + vyy[npml] + vyz[npml];

	    xdum = vzx[npml];
	    ydum = vzy[npml];
	    zdum = vzz[npml];
	    vzx[npml] = PMLdump (bz, dt, ds, dumpx,
			  xdum, txz0[imp-1][j][k], txz0[imp][j][k] );
	    vzy[npml] = PMLdump (bz, dt, ds, dumpy,
			  ydum, tyz0[imp][j-1][k], tyz0[imp][j][k] );
	    vzz[npml] = PMLdump (bz, dt, ds, dumpz2, 
			  zdum, tzz0[imp][j][k], tzz0[imp][j][k+1] );
	    vz0[imp][j][k] = vzx[npml] + vzy[npml] + vzz[npml];

	  } else if ( NDIM == 4 && k != 0){
            vx0[imp][j][k] += bx*fx[imp][j][k]*dt/ds
	      + staggardv4(bx, dt, ds,
	        txx0[imp-1][j][k], txx0[imp][j][k], txx0[imp+1][j][k], 
		txx0[imp+2][j][k],
		txy0[imp][j-2][k], txy0[imp][j-1][k], txy0[imp][j][k],
		txy0[imp][j+1][k],
		txz0[imp][j][k-2], txz0[imp][j][k-1], txz0[imp][j][k],
		txz0[imp][j][k+1] );
            vy0[imp][j][k] += by*fy[imp][j][k]*dt/ds
	      + staggardv4( by, dt, ds,
		txy0[imp-2][j][k], txy0[imp-1][j][k], txy0[imp][j][k],
		txy0[imp+1][j][k],
	        tyy0[imp][j-1][k], tyy0[imp][j][k], tyy0[imp][j+1][k],
		tyy0[imp][j+2][k],
		tyz0[imp][j][k-2], tyz0[imp][j][k-1], tyz0[imp][j][k],
		tyz0[imp][j][k+1] );
            vz0[imp][j][k] += bz*fz[imp][j][k]*dt/ds
              + staggardv4( bz, dt, ds,
		txz0[imp-2][j][k], txz0[imp-1][j][k], txz0[imp][j][k],
		txz0[imp+1][j][k],
		tyz0[imp][j-2][k], tyz0[imp][j-1][k], tyz0[imp][j][k],
		tyz0[imp][j+1][k],
		tzz0[imp][j][k-1], tzz0[imp][j][k], tzz0[imp][j][k+1],
		tzz0[imp][j][k+2] );
	    vx0[imp][j][k] = vx0[imp][j][k]*amp[imp][j][k];
	    vy0[imp][j][k] = vy0[imp][j][k]*amp[imp][j][k];
	    vz0[imp][j][k] = vz0[imp][j][k]*amp[imp][j][k];
	  } else {				/* mode normale (2me) */
            vx0[imp][j][k] += bx*fx[imp][j][k]*dt/ds 
   	      + staggardv2( bx, dt, ds,
		txx0[imp][j][k], txx0[imp+1][j][k],
		txy0[imp][j-1][k], txy0[imp][j][k],
		txz0[imp][j][k-1], txz0[imp][j][k] );
            vy0[imp][j][k] += by*fy[imp][j][k]*dt/ds
	      + staggardv2( by, dt, ds,
		txy0[imp-1][j][k], txy0[imp][j][k],
		tyy0[imp][j][k], tyy0[imp][j+1][k],
		tyz0[imp][j][k-1], tyz0[imp][j][k] );
            if( k == 0 ) fz[imp][j][k] = 0.0;
            vz0[imp][j][k] += bz*fz[imp][j][k]*dt/ds 
	      + staggardv2( bz, dt, ds,
		txz0[imp-1][j][k], txz0[imp][j][k],
		tyz0[imp][j-1][k], tyz0[imp][j][k],
		tzz0[imp][j][k], tzz0[imp][j][k+1] );
	    vx0[imp][j][k] = vx0[imp][j][k]*amp[imp][j][k];
            vy0[imp][j][k] = vy0[imp][j][k]*amp[imp][j][k];
            vz0[imp][j][k] = vz0[imp][j][k]*amp[imp][j][k];
          } /* fin de champ */

        } /* fin de k */
/*	printf("%d %d %d %d \n", my_rank, imp, i, j); */

      } /* fin de j */
    } /* fin d'i */
/*  printf("%d\n", my_rank);  */

    MPI_Barrier(MPI_COMM_WORLD); 
			/* communication pour synclonize */ 
				/* positive direction */
    if( my_rank != (np-1) ){
      i = 0;
      isend = 0;
      for ( imp = mpmx-1; imp <= mpmx; imp++ ){
        for( j = YMIN-delta; j <= YMAX+delta; j++ ){
          for( k = ZMIN-delta; k <= 0; k++ ){
            sxbuf[i] = vx0[imp][j][k];
            sybuf[i] = vy0[imp][j][k];
            szbuf[i] = vz0[imp][j][k];
            i = i + 1;
	    if( sxbuf[i] != 0 ) isend = 1;
	  }
        }
      }
    }

    if( my_rank != (np-1) )
      MPI_Isend(sxbuf, nmax1, MPI_DOUBLE, my_rank+1,
        1, MPI_COMM_WORLD, &req);
    if( my_rank != 0 )
      MPI_Recv(rxbuf, nmax1, MPI_DOUBLE, my_rank-1,
        1, MPI_COMM_WORLD, &status);
    if( my_rank != (np-1) )
      MPI_Wait(&req, &status);
    if( my_rank != (np-1) )
      MPI_Isend(sybuf, nmax1, MPI_DOUBLE, my_rank+1,
        2, MPI_COMM_WORLD, &req);
    if( my_rank != 0 )
      MPI_Recv(rybuf, nmax1, MPI_DOUBLE, my_rank-1,
        2, MPI_COMM_WORLD, &status);
    if( my_rank != (np-1) )
      MPI_Wait(&req, &status);

    if( my_rank != (np-1) )
      MPI_Isend(szbuf, nmax1, MPI_DOUBLE, my_rank+1,
        3, MPI_COMM_WORLD, &req);
    if( my_rank != 0 )
      MPI_Recv(rzbuf, nmax1, MPI_DOUBLE, my_rank-1,
        3, MPI_COMM_WORLD, &status);
    if( my_rank != (np-1) )
      MPI_Wait(&req, &status); 


    if( my_rank != 0 ){
      i = 0;
      irecv = 0;
      for( imp = -1; imp <= 0; imp++ ){
        for( j = YMIN-delta; j <= YMAX+delta; j++ ){
          for( k = ZMIN-delta; k <= 0; k++ ){
            vx0[imp][j][k] = rxbuf[i];
            vy0[imp][j][k] = rybuf[i];
            vz0[imp][j][k] = rzbuf[i];
	    if(rxbuf[i] != 0 ) irecv = 1;
            i = i + 1;
          }
        }
      }
    }

  MPI_Barrier(MPI_COMM_WORLD); 
				/* negative direction */
    if( my_rank != 0 ){
      i = 0;
      for ( imp = 1; imp <= 2; imp++ ){
        for( j = YMIN-delta; j <= YMAX+delta; j++ ){
          for( k = ZMIN-delta; k <= 0; k++ ){
            sxbuf[i] = vx0[imp][j][k];
            sybuf[i] = vy0[imp][j][k];
            szbuf[i] = vz0[imp][j][k];
            i = i + 1;
          }
        }
      }
    }
    if( my_rank != 0 )
      MPI_Isend(sxbuf, nmax1, MPI_DOUBLE, my_rank-1,
        4, MPI_COMM_WORLD, &req);
    if( my_rank != (np-1) )
      MPI_Recv(rxbuf, nmax1, MPI_DOUBLE, my_rank+1,
        4, MPI_COMM_WORLD, &status);
    if( my_rank != 0 )
      MPI_Wait(&req, &status);
    if( my_rank != 0 )
      MPI_Isend(sybuf, nmax1, MPI_DOUBLE, my_rank-1,
        5, MPI_COMM_WORLD, &req);
    if( my_rank != (np-1) )
      MPI_Recv(rybuf, nmax1, MPI_DOUBLE, my_rank+1,
        5, MPI_COMM_WORLD, &status);
    if( my_rank != 0 )
      MPI_Wait(&req, &status);

    if( my_rank != 0 )
      MPI_Isend(szbuf, nmax1, MPI_DOUBLE, my_rank-1,
        6, MPI_COMM_WORLD, &req);
    if( my_rank != (np-1) )
      MPI_Recv(rzbuf, nmax1, MPI_DOUBLE, my_rank+1,
        6, MPI_COMM_WORLD, &status);
    if( my_rank != 0 )
      MPI_Wait(&req, &status);

    if( my_rank != (np-1) ){
      i = 0;
      for( imp = mpmx+1; imp <= mpmx+2; imp++ ){
        for( j = YMIN-delta; j <= YMAX+delta; j++ ){
          for( k = ZMIN-delta; k <= 0; k++ ){
            vx0[imp][j][k] = rxbuf[i];
            vy0[imp][j][k] = rybuf[i];
            vz0[imp][j][k] = rzbuf[i];
            i = i + 1;
          }
        }
      }
    }

/*  printf("%d %d\n", my_rank, l);  */
    MPI_Barrier(MPI_COMM_WORLD); 
    
				/* increment of seismic moment */ 
    it = (int) (time/dtbiem); 
    if ( it < IDUR ){ 
      for ( is = 0; is < ISRC; is++ ){
	if (insrc[is] == 1 ){
          mo = vel[is][it];
          pxx = radxx(strike[is], dip[is], rake[is]);
          pyy = radyy(strike[is], dip[is], rake[is]);
          pzz = radzz(strike[is], dip[is], rake[is]);
          pxy = radxy(strike[is], dip[is], rake[is]);
          pyz = radyz(strike[is], dip[is], rake[is]);
          pxz = radxz(strike[is], dip[is], rake[is]);

	  i = ixhypo[is];
	  j = iyhypo[is];
	  k = izhypo[is];
          icpu = i2icpu( mpmx, i, XMIN);
          imp = i2imp( mpmx, i, XMIN);
          if ( my_rank == icpu){
	        sxx[imp][j][k] += vel[is][it]*pxx;
	        syy[imp][j][k] += vel[is][it]*pyy;
	        szz[imp][j][k] += vel[is][it]*pzz;
	        sxy[imp][j][k] += vel[is][it]*pxy;
	        syz[imp][j][k] += vel[is][it]*pyz;
	        sxz[imp][j][k] += vel[is][it]*pxz;
          }
        } /* fin d'insrc */
      } /* fin d'is (each source) */
    } /* fin d'if */
/*  printf ("prepared source %d %d \n", my_rank, l);   */

    MPI_Barrier(MPI_COMM_WORLD); 
					/* pour contraint */
					/* t = l + 1 */		
    for(imp = 1; imp <= mpmx; imp++){
      i = imp2i( my_rank, mpmx, imp, XMIN);
#pragma omp parallel for default (shared) private (j, k, muxy, muyz, muxz, npml, dumpx, dumpx2, dumpy, dumpy2, dumpz, dumpz2, b1, b2, xdum, ydum, zdum) 
      for(j = YMIN-delta; j <= YMAX+delta; j++){	
        for(k = ZMIN-delta; k <= 0; k++){

          muxy = 0.25*(1.0/mu[imp][j][k] + 1.0/mu[imp+1][j][k]
			+1.0/mu[imp][j+1][k] + 1.0/mu[imp+1][j+1][k]);
          muyz = 0.25*(1.0/mu[imp][j][k] + 1.0/mu[imp][j+1][k]
                        +1.0/mu[imp][j][k+1] + 1.0/mu[imp][j+1][k+1]);
          muxz = 0.25*(1.0/mu[imp][j][k] + 1.0/mu[imp+1][j][k]
			+1.0/mu[imp][j][k+1] + 1.0/mu[imp+1][j][k+1]);
	  muxy = 1.0/muxy;
	  muyz = 1.0/muyz;
	  muxz = 1.0/muxz;
	  if ( idpmlt[imp][j][k] >= 1){ 	/* PML boundary */
	    npml = idpmlt[imp][j][k];
            dumpfactor(i, XMIN, XMAX, dump0*vp[imp][j][k], &dumpx, &dumpx2 );
            dumpfactor(j, YMIN, YMAX, dump0*vp[imp][j][k], &dumpy, &dumpy2 );
            dumpfactor(k, ZMIN, 10, dump0*vp[imp][j][k], &dumpz, &dumpz2 );

	    if ( k == 0 ) {
	      b1 = 4*mu[imp][j][k]*(lam[imp][j][k]+mu[imp][j][k]);
	      b1 = b1/(lam[imp][j][k] + 2*mu[imp][j][k]) ;
	      b2 = 2*mu[imp][j][k]*lam[imp][j][k]/(lam[imp][j][k] 
		 + 2*mu[imp][j][k]);

	      xdum = txxx[npml];
	      ydum = txxy[npml];
              txxx[npml] = PMLdump(b1, dt, ds, dumpx, 
		    		xdum, vx0[imp-1][j][k], vx0[imp][j][k] );
              txxy[npml] = PMLdump(b2, dt, ds, dumpy,
				ydum, vy0[imp][j-1][k], vy0[imp][j][k] );
              txx0[imp][j][k] = txxx[npml] + txxy[npml] ;

	      xdum = tyyx[npml];
	      ydum = tyyy[npml];
	      tyyx[npml] = PMLdump(b2, dt, ds, dumpx, 
				xdum, vx0[imp-1][j][k], vx0[imp][j][k] );
              tyyy[npml] = PMLdump(b1, dt, ds, dumpy,
				ydum, vy0[imp][j-1][k], vy0[imp][j][k] );
              tyy0[imp][j][k] = tyyx[npml] + tyyy[npml];

	      tzz0[imp][j][k] = 0.;
	    } else {
	      b1 = lam[imp][j][k] + 2*mu[imp][j][k];
	      b2 = lam[imp][j][k];

	      xdum = txxx[npml];
	      ydum = txxy[npml];
	      zdum = txxz[npml];
	      txxx[npml] = PMLdump(b1, dt, ds, dumpx,
			xdum, vx0[imp-1][j][k], vx0[imp][j][k] );
	      txxy[npml] = PMLdump(b2, dt, ds, dumpy,
			ydum, vy0[imp][j-1][k], vy0[imp][j][k] );
	      txxz[npml] = PMLdump(b2, dt, ds, dumpz, 
			zdum, vz0[imp][j][k-1], vz0[imp][j][k] );
	      txx0[imp][j][k] = txxx[npml] + txxy[npml] + txxz[npml];

	      xdum = tyyx[npml];
	      ydum = tyyy[npml];
	      zdum = tyyz[npml];
  	      tyyx[npml] = PMLdump(b2, dt, ds, dumpx, 
			xdum, vx0[imp-1][j][k], vx0[imp][j][k] );
	      tyyy[npml] = PMLdump(b1, dt, ds, dumpy,
			ydum, vy0[imp][j-1][k], vy0[imp][j][k] );
	      tyyz[npml] = PMLdump(b2, dt, ds, dumpz,
			zdum, vz0[imp][j][k-1], vz0[imp][j][k] );
              tyy0[imp][j][k] = tyyx[npml] + tyyy[npml] + tyyz[npml];

	      xdum = tzzx[npml];
	      ydum = tzzy[npml];
	      zdum = tzzz[npml];
	      tzzx[npml] = PMLdump(b2, dt, ds, dumpx, 
			xdum, vx0[imp-1][j][k], vx0[imp][j][k] );
	      tzzy[npml] = PMLdump(b2, dt, ds, dumpy,
			ydum, vy0[imp][j-1][k], vy0[imp][j][k] );
	      tzzz[npml] = PMLdump(b1, dt, ds, dumpz, 
			zdum, vz0[imp][j][k-1], vz0[imp][j][k] );
              tzz0[imp][j][k] = tzzx[npml] + tzzy[npml] + tzzz[npml];

	      if ( j != YMAX+delta ){
		ydum = tyzy[npml];
		zdum = tyzz[npml];
  	        tyzy[npml] = PMLdump(muyz, dt, ds, dumpy2,
			ydum, vz0[imp][j][k], vz0[imp][j+1][k] );
	        tyzz[npml] = PMLdump(muyz, dt, ds, dumpz2, 
			zdum, vy0[imp][j][k], vy0[imp][j][k+1] );
	        tyz0[imp][j][k] = tyzy[npml] + tyzz[npml] ;
	      }

	      if ( i != XMAX+delta ){
	        xdum = txzx[npml];
	        zdum = txzz[npml];
	        txzx[npml] = PMLdump(muxz, dt, ds, dumpx2,
			xdum, vz0[imp][j][k], vz0[imp+1][j][k] );
	        txzz[npml] = PMLdump(muxz, dt, ds, dumpz2,
			zdum, vx0[imp][j][k], vx0[imp][j][k+1] );
	        txz0[imp][j][k] = txzx[npml] + txzz[npml];
	      }
	    }
	    if ( i != XMAX+delta && j != YMAX+delta ){
	      xdum = txyx[npml];
	      ydum = txyy[npml];
              txyx[npml] = PMLdump( muxy, dt, ds, dumpx2, 
			xdum, vy0[imp][j][k], vy0[imp+1][j][k] );
              txyy[npml] = PMLdump( muxy, dt, ds, dumpy2,
			ydum, vx0[imp][j][k], vx0[imp][j+1][k] );
              txy0[imp][j][k] = txyx[npml] + txyy[npml] ;
	    }

 	  } else if ( NDIM == 4 && k <= -2 ){
            txx0[imp][j][k] += 
		staggards4 (lam[imp][j][k], mu[imp][j][k], dt, ds,
		  vx0[imp-2][j][k], vx0[imp-1][j][k], vx0[imp][j][k],
		  vx0[imp+1][j][k],
		  vy0[imp][j-2][k], vy0[imp][j-1][k], vy0[imp][j][k],
		  vy0[imp][j+1][k],
		  vz0[imp][j][k-2], vz0[imp][j][k-1], vz0[imp][j][k],
		  vz0[imp][j][k+1] )
		  - sxx[imp][j][k]*dt/(ds*ds*ds) ; 
            tyy0[imp][j][k] += 
		staggards4 (lam[imp][j][k], mu[imp][j][k], dt, ds,
		  vy0[imp][j-2][k], vy0[imp][j-1][k], vy0[imp][j][k],
                  vy0[imp][j+1][k],
		  vx0[imp-2][j][k], vx0[imp-1][j][k], vx0[imp][j][k],
                  vx0[imp+1][j][k],
		  vz0[imp][j][k-2], vz0[imp][j][k-1], vz0[imp][j][k],
                  vz0[imp][j][k+1] )
		  - syy[imp][j][k]*dt/(ds*ds*ds) ; 
            tzz0[imp][j][k] += 
		  staggards4 (lam[imp][j][k], mu[imp][j][k], dt, ds,
		  vz0[imp][j][k-2], vz0[imp][j][k-1], vz0[imp][j][k],       
                  vz0[imp][j][k+1],
		  vx0[imp-2][j][k], vx0[imp-1][j][k], vx0[imp][j][k],
                  vx0[imp+1][j][k],
                  vy0[imp][j-2][k], vy0[imp][j-1][k], vy0[imp][j][k],
                  vy0[imp][j+1][k] )
		  - szz[imp][j][k]*dt/(ds*ds*ds) ; 
            txy0[imp][j][k] += staggardt4( muxy, dt, ds,
		  vx0[imp][j-1][k], vx0[imp][j][k], vx0[imp][j+1][k],
		  vx0[imp][j+2][k],
		  vy0[imp-1][j][k], vy0[imp][j][k], vy0[imp+1][j][k],
		  vy0[imp+2][j][k] )
		  - sxy[imp][j][k]*dt/(ds*ds*ds) ; 
            tyz0[imp][j][k] += staggardt4( muyz, dt, ds,
		  vy0[imp][j][k-1], vy0[imp][j][k], vy0[imp][j][k+1],
		  vy0[imp][j][k+2],
		  vz0[imp][j-1][k], vz0[imp][j][k], vz0[imp][j+1][k],
	    	  vz0[imp][j+2][k] )
		  - syz[imp][j][k]*dt/(ds*ds*ds) ; 
            txz0[imp][j][k] += staggardt4( muxz, dt, ds,
		  vx0[imp][j][k-1], vx0[imp][j][k], vx0[imp][j][k+1],
		  vx0[imp][j][k+2],
		  vz0[imp-1][j][k], vz0[imp][j][k], vz0[imp+1][j][k],
		  vz0[imp+2][j][k] )
		  - sxz[imp][j][k]*dt/(ds*ds*ds) ;
	    txx0[imp][j][k] = txx0[imp][j][k]*amp[imp][j][k];
	    tyy0[imp][j][k] = tyy0[imp][j][k]*amp[imp][j][k];
	    tzz0[imp][j][k] = tzz0[imp][j][k]*amp[imp][j][k];
	    txy0[imp][j][k] = txy0[imp][j][k]*amp[imp][j][k];
	    txz0[imp][j][k] = txz0[imp][j][k]*amp[imp][j][k];
	    tyz0[imp][j][k] = tyz0[imp][j][k]*amp[imp][j][k];
	  } else {			/* champ normale 2me */
            if ( k == 0 ){              /* surface libre */
	      b1 = 4*mu[imp][j][k]*(lam[imp][j][k]+mu[imp][j][k]);
	      b1 = b1/(lam[imp][j][k] + 2*mu[imp][j][k]);
	      b2 = 2*mu[imp][j][k]*lam[imp][j][k];
	      b2 = b2/(lam[imp][j][k] + 2*mu[imp][j][k]);
              txx0[imp][j][0] +=
                b1*( vx0[imp][j][0] - vx0[imp-1][j][0] )*dt/ds
                + b2*( vy0[imp][j][0] - vy0[imp][j-1][0] )*dt/ds;
              tyy0[imp][j][0] +=
                b2*( vx0[imp][j][0] - vx0[imp-1][j][0] )*dt/ds
                + b1*( vy0[imp][j][0] - vy0[imp][j-1][0] )*dt/ds;
              tzz0[imp][j][0] = 0;
            } else {
              txx0[imp][j][k] += staggards2 ( lam[imp][j][k], mu[imp][j][k], dt, ds,
		vx0[imp-1][j][k], vx0[imp][j][k],
	   	vy0[imp][j-1][k], vy0[imp][j][k],
		vz0[imp][j][k-1], vz0[imp][j][k] ) - sxx[imp][j][k]*dt/(ds*ds*ds);
              tyy0[imp][j][k] += staggards2 ( lam[imp][j][k], mu[imp][j][k], dt, ds,
		vy0[imp][j-1][k], vy0[imp][j][k],
		vx0[imp-1][j][k], vx0[imp][j][k],
		vz0[imp][j][k-1], vz0[imp][j][k] ) - syy[imp][j][k]*dt/(ds*ds*ds);
              tzz0[imp][j][k] += staggards2 ( lam[imp][j][k], mu[imp][j][k], dt, ds,
		vz0[imp][j][k-1], vz0[imp][j][k],
		vx0[imp-1][j][k], vx0[imp][j][k],
		vy0[imp][j-1][k], vy0[imp][j][k] ) - szz[imp][j][k]*dt/(ds*ds*ds);
	      if( j != YMAX+delta ) 
		tyz0[imp][j][k] += staggardt2( muyz, dt, ds,
		  vy0[imp][j][k], vy0[imp][j][k+1],
		  vz0[imp][j][k], vz0[imp][j+1][k] ) - syz[imp][j][k]*dt/(ds*ds*ds);
	      if( i != XMAX+delta )
                txz0[imp][j][k] += staggardt2( muxz, dt, ds,
		  vx0[imp][j][k], vx0[imp][j][k+1],
		  vz0[imp][j][k], vz0[imp+1][j][k] ) - sxz[imp][j][k]*dt/(ds*ds*ds);	
	    }
	    if( i != XMAX+delta && j != YMAX+delta )
              txy0[imp][j][k] += staggardt2 (muxy, dt, ds,
		  vx0[imp][j][k], vx0[imp][j+1][k],
		  vy0[imp][j][k], vy0[imp+1][j][k] ) - syz[imp][j][k]*dt/(ds*ds*ds);

	    txx0[imp][j][k] = txx0[imp][j][k]*amp[imp][j][k];
	    tyy0[imp][j][k] = tyy0[imp][j][k]*amp[imp][j][k];
	    tzz0[imp][j][k] = tzz0[imp][j][k]*amp[imp][j][k];
	    txy0[imp][j][k] = txy0[imp][j][k]*amp[imp][j][k];
	    txz0[imp][j][k] = txz0[imp][j][k]*amp[imp][j][k];
	    tyz0[imp][j][k] = tyz0[imp][j][k]*amp[imp][j][k];
	} /* fin de champ */

        if ( k == -1 ) {
          tyz0[imp][j][0] = -tyz0[imp][j][-1];
          txz0[imp][j][0] = -txz0[imp][j][-1];
	  tzz0[imp][j][1] = -tzz0[imp][j][-1];
        }
	if ( k == -2 ) {
	  txz0[imp][j][1] = -txz0[imp][j][-2];
	  tyz0[imp][j][1] = -tyz0[imp][j][-2];
	}

	if ( imp >= 1 && imp <= mpmx && j >= YMIN && j <= YMAX && 
	     k >= ZMIN && k <= 0 ){
	sxx[imp][j][k] = 0.0;
	syy[imp][j][k] = 0.0;
	szz[imp][j][k] = 0.0;
	sxy[imp][j][k] = 0.0;
	syz[imp][j][k] = 0.0;
	sxz[imp][j][k] = 0.0;
        }
      } /* fin de k */
    } /* fin de j */
  } /* fin d'i */
  MPI_Barrier(MPI_COMM_WORLD); 

                        /* communication pour synclonize */ 
  if( my_rank != (np-1) ){
    i = 0;
    for ( imp = mpmx-1; imp <= mpmx; imp++ ){
      for( j = YMIN-delta; j <= YMAX+delta; j++ ){
        for( k = ZMIN-delta; k <= 1; k++ ){
          sxxbuf[i] = txx0[imp][j][k];
          syybuf[i] = tyy0[imp][j][k];
          szzbuf[i] = tzz0[imp][j][k];
	  sxybuf[i] = txy0[imp][j][k];
	  syzbuf[i] = tyz0[imp][j][k];
	  sxzbuf[i] = txz0[imp][j][k];
          i = i + 1;
        }
      }
    }
  }
  if( my_rank != (np-1) )
      MPI_Isend(sxxbuf, nmax2, MPI_DOUBLE, my_rank+1,
        1, MPI_COMM_WORLD, &req);
  if( my_rank != 0 )
      MPI_Recv(rxxbuf, nmax2, MPI_DOUBLE, my_rank-1,
        1, MPI_COMM_WORLD, &status);
  if( my_rank != (np-1) )
      MPI_Wait(&req, &status);

  if( my_rank != (np-1) )
      MPI_Isend(syybuf, nmax2, MPI_DOUBLE, my_rank+1,
        2, MPI_COMM_WORLD, &req);
  if( my_rank != 0 )
      MPI_Recv(ryybuf, nmax2, MPI_DOUBLE, my_rank-1,
        2, MPI_COMM_WORLD, &status);
  if( my_rank != (np-1) )
      MPI_Wait(&req, &status);

  if( my_rank != (np-1) )
      MPI_Isend(szzbuf, nmax2, MPI_DOUBLE, my_rank+1,
        3, MPI_COMM_WORLD, &req);
  if( my_rank != 0 )
      MPI_Recv(rzzbuf, nmax2, MPI_DOUBLE, my_rank-1,
        3, MPI_COMM_WORLD, &status);
  if( my_rank != (np-1) )
      MPI_Wait(&req, &status);

  if( my_rank != (np-1) )
      MPI_Isend(sxybuf, nmax2, MPI_DOUBLE, my_rank+1,
        4, MPI_COMM_WORLD, &req);
  if( my_rank != 0 )
      MPI_Recv(rxybuf, nmax2, MPI_DOUBLE, my_rank-1,
        4, MPI_COMM_WORLD, &status);
  if( my_rank != (np-1) )
      MPI_Wait(&req, &status);

  if( my_rank != (np-1) )
      MPI_Isend(syzbuf, nmax2, MPI_DOUBLE, my_rank+1,
        5, MPI_COMM_WORLD, &req);
  if( my_rank != 0 )
      MPI_Recv(ryzbuf, nmax2, MPI_DOUBLE, my_rank-1,
        5, MPI_COMM_WORLD, &status);
  if( my_rank != (np-1) )
      MPI_Wait(&req, &status);

  if( my_rank != (np-1) )
      MPI_Isend(sxzbuf, nmax2, MPI_DOUBLE, my_rank+1,
        6, MPI_COMM_WORLD, &req);
  if( my_rank != 0 )
      MPI_Recv(rxzbuf, nmax2, MPI_DOUBLE, my_rank-1,
        6, MPI_COMM_WORLD, &status);
  if( my_rank != (np-1) )
      MPI_Wait(&req, &status);

  if( my_rank != 0 ){
    i = 0;
    for( imp = -1; imp <= 0; imp++ ){
      for( j = YMIN-delta; j <= YMAX+delta; j++ ){
        for( k = ZMIN-delta; k <= 1; k++ ){
            txx0[imp][j][k] = rxxbuf[i];
            tyy0[imp][j][k] = ryybuf[i];
            tzz0[imp][j][k] = rzzbuf[i];
            txy0[imp][j][k] = rxybuf[i];
            tyz0[imp][j][k] = ryzbuf[i];
            txz0[imp][j][k] = rxzbuf[i];
            i = i + 1;
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
				/* negative direction */
  if( my_rank != 0 ){
    i = 0;
    for ( imp = 1; imp <= 2; imp++ ){
      for( j = YMIN-delta; j <= YMAX+delta; j++ ){
        for( k = ZMIN-delta; k <= 1; k++ ){
          sxxbuf[i] = txx0[imp][j][k];
          syybuf[i] = tyy0[imp][j][k];
          szzbuf[i] = tzz0[imp][j][k];
          sxybuf[i] = txy0[imp][j][k];
          syzbuf[i] = tyz0[imp][j][k];
          sxzbuf[i] = txz0[imp][j][k];
          i = i + 1;
        }
      }
    }
  }
  if( my_rank != 0 )
      MPI_Isend(sxxbuf, nmax2, MPI_DOUBLE, my_rank-1,
        7, MPI_COMM_WORLD, &req);
  if( my_rank != (np-1) )
      MPI_Recv(rxxbuf, nmax2, MPI_DOUBLE, my_rank+1,
        7, MPI_COMM_WORLD, &status);
  if( my_rank != 0 )
      MPI_Wait(&req, &status);

  if( my_rank != 0 )
      MPI_Isend(syybuf, nmax2, MPI_DOUBLE, my_rank-1,
        8, MPI_COMM_WORLD, &req);
  if( my_rank != (np-1) )
      MPI_Recv(ryybuf, nmax2, MPI_DOUBLE, my_rank+1,
        8, MPI_COMM_WORLD, &status);
  if( my_rank != 0 )
      MPI_Wait(&req, &status);

  if( my_rank != 0 )
      MPI_Isend(szzbuf, nmax2, MPI_DOUBLE, my_rank-1,
        9, MPI_COMM_WORLD, &req);
  if( my_rank != (np-1) )
      MPI_Recv(rzzbuf, nmax2, MPI_DOUBLE, my_rank+1,
        9, MPI_COMM_WORLD, &status);
  if( my_rank != 0 )
      MPI_Wait(&req, &status);

  if( my_rank != 0 )
      MPI_Isend(sxybuf, nmax2, MPI_DOUBLE, my_rank-1,
        10, MPI_COMM_WORLD, &req);
  if( my_rank != (np-1) )
      MPI_Recv(rxybuf, nmax2, MPI_DOUBLE, my_rank+1,
        10, MPI_COMM_WORLD, &status);
  if( my_rank != 0 )
      MPI_Wait(&req, &status);

  if( my_rank != 0 )
      MPI_Isend(syzbuf, nmax2, MPI_DOUBLE, my_rank-1,
        11, MPI_COMM_WORLD, &req);
  if( my_rank != (np-1) )
      MPI_Recv(ryzbuf, nmax2, MPI_DOUBLE, my_rank+1,
        11, MPI_COMM_WORLD, &status);
  if( my_rank != 0 )
      MPI_Wait(&req, &status);

  if( my_rank != 0 )
      MPI_Isend(sxzbuf, nmax2, MPI_DOUBLE, my_rank-1,
        12, MPI_COMM_WORLD, &req);
  if( my_rank != (np-1) )
      MPI_Recv(rxzbuf, nmax2, MPI_DOUBLE, my_rank+1,
        12, MPI_COMM_WORLD, &status);
  if( my_rank != 0 )
      MPI_Wait(&req, &status);

  if( my_rank != (np-1) ){
    i = 0;
    for( imp = mpmx+1; imp <= mpmx+2; imp++ ){
      for( j = YMIN-delta; j <= YMAX+delta; j++ ){
        for( k = ZMIN-delta; k <= 1; k++ ){
            txx0[imp][j][k] = rxxbuf[i];
            tyy0[imp][j][k] = ryybuf[i];
            tzz0[imp][j][k] = rzzbuf[i];
            txy0[imp][j][k] = rxybuf[i];
            tyz0[imp][j][k] = ryzbuf[i];
            txz0[imp][j][k] = rxzbuf[i];
            i = i + 1;
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
                                        /* deplacement */
  for ( imp = 1; imp <= mpmx; imp++ ){
    for ( j = YMIN; j <= YMAX; j++ ) {
        wx[imp][j] += vx0[imp][j][0]*dt;
        wy[imp][j] += vy0[imp][j][0]*dt;
        wz[imp][j] += vz0[imp][j][-1]*dt;

        vdum = sqrt(vx0[imp][j][0]*vx0[imp][j][0] +
                    vy0[imp][j][0]*vy0[imp][j][0] +
                    vz0[imp][j][-1]*vz0[imp][j][-1] );
	if ( vdum > vmax[imp][j] ) {
	  vmax[imp][j] = vdum;
          vazm[imp][j] = 180.0*atan2(vx0[imp][j][0],vy0[imp][j][0])/pi;
          vdum = sqrt(vx0[imp][j][0]*vx0[imp][j][0] +
                    vy0[imp][j][0]*vy0[imp][j][0] );
          vinc[imp][j] = 180.0*atan2(vdum, vz0[imp][j][-1])/pi ;
	} 
    }
    for ( k = ZMIN; k <= 0; k++ ) {
        wx2[imp][k] += vx0[imp][iyhypo0][k]*dt;
        wy2[imp][k] += vy0[imp][iyhypo0][k]*dt;
	if ( k != 0 ){
          wz2[imp][k] += vz0[imp][iyhypo0][k]*dt;
	} else {
          wz2[imp][k] += vz0[imp][iyhypo0][-1]*dt;
	} 
    }
  }

  if ( my_rank == 0 && ( l % 100 == 0 ) )
 	 fprintf (stdout, "fini time %i %f\n", l, time);
  
					/* seismograms */
  for ( ir = 0; ir < IOBS; ir++ ){
    if( ista[ir] == 1 ){
      i = double2int(xobs[ir]/ds);
      j = double2int(yobs[ir]/ds);
      k = double2int(zobs[ir]/ds);

      icpu = i2icpu( mpmx, i, XMIN);
      imp = i2imp( mpmx, i, XMIN);
      if ( my_rank == icpu ){
        seisx[ir][l-1] = vx0[imp][j][k];
        seisy[ir][l-1] = vy0[imp][j][k];
        seisz[ir][l-1] = vz0[imp][j][k-1];
      }
    }
  }
					/* output */ 
  strcpy(flname, outdir);
  strcpy(flname2, outdir);
  strcpy(flname4, outdir);
                                        
  strcat(flname, char1);
  strcat(flname2, char2);
  strcat(flname4, char5);

  sprintf(number, "%5.5d", l);
  strcat(flname, number);
  strcat(flname2, number);
  strcat(flname4, number);

  sprintf(number, "%3.3d", my_rank);
  strcat(flname, number);
  strcat(flname2, number);
  strcat(flname4, number);

  if ( l == -2000 ){
    fp4 = fopen(flname4, "w"); 
    for ( imp = 1; imp <= mpmx; imp++ ){
      i = imp2i( my_rank, mpmx, imp, XMIN);
      if ( i >= XMIN && i <= XMAX ){
        for ( j = YMIN; j <= YMAX; j++ ){
	  for ( k = ZMIN; k <= 0; k++ ){
            fprintf(fp4, "%7.2f %7.2f %7.2f %8.3e %8.3e %8.3e\n",
                ds*i/1000, ds*j/1000, ds*k/1000, 
                vx0[imp][j][k], vy0[imp][j][k], vz0[imp][j][k]);
          }
        }
      }
    }
    fclose(fp4);
  }

  if ( (l % 500) == 0 ){
    fp1 = fopen(flname, "w"); 
    fp2 = fopen(flname2, "w");

    for ( imp = 1; imp <= mpmx; imp++ ){
      i = imp2i( my_rank, mpmx, imp, XMIN);
      if ( i >= XMIN && i <= XMAX ){
        for ( j = YMIN; j <= YMAX; j++ ){
	  if( ( ((int)(ds*i))%500) == 0 && ( ((int)(ds*j))%500) == 0 ){

            fprintf(fp1, "%7.2f %7.2f %8.3e %8.3e %8.3e %8.3f %8.3f %8.3f\n", 
		ds*i/1000, ds*j/1000, 
                vx0[imp][j][0], vy0[imp][j][0], vz0[imp][j][-1], 
		wx[imp][j], wy[imp][j], wz[imp][j]);
	  }
	}
      }

      for ( k = ZMIN; k <= 0; k++ ){
	if( ( ((int)(ds*i))%500) == 0 && ( ((int)(ds*k))%500) == 0 ){
	  if ( k != 0 ){
	    zdum = vz0[imp][iyhypo0][k];
	  } else {
	    zdum = vz0[imp][iyhypo0][-1];
	  }
          if ( i >= XMIN && i <= XMAX && k >= ZMIN && k <= -0)
          fprintf(fp2, "%7.2f %7.2f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
		ds*i/1000, ds*k/1000,
            vx0[imp][iyhypo0][k], vy0[imp][iyhypo0][k], zdum, 
                wx2[imp][k], wy2[imp][k], wz2[imp][k]);

        }
      }
    }
    fclose(fp1); 
    fclose(fp2);
  }

  if ( (l % 200) == 0 ){
    for ( ir = 0; ir < IOBS; ir++ ){
      if ( ista[ir] == 1 ){
        i = double2int(xobs[ir]/ds);
        j = double2int(yobs[ir]/ds);
        k = double2int(zobs[ir]/ds);
        icpu = i2icpu( mpmx, i, XMIN);
        imp = i2imp( mpmx, i, XMIN);
        if ( my_rank == icpu ){
          strcpy(flname5, outdir);
          strcat(flname5, char3);
          sprintf(number, "%4.4d", ir+1);
          strcat(flname5, number);
          strcat(flname5, ".dat");

          if ( l == TMAX) printf("%d %d %s\n", ir, icpu, flname5);
          fp5 = fopen(flname5, "w");
          fprintf(fp5, "%d %f %f %f %d %d %d %f\n", 
            ir+1, xobs[ir], yobs[ir], zobs[ir], i, j, k, vs[imp][j][k]);
          fprintf(fp5, "%d %f\n", l, dt);
          for (l1 = 0; l1 < l; l1++)
            fprintf(fp5, "%e %e %e\n", 
		seisx[ir][l1], seisy[ir][l1], seisz[ir][l1]);
          fclose(fp5);
        }
      }
    }

  } /* output */


}						/* fin de programe */
                        /* output of time series at stations */
strcpy(flname, outdir);
strcat(flname, char4);
sprintf(number, "%2.2d", my_rank);
strcat(flname, number);
fp1 = fopen(flname, "w");
for ( imp = 1; imp <= mpmx; imp++ ){
  i = imp2i( my_rank, mpmx, imp, XMIN);
  if ( i >= XMIN && i <= XMAX )
    for ( j = YMIN; j <= YMAX; j++ )
      if( ( ((int)(ds*i))%800) == 0 && ( ((int)(ds*j))%800) == 0 )
          fprintf(fp1, "%7.2f %7.2f %8.3f %8.3f %8.3f\n",
                ds*i/1000, ds*j/1000,
                vmax[imp][j], vazm[imp][j], vinc[imp][j] );
}
fclose(fp1);

MPI_Finalize(); 
return 0;
}
						/* functions */

int double2int(double x)
{
    if ( fabs(x- (int) x) < 0.5 ){
      return (int) x;
    } else {
      if ( x > 0 ){
	return ((int) x + 1);
      } else {
	return ((int) x - 1);
      }
    }
}

int imp2i(int my_rank, int mpmx, int imp, int XMIN)
{
  return XMIN - delta + my_rank*mpmx + imp - 1;
}

int i2imp(int mpmx, int i, int XMIN)
{
  if ( ( (i + 1 - XMIN + delta)%mpmx ) == 0 ){
    return mpmx;
  } else {
    return (i + 1 - XMIN + delta) % mpmx;
  }
}

int i2icpu(int mpmx, int i, int XMIN) 
{
  if ( ( (i + 1 - XMIN + delta)%mpmx ) == 0 ){
    return (i + 1 - XMIN + delta) / mpmx - 1;
  } else {
    return (i + 1 - XMIN + delta) / mpmx;
  }
}

double radxx(double strike, double dip, double rake)
{
  return  cos(rake)*sin(dip)*sin(2.*strike)
        - sin(rake)*sin(2.*dip)*cos(strike)*cos(strike) ;
}
double radyy(double strike, double dip, double rake)
{
  return  - ( cos(rake)*sin(dip)*sin(2.*strike)
        + sin(rake)*sin(2.*dip)*sin(strike)*sin(strike) );
}
double radzz(double strike, double dip, double rake)
{
  return sin(rake)*sin(2.*dip);
}
double radxy(double strike, double dip, double rake)
{
  return cos(rake)*sin(dip)*cos(2.*strike)
        + 0.5*sin(rake)*sin(2.*dip)*sin(2.*strike);
}
double radyz(double strike, double dip, double rake)
{
  return cos(rake)*cos(dip)*cos(strike) + sin(rake)*cos(2.*dip)*sin(strike);
}
double radxz(double strike, double dip, double rake)
{
  return cos(rake)*cos(dip)*sin(strike) - sin(rake)*cos(2.*dip)*cos(strike);
}

double staggardv4 (double b, double dt, double dx,
	double x1, double x2, double x3, double x4,
	double y1, double y2, double y3, double y4,
	double z1, double z2, double z3, double z4)
{
  return 9*b*dt*( (x3 - x2 ) + (y3 - y2) + (z3 - z2) )/(8*dx)
              - b*dt*( (x4 - x1) + (y4 - y1) + (z4 - z1) )/(24*dx);
}

double staggardv2 (double b, double dt, double dx,
	double x1, double x2,
	double y1, double y2, 
	double z1, double z2 )
{
  return b*dt*( (x2 - x1) + (y2 - y1) + (z2 - z1) )/dx;
}

double staggards4 (double lam, double mu, double dt, double dx,
	double x1, double x2, double x3, double x4,
	double y1, double y2, double y3, double y4,
	double z1, double z2, double z3, double z4 )
{
  return 9*dt*( (lam+2*mu)*(x3 - x2) + lam*(y3 - y2) + lam*(z3 - z2) )/(8*dx)
         - dt*( (lam+2*mu)*(x4 - x1) + lam*(y4 - y1) + lam*(z4 - z1) )/(24*dx);
}

double staggards2 (double lam, double mu, double dt, double dx,
	double x1, double x2,
	double y1, double y2,
	double z1, double z2 )
{
  return dt*( (lam+2*mu)*(x2 - x1) + lam*(y2 - y1) + lam*(z2 - z1) )/dx;
}

double staggardt4 (double mu, double dt, double dx,
	double x1, double x2, double x3, double x4,
	double y1, double y2, double y3, double y4 )
{
  return 9*dt*mu*( (x3 - x2) + (y3 - y2) )/(8*dx)
          - dt*mu*( (x4 - x1) + (y4 - y1) )/(24*dx);
}

double staggardt2 (double mu, double dt, double dx,
	double x1, double x2,
	double y1, double y2 )
{
  return dt*mu*( (x2 - x1) + (y2 - y1) )/dx;
}

double PMLdump(double b, double dt, double dx, double dump,
	double v, double x1, double x2 )
{
  return ( (2 - dt*dump)*v + b*( x2 - x1 )*(2*dt/dx) )/(2 + dt*dump);
}

void dumpfactor(int i, int XMIN, int XMAX, double dump, 
	double *dumpx, double *dumpx2 )
{
            if ( i < XMIN ){
              *dumpx = dump* (i-XMIN)*(i-XMIN)/(delta*delta);
              *dumpx2 = dump* (i+0.5-XMIN)*(i+0.5-XMIN)/(delta*delta);
            } else if ( i >= XMAX ){
              *dumpx = dump* (i-XMAX)*(i-XMAX)/(delta*delta);
              *dumpx2 = dump* (i+0.5-XMAX)*(i+0.5-XMAX)/(delta*delta);
            } else {
              *dumpx = 0.0;
              *dumpx2 = 0.0;
            }
	    return;
}

double hardrock(double z){
  if( z <= 0.75 ){
    return 3260.;
  } else if ( z <= 2.70 ){
    return 3324.*pow(z, 0.067);
  } else if ( z <= 8.0 ){
    return 3447.*pow(z, 0.0209);
  } else {
    return 0;
  }
}

double softrock(double z){
  if ( z <= 0.001 ){
    return 245.;
  } else if (z <= 0.03 ){
    return 2206.*pow(z, 0.272);
  } else if (z <= 0.19 ){
    return 3542.*pow(z, 0.407);
  } else if (z <= 4. ){
    return 2505.*pow(z, 0.199);
  } else if (z <= 8.){
    return 2927.*pow(z, 0.086);
  } else {
    return 0;
  }
}
