#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mat.h"
#include "matrix.h"
#include "complex.h"
#include "hdf5.h"
#include <gsl/gsl_cblas.h>
#include "nrutil.h"

#define FILE "dene.h5"

#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )

#define dist(x1, y1, z1, x2, y2, z2) sqrtf((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))

#define NX 30
#define NY 30
#define NZ 30
#define NMAX 100

/*

libraries in /home/daveh/MATLAB/R2015b/bin/glnxa64

*/

double *pr, *pim; // pointer for reading matlab files

void wtnc(fcomplex a[], unsigned long nn[], int ndim, int isign,
	  void (*wtstep)(float [], unsigned long, int));
extern void pwt(float a[], unsigned long n, int isign);

void MatrixComplexInverse(fcomplex *invA, fcomplex *A, int n);
//void sincosf(float,float *,float *);

int diagnose(const char *file) {
  MATFile *pmat;
  const char **dir;
  const char *name;
  int	  ndir, ndim;
  int	  i,j;
  mxArray *pa;
  mwSize num;

  printf("Reading file %s...\n\n", file);

  /*
   * Open file to get directory
   */
  pmat = matOpen(file, "r");
  if (pmat == NULL) {
    printf("Error opening file %s\n", file);
    return(1);
  }
  
  /*
   * get directory of MAT-file
   */
  dir = (const char **)matGetDir(pmat, &ndir);
  if (dir == NULL) {
    printf("Error reading directory of file %s\n", file);
    return(1);
  } else {
    printf("Directory of %s:\n", file);
    for (i=0; i < ndir; i++)
      printf("%s\n",dir[i]);
  }
  mxFree(dir);

  /* In order to use matGetNextXXX correctly, reopen file to read in headers. */
  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(1);
  }
  pmat = matOpen(file, "r");
  if (pmat == NULL) {
    printf("Error reopening file %s\n", file);
    return(1);
  }

  /* Read in each array. */
  printf("\nReading in the actual array contents:\n");
  for (i=0; i<ndir; i++) {
      pa = matGetNextVariable(pmat, &name);
      if (pa == NULL) {
	  printf("Error reading in file %s\n", file);
	  return(1);
      } 
      /*
       * Diagnose array pa
       */

      ndim = mxGetNumberOfDimensions(pa);

      printf("According to its contents, array %s has %d dimensions\n",
	     name, ndim);
      if (mxIsFromGlobalWS(pa))
	printf("  and was a global variable when saved\n");
      else
	printf("  and was a local variable when saved\n");      
      
      num = mxGetNumberOfElements(pa);       

      pr = mxGetPr(pa); 
      if(mxIsComplex(pa)){
	printf("  and is a complex variable\n");
	pim = mxGetPi(pa);
      }      

      // mxDestroyArray(pa); keep around
  }

  if (matClose(pmat) != 0) {
      printf("Error closing file %s\n",file);
      return(1);
  }
  printf("Done\n");
  return(num);
}



void main(void){
  hid_t       file_id; /* identifiers */
  hid_t       data_id0,data_id1,data_id2,data_id3;
  hid_t       dataspace_id1, dataspace_id2;
  hsize_t     dims[3];
  herr_t      status;
  float       *xx,*yy,*zz;

  int i, j, k, l, m, n, p=0;
  int nfreq, ntag, nrx, ndat, nmod=NX*NY*NZ;
  float wn, pi, c, pse, dist;
  float *freqs,*rxx,*rxy,*rxz,*txx,*txy,*txz;
  float *image, x, y, z, *d1, *d2, gamma;
  fcomplex *sparams, *sparams2, poff, sum, fac;
  fcomplex *rcal, *tcal;
  fcomplex *s0, *r0, *q0, *p0;
  float tk,tkold,fac1,fac2,anorm;
  fcomplex *g, *gg, *gp, *pe, *ata, *atai, zero, one, mone; 
  fcomplex *xk,*xp, calpha;
  float tnew,told,alpha,rad,rad2,arg,corr,cmax;
  int inc, imax, imaxp;
  int *lib, im;
  long pointer, pointer2;
  unsigned long nn[3]={NX,NY,NZ};
  float sinx, cosx;

  float x0=-0.5, x1=0.5, y0=-0.5, y1=0.5, z0=-0.25, z1=0.25;

  pi = 4.0*atanf(1.0);
  c = 2.998e8;

  image=(float*) calloc(NX*NY*NZ,sizeof(float));  

  // begin by reading matlab spec files

  nfreq = diagnose("freq.mat");
  freqs = (float*)malloc(sizeof(float)*nfreq);
  for(i=0;i<nfreq;i++)
    freqs[i]=pr[i];
    
  nrx = diagnose("rxPosition.mat")/3;
  rxx = (float*)malloc(sizeof(float)*nrx);
  rxy = (float*)malloc(sizeof(float)*nrx);
  rxz = (float*)malloc(sizeof(float)*nrx);
  for(i=0;i<nrx;i++){
    rxx[i]=pr[i];
    rxy[i]=pr[i+nrx];
    rxz[i]=pr[i+2*nrx];   
  }

  ntag = diagnose("tagPosition.mat")/3;
  txx = (float*)malloc(sizeof(float)*ntag);
  txy = (float*)malloc(sizeof(float)*ntag);
  txz = (float*)malloc(sizeof(float)*ntag);
  for(i=0;i<ntag;i++){
    txx[i]=pr[i];
    txy[i]=pr[i+ntag];
    txz[i]=pr[i+2*ntag];    
  }

  d1=(float*)malloc(sizeof(float)*NX*NY*NZ*ntag);
  d2=(float*)malloc(sizeof(float)*NX*NY*NZ*nrx);

  // set up cals

  rcal=(fcomplex*) malloc(sizeof(fcomplex)*nfreq*nrx);
  tcal=(fcomplex*) malloc(sizeof(fcomplex)*nfreq*ntag);

  // prefil distances

  printf("Prefilling ...\n");

  for(l=0;l<NX;l++){
    x=x0+(x1-x0)*(float)l/(float)(NX-1);
    for(m=0;m<NY;m++){
      y=y0+(y1-y0)*(float)m/(float)(NY-1);
      for(n=0;n<NZ;n++){
	z=z0+(z1-z0)*(float)n/(float)(NZ-1);
	for(i=0;i<ntag;i++)
	  d1[l+NX*(m+NY*(n+NZ*i))]=dist(x,y,z,txx[i],txy[i],txz[i]);
	for(k=0;k<nrx;k++){
	  d2[l+NX*(m+NY*(n+NZ*k))]=dist(x,y,z,rxx[k],rxy[k],rxz[k]);
	}
      }
    }
  }

  // read cal data

  ndat = diagnose("sParamNoObj.mat"); // should be ntag*nrx*nfreq
  sparams = (fcomplex*)malloc(sizeof(fcomplex)*ndat);
  for(i=0;i<ndat;i++){
    //printf("%i %i %f %f \n",i,ndat,pr[i],pim[i]);
    sparams[i]=Complex(pr[i],pim[i]);    
  }

  // calibrate receivers; cal for each freq, rx
  

  for(j=0;j<nfreq;j++){
    wn = 2.0*pi*freqs[j]/c;
    rcal[0+nrx*j]=Complex(1.0,0.0); // regard rx 0 as the standard
    for(k=1;k<nrx;k++){ 
      sum=Complex(0.0,0.0);
      for(i=0;i<ntag;i++){ // average over tags
	// anticipated relative phase
	pse = -wn*(dist(rxx[k],rxy[k],rxz[k],txx[i],txy[i],txz[i])
		   - dist(rxx[0],rxy[0],rxz[0],txx[i],txy[i],txz[i]));
	fac = Complex(cos(pse),sin(pse));
	// measured relative phase
	poff = Cmul(sparams[i+ntag*(k+nrx*j)],
		    Conjg(sparams[i+ntag*(0+nrx*j)]));
	sum=Cadd(sum,Cmul(fac,Conjg(poff))); // avg over tags
      }
      rcal[k+nrx*j]=RCdiv(sum,sqrt(Cmod(sum))); 
      //printf("%i %i %f \n",j,k,(180.0/pi)*atan2f(rcal[k+nrx*j].i,rcal[k+nrx*j].r));
    }
  }
  
  // calibrate tags, making use of receiver cals; cal for each freq, tag
    
  for(j=0;j<nfreq;j++){
      wn = 2.0*pi*freqs[j]/c;
      for(i=0;i<ntag;i++){
	sum=Complex(0.0,0.0);
	for(k=0;k<nrx;k++){ 
	  // anticipated accumulated phase
	  pse = -wn*dist(rxx[k],rxy[k],rxz[k],txx[i],txy[i],txz[i]);
	  fac = Complex(cos(pse),sin(pse));
	  // measured accumulated phase w/ rx correction
	  poff = Cmul(sparams[i+ntag*(k+nrx*j)],rcal[k+nrx*j]);		     
	  sum=Cadd(sum,Cmul(fac,Conjg(poff))); // avg over receivers
	}
	tcal[i+ntag*j]=RCdiv(sum,sqrt(Cmod(sum))); 
	//printf("%i %i %f \n",j,i,(180.0/pi)*atan2f(tcal[i+ntag*j].i,rcal[i+ntag*j].r));
      }
  }
  
  // read actual data, correct rx and tag phases
  
  ndat = diagnose("sParamObj.mat"); // should be ntag*nrx*nfreq
  sparams2 = (fcomplex*)malloc(sizeof(fcomplex)*ndat);
  for(i=0;i<ndat;i++){
    sparams2[i]=Complex(pr[i],pim[i]); 
    //printf("%i %i %f %f \n",i,ndat,pr[i],pim[i]);
  }
  for(i=0;i<ntag;i++)
    for(j=0;j<nfreq;j++)
      for(k=0;k<nrx;k++){
	
	
	sparams2[i+ntag*(k+nrx*j)]=Csub(sparams2[i+ntag*(k+nrx*j)],
					sparams[i+ntag*(k+nrx*j)]);
	

	fac = Cmul(rcal[k+nrx*j],tcal[i+ntag*j]);
	sparams2[i+ntag*(k+nrx*j)]=Cmul(sparams2[i+ntag*(k+nrx*j)],
					fac);	
      }

  printf("Done with cal ...\n");

  // populate g matrix

  // r0 is residual vector, q0 is largest residual vector  
  r0 = (fcomplex*) malloc(sizeof(fcomplex)*ndat);
  q0 = (fcomplex*) malloc(sizeof(fcomplex)*ndat);
  p0 = (fcomplex*) malloc(sizeof(fcomplex)*nmod);
  s0 = (fcomplex*) malloc(sizeof(fcomplex)*nmod);

  lib = (int*) malloc(sizeof(int)*nmod);

  xk = (fcomplex*) malloc(sizeof(fcomplex)*nmod);
  xp = (fcomplex*) malloc(sizeof(fcomplex)*nmod);

  // g is a 1D vector made of stacked vectors in Col major format
  ata = (fcomplex*)malloc(sizeof(fcomplex)*NMAX*NMAX); // huge
  atai = (fcomplex*)malloc(sizeof(fcomplex)*NMAX*NMAX); // huge
  g = (fcomplex*)malloc(sizeof(fcomplex)*ndat*nmod); // huge
  //gg = (fcomplex*)malloc(sizeof(fcomplex)*ndat*nmod); // huge
  gp = (fcomplex*)malloc(sizeof(fcomplex)*ndat*NMAX); // huge
  pe = (fcomplex*)malloc(sizeof(fcomplex)*ndat*NMAX); // huge
  if(g==NULL || gp==NULL || ata==NULL || atai==NULL || pe==NULL){
    printf("out of memory \n");
    exit(1);
  }
  
  for(l=0;l<NX;l++)
    for(m=0;m<NY;m++)
      for(n=0;n<NZ;n++){	
	for(j=0;j<nfreq;j++){
	  wn = 2.0*pi*freqs[j]/c; 
	  for(i=0;i<ntag;i++){
	    for(k=0;k<nrx;k++){
	      pse=-wn*(d1[l+NX*(m+NY*(n+NZ*i))]+d2[l+NX*(m+NY*(n+NZ*k))]);
	      sincosf(pse,&sinx,&cosx);
	      pointer=(long)(i+ntag*(k+nrx*j))+(long)ndat*(long)(l+NX*(m+NY*n));
	      //printf("%i %i %li \n",ndat,nmod,pointer);
	      pointer2=(long)(l+NX*(m+NY*n))+(long)nmod*(long)(i+ntag*(k+nrx*j));
	      g[pointer]=Complex(cosx,sinx);
	      //gg[pointer2]=g[pointer];
	    }
	  }	  
	}
	//printf("%i %i %i \n",l,m,n);
      }
  printf("Done initializing G\n");
      
  /*
  // convert to wavelet space 
  pwtset(4);
  for(j=0;j<nfreq;j++)
    for(i=0;i<ntag;i++)
      for(k=0;k<nrx;k++){
	pointer2=(long)nmod*(long)(i+ntag*(k+nrx*j));
	wtnc(&gg[pointer2]-1,nn-1,3,1,&pwt);
      }
  
  for(l=0;l<NX;l++)
    for(m=0;m<NY;m++)
      for(n=0;n<NZ;n++)
	for(j=0;j<nfreq;j++) 
	  for(i=0;i<ntag;i++)
	    for(k=0;k<nrx;k++){
	      pointer=(long)(i+ntag*(k+nrx*j))+(long)ndat*(long)(l+NX*(m+NY*n));
	      pointer2=(long)(l+NX*(m+NY*n))+(long)nmod*(long)(i+ntag*(k+nrx*j));
	      g[pointer]=gg[pointer2];
	    }
  printf("Done transforming \n");
  */

  // prepare iteration ... note that image is complex here

  mone = Complex(-1.0,0.0);
  one = Complex(1.0,0.0);
  zero = Complex(0.0,0.0);
  inc = 1;
      
  alpha = 1.0e4;

  // initialize residual ... could be anything
  cblas_ccopy(ndat,sparams2,inc,r0,inc);

  for(j=0;j<NMAX;j++){ // arbitrary, but expense grows quickly with no.

    // find column of library most correlated with residual
    // p0 = Gt r
    cblas_cgemv(CblasColMajor,CblasConjTrans,
		ndat,nmod,&one,g,ndat,r0,inc,&zero,p0,inc);

    cmax=0.0;    
    imax=0;
    
    for(i=0;i<nmod;i++){ // test each term in library for correlation
      im = 0;
      for(j=0;j<p;j++){
	if(i==lib[j])
	  im = 1;
      }
      if(!im){ // make sure term hasn't already been incorporated
	corr=Cmod(p0[i]); // already normalized
	if(corr>cmax){
	  cmax=corr;
	  imax=i;
	}
      }
    }
    lib[p]=imax;

    // copy component of design matrix g into subarray g'

    for(i=0;i<ndat;i++){
      pointer=(long)i+(long)ndat*(long)imax;
      pointer2=(long)i+(long)ndat*(long)p;
      gp[pointer2]=g[pointer]; // can speed this up
      //printf("%i %f %f \n",i,gp[pointer2].r,gp[pointer2].i);
    }
    p++; // note increasing subarray size
    printf("%i %i %f \n",p,imax,cmax);

    // ata = g't g'
    cblas_cgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,p,p,ndat,&one,
    		gp,ndat,gp,ndat,&zero,ata,p);        

    // regularize
    for(i=0;i<p;i++){      
      ata[i*(1+p)]=Cadd(ata[i*(1+p)],Complex(alpha,0.0));
    }

    // take inverse 
    MatrixComplexInverse(atai,ata,p); 

    // form pseudoinverse (could also used underdetermined form, regularize)
    // pe = ata ^-1 at
    cblas_cgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,p,ndat,p,&one,
    		atai,p,gp,ndat,&zero,pe,p);
    
    // update model estimate
    cblas_ccopy(ndat,sparams2,inc,r0,inc);
    cblas_cgemv(CblasColMajor,CblasNoTrans,
		p,ndat,&one,pe,p,r0,inc,&zero,xp,inc);

    // put columns in right place
    for(i=0;i<p;i++)
      xk[lib[i]]=xp[i];
    
    // update residual
    // r = d - gm'
    cblas_cgemv(CblasColMajor,CblasNoTrans,
		ndat,nmod,&mone,g,ndat,xk,inc,&one,r0,inc);   

  }

  // transform image back from wavelet space     
  //wtnc(xk-1,nn-1,3,-1,&pwt);


  // now form image

  printf("Imaging ...\n");

  free(g);

  for(l=0;l<NX;l++)
    for(m=0;m<NY;m++)
      for(n=0;n<NZ;n++)
	image[n+NZ*(m+NY*l)]=Cmod(xk[l+NX*(m+NY*n)]);  
  
  // hdf5 stuff

  /* Create a new file using default properties. */
  file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  xx = (float*)malloc(sizeof(float)*(NX+1)*(NY+1)*(NZ+1));
  yy = (float*)malloc(sizeof(float)*(NX+1)*(NY+1)*(NZ+1));
  zz = (float*)malloc(sizeof(float)*(NX+1)*(NY+1)*(NZ+1));

  for(l=0;l<NX+1;l++){
    x=x0+(x1-x0)*(float)l/(float)(NX);
    for(m=0;m<NY+1;m++){
      y=y0+(y1-y0)*(float)m/(float)(NY);
      for(n=0;n<NZ+1;n++){
	z=z0+(z1-z0)*(float)n/(float)(NZ);

	xx[n+(NZ+1)*(m+(NY+1)*l)]=x;
	yy[n+(NZ+1)*(m+(NY+1)*l)]=y;
	zz[n+(NZ+1)*(m+(NY+1)*l)]=z;

      }
    }
  }

  /* Create the data space for the dataset. */
  dims[0] = NX; 
  dims[1] = NY; 
  dims[2] = NZ;
  dataspace_id1 = H5Screate_simple(3, dims, NULL);
  dims[0] = NX+1; 
  dims[1] = NY+1; 
  dims[2] = NZ+1;
  dataspace_id2 = H5Screate_simple(3, dims, NULL);
	      
  /* Create the dataset. */
  data_id0 = H5Dcreate2(file_id, "/x", H5T_NATIVE_FLOAT, dataspace_id2, 
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  data_id1 = H5Dcreate2(file_id, "/y", H5T_NATIVE_FLOAT, dataspace_id2, 
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  data_id2 = H5Dcreate2(file_id, "/z", H5T_NATIVE_FLOAT, dataspace_id2, 
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  data_id3 = H5Dcreate2(file_id, "/dene_0001", H5T_NATIVE_FLOAT, dataspace_id1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  /* write dataset */
  status = H5Dwrite(data_id0, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,xx);

  status = H5Dwrite(data_id1, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,yy);

  status = H5Dwrite(data_id2, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,zz);

  status = H5Dwrite(data_id3, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,image);


  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(data_id0);
  status = H5Dclose(data_id1);
  status = H5Dclose(data_id2);
  status = H5Dclose(data_id3);

  /* Terminate access to the data space. */ 
  status = H5Sclose(dataspace_id1);
  status = H5Sclose(dataspace_id2);		      

  /* Close the file. */
  status = H5Fclose(file_id);

}

//........................................................................................
void cgeTranspose(fcomplex *Transposed, fcomplex *M ,int n)
{

int i,j;
for(i=0;i<n;i++)

for(j=0;j<n;j++) Transposed[i+n*j] = M[i*n+j];
}

//.........................................................................................
void MatrixComplexInverse(fcomplex *invA, fcomplex *A, int n)
{

int LWORK=10*n;

int *permutations;

fcomplex *WORK, *tempA;

tempA = (fcomplex*) malloc( n*n*sizeof(fcomplex) );

permutations = (int*) malloc( 2*n*sizeof(int) );

WORK = (fcomplex *)malloc(LWORK*sizeof(fcomplex));

int INFO;

cgeTranspose(tempA,A,n);

cgetrf_( &n, &n, tempA , &n, permutations , &INFO );

if (INFO != 0) {
 printf("ComplexMatrixInverse: Error at zgetrf  \n"); exit(0);
 }

cgetri_( &n, tempA , &n, permutations , WORK, &LWORK, &INFO );

if (INFO != 0) {
 printf("ComplexMatrixInverse: Error at zgetri  \n"); exit(0);
 }

cgeTranspose(invA,tempA,n);

free(WORK);
free(tempA);
free(permutations);

}
