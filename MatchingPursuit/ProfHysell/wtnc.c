#define NRANSI
#include "nrutil.h"
// complex matrix a
#include "complex.h"

void wtnc(fcomplex a[], unsigned long nn[], int ndim, int isign,
	void (*wtstep)(float [], unsigned long, int))
{
	unsigned long i1,i2,i3,k,n,nnew,nprev=1,nt,ntot=1;
	int idim;
	float *wksp,*wkspr,*wkspi;

	for (idim=1;idim<=ndim;idim++) ntot *= nn[idim];
	wkspr=vector(1,ntot);
	wkspi=vector(1,ntot);

	for (idim=1;idim<=ndim;idim++) {
		n=nn[idim];
		nnew=n*nprev;
		if (n > 4) {
		  for (i2=0;i2<ntot;i2+=nnew) {
		    for (i1=1;i1<=nprev;i1++) {
		      for (i3=i1+i2,k=1;k<=n;k++,i3+=nprev){
			wkspr[k]=a[i3].r;
			wkspi[k]=a[i3].i;
		      }
		      if (isign >= 0) {
			for(nt=n;nt>=4;nt >>= 1){
			  (*wtstep)(wkspr,nt,isign);
			  (*wtstep)(wkspi,nt,isign);
			}
		      } else {
			for(nt=4;nt<=n;nt <<= 1){
			  (*wtstep)(wkspr,nt,isign);
			  (*wtstep)(wkspi,nt,isign);
			}
		      }
		      
		      for (i3=i1+i2,k=1;k<=n;k++,i3+=nprev) 
			a[i3]=Complex(wkspr[k],wkspi[k]);
		    }
		  }
		}
		nprev=nnew;
	}
	free_vector(wkspr,1,ntot);
	free_vector(wkspi,1,ntot);

}

#undef NRANSI
