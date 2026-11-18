/* NUMERIA ES8-P242 independent outward-rounded LDL pivot replay.
 * Reads dyadic midpoint matrices directly from P237 JSON after extraction.
 * NO Arb, python-flint, or use of printed P237 center pivots.
 * The public libmpfr ABI is declared here because mpfr.h is absent in this environment.
 * Link against system libmpfr.so.6 + libgmp. Every arithmetic operation is directed.
 * STATUS: finite symmetric matrix only. Not a Weil/RH theorem.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

typedef long mpfr_prec_t;
typedef long mpfr_exp_t;
typedef int mpfr_sign_t;
typedef enum {MPFR_RNDN=0,MPFR_RNDZ=1,MPFR_RNDU=2,MPFR_RNDD=3,MPFR_RNDA=4} mpfr_rnd_t;
typedef struct {mpfr_prec_t _mpfr_prec; mpfr_sign_t _mpfr_sign; mpfr_exp_t _mpfr_exp; mp_limb_t *_mpfr_d;} __mpfr_struct;
typedef __mpfr_struct mpfr_t[1];
extern void mpfr_init2(mpfr_t,long);
extern void mpfr_clear(mpfr_t);
extern int mpfr_set_str(mpfr_t,const char*,int,mpfr_rnd_t);
extern int mpfr_mul_2si(mpfr_t,mpfr_t,long,mpfr_rnd_t);
extern int mpfr_add(mpfr_t,mpfr_t,mpfr_t,mpfr_rnd_t);
extern int mpfr_sub(mpfr_t,mpfr_t,mpfr_t,mpfr_rnd_t);
extern int mpfr_mul(mpfr_t,mpfr_t,mpfr_t,mpfr_rnd_t);
extern int mpfr_div(mpfr_t,mpfr_t,mpfr_t,mpfr_rnd_t);
extern int mpfr_cmp(mpfr_t,mpfr_t);
extern int mpfr_cmp_si(mpfr_t,long);
extern int mpfr_number_p(mpfr_t);
extern double mpfr_get_d(mpfr_t,mpfr_rnd_t);
extern int mpfr_snprintf(char*,size_t,const char*,...);

typedef struct{mpfr_t l,h;} I;
static mpfr_prec_t prec=512;
static void init(I *x){mpfr_init2(x->l,prec);mpfr_init2(x->h,prec);}
static void clear(I *x){mpfr_clear(x->l);mpfr_clear(x->h);}
/* For the dyadic exact center x=m*2^e, directed parse and binary scaling. */
static void parse(I*x,const char*m,long e){
 if(mpfr_set_str(x->l,m,10,MPFR_RNDD)||mpfr_set_str(x->h,m,10,MPFR_RNDU)){
  fprintf(stderr,"bad mantissa\n");exit(2);
 }
 mpfr_mul_2si(x->l,x->l,e,MPFR_RNDD);mpfr_mul_2si(x->h,x->h,e,MPFR_RNDU);
}
extern int mpfr_set(mpfr_t,mpfr_t,mpfr_rnd_t);
extern long mpfr_get_z_2exp(mpz_t,mpfr_t);
extern int mpfr_set_si(mpfr_t,long,mpfr_rnd_t);
extern int mpfr_abs(mpfr_t,mpfr_t,mpfr_rnd_t);

static void quot_pos(I*o,I*a,I*b){
 if(mpfr_cmp_si(b->l,0)<=0){fprintf(stderr,"nonpositive denominator\n");exit(4);}
 if(mpfr_cmp_si(a->l,0)>=0){
   mpfr_div(o->l,a->l,b->h,MPFR_RNDD);
   mpfr_div(o->h,a->h,b->l,MPFR_RNDU);
 }else if(mpfr_cmp_si(a->h,0)<=0){
   mpfr_div(o->l,a->l,b->l,MPFR_RNDD);
   mpfr_div(o->h,a->h,b->h,MPFR_RNDU);
 }else{
   mpfr_div(o->l,a->l,b->l,MPFR_RNDD);
   mpfr_div(o->h,a->h,b->l,MPFR_RNDU);
 }
}
static mpfr_t scratch;
static void mul(I*o,I*a,I*b){
 /* Sign-specialized outward interval multiplication (one lower, one upper).
    Four-corner fallback only if at least one factor straddles zero. */
 int a_pos=mpfr_cmp_si(a->l,0)>=0, a_neg=mpfr_cmp_si(a->h,0)<=0;
 int b_pos=mpfr_cmp_si(b->l,0)>=0, b_neg=mpfr_cmp_si(b->h,0)<=0;
 if(a_pos&&b_pos){
  mpfr_mul(o->l,a->l,b->l,MPFR_RNDD);mpfr_mul(o->h,a->h,b->h,MPFR_RNDU);return;
 }
 if(a_neg&&b_neg){
  mpfr_mul(o->l,a->h,b->h,MPFR_RNDD);mpfr_mul(o->h,a->l,b->l,MPFR_RNDU);return;
 }
 if(a_pos&&b_neg){
  mpfr_mul(o->l,a->h,b->l,MPFR_RNDD);mpfr_mul(o->h,a->l,b->h,MPFR_RNDU);return;
 }
 if(a_neg&&b_pos){
  mpfr_mul(o->l,a->l,b->h,MPFR_RNDD);mpfr_mul(o->h,a->h,b->l,MPFR_RNDU);return;
 }
 /* General case, all four endpoint products. */
 mpfr_mul(o->l,a->l,b->l,MPFR_RNDD);
 mpfr_mul(scratch,a->l,b->h,MPFR_RNDD);if(mpfr_cmp(scratch,o->l)<0)mpfr_set(o->l,scratch,MPFR_RNDN);
 mpfr_mul(scratch,a->h,b->l,MPFR_RNDD);if(mpfr_cmp(scratch,o->l)<0)mpfr_set(o->l,scratch,MPFR_RNDN);
 mpfr_mul(scratch,a->h,b->h,MPFR_RNDD);if(mpfr_cmp(scratch,o->l)<0)mpfr_set(o->l,scratch,MPFR_RNDN);
 mpfr_mul(o->h,a->l,b->l,MPFR_RNDU);
 mpfr_mul(scratch,a->l,b->h,MPFR_RNDU);if(mpfr_cmp(scratch,o->h)>0)mpfr_set(o->h,scratch,MPFR_RNDN);
 mpfr_mul(scratch,a->h,b->l,MPFR_RNDU);if(mpfr_cmp(scratch,o->h)>0)mpfr_set(o->h,scratch,MPFR_RNDN);
 mpfr_mul(scratch,a->h,b->h,MPFR_RNDU);if(mpfr_cmp(scratch,o->h)>0)mpfr_set(o->h,scratch,MPFR_RNDN);
}
static void sub_inplace(I*x,I*y){
 mpfr_sub(x->l,x->l,y->h,MPFR_RNDD);
 mpfr_sub(x->h,x->h,y->l,MPFR_RNDU);
 if(mpfr_cmp(x->l,x->h)>0){fprintf(stderr,"invalid interval order\n");exit(5);}
}
static void emit_mpfr(FILE*f,mpfr_t x){
 mpz_t m;mpz_init(m);long e=mpfr_get_z_2exp(m,x);
 char *s=mpz_get_str(NULL,10,m);
 fprintf(f,"%s %ld",s,e);
 /* mpz_get_str uses malloc() by default from GMP allocator */
 free(s);mpz_clear(m);
}
#define AT(j,i) (S[(size_t)(j)*((j)+1)/2+(i)])
int main(int argc,char**argv){
 if(argc!=5){fprintf(stderr,"usage: %s matrix.txt pivots.txt bits n\n",argv[0]);return 2;}
 prec=atol(argv[3]);int n=atoi(argv[4]);if(prec<128||prec>4096||n<1||n>1000)return 2;
 FILE*f=fopen(argv[1],"r"),*out=fopen(argv[2],"w");
 if(!f||!out){perror("open");return 2;}
 I*S=malloc((size_t)n*(n+1)/2*sizeof(I));
 I*col=malloc((size_t)n*sizeof(I));
 if(!S||!col)return 2;
 for(int i=0;i<n;i++)init(&col[i]);
 for(int j=0;j<n;j++)for(int i=0;i<=j;i++)init(&AT(j,i));
 char*m=malloc(3000);int j,i;long e;
 for(int jj=0;jj<n;jj++)for(int ii=0;ii<=jj;ii++){
   if(fscanf(f,"%d %d %2999s %ld",&i,&j,m,&e)!=4 || i!=ii||j!=jj){
     fprintf(stderr,"bad input at %d %d (got %d %d)\n",ii,jj,i,j);return 2;
   }
   parse(&AT(jj,ii),m,e);
 }
 char trailing;
 if(fscanf(f," %c",&trailing)==1){fprintf(stderr,"unexpected trailing matrix data\n");return 2;}
 free(m);fclose(f);
 I product;init(&product);mpfr_init2(scratch,prec);
 fprintf(out,"# index low_mantissa low_exp high_mantissa high_exp\n");
 for(int k=0;k<n;k++){
   I*p=&AT(k,k);
   if(!mpfr_number_p(p->l)||!mpfr_number_p(p->h)||mpfr_cmp_si(p->l,0)<=0){
     fprintf(stderr,"INCONCLUSIVE pivot %d approx lower %.8g upper %.8g\n",k,
          mpfr_get_d(p->l,MPFR_RNDD),mpfr_get_d(p->h,MPFR_RNDU));
     return 3;
   }
   fprintf(out,"%d ",k);emit_mpfr(out,p->l);fprintf(out," ");emit_mpfr(out,p->h);fprintf(out,"\n");
   for(int i=k+1;i<n;i++)quot_pos(&col[i],&AT(i,k),p);
   for(int j=k+1;j<n;j++){
     for(int i=k+1;i<=j;i++){
       mul(&product,&col[i],&AT(j,k));
       sub_inplace(&AT(j,i),&product);
     }
   }
   /* After the Schur update, save column k of unit-lower factor L. */
   for(int i=k+1;i<n;i++){
     mpfr_set(AT(i,k).l,col[i].l,MPFR_RNDN);
     mpfr_set(AT(i,k).h,col[i].h,MPFR_RNDN);
   }
   if((k+1)%25==0 || k+1==n){
      fprintf(stderr,"MPFR directed pivots certified: %d/%d (min at step %.6e)\n",k+1,n,mpfr_get_d(p->l,MPFR_RNDD));fflush(stderr);
   }
 }
 /* Independent inverse-Frobenius enclosure via interval triangular solves.
    Solve L D L^T x = e_j for each column j. Exact true L,D lie
    within the factor interval matrices retained above. */
 I*y=malloc((size_t)n*sizeof(I)),*z=malloc((size_t)n*sizeof(I)),*x=malloc((size_t)n*sizeof(I));
 if(!y||!z||!x)return 5;
 for(int i=0;i<n;i++){init(&y[i]);init(&z[i]);init(&x[i]);}
 mpfr_t nu,ab1,ab2,term;mpfr_init2(nu,prec);mpfr_init2(ab1,prec);
 mpfr_init2(ab2,prec);mpfr_init2(term,prec);
 mpfr_set_si(nu,0,MPFR_RNDU);
 for(int j=0;j<n;j++){
   for(int i=0;i<n;i++){
     mpfr_set_si(y[i].l,i==j?1:0,MPFR_RNDN);
     mpfr_set_si(y[i].h,i==j?1:0,MPFR_RNDN);
     for(int k=j;k<i;k++){
       mul(&product,&AT(i,k),&y[k]);sub_inplace(&y[i],&product);
     }
   }
   for(int i=0;i<n;i++)quot_pos(&z[i],&y[i],&AT(i,i));
   for(int i=n-1;i>=0;i--){
     mpfr_set(x[i].l,z[i].l,MPFR_RNDN);
     mpfr_set(x[i].h,z[i].h,MPFR_RNDN);
     for(int k=i+1;k<n;k++){
       mul(&product,&AT(k,i),&x[k]);sub_inplace(&x[i],&product);
     }
   }
   for(int i=0;i<n;i++){
     mpfr_abs(ab1,x[i].l,MPFR_RNDU);mpfr_abs(ab2,x[i].h,MPFR_RNDU);
     if(mpfr_cmp(ab2,ab1)>0)mpfr_set(ab1,ab2,MPFR_RNDN);
     mpfr_mul(term,ab1,ab1,MPFR_RNDU);
     mpfr_add(nu,nu,term,MPFR_RNDU);
   }
   if((j+1)%75==0||j+1==n){
      fprintf(stderr,"MPFR inverse enclosure computed column: %d/%d\n",j+1,n);fflush(stderr);
   }
 }
 if(!mpfr_number_p(nu)||mpfr_cmp_si(nu,0)<=0)return 5;
 fprintf(out,"# inverse_frobenius_squared_upper ");emit_mpfr(out,nu);fprintf(out,"\n");
 fprintf(stderr,"INVERSE_FROB_SQUARED_UPPER approx %.12g\n",mpfr_get_d(nu,MPFR_RNDU));
 mpfr_clear(nu);mpfr_clear(ab1);mpfr_clear(ab2);mpfr_clear(term);
 for(int i=0;i<n;i++){clear(&y[i]);clear(&z[i]);clear(&x[i]);}
 free(y);free(z);free(x);
 fclose(out);
 for(int k=0;k<n;k++){clear(&col[k]);}
 free(col);
 for(int j=0;j<n;j++)for(int i=0;i<=j;i++){clear(&AT(j,i));}
 free(S);
 clear(&product);mpfr_clear(scratch);fprintf(stderr,"CERTIFIED %d/%d\n",n,n);
 return 0;
}
