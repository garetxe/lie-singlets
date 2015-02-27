#include  "lie.h"

poly* Adams(index n,poly* p)
{ 
if (n==1) return p; /* avoid work in this trivial case */
  { index i,j, r=Lierank(grp); poly* dom_ch=Domchar_p(p);
    for (i=0; i<dom_ch->nrows; ++i)
      for (j=0; j<r; j++)
	dom_ch->elm[i][j] *= n;
    
    { poly* result=Vdecomp(dom_ch);
      freepol(dom_ch);
      return result; }
  }
}

poly* SAtensor(boolean alt,index m,poly* p)
{ index n,r=Lierank(grp); poly** adams,** q,* result;
  if (m==0) return poly_one(r);  else if (m==1) return p;

  adams=alloc_array(poly*,m+1); 
  for (n=1; n<=m; ++n) adams[n]=Adams(n,p);
  q=alloc_array(poly*,m+1);
  q[0]=poly_one(r);
  for (n=1; n<=m; ++n)
  { 
    { index i; q[n]=Tensor(p,q[n-1]); /* the initial term of the summation */
      for (i=2; i<=n; ++i) q[n] =
        Add_pol_pol(q[n],Tensor(adams[i],q[n-i]),alt&&i%2==0);
    }
    
    { index i; bigint* big_n=entry2bigint(n);  setshared(big_n);
      for (i=0; i<q[n]->nrows; ++i)
      { bigint** cc= &q[n]->coef[i]
             ,* c= (clrshared(*cc),isshared(*cc)) ? copybigint(*cc,NULL) : *cc;
        *cc=divq(c,big_n); setshared(*cc);
        
        { if (c->size != 0)
            error("Internal error (SAtensor): remainder from %ld.\n" ,(long)n);
          freemem(c);
        }
      }
      clrshared(big_n); freemem(big_n);
    }
  }
  result=q[m];
{ for (n=1; n<=m; ++n) freepol(adams[n]); } freearr(adams);
{ for (n=0; n<m; ++n)  freepol(q[n]); } freearr(q);
 return result;
}

#define CACHE_DEBUG 1

/* This assumes that the rank of the group and the entries of the
   partition will not be larger than 99... */
static char *fname_from_partition(entry *partition, int n)
{
  char *fname = malloc(3*n+10+4);
  char *p = fname;
  int k;

  if (fname == NULL)
    return NULL;

  p += sprintf(p, "[%02d]", Lierank(grp));

  for (k=0;k<n && partition[k]>0; k++) {
    if (k>0) {
      *p = ',';
      p++;
    }
    p += sprintf(p, "%d", partition[k]);
  }

  sprintf(p, ".cache");

  return fname;
}

static poly *load_poly_from_cache(entry *partition, int n)
{
  char *fname = fname_from_partition(partition, n);
  FILE *fp = fopen(fname, "r");
  index nrows, ncols;
  poly *p;
  int i, j;

  if (fp == NULL) {
    free(fname);
    return NULL;
  }

  if (fread(&nrows, sizeof(nrows), 1, fp) != 1) {
    perror("nrows");
    return NULL;
  }

  if (fread(&ncols, sizeof(ncols), 1, fp) != 1) {
    perror("ncols");
    return NULL;
  }

  p = mkpoly(nrows, ncols);

  if (fread(&p->rowsize, sizeof(p->rowsize), 1, fp) != 1) {
    perror("rowsize");
    return NULL;
  }  

  for (i=0; i<nrows; i++) {
    for (j=0; j<ncols; j++) {
      if (fread(&p->elm[i][j], sizeof(entry), 1, fp) != 1) {
	perror("entry");
	return NULL;
      }
    }
  }

  for (i=0; i<nrows; i++) {
    int len;
    string s;
    bigint *a, *b;

    if (fread(&len, sizeof(len), 1, fp) != 1) {
      perror("len");
      return NULL;
    }

    s = calloc(len+1, 1);
    if (s == NULL) {
      return NULL;
    }

    if (fread(s, len, 1, fp) != 1) {
      perror("coeff");
      return NULL;
    }

    /* str2bigint does not understand negative numbers, make its life
       easier. */
    if (s[0] == '-') {
      a = entry2bigint(-1);
      b = str2bigint(s+1);
    } else {
      a = entry2bigint(1);
      b = str2bigint(s);
    }

    p->coef[i] = mult(a, b);
    setshared(p->coef[i]);

    free(s);
  }

#if 0 && defined(CACHE_DEBUG)
  fprintf(stderr, "[CD]\tLoaded <%s>\n", fname);
#endif

  fclose(fp);
  free(fname);

  return p;
}

static void save_poly_to_cache(entry *partition, int n, poly *p)
{
  char *fname = fname_from_partition(partition, n);
  FILE *fp = fopen(fname, "w");
  index i, j;

#ifdef CACHE_DEBUG
  fprintf(stderr, "[CD]\t\t\tSaving <%s>\n", fname);
#endif

  fwrite(&p->nrows, sizeof(p->nrows), 1, fp);
  fwrite(&p->ncols, sizeof(p->ncols), 1, fp);
  fwrite(&p->rowsize, sizeof(p->rowsize), 1, fp);

  for (i=0; i<p->nrows; i++) {
    for (j=0; j<p->ncols; j++) {
      fwrite(&p->elm[i][j], sizeof(entry), 1, fp);
    }
  }

  for (i=0; i<p->nrows; i++) {
    string s = bigint2str(p->coef[i]);
    int len = strlen(s);

    fwrite(&len, sizeof(int), 1, fp);
    fwrite(s, len, 1, fp);

    freem(s);
  }

  fflush(fp);
  fsync(fileno(fp));
  fclose(fp);

  free(fname);
}

/* Given a bigint N, construct the polynomial NX[0,0,...] */
static poly *N_singlets(bigint *N)
{
  poly *one = poly_one(Lierank(grp));
  poly *result;

  setshared(one);
  result = Mul_bin_pol(N, one);
  clrshared(one);
  
  return result;
}

/* Given a polynomial, find the number of singlets it contains. Notice
   that the returned bigint is not newly allocated, but it is rather
   the coeff in the polynomial, so there is no need to free it. */
static bigint *singlets_in_poly(poly *p)
{
  index i, j;

  for (i=0; i < p->nrows; i++) {
    for (j=0; j < p->ncols; j++) {
      if (p->elm[i][j] != 0)
	break;
    }

    if (j == p->ncols) {
      /* We found a singlet */
      return p->coef[i];
    }
  }

  /* If we reach this point we found no singlet, so return the bigint
     called null, representing 0. */
  return null;
}

/* Compute the product of the given reps (in Adams), indexed by kappa,
   a vector of length k. This tries to load from the cache as much of
   the computation as possible. 
*/
static poly *Adams_products_full(entry *kappa, index k, poly **adams,
				 index i, index maxi, index lmu)
{
  poly *prod, *a, *b;

  /* If we have a single term life is easy */
  if (k == 1) {
    return adams[kappa[0]];
  }

  /* Try loading everything from the cache */
  prod = load_poly_from_cache(kappa, k);

  if (prod) {
    return prod;
  }

#ifdef CACHE_DEBUG
  /* Print some useful info to know where we are in the computation */
  {
    int j;

    fprintf(stderr, "[CD]\t");

    if (k == lmu) {
      fprintf(stderr, "[%d / %d]\t", i+1, maxi);
    } else {
      fprintf(stderr, "\t\t");
    }

    for (j=0; j<k; j++) {
      fprintf(stderr, "%d ", kappa[j]);
    }
    fprintf(stderr, "\n");
  }
#endif

  /* We haven't computed this before, compute it now. */
  
  /* The most efficient strategy depends on whether we only are
     computing the singlets or we are computing everything. If we are
     computing everything we give priority to computing */
  b = adams[kappa[k-1]];
  if (k > 2) {
    a = Adams_products_full(kappa, k-1, adams, i, maxi, lmu);
  } else {
    /* When we have two terms we can compute directly */
    a = adams[kappa[0]];
  }
  prod = Tensor(a,b);

  freepol(a);
  freepol(b);

  /* Save the result for future computations */
  save_poly_to_cache(kappa, k, prod);

  return prod;
}

/* True zero iff the given vectors are each others transpose. */
static int are_transpose(entry *a, index na,
			 entry *b, index nb)
{
  int i;

  if (na != nb) {
    return 0;
  }

  for (i = 0; i < na; i++) {
    if (a[i] != b[nb-i-1]) {
      return 0;
    }
  }

  return 1;
}

/* Compute the given product of Adams tensors. */
static poly *Adams_products(entry *kappa, index k, poly **adams,
			    index _i, index maxi, index lmu,
			    int only_singlets)
{
  poly *first;
  poly *rest;
  bigint *nsinglets;
  poly *result;
  index i, j;

  /* If we are interested in the full result we just let
     Adams_products_full do the real work. */
  if (!only_singlets) {
    return Adams_products_full(kappa, k, adams,
			       _i, maxi, lmu);
  }

  /* We just want the singlets. */
  if (k == 1) {
    /* We have only been requested the singlet part, return that. */
    return N_singlets(singlets_in_poly(adams[kappa[0]]));
  }

  /* We divide the computation into the first term (which is generally
     the largest) time the rest, so we avoid the most expensive tensor
     product. */
  first = adams[kappa[0]];
  rest = Adams_products_full(kappa+1, k-1, adams, _i, maxi, lmu-1);

  nsinglets = null;
  for (i=0; i<first->nrows; i++) {
    for (j=0; j<rest->nrows; j++) {
      if (are_transpose(first->elm[i], first->ncols,
			rest->elm[j], rest->ncols)) {
	/* Found a singlet, add the product of the respective
	   coefficients to the number of singlets. */
	bigint *prod, *tmp;

	setshared(first->coef[i]);
	setshared(rest->coef[j]);
	prod = mult(first->coef[i], rest->coef[j]);
	clrshared(rest->coef[j]);
	clrshared(first->coef[i]);
	nsinglets = add(nsinglets, prod);
      }
    }
  }

  result = N_singlets(nsinglets);

  /* zzz Not sure what's wrong with the memory management... */
  /*
    freep(rest);
  */

  return result;
}


poly* Plethysm(entry* lambda,index l,index n,poly* p)
{
  int only_singlets = getenv("LIE_ONLY_SINGLETS") != NULL;

  if (n==0) return poly_one(Lierank(grp));  else if (n==1) return p;

  { index i,j;
    poly* sum= poly_null(Lierank(grp)),**adams=alloc_array(poly*,n+1);
    poly* chi_lambda=MN_char(lambda,l);
    

    for (i=1; i<=n; ++i) { adams[i]=Adams(i,p); setshared(adams[i]); }
    
#ifdef CACHE_DEBUG
    fprintf(stderr, "[CD] Computing Adams [%d, %d, %d]\n", chi_lambda->nrows, l, n);
#endif

    for (i=0;i<chi_lambda->nrows;i++) {
	entry* mu=chi_lambda->elm[i]; poly* prod;
	int lmu; /* length of mu */

	for (lmu = 0; lmu < n && mu[lmu] > 0; lmu++);
  
	prod = Adams_products(mu, lmu, adams, i, chi_lambda->nrows, lmu,
			      only_singlets);

	sum= Addmul_pol_pol_bin(sum,prod,mult(chi_lambda->coef[i],Classord(mu,n)));
    }

    freemem(chi_lambda);

    setshared(p); /* protect |p|; it coincides with |adams[1]| */
    for (i=1; i<=n; ++i)
      {
	clrshared(adams[i]);
	freepol(adams[i]);
      }

    freearr(adams);
  clrshared(p);


    { bigint* fac_n=fac(n);  setshared(fac_n); /* used repeatedly */
      for (i=0; i<sum->nrows; ++i)
      { bigint** cc= &sum->coef[i]
             ,* c= (clrshared(*cc),isshared(*cc)) ? copybigint(*cc,NULL) : *cc;
        *cc=divq(c,fac_n); setshared(*cc);
        
	/*
	  zzz There is something wrong with the memory management in
	  the singlets only case...
	*/
	if (!only_singlets) {
	  if (c->size!=0) error("Internal error (plethysm).\n");
	  else freemem(c);
	}
      }
      clrshared(fac_n); freemem(fac_n);
    }

    return sum;
  }
}

