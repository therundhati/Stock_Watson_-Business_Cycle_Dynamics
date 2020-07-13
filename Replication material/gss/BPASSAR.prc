

/* -- BandPass Procedures -- */
/* bpassar.prc -- Bandpass with AR padding for ends */
proc(2) = bpassar(x,updc,lpdc,n,nar,tcode,seflag);
/* -- input:
      x = series to be filtered
      updc = Period corresponding to upper cutoff frequency
      lpdc = Period corresponding to lower cutoff frequency
      n = number of terms in moving average filter
      nar = order of AR used for padding
      tcode = transformation code for AR model
              0 -- no transformation
              1 -- first difference
      seflag == flag for standard errors (associated with endpoint problem)
              0 -- no SEs computed
              1 -- SEs computed

      Return is filtered value of series and standard errors
*/
local xtran, xpad, xpadf, xbp, sebp;
@ -- Pad Series -- @
if tcode .== 0; xtran=x; endif;
if tcode .== 1; xtran=x[2:rows(x)]-x[1:rows(x)-1]; endif;
xpad=pad(x,xtran,tcode,n,nar);

@ -- Bandpass Padded Series -- @
xpadf = bpass(xpad,updc,lpdc,n);
xbp=xpadf[n+1:rows(xpadf)-n];
if seflag .== 0; sebp=miss(zeros(rows(xbp),1),0); endif;
if seflag .== 1;
 sebp=sebpass(x,xtran,tcode,n,nar,updc,lpdc);
endif;

retp(xbp,sebp);
endp;

/* bpass1sd.prc -- Bandpass, 1 sided ar padding -- */
proc(1) = bpass1sd(x,updc,lpdc,n,nar,tcode,narmin);
/* -- input:
      x = series to be filtered
      updc = Period corresponding to upper cutoff frequency
      lpdc = Period corresponding to lower cutoff frequency
      n = number of terms in moving average filter
      nar = order of AR used for padding
      tcode = transformation code for AR model
              0 -- no transformation
              1 -- first difference
      narmin = minimum number of terms for AR estimates

      
*/
local xtran, xpad, xpadf, xbp, sebp, xtmp, i;

xbp=miss(zeros(rows(x),1),0);

if ismiss(x) .== 1;
 "X contains missing values in bpass1sd -- processing stops";
 stop;
endif;


i=2; do while i <= rows(x);
 xtmp=x[1:i];
 @ -- Pad Series -- @
 if tcode .== 0; xtran=xtmp; endif;
 if tcode .== 1; xtran=xtmp[2:rows(xtmp)]-xtmp[1:rows(xtmp)-1]; endif;  
 if rows(xtran) .>= narmin;
   xpad=pad(xtmp,xtran,tcode,n,nar);
   @ -- Bandpass Padded Series -- @
   xpadf = bpass(xpad,updc,lpdc,n);
   xbp[i]=xpadf[rows(xpadf)-n];
 endif;
i=i+1; endo;
retp(xbp);
endp;

/*  Sebpass.gss
    Compute SEs for Bandpass data (using AR padding) -- 
    SEs associated with "endpoint" problems associated with
    estimating pre and postsample variables
*/

@ -- Compute BP filter weights -- @

proc(1) = sebpass(x,xtran,tcode,n,nar,updc,lpdc);

local avec, w, xma, c, i, j, wc, se, v;
if (rows(x) < n);
 "Rows of X are less than n is SEpass";
 "Cannot compute SEs";
 se=miss(zeros(rows(x),1),0);
 retp(se);
endif;

avec=bpweight(updc,lpdc,n);       @ Bandpass weights @  
@ -- Compute 1 sided weights -- @
w=avec[n+2:rows(avec)];    
xma=ar2ma(x,xtran,tcode,n,nar);   @ MA Rep computed from AR -- truncated at n 
                                    Normalized with SD of error = 1 @
c=zeros(n,n);
i=1; do while i <= n;
 c[i:n,i]=xma[1:n+1-i];
i=i+1; endo;

se=zeros(rows(x),1);
i=1; do while i <= rows(x);
 @ Find Distance from end @
 j=rows(x)-i+1;
 j=minc(j|i);
 if j <= n;
  wc=w'c[.,j:n];
  v=wc*wc';
  se[i]=sqrt(v);
 endif;
i=i+1; endo;

retp(se);
endp;

 
/* AR2MA -- Compute MA Rep from Estimated AR
            -- truncate after n terms 
            -- Allow Levels or First Differences 
             MA Rep computed from AR -- truncated at n 
             Normalized with SD of error = 1
*/
            
proc(1)=ar2ma(y,yf,tcode,n,nar);

/* -- Pad Data series y out using AR Forecasts and Backcasts
      y -- series to be padded
     yf -- data to use in AR (levels or diffs)
     fcode -- 1 if yf is first dif of y
              otherwise yf=y
     n -- number of terms for MA rep
     nar -- order of autoregression to use
*/
local w, x, i, beta, e, ve, se, comp, z, xma;

@ -- Compute MA Rep for Transformed Series -- @
xma=miss(zeros(n,1),1);
@ Pad out future @
w=yf[nar+1:rows(yf)];
x=ones(rows(w),1);
i=1; do while i<=nar;
 x=x~yf[nar+1-i:rows(yf)-i];
i=i+1; endo;
beta=invpd(x'x)*(x'w);
e=w-x*beta;
ve=e'e/(rows(x)-cols(x));
se=sqrt(ve);              @ Standard Error of Regression @
beta=beta[2:rows(beta)];
comp=zeros(nar,nar);      @ Companion Matrix @
comp[1,1:nar]=beta';
if nar .> 1;
 comp[2:nar,1:nar-1]=eye(nar-1);
endif;
z=zeros(nar,1);
z[1]=1;
xma[1]=1;
i=2; do while i <= n;
 z=comp*z;
 xma[i]=z[1];
i=i+1; endo;
xma=xma*se;
if tcode .== 1;
 xma=cumsumc(xma);
endif;
retp(xma);
endp;



/*
    bpweight.prc
    GAUSS Translation of King's bpass.m
    1/26/94, mww
    computes bandpass filter weights using upper and lower
    cutoff periods.

*/
proc(1) = bpweight(updc,lpdc,n);

/* -- input:
      x = series to be filtered
      updc = Period corresponding to upper cutoff frequency
      lpdc = Period corresponding to lower cutoff frequency
      n = number of terms in moving average filter

      Note: updc=2 makes this a high pass filter
      (since period = 2 implies om=pi).

      if lpdc .> rows(x), frequency is set to zero
*/
local avec, omubar, omlbar, step, np, om, omp, kk, lam, akvec, j, i;

/*
"Band Pass Filtering of Time Series: ";
"Components between Periods of ";;  updc~lpdc;; "Time Units";;
*/

@ Implied Frequencies @
omubar=2*pi/updc;
omlbar=2*pi/lpdc;
@ if (lpdc .>= rows(x)); omlbar=0; endif; @

/*
 To construct a low pass filter, with a cutoff frequency of "ombar",
 we note that the transfer function of the approximating filter
 is given by:

  alpha(om) = a0 + a1 cos(om) + ... aK cos(K om)

 and the ak's are given by:

 a0 = ombar/(pi)

 ak = sin(k ombar)/(k pi)

 where ombar is the cutoff frequency.

  We employ the fact that a bandpass filter is the difference between two
  low pass filters,
    bp(L) = bu(L) - bl(L)
  with bu(L) being the filter with the high cutoff point and bl(L) being
  that with the low cutoff point.
*/

@ Define the vector of K's to be studied @

@ Set the grid for om @
step=.01; np=(2/step)+1;
om=seqa(-1,step,np);

om=pi*om;
omp=om/pi;

@ Initialize output matrix @

@  Loop over specified K's @

 @  Construct Filter Weights @

 akvec=zeros(n+1,1);

 akvec[1]=(omubar-omlbar)/(pi);

 kk=1; do while kk <= n;
  akvec[kk+1]=(sin(kk*omubar)-sin(kk*omlbar))/(kk*pi);
 kk=kk+1; endo;
/*
  Impose constraint that transfer is
    (i)  0 at om = 0 if oml .> 0;
    (ii) 1 at om = 0 if oml .== 0;
  This amounts to requiring that weights sum to zero.
  Initial sum of weights:
*/
lam=akvec[1]+2*sumc(akvec[2:n+1]);
/* amount to add to each weight to get sum to add to zero */
if (omlbar .> .00000001);
   lam = -lam/(2*(n+1));
 else;
   lam = (1-lam)/(2*(n+1));
endif;
akvec=akvec+lam;
akvec[1]=akvec[1]+lam;

@ Set vector of weights @

avec=zeros(2*n+1,1);
avec[n+1]=akvec[1];
i=1; do while i <= n;
 avec[n+1-i]=akvec[i+1];
 avec[n+1+i]=akvec[i+1];
i=i+1; endo;

retp(avec);
endp;

/*
    bpass.prc
    GAUSS Translation of King's bpass.m
    1/26/94, mww
    program to compute bandpass filtered series using upper and lower
    cutoff periods.

*/
proc(1) = bpass(x,updc,lpdc,n);

/* -- input:
      x = series to be filtered
      updc = Period corresponding to upper cutoff frequency
      lpdc = Period corresponding to lower cutoff frequency
      n = number of terms in moving average filter

      Note: updc=2 makes this a high pass filter
      (since period = 2 implies om=pi).

      if lpdc .> rows(x), frequency is set to zero
*/
local avec, xf, j, t;
t=rows(x);    @ number of observations @
avec=bpweight(updc,lpdc,n);  @ Bandpass weights @

xf=miss(zeros(t,1),0);
j=n+1; do while j <= t-n;
 xf[j]=x[j-n:j+n]'*avec;
j=j+1; endo;

retp(xf);
endp;

proc(1)=pad(y,yf,fcode,n,nar);

/* -- Pad Data series y out using AR Forecasts and Backcasts
      y -- series to be padded
     yf -- data to use in AR (levels or diffs)
     fcode -- 1 if yf is first dif of y
              otherwise yf=y
     n -- number of terms to pad forward and backward
     nar -- order of autoregression to use
*/
local w, x, i, beta, v, forc, ypad;

@ Pad out future @
w=yf[nar+1:rows(yf)];
x=ones(rows(w),1);
i=1; do while i<=nar;
 x=x~yf[nar+1-i:rows(yf)-i];
i=i+1; endo;
beta=invpd(x'x)*(x'w);
v=rev(yf[rows(yf)-nar+1:rows(yf)]);
forc=zeros(n,1);
i=1; do while i <= n;
 forc[i]=beta'(1|v);
 v[2:rows(v)]=v[1:rows(v)-1];
 v[1]=forc[i];
i=i+1; endo;
if fcode .== 1;
 forc=cumsumc(forc) + ones(n,1)*y[rows(y),1];
endif;
ypad=y|forc;

@ Pad out past, by reversing series @
yf=rev(yf);
if fcode .== 1;
 yf=-yf;
endif;
w=yf[nar+1:rows(yf)];
x=ones(rows(w),1);
i=1; do while i<=nar;
 x=x~yf[nar+1-i:rows(yf)-i];
i=i+1; endo;
beta=invpd(x'x)*(x'w);
v=rev(yf[rows(yf)-nar+1:rows(yf)]);
forc=zeros(n,1);
i=1; do while i <= n;
 forc[i]=beta'(1|v);
 v[2:rows(v)]=v[1:rows(v)-1];
 v[1]=forc[i];
i=i+1; endo;
if fcode .== 1;
 forc=cumsumc(forc) + ones(n,1)*y[1,1];
endif;
forc=rev(forc);
ypad=forc|ypad;

retp(ypad);
endp;

/* ----------------------------------------------- */
Proc(1) = bp_filt(x,xfcst,filt);

/*
Construct BP filtered Series -- future padded
 past -- truncated 
 
Input:

    X -- T x N matrix of series to be filtered using a symmetric
         two sided filter filt
    XFCST -- BPN x N matrix of forecasted values of X
    Filt -- 2*BPN + 1 vector of filter weights
    
Output 

    Xfilt -- Txn -- filtered value of X
                    (First BPN elements are missing in this version)
 
   
*/

@ -- Construct Date T State Vector -- @
local xfilt, bpn, y, i;

xfilt=miss(zeros(rows(x),cols(x)),0);
bpn=(rows(filt)-1+.00001)/2;
bpn=trunc(bpn);

y=x|xfcst;

i=bpn+1; do while i <= rows(x);
 xfilt[i,.]=filt'y[i-bpn:i+bpn,.];
i=i+1; endo;

retp(xfilt);
endp;

/* ------------------------------------------------- */
Proc(1) = bp_vcv(c,g,a,filt1);

/*
Construct Variance Matrix of Filtered Estimates (date by date)
Using Companion Form used to construct foreasts

Model is y(t) = a x(t)
         x(t) = c x(t-1) + g*e(t)
         var(e(t))=I
         
Input:
	C, G, A -- Matrices of State Space form of Model
	Filt1 -- Forward part of filter
       
       
Output:
	V: N x BPN matrix of Variances:
	   The i'th column, say V(i) is the vectorized version of
	   the Variance for the Filtered Series at period
	   T-BPN+i        
*/

@ -- Step 1 -- form powers of c times G -- @
local bpn, theta, gam, i, lam, psi, v1, tmp, v;

bpn=rows(filt1);
theta=g;
gam=vec(theta);
i=2; do while i <= bpn;
 theta=c*theta;
 gam=gam~vec(theta);
i=i+1; endo;
lam=zeros(bpn,bpn);
i=1; do while i <= bpn;
 lam[1:bpn+1-i,i]=filt1[i:bpn,1];
i=i+1; endo;
psi=gam*lam;
v1=zeros(rows(a)*rows(a),bpn);
i=1; do while i <= cols(psi);
  tmp=reshape(psi[.,i],cols(g),rows(c));
  tmp=tmp';
  tmp=a*tmp;
  tmp=tmp*tmp';
  v1[.,i]=vec(tmp);
i=i+1; endo;
v1=v1';
v1=rev(v1);
v=cumsumc(v1);
retp(v);

endp;
