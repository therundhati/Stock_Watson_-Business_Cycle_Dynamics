proc(1) = poist(y,tau,imean,itst);

/* -- poist.prc, mww, 5/29/96 (revised)
      Calculate POI test statistic for stochastic trend model:

      H0: SD(a) = 0
          vs.
      HA: SD(a) = tau

      y(t) = b(t) + e(t)
      b(t-1) = b(t-1) + a(t)

      the value of b(0) is an unknown nuissance parameter
      the SD of e(t) is 1

imean = 0   =>  Construct Exact GLS Estimate
      = 1       Construct Asymptotic Approximation

itst = 0   =>  Construct Exact POI Statistic (given estimate of mean)
     = 1       Construct Asymptotic Approximation

*/

local e0, ss0, vague, b1, p1, r, q, t, balt, n, h, nt, ht, e1, ss1;
local b2,p2, b1t, b2t, p1t, p2t, b3, p3, g;
local theta, s, hs, ee, yf, poitst;

theta=1-tau;   @  Approximate MA(1) Coefficient for FD of y @
g=tau*rows(y); @  Local to unity parameter @

@ -- Residuals under the null -- @
e0=y-meanc(y);
ss0=e0'e0;


@ -- Estimate of b(0) under alternative -- @

if imean .== 0;
 vague=1.0E+8;   @ -- Vague Prior -- @
 b1=0;
 p1=vague;
 r=1;
 q=tau.^2;

 b1t=zeros(rows(y),1);
 b2t=zeros(rows(y),1);
 p1t=zeros(rows(y),1);
 p2t=zeros(rows(y),1);

 t=1; do while t <= rows(y);
  b1t[t]=b1; p1t[t]=p1;
  {b1,p1,b2,p2,n,h}=rwkfilt(y[t],b1,p1,r,q);
  b2t[t]=b2; p2t[t]=p2;
 t=t+1; endo;

 b3=b1; p3=p1;
 t=rows(y); do while t >= 1;
  b2=b2t[t];
  p2=p2t[t];
  b1=b1t[t];
  p1=p1t[t];
  {b3,p3}=rwksmth(b1,b2,b3,p1,p2,p3);
 t=t-1; endo;
 balt=b3;
endif;
if imean .== 1;
 s=seqa(1,1,rows(y))/rows(y);
 hs=exp(-g*s) + exp(-g*(2*ones(rows(y),1)-s));
 hs=hs*g/(rows(y)*(1-exp(-2*g)));
 balt=hs'y;
endif;

if itst .== 0;
 r=1;
 q=tau.^2;
 @ -- Find Filtered Residuals under alternative -- @
 nt=zeros(rows(y),1);
 ht=zeros(rows(y),1);
 b1=balt;
 p1=0;
 t=1; do while t <= rows(y);
  {b1,p1,b2,p2,n,h}=rwkfilt(y[t],b1,p1,r,q);
  nt[t]=n;
  ht[t]=h;
 t=t+1; endo;
 e1=nt./sqrt(ht);
endif;

if itst .== 1;
 y=y-balt;
 ee=0;
 yf=recserar(y,ee,theta);
 yf=0|yf[1:rows(y)-1];
 yf=g*yf/rows(y);
 e1=y-yf;
 e1=e1*sqrt(1-tau);   @ This improves small sample approximation @
endif;

ss1=e1'e1;

poitst=rows(y)*(1-(ss1/ss0));
retp(poitst);

endp;
