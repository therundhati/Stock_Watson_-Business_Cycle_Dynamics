/* yulewalk.prc ... Compute Autocovariances from AR coefficients
   
   compute autocovariances for 
   y(t) = phi(1)*y(t-1) + ...+ phi(p)y(t-p) + e(t)
   
   input:
        phi = AR coefficient vector (phi(1)|phi(2)|...|phi(p))
        vare = variance of e
        nacv = number of autocovariances desired
   output:
        acv = nacv+1 x 1 vector with variance and NACV autocovariances
*/

proc(1) = yulewalk(phi,vare,nacv);

local nar, tmp, i, tmp2, f, fi, acv;

nar=rows(phi);
tmp=zeros(nar+1,2*nar+1);
for i (0,nar,1);
 tmp[i+1,nar+2-i:2*nar+1-i]=phi';
endfor;
tmp2=zeros(nar+1,nar+1);
tmp2[.,1]=tmp[.,nar+1];
for i (1,nar,1);
 tmp2[.,1+i] = tmp[.,nar+1-i]+tmp[.,nar+1+i];
endfor;
f=eye(nar+1)-tmp2;
fi=inv(f);
acv=vare*fi[.,1];
if nacv .> nar;
 i = nar; do while i <= nacv;
  tmp=(rev(acv[rows(acv)-nar+1:rows(acv)]))'phi;
  acv=acv|tmp;
 i=i+1; endo;
endif;
acv=acv[1:1+nacv];

retp(acv);
endp;
