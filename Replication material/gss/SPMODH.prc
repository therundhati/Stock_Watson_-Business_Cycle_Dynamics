proc(1) = spmodh(h,f,g,w);
/* Procedure for Calculating Spectral Density Matrix for a VAR of the
   form:
     x(t)=h*s(t)
     s(t)=f*s(t-1) + g*e(t)

     where var(e(t))=I

   Output:  ss = spectrum of x at frequency w.
   Note: Spectrum is not divided by 2*pi

*/
local im, z, sm, smi, ss;


 let im = 0+1i;
 z=exp(-w*im);
 sm=eye(rows(f));
 sm=sm-z*f;
 smi=inv(sm);
 ss=h*smi*g*g'*smi'*h';

retp(ss);
endp;
