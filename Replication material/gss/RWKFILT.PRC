@ rwkfilt.prc, 4/6/96, mww
  This is the kalman Filter for the model:

  y(t) = b(t) + e(t)   (Scalar)
  b(t) = b(t-1) + a(t) (Scalar)

  where var(e)=r and var(a)=q;

  Input:  y(t)
          b1 = b(t-1/t-1)
          p1 = p(t-1/t-1)
          r
          q

  output: b1 = b(t/t)
          p1 = p(t/t)
          n  = y(t)-y(t/t-1)
          h  = var(n)

@
proc(6) = rwkfilt(y,b1,p1,r,q);
local b2, p2, k, h, n;

 b2=b1;
 p2=p1+q;
 n=y-b2;
 h=p2+r;
 k=p2/h;
 b1=b2+k*n;
 p1=p2-k*p2;

retp(b1,p1,b2,p2,n,h);
endp;

