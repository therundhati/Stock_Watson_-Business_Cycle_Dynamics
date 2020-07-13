@ rwksmth.prc, 4/8/96, mww
  This is the kalman smmother for the model:

  y(t) = x(t)*b(t) + e(t)   (Scalar)
  b(t) = b(t-1) + a(t) (Scalar)

  where var(e)=r and var(a)=q;

  Input:  b1 = b(t/t)
          b2 = b(t+1/t)
          b3 = b(t+1/T);
          p1 = p(t/t)
          p2 = p(t+1/t)
          p3 = p(t+1/T);

          b1 = b(t-1/t-1)
          p1 = p(t-1/t-1)

  output: b3 = b(t/T)
          p3 = p(t/T)

@
proc(2) = rwksmth(b1,b2,b3,p1,p2,p3);

@ Kalman Smoother  @
 local p2i, as;
 p2i=1/p2;
 as=p1*p2i;
 p3=p1+as*(p3-p2)*as';
 b3=b1+as*(b3-b2);
retp(b3,p3);
endp;
