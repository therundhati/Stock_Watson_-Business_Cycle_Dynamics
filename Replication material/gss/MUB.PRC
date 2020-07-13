proc(1) = mub(tstat,lam,tabt);

/* -- Procedure for Determining MUB estimator in look-up table -- */

local d, db, ac, r1, dd, dn, f, lamint;
if tstat .<= tabt[1];
 retp(lam[1]);
elseif tstat .>= tabt[rows(tabt)];
 retp(lam[rows(lam)]);
else;

 d=tstat*ones(cols(tabt),1);
 db=(d .<= tabt) - (d .> tabt);
 ac=cumsumc(db);
 r1=minindc(ac);

 dd=tabt[r1+1]-tabt[r1];
 dn=tstat-tabt[r1];
 f=dn/dd;
 lamint=lam[r1]+f*(lam[r1+1]-lam[r1]);
 retp(lamint);
endif;
endp;
