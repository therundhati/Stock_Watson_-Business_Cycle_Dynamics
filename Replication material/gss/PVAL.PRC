proc(1) = pval(ts,tabb,p);

/* Calculate p-value of test
   ts = value of test statistic
   tabb = CDF of statistic (TS Values)
   p = CDF of statistics (Probs)
*/
local d, di, pp;

d=tabb-ts*ones(1,cols(tabb));
d=abs(d');
di=minindc(d);
pp=p[di];
retp(1-pp);
endp;



