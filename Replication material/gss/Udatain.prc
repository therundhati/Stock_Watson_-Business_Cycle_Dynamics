proc(1)=udatain(sy,missc,tobs);
/* udatain.prc, 5/24/94, mww
   This proc reads in data from file sy.
   It eliminates data in y when data is missing.
   Input:
         sy = string with file name for series y
         missc = missing value code
         tobs = total number of obs in series
   Output
          y = series y output
*/
local y, iy, ii;
load y[tobs,1]=^sy;
y=miss(y,missc);    @ -- Convert to GAUSS Missing Value Code -- @
retp(y);
endp;
