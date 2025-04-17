function integral_val = compute_integral(qval,rdata,gofrdata,maxR,lorch_wind)

  %Using Lorch type window
  if (lorch_wind == 1)
    extra_fac    = (pi*rdata)/maxR;
    wofr         = sin(extra_fac);
  elseif (lorch_wind == 2)
    extra_fac    = (pi*rdata)/maxR;
    exfac_inv    = 1./extra_fac;
    wofr         = 3*(exfac_inv.^3).*(sin(extra_fac)-extra_fac.*cos(extra_fac));
  else
    wofr         = 1;
  end
  yval         = rdata.*wofr.*sin(qval.*rdata).*(gofrdata - 1.0);
  main_int     = trapz(rdata,yval);
  %main_int     = sum(yval)*(rdata(2)-rdata(1));
  integral_val = main_int;




##  if (lorch_wind == 1)
##    extra_fac    = (pi*rdata)/maxR;
##    yval         = sin(extra_fac).*sin(qval.*rdata).*(gofrdata - 1.0);
##    prefac       = 4*maxR/qval; % pi*r cancels out
##    main_int     = trapz(rdata,yval);
##    integral_val = prefac*main_int;
##  elseif (lorch_wind == 2)
##    extra_fac    = (pi*rdata)/maxR;
##    exfac_inv    = 1./extra_fac;
##    yval         = 3*(exfac_inv.^3).*(sin(extra_fac)-extra_fac).*rdata.*sin(qval.*rdata).*(gofrdata - 1.0);
##    prefac       = 4*pi/qval;
##    main_int     = trapz(rdata,yval);
##    integral_val = prefac*main_int;
##  else
##    yval         = rdata.*sin(qval*rdata).*(gofrdata - 1.0);
##    prefac       = 4*pi/qval;
##    main_int     = trapz(rdata,yval);
##    integral_val = prefac*main_int;
##  end




