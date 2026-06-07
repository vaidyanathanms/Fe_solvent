function [coeffs,errvals] = compute_fitmsd(time_arr,diff_arr,fit_start,fit_end)
  % diffusivity
  xfitdata = time_arr(fit_start:fit_end,1);
  yfitdata = diff_arr(fit_start:fit_end,1);
  [coeffs,s] = polyfit(xfitdata, yfitdata, 1);
  errvals = sqrt (diag (s.C)/s.df)*s.normr;


