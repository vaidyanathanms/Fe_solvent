function normv = compute_norm_vals(lqvec,atfrac,bqfacs,idvals)

  normv = zeros(lqvec,1);

  for icnt = 1:length(idvals)

    normv  = normv + atfrac(icnt)*bqfacs(:,icnt);

  end
  normval = normv.*normv;

