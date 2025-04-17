function calc_vec = compute_sofq_from_coeff(coeffs,qvec)
  %Ref: https://lampz.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
  %Ref: http://it.iucr.org/Cb/ch6o1v0001/

  if length(coeffs) != 9
    error('Wrong number of ooefficients\n')
  endif
  fac = 16*pi*pi;
  % Account for the conversion between sin(\theta)/\lambda to q; q = 4*pi*sin(theta)/\lambda
  calc_vec = zeros(length(qvec),1);
  for qncnt = 1:length(qvec)
    qval = qvec(qncnt);
    for i = 1:2:7
      calc_vec(qncnt,1) = calc_vec(qncnt,1) + coeffs(i)*exp(-((coeffs(i+1)*qval^2)/fac));
    end
    calc_vec(qncnt,1) = calc_vec(qncnt,1) + coeffs(9);
  end


