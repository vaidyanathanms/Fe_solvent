function partialSq = compute_partial_sofq(lqvec,idvals,Sofqall)
  partialSq = zeros(lqvec,1);
  for k11 = 1:length(idvals)
    for k12 = 1:length(idvals)
      for k13 = 1:lqvec
        partialSq(k13,1) = partialSq(k13,1) + Sofqall(idvals(k11),idvals(k12),k13);
      end
    end
  end

