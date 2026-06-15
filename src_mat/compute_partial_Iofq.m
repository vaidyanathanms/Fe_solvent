function partialIq = compute_partial_Iofq(lqvec,idvalsA,idvalsB,Iofqall)
  partialIq = zeros(lqvec,1);
  for k11 = 1:length(idvalsA)
    for k12 = 1:length(idvalsB)
      for k13 = 1:lqvec
        partialIq(k13,1) = partialIq(k13,1) + Iofqall(idvalsA(k11),idvalsB(k12),k13);
      end
    end
  end

