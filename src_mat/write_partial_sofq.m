function fwrite = write_partial_sofq(qvec,sofq_in,fwin,idins,alltypes)
  fprintf(fwin,'%s\t %s\t','q','Sofq_ids:');
  for i1 = 1:length(idins)
    fprintf(fwin,'%s\t',alltypes{idins(i1)});
  end
  fprintf(fwin,'\n')

  for k1 = 1:length(qvec)
    fprintf(fwin,'%g\t%g\n',qvec(k1),sofq_in(k1));
  end

  fwrite = 1;


