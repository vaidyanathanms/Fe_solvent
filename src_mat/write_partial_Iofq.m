function fwrite = write_partial_Iofq(qvec,Iofq_in,fwin,idinsA,idinsB,alltypes)

  fprintf(fwin,'%s\t %s\t','q','Sofq_ids_TypeA:');
  for i1 = 1:length(idinsA)
    fprintf(fwin,'%s\t',alltypes{idinsA(i1)});
  end
  fprintf(fwin,'%s\t','Sofq_ids_TypeB:');
  for i1 = 1:length(idinsB)
    fprintf(fwin,'%s\t',alltypes{idinsB(i1)});
  end


  fprintf(fwin,'\n')

  for k1 = 1:length(qvec)
    fprintf(fwin,'%g\t%g\n',qvec(k1),Iofq_in(k1));
  end

  fwrite = 1;


