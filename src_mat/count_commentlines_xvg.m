function nxvg_comments = count_commentlines_xvg(filename)
  nxvg_comments = 0;
  fxvg = fopen(filename);

  while ~feof(fxvg)
    line = fgetl(fxvg);
     if ischar(line) && !startsWith(line, "#") && !startsWith(line, "@")
       break
     else
       nxvg_comments++;
     end
   end
   fclose(fxvg);

