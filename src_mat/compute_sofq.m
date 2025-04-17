% Structure factor from radial distribution functions
% IUCR values: https://it.iucr.org/Cb/ch6o1v0001/

clc;
clear;
close all;


%% Input data
fe_charge    = 2; %or 3
lorch_window = 1; % 0 - No Lorch Window; 1 - Lorch Window; 2 - Lorch Window Amalie
norm_den     = 'add'; % 'denom','add','none'
molality     = cellstr({"0.1","0.2","0.4","1.0"}); % String type
nmolec_wat   = [5551; 5551; 5551; 5551]; % number of water molecules
nmolec_fe    = [10; 20; 40; 100]; % number of iron molecules
nmolec_tfsi  = fe_charge*nmolec_fe; % number of tfsi molecules
% lx: 1st row - intentionally -1, 2nd row - Fe+2; 3rd row - Fe+3
avg_boxlen   = [-1, -1, -1, -1;
                55.7297, 56.1992, 57.1109, 59.6802;
                55.9235, 56.5808, 57.8513, 61.3732]; %in angstroms

Rcut_arr     = avg_boxlen/2;
clr_arr      = {'g','b','r','m','c',"#ffa500",'k',"#f4a460"};

% Atom order in RDF files - SHOULD BE IN THE SAME ORDER AS RDF DATA
rcol  = 1;
type_str = {'Fe'; 'NTFSI'; 'OTFSI'; 'FTFSI'; 'CTFSI'; 'STFSI'; 'OWat'; 'HWat'};
ntypes   = length(type_str);

%% Compute factors using the fitting formula
qmin    = 4*pi/min(avg_boxlen(fe_charge,:));
qmax    = 1.3;
delq    = 0.005;
qvec    = qmin:delq:qmax;
fprintf('qmin: %g; qmax: %g, delq: %g\n',qmin,qmax,delq)

Fe2_coeff = [11.04240, 4.65380,	7.374000,	0.30530, 4.134600, 12.05460, 0.439900, 31.2809,	1.009700];
N_coeff   = [12.21260, 0.00570,	3.132200,	9.89330, 2.012500, 28.99750, 1.166300, 0.58260, -11.5290];
O_coeff   = [3.048500, 13.2771,	2.286800,	5.70110, 1.546300, 0.323900, 0.867000, 32.9089,	0.250800];
F_coeff   = [3.539200, 10.2825,	2.641200,	4.29440, 1.517000, 0.261500, 1.024300, 26.1476,	0.277600];
C_coeff   = [2.310000, 20.8439,	1.020000,	10.2075, 1.588600, 0.568700, 0.865000, 51.6512,	0.215600];
S_coeff   = [6.905300, 1.46790,	5.203400,	22.2151, 1.437900, 0.253600, 1.586300, 56.1720,	0.866900];
OW_coeff  = [3.048500, 13.2771,	2.286800,	5.70110, 1.546300, 0.323900, 0.867000, 32.9089,	0.250800]; % same as O_coeff
H_coeff   = [0.489918, 20.6593, 0.262003, 7.74039, 0.196767, 49.55190, 0.049879, 2.20159, 0.001305];
Fe3_coeff = [11.17640, 4.61470,	7.386300,	0.30050, 3.394800, 11.67290, 0.072400, 38.5566,	0.970700];

if fe_charge == 2
  Fe_bfac = compute_sofq_from_coeff(Fe2_coeff,qvec);
else
  Fe_bfac = compute_sofq_from_coeff(Fe3_coeff,qvec);
end


N_bfac   = compute_sofq_from_coeff(N_coeff,qvec);
O_bfac   = compute_sofq_from_coeff(O_coeff,qvec);
F_bfac   = compute_sofq_from_coeff(F_coeff,qvec);
C_bfac   = compute_sofq_from_coeff(C_coeff,qvec);
S_bfac   = compute_sofq_from_coeff(S_coeff,qvec);
OW_bfac  = compute_sofq_from_coeff(O_coeff,qvec); % same as O_bfac
H_bfac   = compute_sofq_from_coeff(H_coeff,qvec);

all_bqfacs = [Fe_bfac, N_bfac, O_bfac, F_bfac, C_bfac, S_bfac, OW_bfac, H_bfac];

if length(all_bqfacs(1,:)) != ntypes
  error("Wrong number of coefficients in all_bqfacs \n")
end


tot_sofq = zeros(length(qvec),length(molality(:)));  % Zero total structure factor

%% Main calculations
for molcnt = 1:length(molality(:)) %

  printf('Analyzing  Fe-%d at molality = %s M with lorch window type: %d and norm_type: %s\n',...
  fe_charge,char(molality{molcnt}),lorch_window,norm_den)

  % Fraction Calculations
  % Fe atoms
  natoms_Fe    = nmolec_fe(molcnt);

  % TFSI atoms
  natoms_N     = nmolec_tfsi(molcnt);
  natoms_OTFSI = nmolec_tfsi(molcnt)*4;
  natoms_F     = nmolec_tfsi(molcnt)*6;
  natoms_C     = nmolec_tfsi(molcnt)*2;
  natoms_S     = nmolec_tfsi(molcnt)*2;

  % Water atoms
  natoms_OW    = nmolec_wat(molcnt);
  natoms_H     = nmolec_wat(molcnt)*2;

  % Total atoms
  natoms_O     = natoms_OW + natoms_OTFSI;
  natoms_tfsi  = natoms_N  + natoms_OTFSI + natoms_F + natoms_C + natoms_S;
  natoms_wat   = natoms_OW + natoms_H;
  ntot_atoms   = natoms_wat + natoms_Fe + natoms_tfsi;


  if (natoms_H + natoms_C + natoms_N + natoms_O + natoms_F + natoms_S + natoms_Fe) != ntot_atoms
    printf('%d \t %d \n',natoms_H + natoms_C + natoms_N + natoms_O + natoms_F + natoms_S + natoms_Fe, ntot_atoms);
    error('Wrong total number of atoms \n')
  end

  at_frac = [natoms_Fe, natoms_N, natoms_OTFSI, natoms_F, natoms_C, natoms_S, natoms_OW, natoms_H];
  at_frac = at_frac/ntot_atoms;
  natoms_dens  = ntot_atoms/(avg_boxlen(fe_charge,molcnt)^3);
  Rcut_lorch = Rcut_arr(fe_charge,molcnt);

  % Compute denominator (CHECK the DERIVATION)
  normval = zeros(length(qvec),1);
  for dcnt = 1:ntypes
    normval  = normval + at_frac(dcnt)*all_bqfacs(:,dcnt);
  end
  normval = normval.*normval;

  % Zero sofq data
  Sofq = zeros(ntypes,ntypes,length(qvec));

  % Load RDF files
  for nt1 = 1:ntypes
    fname    = sprintf('./../../FeTFSI/results_all/rdf_fetfsi%d/rdfall_%s/rdf_%s_all.xvg',fe_charge,char(molality{molcnt}),char(type_str{nt1}));
    rdf_arr  = dlmread(fname,'',31,0); %skip 31 lines

    % Heart of the calculations
    for qcnt = 1:length(qvec)
      qval   = qvec(qcnt);

      for nt2 = 1:ntypes

        qprefac  = 4*pi*natoms_dens*at_frac(nt1)*at_frac(nt2)*all_bqfacs(qcnt,nt1)*all_bqfacs(qcnt,nt2)/qval;
        rarr     = 10*rdf_arr(:,1); % Factor 10 for conversion to angstroms
        integral_val = compute_integral(qval,rarr(2:length(rarr)),rdf_arr(2:length(rarr),nt2+1),Rcut_lorch,lorch_window);
        Sofq(nt1,nt2,qcnt) = qprefac*integral_val;

      end

    end

  end

  % Compute total structure factor
  for i1 = 1:ntypes
    for j1 = 1:ntypes
      for k1 = 1:length(qvec)
        tot_sofq(k1,molcnt) = tot_sofq(k1,molcnt) + Sofq(i1,j1,k1);
      end
    end
  end

  if strcmp(norm_den,'denom')
    tot_sofq(:,molcnt) = tot_sofq(:,molcnt)./normval;
  elseif strcmp(norm_den,'add')
    tot_sofq(:,molcnt) = tot_sofq(:,molcnt) + normval;
  elseif strcmp(norm_den,'none')
    tot_sofq(:,molcnt) = tot_sofq(:,molcnt);
  else
    error('Unknown norm_den keyword\n')
  endif

  %% Save to file
  fid = fopen(sprintf('./../../FeTFSI/sofq_results/sofq_fe_%d_mol_%s_lorwind_%d_norm_%s.dat',fe_charge,char(molality{molcnt}),lorch_window,norm_den),'w');
  fprintf(fid,'%s\t','q');
  for i1 = 1:ntypes
    for j1 = 1:ntypes
      fprintf(fid,'%s-%s\t',type_str{i1},type_str{j1});
    end
  end
  fprintf(fid,'TotSofq\n')

  for k1 = 1:length(qvec)
    fprintf(fid,'%g\t',qvec(k1));
    for i1 = 1:ntypes
      for j1 = 1:ntypes
        fprintf(fid,'%g\t',Sofq(i1,j1,k1));
      end
    end
    fprintf(fid,'%g\n',tot_sofq(k1,molcnt))
  end
  fclose(fid);

  %% Analyze water-water
  inpids = [7;8];
  fpid = fopen(sprintf('./../../FeTFSI/sofq_results/psofq_WW_%d_mol_%s_lorwind_%d_norm_%s.dat',fe_charge,char(molality{molcnt}),lorch_window,norm_den),'w');
  normvalWW = compute_norm_vals(length(qvec),at_frac,all_bqfacs,inpids);
  sofq_WW   = compute_partial_sofq(length(qvec),inpids,Sofq);
  outw = write_partial_sofq(qvec,sofq_WW,fpid,inpids,type_str);
  fclose(fpid);

  %% Analyze water-TFSI
  inpids = [2;3;4;5;6;7;8];
  fpid = fopen(sprintf('./../../FeTFSI/sofq_results/psofq_WA_%d_mol_%s_lorwind_%d_norm_%s.dat',fe_charge,char(molality{molcnt}),lorch_window,norm_den),'w');
  normvalWA = compute_norm_vals(length(qvec),at_frac,all_bqfacs,inpids);
  sofq_WA   = compute_partial_sofq(length(qvec),inpids,Sofq);
  outw = write_partial_sofq(qvec,sofq_WA,fpid,inpids,type_str);
  fclose(fpid)

  %% Analyze Water-Fe
  inpids = [1;7;8];
  fpid = fopen(sprintf('./../../FeTFSI/sofq_results/psofq_WF_%d_mol_%s_lorwind_%d_norm_%s.dat',fe_charge,char(molality{molcnt}),lorch_window,norm_den),'w');
  normvalWF = compute_norm_vals(length(qvec),at_frac,all_bqfacs,inpids);
  sofq_WF   = compute_partial_sofq(length(qvec),inpids,Sofq);
  outw = write_partial_sofq(qvec,sofq_WF,fpid,inpids,type_str);
  fclose(fpid);

  %% Analyze TFSI-TFSI
  inpids = [2;3;4;5;6];
  fpid = fopen(sprintf('./../../FeTFSI/sofq_results/psofq_AA_%d_mol_%s_lorwind_%d_norm_%s.dat',fe_charge,char(molality{molcnt}),lorch_window,norm_den),'w');
  normvalAA = compute_norm_vals(length(qvec),at_frac,all_bqfacs,inpids);
  sofq_AA   = compute_partial_sofq(length(qvec),inpids,Sofq);
  outw = write_partial_sofq(qvec,sofq_AA,fpid,inpids,type_str);
  fclose(fpid);

  %% Analyze TFSI-Fe
  inpids = [1;2;3;4;5;6];
  fpid = fopen(sprintf('./../../FeTFSI/sofq_results/psofq_AF_%d_mol_%s_lorwind_%d_norm_%s.dat',fe_charge,char(molality{molcnt}),lorch_window,norm_den),'w');
  normvalAF = compute_norm_vals(length(qvec),at_frac,all_bqfacs,inpids);
  sofq_AF   = compute_partial_sofq(length(qvec),inpids,Sofq);
  outw = write_partial_sofq(qvec,sofq_AF,fpid,inpids,type_str);
  fclose(fpid);

  %% Analyze Fe-Fe
  inpids = [1];
  fpid = fopen(sprintf('./../../FeTFSI/sofq_results/psofq_FF_%d_mol_%s_lorwind_%d_norm_%s.dat',fe_charge,char(molality{molcnt}),lorch_window,norm_den),'w');
  normvalFF = compute_norm_vals(length(qvec),at_frac,all_bqfacs,inpids);
  sofq_FF   = compute_partial_sofq(length(qvec),inpids,Sofq);
  outw = write_partial_sofq(qvec,sofq_FF,fpid,inpids,type_str);
  fclose(fpid);

  %% Plot partial structure factors
  h = figure;
  hold on
  box on
  grid on
  leg_arr = {};
  plot(qvec,sofq_WW,'color',clr_arr{1},'LineWidth',2)
  leg_arr{1} = 'Water-Water';
  plot(qvec,sofq_WA,'color',clr_arr{2},'LineWidth',2)
  leg_arr{2} = 'Water-TFSI';
  plot(qvec,sofq_WF,'color',clr_arr{3},'LineWidth',2)
  leg_arr{3} = 'Water-Fe';
  plot(qvec,sofq_AA,'color',clr_arr{4},'LineWidth',2)
  leg_arr{4} = 'TFSI-TFSI';
  plot(qvec,sofq_AF,'color',clr_arr{5},'LineWidth',2)
  leg_arr{5} = 'TFSI-Fe';
  plot(qvec,sofq_FF,'color',clr_arr{6},'LineWidth',2)
  leg_arr{6} = 'Fe-Fe';
  xlim([0.9*qmin 1.1*qmax])
  %ylim([0.9*min(min(tot_sofq)) 1.1*max(max(tot_sofq))])
  xlabel('q','FontSize',20)
  ylabel('I(q)','FontSize',20)
  set(gca, 'FontSize', 12)
  set(gca, 'xscale', 'log')
  %set(gca, 'yscale', 'log')
  set(gca, 'xtick', 0.1:0.3:1);
  set(gca,'XtickLabel', [0.1 0.4 0.7 1.0])
  legend(leg_arr,'location','bestoutside')
  saveas(h,sprintf('./../../FeTFSI/figures/partialsq_%d_mol_%s_%d_lorwind_%d',fe_charge,char(molality{molcnt}),fe_charge,lorch_window),'png')
  %close(h)

end



clr_arr = {'g','b','r','m'};
h = figure;
hold on
box on
grid on
for i = 1:length(molality(:))
  plot(qvec,tot_sofq(:,i),'color',clr_arr{i},'LineWidth',2)
  leg_arr{i} = molality{i};
end
xlim([0.9*qmin 1.1*qmax])
ylim([0.9*min(min(tot_sofq)) 1.1*max(max(tot_sofq))])
xlabel('q','FontSize',20)
ylabel('I(q)','FontSize',20)
set(gca, 'FontSize', 12)
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
set(gca, 'xtick', 0.1:0.3:1);
set(gca,'XtickLabel', [0.1 0.4 0.7 1.0])
legend(leg_arr,'location','bestoutside')
saveas(h,sprintf('./../../FeTFSI/figures/Fe_%d_lorwind_%d',fe_charge,lorch_window),'png')

