%% To fit diffusivity data
clear;
close all;
clc;
pkg load optim;
clr_arr      = {'g','b','r','m'};

%% Inputs
molality = {'0.1','0.2','0.4','1.0'};
trial_num = 4;
fecharge  = 3; % 2 or 3

% Directory names
maindirname = sprintf('../../FeTFSI/results_all_trial%d/diff_fetfsi%d',trial_num,fecharge);
outdirname  = sprintf('../../FeTFSI/analyzed_results/diff_all');
figdirname  = sprintf('../../FeTFSI/figures');

% Write output to file
if ~exist(outdirname, 'dir')
    mkdir(outdirname);
    disp(['Directory "', outdirname, '" created.']);
else
    disp(['Directory "', outdirname, '" already exists.']);
end
fout = fopen(sprintf('%s/computed_diffusivity_Fe%d.dat',outdirname,fecharge),'w');
fprintf(fout,'%s\t %s\t %s\n','Molality','D_TFSI(x10^(8)m^(2)/s)', 'D_TFSI_err(x10^(8)m^(2)/s)')

% Load Diffusivity files
h=figure;
hold on
box on
grid on
leg_arr={};
xlabel('Time (ps)','FontSize',20)
ylabel('MSD (Å^2)','FontSize',20,'interpreter','Tex')
set(gca, 'FontSize', 12)

diffusivity_vals = zeros(length(molality),1);
differr_vals = zeros(length(molality),1);

for fcnt = 1:length(molality)

  printf('Analyzing molality: %s\n', molality{fcnt});
  diff_fname = sprintf('%s/iondiff_com_%s.lammpstrj',maindirname,molality{fcnt});
  if exist(diff_fname, 'file') ~= 2
    disp(sprintf('WARNING: %s does not exist.', diff_fname));
    continue
  end

  % Trajfile
  diff_arr  = dlmread(diff_fname,'',1,0); %skip header line

  % File parameters
  lendata = length(diff_arr(:,1));
  fit_start = floor(0.1*lendata); fit_end = floor(0.8*lendata);

  xdata  = diff_arr(:,1);
  ydata  = 1/6*(diff_arr(:,2) + diff_arr(:,3) + diff_arr(:,4));
  [coeffs,errvals] = compute_fitmsd(xdata,ydata,fit_start,fit_end);
  slope = coeffs(1); intercept = coeffs(2);

  xfitdata = xdata(fit_start:fit_end);
  yfitdata = slope*xfitdata + intercept;
  residuals = ydata(fit_start:fit_end,1) - yfitdata;


  plot(xdata-xdata(1,1),ydata,'color',clr_arr{fcnt},'LineWidth',2,'LineStyle','-')
  plot(xfitdata-xdata(1,1),yfitdata,'color',clr_arr{fcnt},'LineWidth',2,'LineStyle','--')
  leg_arr{2*fcnt-1} = ['Data f = ' molality{fcnt}];
  leg_arr{2*fcnt}   = ['Fitdata f = ' molality{fcnt}];

  % Write to file
  fprintf(fout,'%s\t %g\t %g\n', molality{fcnt},slope,errvals(1))

  diffusivity_vals(fcnt,1) = slope;
  differr_vals(fcnt,1) = errvals(1);

end
fclose(fout);
legend(leg_arr,'location','bestoutside');
saveas(h,sprintf('%s/msd_Fe%d',figdirname,fecharge),'png')


