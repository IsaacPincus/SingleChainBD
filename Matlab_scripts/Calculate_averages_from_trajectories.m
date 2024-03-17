clear variables

% if monarch
addpath('/home/ipincus/wm73/ipincus/Single-Chain/Code/single-chain-development/Matlab_scripts')

% Processor names, padded with zeros
procs_nums = 0:999;
procs = string(procs_nums);
procs(1:10) = "00" + procs(1:10);
procs(11:100) = "0" + procs(11:100);

% number of processors used
N_proccessors = 24;

% Options
VR = 0;
Rodlike = 0;

dt_string = '01';

% timestep width to use
dt_num = str2num(dt_string);

time = [];
config = [];
grad = [];
cofm = [];
time_VR = [];
config_VR = [];
grad_VR = [];

%%% Reading files

sim_params = read_inputc_SC();

sigma = sim_params.Q0s;
sqrtb = sim_params.sqrtb;
sr = sim_params.gdots;
dt = sim_params.dtsne(dt_num);
tol = sim_params.tol(dt_num);
NBeads = sim_params.NBeads;
C = sim_params.BendingStiffness;

step = 1;
tic
for i = 1:N_proccessors
    proc = procs(i);

    file_noVR = strcat('net_dt',dt_string,'_proc', proc, '.nc');
    time_proc = ncread(file_noVR, 'Time');
    config_proc = ncread(file_noVR, 'configuration');
    grad_proc = ncread(file_noVR, 'Gradient');
    cofm_proc = ncread(file_noVR, 'cofm');
    time = [time; time_proc(:,1:step:end)];
    config = [config; config_proc(:,1:step:end,:,:)];
    grad = [grad; grad_proc(:,1:step:end,:,:)];
    cofm = [cofm; cofm_proc(:,1:step:end,:)];

    if VR==1
        file_VR = strcat('net_VR_dt',dt_string,'_proc', proc, '.nc');
        time_VR_proc = ncread(file_VR, 'Time');
        config_VR_proc = ncread(file_VR, 'configuration');
        grad_VR_proc = ncread(file_VR, 'Gradient');
        time_VR = [time_VR; time_VR_proc];
        config_VR = [config_VR; config_VR_proc];
        grad_VR = [grad_VR; grad_VR_proc];
    end

end
toc

% 
% time = time(:,1:step:end);
% config = config(:,1:step:end,:,:);
% grad = grad(:,1:step:end,:,:);
% cofm = cofm(:,1:step:end,:);

if Rodlike
    sr = sr*(4*sigma^2);
    config = config/sigma;
    grad = grad*sigma;
    time = time/(4*sigma^2);
    dQ = sqrtb/sigma;
%     dt = dt/(4*sigma^2);
    tol = tol/sigma;
    config_VR = config_VR/sigma;
    grad_VR = grad_VR*sigma;
    time_VR = time_VR/(4*sigma^2);
end

t = time(1,:);
grad_size = size(grad);
N_trajectories = grad_size(1);
N_samples = grad_size(2);
N_dims = grad_size(3);
N_beads = grad_size(4);

save('matlab_data.mat', 't')

%%
% 
% f1 = figure();
% axes1 = gca;
% % axes1.YScale = 'log';
% % xlabel('$t / \lambda_R$', 'Interpreter', 'latex', 'FontSize', 24)
% % ylabel('$<\hat{Re}(0)\cdot \hat{Re}(t)>$', 'Interpreter', 'latex', 'FontSize', 32)
% fsize=20;
% % pbaspect([1. 1. 1]);
% hold(axes1,'on');
% box(axes1,'on');
% set(axes1,'FontSize',14,'LineWidth',2,'TickLength',[0.015 0.025]);
% hold on
% cofm = cofm - cofm(:,1,:);
% msd_trajs = vecnorm(cofm, 2, 3).^2;
% msd_avg = mean(msd_trajs(:,:));
% msd_err = std(msd_trajs(:,:))/sqrt(length(msd_trajs));
% errorbar(t, msd_avg, msd_err, 'ro')
% saveas(f1, 'msdh02.png')
% writematrix([t;msd_avg;msd_err], 'msdh02.txt');
save('matlab_data.mat', 'msd_avg', 'msd_err', '-append')

%% Stress tensor calculations

tau_flow = pagemtimes( permute(config(:,:,:,:,1), [3,5,4,1,2]),...
                       permute(grad(:,:,:,:,1), [5,3,4,1,2]));

if VR == 1
%     tau_EQ = mmx('mult', permute(config_VR(:,:,:,:,1), [3,5,4,1,2]),...
%                          permute(grad_VR(:,:,:,:,1), [5,3,4,1,2]));
    tau_EQ = pagemtimes( permute(config_VR(:,:,:,:,1), [3,5,4,1,2]),...
                         permute(grad_VR(:,:,:,:,1), [5,3,4,1,2]));
end

tau_total = tau_flow-tau_EQ;

tau12_flow = squeeze(sum(tau_flow(1,2,:,:,:),3));
tau12_EQ = squeeze(sum(tau_EQ(1,2,:,:,:),3));
tau12_total = tau12_flow-tau12_EQ;

eta_trajs = -squeeze(sum(tau_total(1,2,:,:,:),3))/sr;
eta = squeeze(mean(eta_trajs,1));
eta_err = std(eta_trajs,1,1)/sqrt(N_trajectories);
save('matlab_data.mat', 'eta', 'eta_err', '-append')

%% Rg and Re

gyration_tensor = permute(...
                   pagemtimes( permute(config, [3,4,1,2]), ...
                               permute(config, [4,3,1,2])),...
                               [3,4,1,2]);
                               
gyration_tensor = gyration_tensor/N_beads;
% particle_distances_from_origin_squared = sum(config.^2,3);
% Rg_squared = squeeze(sum(particle_distances_from_origin_squared,4))/N_beads;
Rg_squared = squeeze(gyration_tensor(:,:,1,1) + ...
                     gyration_tensor(:,:,2,2) + ...
                     gyration_tensor(:,:,3,3));
Rg = sqrt(Rg_squared);

Rgx = squeeze(gyration_tensor(:,:,1,1));
Rgy = squeeze(gyration_tensor(:,:,2,2));
Rgz = squeeze(gyration_tensor(:,:,3,3));

Rg_mean = squeeze(mean(Rg));
Rg_err = squeeze(std(Rg)/sqrt(N_trajectories));

Rgx_mean = squeeze(mean(Rgx));
Rgx_err = squeeze(std(Rgx)/sqrt(N_trajectories));
Rgy_mean = squeeze(mean(Rgy));
Rgy_err = squeeze(std(Rgy)/sqrt(N_trajectories));
Rgz_mean = squeeze(mean(Rgz));
Rgz_err = squeeze(std(Rgz)/sqrt(N_trajectories));

% f2 = figure();
% axes1 = gca;
% axes1.YScale = 'log';
% xlabel('$t / \lambda_R$', 'Interpreter', 'latex', 'FontSize', 24)
% ylabel('${\langle R_g^2\rangle }^{1/2}$', 'Interpreter', 'latex', 'FontSize', 32)
% fsize=20;
% % pbaspect([1. 1. 1]);
% hold(axes1,'on');
% box(axes1,'on');
% set(axes1,'FontSize',14,'LineWidth',2,'TickLength',[0.015 0.025]);
% hold on
% e1 = errorbar(time(1,1:end), Rg_mean(1:end), Rg_err(1:end));
% saveas(f2, 'Rgh02.png')
% writematrix([time(1,1:end);Rg_mean_trajs(1:end);Rg_err_trajs(1:end)], 'rgh02.txt');
save('matlab_data.mat', 'Rg_mean', 'Rg_err', ...
    'Rgx_mean', 'Rgx_err', ...
    'Rgy_mean', 'Rgy_err', ...
    'Rgz_mean', 'Rgz_err', ...
    '-append')

%% Included angle

connectors = diff(config, 1, 4);
Ql = squeeze(sqrt(connectors(:,:,1,:).^2 + connectors(:,:,2,:).^2 + connectors(:,:,3,:).^2));

Ql_time = squeeze(Ql(:,100:end,:));

xi = acos(squeeze(dot(connectors(:,:,:,2:end), connectors(:,:,:,1:end-1), 3))./...
            (Ql(:,:,2:end).*Ql(:,:,1:end-1)));

%% relaxation time from steady state

% sample_cut = 
samp1 = 1;

Re_dist = vecnorm(config(:,:,:,end)-config(:,:,:,1),2,3);
% re_vec = squeeze((config(:,:,:,end) - config(:,:,:,1)))./Re_dist;
% re_vec = squeeze((config(:,:,:,end) - config(:,:,:,1)));
re_vec = squeeze((config(:,samp1:end,:,end) - config(:,samp1:end,:,1)));
re_vec = re_vec./vecnorm(re_vec,2,3);

ut = zeros(N_samples, 1);
ut_err = zeros(N_samples,1);
% ut_alt = zeros(N_samples, 1);
% for tid=1:N_samples
%     ut(1) = ut(1) + sum(dot(re_vec(:,tid,:),re_vec(:,tid,:), 3));
% end
% ut(1) = ut(1)/(N_trajectories*N_samples);
temp = sum(dot(re_vec(:,:,:),re_vec(:,:,:), 3));
ut(1) = sum(temp, 'omitnan')/(N_trajectories*N_samples);
ut_err(1) = std(temp, 'omitnan')/sqrt(N_trajectories*N_samples);
for tid=2:N_samples
    clear temp
    temp = sum(dot(re_vec(:,tid:end,:),re_vec(:,1:end-tid+1,:), 3));
    ut(tid) = sum(temp, 'omitnan')/(N_trajectories*(N_samples-tid));
    ut_err(tid) = std(temp, 'omitnan')/sqrt(N_trajectories*(N_samples-tid));
end

Re_autocorr = ut;
Re_autocorr_err = ut_err;

% figure();
% axes1 = gca;
% axes1.YScale = 'log';
% xlabel('$t / \lambda_R$', 'Interpreter', 'latex', 'FontSize', 24)
% ylabel('$<\hat{Re}(0)\cdot \hat{Re}(t)>$', 'Interpreter', 'latex', 'FontSize', 32)
% fsize=20;
% % pbaspect([1. 1. 1]);
% hold(axes1,'on');
% box(axes1,'on');
% set(axes1,'FontSize',14,'LineWidth',2,'TickLength',[0.015 0.025]);
% hold on
% % plot(t, ut, 'ro');
% errorbar(t, ut, ut_err, 'ro', 'DisplayName', '15kbp DNA');
% hold off

save('matlab_data.mat', 'Re_autocorr', 'Re_autocorr_err', '-append')

%% Stress relaxation

ut = zeros(N_samples, 1);
ut_err = zeros(N_samples,1);
% ut_alt = zeros(N_samples, 1);
% for tid=1:N_samples
%     ut(1) = ut(1) + sum(dot(re_vec(:,tid,:),re_vec(:,tid,:), 3));
% end
% ut(1) = ut(1)/(N_trajectories*N_samples);
temp = -sum(tau12_flow);
ut(1) = sum(temp, 'omitnan')/(N_trajectories*N_samples);
ut_err(1) = std(temp/N_samples, 'omitnan')/sqrt(N_trajectories);
for tid=2:N_samples
    clear temp
%     temp = sum(dot(re_vec(:,tid:end,:),re_vec(:,1:end-tid+1,:), 3));
    temp = sum(tau12_flow(:,tid:end-tid+1));
    ut(tid) = sum(temp, 'omitnan')/(N_trajectories*(N_samples-tid));
    ut_err(tid) = std(temp/(N_samples-tid), 'omitnan')/sqrt(N_trajectories);
end

% figure();
% axes1 = gca;
% axes1.YScale = 'log';
% xlabel('$t / \lambda_R$', 'Interpreter', 'latex', 'FontSize', 24)
% ylabel('$<\hat{Re}(0)\cdot \hat{Re}(t)>$', 'Interpreter', 'latex', 'FontSize', 32)
% fsize=20;
% % pbaspect([1. 1. 1]);
% hold(axes1,'on');
% box(axes1,'on');
% set(axes1,'FontSize',14,'LineWidth',2,'TickLength',[0.015 0.025]);
% hold on
% % plot(t, ut, 'ro');
% errorbar(t, ut, ut_err, 'ro', 'DisplayName', '15kbp DNA');
% hold off

stress_autocorr = ut;
stress_autocorr_err = ut_err;

save('matlab_data.mat', 'stress_autocorr', 'stress_autocorr_err', '-append')

%% Plot of select trajectories

% Center of mass is meant to be set to 0, should be tiny
time_cut = 1;
Rcm = mean(config, 4);
connectors = diff(config, 1, 4);
Ql = squeeze(sqrt(connectors(:,:,1,:).^2 + connectors(:,:,2,:).^2 + connectors(:,:,3,:).^2));
Ql_time = squeeze(Ql(:,time_cut:end,:));

Q = linspace(min(Ql(:)), max(Ql(:)), 100);

% figure('Position', [500 200 600 550]);
% hold on
% axes1 = gca;
% fsize=20;
% pbaspect([1. 1. 1]);
% hold(axes1,'on');
% box(axes1,'on');
% set(axes1,'FontSize',16,'LineWidth',2,'TickLength',[0.015 0.025]);
% title(strcat('${\Delta t}^*_H = $', num2str(dt), ', ${tol}^*_R = $', num2str(tol)), ...
%     'Interpreter', 'Latex');
% hold on
% %psi Plot
% plot(Q, psiQ, 'DisplayName', ['$\delta Q = $', num2str(sqrtb)]);
% %Histogram
% h1 = histogram(Ql_time(:), Q, 'Normalization', 'pdf');
% 
% vals = h1.Values;
% psi_mids = 0.5*(psiQ(2:end)+psiQ(1:end-1));
% err_psiQ = sum((vals-psi_mids).^2);
% hold off

[Q_pdf, Q_edges] = histcounts(Ql_time(:), Q, 'Normalization', 'pdf');

save('matlab_data.mat', 'Q_pdf', 'Q_edges', '-append')

%%

xi = acos(squeeze(dot(connectors(:,:,:,2:end), connectors(:,:,:,1:end-1), 3))./...
                    (Ql(:,:,2:end).*Ql(:,:,1:end-1)));

theta = 0:pi/nbins:pi;

xi_time = squeeze(xi(:,time_cut:end,:));
xi_psi_eq = sin(theta)/2;
xi_ana = (C/(1-exp(-2*C)))*sin(theta).*exp(-C*(1-cos(theta)));

% figure('Position', [500 200 600 550]);
% hold on
% axes1 = gca;
% fsize=20;
% pbaspect([1. 1. 1]);
% hold(axes1,'on');
% box(axes1,'on');
% set(axes1,'FontSize',16,'LineWidth',2,'TickLength',[0.015 0.025]);
% title(strcat('${\Delta t}^*_H = $', num2str(dt), ', ${tol}^*_R =$ ', num2str(tol)), ...
%     'Interpreter', 'Latex');
% hold on
% 
% %Histogram
% h2 = histogram(xi_time(:), theta, 'Normalization', 'pdf', 'displayname', 'pdf');
% vals = h2.Values;
% psi_mids = 0.5*(xi_psi_eq(2:end)+xi_psi_eq(1:end-1));
% err_psi_xi = sum((vals-psi_mids).^2);
% %psi Plot
% % plot(theta, xi_psi_eq);
% plot(theta, xi_ana, 'displayname', 'analytical');
% 
% dim = [0.6 0.18 0.7 0.7];
% str = ['err = ', num2str(err_psi_xi)];
% annotation('textbox',dim,'String',str,'FitBoxToText',...
%     'on','Interpreter','latex','FontSize',14,'EdgeColor','None');
% [h,icons,plots,legend_text]=legend({},'Location','best','FontSize',20,...
%     'Interpreter','latex','Box','off');
% xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', 24)
% ylabel('$\psi(\xi)$', 'Interpreter', 'latex', 'FontSize', 40)
% xlim([0, pi]);
% 
% %     saveFile = fopen(outputFileName, 'w');
% %     fprintf(saveFile, '%g %g \n', [err_psiQ, err_psi_xi]);
% %     fclose(saveFile);
% 
% f = [err_psiQ, err_psi_xi];

[xi_pdf, xi_edges] = histcounts(xi_time(:), theta, 'Normalization', 'pdf');
costheta = mean(cos(xi_time(:)));

save('matlab_data.mat', 'xi_pdf', 'xi_edges', 'costheta', '-append')

%% link-link correlation
t_cut = 1;
u_vecs = diff(config,1,4)./vecnorm(diff(config,1,4),2,3);
correl = zeros(1,N_beads-2);
diff_count = zeros(1, N_beads-2);
for nu=1:N_beads-2
    for mu=nu+1:N_beads-1
        mu_min_nu = mu-nu;
        correl(mu_min_nu) = correl(mu_min_nu) + mean(mean(...
              dot(u_vecs(:,t_cut:end,:,nu), u_vecs(:,t_cut:end,:,mu), 3), ...
              1), 2);
        diff_count(mu_min_nu) = diff_count(mu_min_nu) + 1;
    end
end
correl = correl./diff_count;

% figure();
% axes1 = gca;
% axes1.YScale = 'log';
% xlabel('$(\nu - \mu)$', 'Interpreter', 'latex', 'FontSize', 24)
% ylabel('${\langle u_\nu \cdot u_\mu \rangle }$', 'Interpreter', 'latex', 'FontSize', 32)
% fsize=20;
% % pbaspect([1. 1. 1]);
% hold(axes1,'on');
% box(axes1,'on');
% set(axes1,'FontSize',14,'LineWidth',2,'TickLength',[0.015 0.025]);
% hold on
% plot(1:N_beads-2, correl)

link_link_correl = correl;
save('matlab_data.mat', 'link_link_correl', '-append')

%% lim R dot u0 q

% time_min = 400;
% time_max = 440;
% 
% for i=1:N_beads
%     Re_vec = config(:,time_min:time_max,:,i) - config(:,time_min:time_max,:,1);
%     u0 = config(:,time_min:time_max,:,2) - config(:,time_min:time_max,:,1);
%     u0 = u0./vecnorm(u0,2,3);
%     q = dot(u0, Re_vec,3);
% 
%     per_len(i) = mean(mean(q));
% end
% 
% figure();
% plot(1:N_beads, per_len)

%% plotting dumbbell dists
% [p_ks, diffs, diffs_errs] = ...
%                     eq_dist_plot_dumbbell(config*sigma, sim_params, 1);

                
%% end to end distribution func

time_min = 1;
time_max = N_samples;
Re_dist = vecnorm(config(:,time_min:time_max,:,end)-config(:,time_min:time_max,:,1),2,3);
Re = Re_dist(:);

nbins = 100;
Re_bins = linspace(min(Re_dist(:)), max(Re_dist(:)), nbins);
Re_mids = (Re_bins(1:end-1) + Re_bins(2:end))/2;
[counts, bins] = histcounts(Re, Re_bins);
count_errors = sqrt(counts);
Re_dist_pdf = counts;
normalisation = trapz(Re_mids,counts);
psiRe_sim = counts/normalisation;
psiRe_sim_err = count_errors/normalisation;

% figure('Position', [500 200 840 670]);
% hold on
% axes1 = gca;
% fsize=20;
% % pbaspect([1. 1. 1]);
% hold(axes1,'on');
% box(axes1,'on');
% set(axes1,'FontSize',16,'LineWidth',2,'TickLength',[0.015 0.025]);
% errorbar(Re_mids, psiRe_sim, psiRe_sim_err, 'bo-', 'Linewidth', 1.5,...
%                 'DisplayName', 'Simulation PDF', 'Markersize', 4,...
%                 'MarkerFaceColor', 'b');
% 
% [h,icons,plots,legend_text]=legend({},'Location','best','FontSize',16,'Interpreter','latex','Box','off');
% xlabel('$Q^*_R$', 'Interpreter', 'latex', 'FontSize', 24)
% ylabel('$\psi_Q$, PDF', 'Interpreter', 'latex', 'FontSize', 32)

save('matlab_data.mat', 'psiRe_sim', 'psiRe_sim_err', '-append')
