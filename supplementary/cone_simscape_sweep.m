function cone_simscape_sweep()
% CONE_SIMSCAPE_SWEEP  Friction-cone half-angle validation across mu.
%
% Builds a FRESH Simscape Multibody model for each friction coefficient mu,
% with mu baked into the Mu Gain block at creation time.  This defeats the
% Simulink compile-cache problem: set_param on a Gain block between sims is
% confirmed (get_param shows it) but the compiled simulation keeps using the
% first value.  Rebuilding from scratch per-mu (5 builds, cheap) is the only
% reliable approach.
%
% CONTACT MODEL (IMPORTANT -- differs from robust_motion_cone_simscape.m):
%   This mu-sweep uses a CLEAN one-sided spring-damper  Fn = max(0, Kn*GapX +
%   Dn*vrelX)  with NO Fmin floor and NO pre-compression offset.  The main
%   validation's Fmin floor (kept there to sustain steady-state contact for the
%   full motion cone) would spuriously push the slider in x (runaway vx ~
%   Fmin/Dbase) and mask the friction cone for mu >= 0.30 -- every mu would
%   saturate near atan(0.20).  Removing the floor lets the pusher's x-penetration
%   set Fn so the friction cone atan(mu) emerges exactly (see build_cone_model
%   below, and simscape_validation_status.md, Lab 2).
%
% Sweep:
%   mu      in {0.20, 0.30, 0.40, 0.50, 0.60}
%   angle   in a dense grid (4 .. 70 deg)
% For each mu we detect the effective cone half-angle from the direction
% plateau and compare it to the theory atan(mu). The result is a
% cone-half-angle vs mu figure (sim vs theory) saved to ral_submission/figures.
%
% Requires: Simscape Multibody (sm_lib), nesl_utility.

fprintf('=== Cone Simscape Sweep (per-mu rebuild) ===\n\n');

%% ------------------------------------------------------------------
%  Shared parameters (must match robust_motion_cone_simscape.m)
%  ------------------------------------------------------------------
sp = struct();
sp.Dbase = 20;        % base damping [N/(m/s)]
sp.Kn    = 10000;     % normal contact stiffness [N/m]
sp.Dn    = 50;        % normal contact damping [N/(m/s)]
sp.Kt    = 5000;      % tangential (stick) stiffness [N/m]
sp.Dt    = 100;       % tangential damping [N/(m/s)]
sp.pusher_speed = 0.05;  % pusher speed [m/s]
sp.sim_time     = 3.0;   % time per run [s]
sp.pre_pen      = 2e-4;  % pre-compression [m]
sp.Fmin         = 5.0;   % min sustained normal force [N]
contact_arm = 0.05;      % pusher contact arm offset [m]

% Sweep grid
mu_vals      = [0.20 0.30 0.40 0.50 0.60];
test_angles  = [4 6 8 10 12 14 16 18 20 22 24 26 28 30 35 40 45 50 55 60 65 70];
n_mu         = numel(mu_vals);
n_angles     = numel(test_angles);

% NOTE: a UNIQUE model name is generated per mu (inside the loop) so that
% Simulink's simulation compile-cache (keyed by model name) cannot reuse the
% first mu's compiled binary.  Reusing one name leaves get_param showing the
% new mu while the compiled sim silently keeps the first value.

%% ------------------------------------------------------------------
%  Sweep: rebuild model per mu, angle sweep on the in-place model
%  ------------------------------------------------------------------
sweep = struct();
sweep.mu          = mu_vals;
sweep.half_sim    = nan(n_mu,1);          % effective cone half-angle from sim [deg]
sweep.half_theory = rad2deg(atan(mu_vals));
sweep.angle       = test_angles;
sweep.vdir        = nan(n_mu, n_angles);   % slider motion direction per (mu,angle)

for k = 1:n_mu
    mu = mu_vals(k);
    mdl = sprintf('coneSweep_mu%02d', k);   % unique name -> fresh compile per mu
    fprintf('--- mu = %.2f  (theory half-angle atan(mu) = %.2f deg) ---\n', mu, rad2deg(atan(mu)));
    % Build a FRESH model for this mu value.  mu is baked into the Mu Gain
    % block at creation; the compiled sim then uses the correct value.
    build_cone_model(mdl, mu, sp, contact_arm);
    fprintf('  [diag] Mu block Gain param = %s (asked %.4f)\n', ...
        get_param([mdl '/Mu'],'Gain'), mu);

    for j = 1:n_angles
        theta = test_angles(j);
        tr = deg2rad(theta);
        vx = sp.pusher_speed * cos(tr);
        vy = sp.pusher_speed * sin(tr);
        % Re-parameterize pusher velocity targets (in-place, no rebuild)
        set_param([mdl '/PJ'], ...
            'PxVelocityTargetValue',sprintf('%.6f',vx), ...
            'PyVelocityTargetValue',sprintf('%.6f',vy));

        try
            out = sim(mdl);
            T = out.Tvec(:); svx = out.Svx(:); svy = out.Svy(:);
            n = length(T);
            if n < 10
                fprintf('  %4.0f deg | too few samples\n', theta);
                continue;
            end
            i_ss = round(0.6*n):n;
            vx_ss = mean(svx(i_ss),'omitnan');
            vy_ss = mean(svy(i_ss),'omitnan');
            vdir  = rad2deg(atan2(vy_ss, max(abs(vx_ss),1e-10)));
            sweep.vdir(k,j) = vdir;
            fprintf('  %4.0f deg | dir = %6.2f deg\n', theta, vdir);
        catch ME
            msg = ME.message;
            if numel(msg) > 80, msg = msg(1:80); end
            fprintf('  %4.0f deg | FAILED: %s\n', theta, msg);
        end
    end

    % Effective cone half-angle = median |vdir| in the deep-sliding region
    deep = abs(test_angles) > 55;
    if any(deep & isfinite(sweep.vdir(k,:)))
        sweep.half_sim(k) = median(abs(sweep.vdir(k,deep & isfinite(sweep.vdir(k,:)))));
    end
    fprintf('  -> effective half-angle (sim) = %.2f deg | theory = %.2f deg\n\n', ...
        sweep.half_sim(k), sweep.half_theory(k));
    % Tear down this mu's model AND its Simulink accelerator cache so the
    % next iteration's unique-named model gets a guaranteed-fresh compile.
    close_system(mdl,0);
    if exist([mdl '.slx'],'file'), delete([mdl '.slx']); end
    if exist([mdl '.slxc'],'dir'), rmdir([mdl '.slxc'],'s'); end
end

%% ------------------------------------------------------------------
%  Machine-readable export (for code/experiments/make_fig7_composite.py)
%  ------------------------------------------------------------------
export = struct('mu', sweep.mu, 'theory_deg', sweep.half_theory, ...
    'sim_deg', sweep.half_sim, 'err_deg', sweep.half_sim(:) - sweep.half_theory(:));
export_path = fullfile(fileparts(mfilename('fullpath')), 'mu_sweep_results.json');
fid = fopen(export_path, 'w');
fprintf(fid, '%s', jsonencode(export));
fclose(fid);
fprintf('Exported: %s\n', export_path);

%% ------------------------------------------------------------------
%  Figure: cone half-angle vs mu (sim vs theory)
%  ------------------------------------------------------------------
hf = figure('Name','Cone half-angle vs mu','Position',[100 100 900 650]);

% (a) Half-angle vs mu
subplot(2,1,1);
plot(sweep.mu, sweep.half_theory, 'k--o','LineWidth',1.5,'MarkerSize',7,'MarkerFaceColor','k'); hold on;
plot(sweep.mu, sweep.half_sim,    'bs-','LineWidth',1.5,'MarkerSize',7,'MarkerFaceColor','b');
xlabel('Friction coefficient \mu');
ylabel('Cone half-angle (deg)');
legend('Theory atan(\mu)','Simscape (effective)', 'Location','northwest');
grid on; box on;
title('(a) Friction-cone half-angle: simulation vs theory','FontWeight','normal');

% (b) Direction response for each mu
subplot(2,1,2);
cm = lines(n_mu);
for k = 1:n_mu
    plot(test_angles, sweep.vdir(k,:), '-o','Color',cm(k,:),'LineWidth',1.2,'MarkerSize',4,'MarkerFaceColor',cm(k,:)); hold on;
end
plot(test_angles, test_angles, '--','Color',[0.6 0.6 0.6],'LineWidth',0.8);
leg = cell(1,n_mu+1);
for k=1:n_mu, leg{k} = sprintf('\\mu=%.2f', sweep.mu(k)); end
leg{end} = 'Push dir';
legend(leg, 'Location','northwest','FontSize',8);
xlabel('Push direction, \theta (deg)');
ylabel('Slider motion direction (deg)');
grid on; box on;
title('(b) Slider direction response across \mu','FontWeight','normal');

sgtitle('Simscape Sweep: Motion-Cone Half-Angle vs Friction (per-mu rebuild)', ...
    'FontWeight','bold','FontSize',12);

fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
fig_dir = fullfile(fig_dir);
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
fig_path = fullfile(fig_dir, 'simscape_cone_halfangle_vs_mu.png');
saveas(hf, fig_path);
fprintf('Saved: %s\n', fig_path);

% Summary table
fprintf('\n========================================\n');
fprintf('  mu    theory(deg)   sim(deg)    err(deg)\n');
for k = 1:n_mu
    e = sweep.half_sim(k) - sweep.half_theory(k);
    fprintf('  %.2f   %8.2f   %8.2f   %8.2f\n', sweep.mu(k), sweep.half_theory(k), sweep.half_sim(k), e);
end
fprintf('========================================\n');
fprintf('Cone Simscape Sweep Complete.\n');
end

%% ====================================================================
%  build_cone_model  - build a FRESH cone model with mu baked into the
%  Mu Gain block at creation time (defeats Simulink compile-cache).
%  ====================================================================
function build_cone_model(mdl, mu, sp, contact_arm)
    if bdIsLoaded(mdl), close_system(mdl,0); end
    if exist([mdl '.slx'],'file'), delete([mdl '.slx']); end
    new_system(mdl);

    add_block('nesl_utility/Solver Configuration', [mdl '/SC']);
    add_block('sm_lib/Frames and Transforms/World Frame', [mdl '/W']);
    add_block('sm_lib/Utilities/Mechanism Configuration', [mdl '/M']);
    add_block('sm_lib/Joints/Planar Joint', [mdl '/SJ']);
    set_param([mdl '/SJ'], ...
        'PxSensePosition','on','PxSenseVelocity','on', ...
        'PySensePosition','on','PySenseVelocity','on', ...
        'PxDampingCoefficient',sprintf('%.1f',sp.Dbase), ...
        'PyDampingCoefficient',sprintf('%.1f',sp.Dbase));
    add_block('sm_lib/Body Elements/Brick Solid', [mdl '/S']);
    set_param([mdl '/S'],'BrickDimensions','[0.5 0.3 0.1]','Density','1200');

    % Pusher
    add_block('sm_lib/Joints/Planar Joint', [mdl '/PJ']);
    set_param([mdl '/PJ'], ...
        'PxSensePosition','on','PxSenseVelocity','on', ...
        'PySensePosition','on','PySenseVelocity','on', ...
        'PxVelocityTargetSpecify','on','PxVelocityTargetValue','0', ...
        'PxVelocityTargetValueUnits','m/s', ...
        'PyVelocityTargetSpecify','on','PyVelocityTargetValue','0', ...
        'PyVelocityTargetValueUnits','m/s');
    add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/PO']);
    set_param([mdl '/PO'],'TranslationCartesianOffset',sprintf('[0 0 %.3f]',contact_arm));
    add_block('sm_lib/Body Elements/Brick Solid', [mdl '/P']);
    set_param([mdl '/P'],'BrickDimensions','[0.04 0.02 0.08]','Density','1200');

    % External force on slider
    add_block('sm_lib/Forces and Torques/External Force and Torque', [mdl '/EFT']);
    set_param([mdl '/EFT'],'EnableForceX','on','EnableForceY','on');

    % Topology
    add_line(mdl,'W/RConn1','SJ/LConn1');
    add_line(mdl,'SJ/RConn1','S/RConn1');
    add_line(mdl,'W/RConn1','PJ/LConn1');
    add_line(mdl,'PJ/RConn1','PO/LConn1');
    add_line(mdl,'PO/RConn1','P/RConn1');
    add_line(mdl,'SC/RConn1','W/RConn1');
    add_line(mdl,'S/RConn1','EFT/RConn1');

    % X sensors
    add_block('nesl_utility/PS-Simulink Converter',[mdl '/Spx']); set_param([mdl '/Spx'],'Unit','m'); add_line(mdl,'SJ/RConn2','Spx/LConn1');
    add_block('nesl_utility/PS-Simulink Converter',[mdl '/Svx']); set_param([mdl '/Svx'],'Unit','m/s'); add_line(mdl,'SJ/RConn3','Svx/LConn1');
    add_block('nesl_utility/PS-Simulink Converter',[mdl '/Ppx']); set_param([mdl '/Ppx'],'Unit','m'); add_line(mdl,'PJ/RConn2','Ppx/LConn1');
    add_block('nesl_utility/PS-Simulink Converter',[mdl '/Pvx']); set_param([mdl '/Pvx'],'Unit','m/s'); add_line(mdl,'PJ/RConn3','Pvx/LConn1');
    % Y sensors
    add_block('nesl_utility/PS-Simulink Converter',[mdl '/Spy']); set_param([mdl '/Spy'],'Unit','m'); add_line(mdl,'SJ/RConn4','Spy/LConn1');
    add_block('nesl_utility/PS-Simulink Converter',[mdl '/Svy']); set_param([mdl '/Svy'],'Unit','m/s'); add_line(mdl,'SJ/RConn5','Svy/LConn1');
    add_block('nesl_utility/PS-Simulink Converter',[mdl '/Ppy']); set_param([mdl '/Ppy'],'Unit','m'); add_line(mdl,'PJ/RConn4','Ppy/LConn1');
    add_block('nesl_utility/PS-Simulink Converter',[mdl '/Pvy']); set_param([mdl '/Pvy'],'Unit','m/s'); add_line(mdl,'PJ/RConn5','Pvy/LConn1');

    % Normal force: Fn = Kn*GapX + Dn*vrelX, ONE-SIDED (>= 0).
    % NO Fmin floor and NO PrePen offset: a constant normal-force floor
    % spuriously pushes the slider in x (runaway vx ~ Fmin/Dbase) and masks
    % the friction cone for mu >= 0.30.  The clean one-sided spring-damper
    % lets the pusher's x-penetration set Fn so the cone atan(mu) emerges.
    add_block('simulink/Math Operations/Add',[mdl '/GapX']); set_param([mdl '/GapX'],'Inputs','+-');
    add_line(mdl,'Ppx/1','GapX/1'); add_line(mdl,'Spx/1','GapX/2');
    add_block('simulink/Math Operations/Gain',[mdl '/Kx']); set_param([mdl '/Kx'],'Gain',sprintf('%.0f',sp.Kn)); add_line(mdl,'GapX/1','Kx/1');
    add_block('simulink/Math Operations/Add',[mdl '/VrelX']); set_param([mdl '/VrelX'],'Inputs','+-');
    add_line(mdl,'Pvx/1','VrelX/1'); add_line(mdl,'Svx/1','VrelX/2');
    add_block('simulink/Math Operations/Gain',[mdl '/Dx']); set_param([mdl '/Dx'],'Gain',sprintf('%.0f',sp.Dn)); add_line(mdl,'VrelX/1','Dx/1');
    add_block('simulink/Math Operations/Add',[mdl '/SumFn'],'Inputs','2');
    add_line(mdl,'Kx/1','SumFn/1'); add_line(mdl,'Dx/1','SumFn/2');
    add_block('simulink/Discontinuities/Saturation',[mdl '/ClampFn']);
    set_param([mdl '/ClampFn'],'LowerLimit','0','UpperLimit','inf');
    add_line(mdl,'SumFn/1','ClampFn/1');

    % Tangential force (elastic stick + damping), Coulomb-limited by mu
    add_block('simulink/Math Operations/Add',[mdl '/GapY']); set_param([mdl '/GapY'],'Inputs','+-');
    add_line(mdl,'Ppy/1','GapY/1'); add_line(mdl,'Spy/1','GapY/2');
    add_block('simulink/Math Operations/Gain',[mdl '/Ky']); set_param([mdl '/Ky'],'Gain',sprintf('%.0f',sp.Kt)); add_line(mdl,'GapY/1','Ky/1');
    add_block('simulink/Math Operations/Add',[mdl '/VrelY']); set_param([mdl '/VrelY'],'Inputs','+-');
    add_line(mdl,'Pvy/1','VrelY/1'); add_line(mdl,'Svy/1','VrelY/2');
    add_block('simulink/Math Operations/Gain',[mdl '/Dy']); set_param([mdl '/Dy'],'Gain',sprintf('%.0f',sp.Dt)); add_line(mdl,'VrelY/1','Dy/1');
    add_block('simulink/Math Operations/Add',[mdl '/SumFtDes']);
    add_line(mdl,'Ky/1','SumFtDes/1'); add_line(mdl,'Dy/1','SumFtDes/2');
    add_block('simulink/Math Operations/Gain',[mdl '/Mu']); set_param([mdl '/Mu'],'Gain',sprintf('%.4f',mu));   % mu baked in at creation
    add_line(mdl,'ClampFn/1','Mu/1');
    add_block('simulink/Math Operations/Unary Minus',[mdl '/NegMu']); add_line(mdl,'Mu/1','NegMu/1');
    add_block('simulink/Discontinuities/Saturation Dynamic',[mdl '/SatFt']);
    add_line(mdl,'SumFtDes/1','SatFt/1'); add_line(mdl,'NegMu/1','SatFt/2'); add_line(mdl,'Mu/1','SatFt/3');

    % Outputs to EFT
    add_block('nesl_utility/Simulink-PS Converter',[mdl '/Fxout']); set_param([mdl '/Fxout'],'Unit','N');
    add_line(mdl,'ClampFn/1','Fxout/1'); add_line(mdl,'Fxout/RConn1','EFT/LConn1');
    add_block('nesl_utility/Simulink-PS Converter',[mdl '/Fyout']); set_param([mdl '/Fyout'],'Unit','N');
    add_line(mdl,'SatFt/1','Fyout/1'); add_line(mdl,'Fyout/RConn1','EFT/LConn2');

    % Data capture
    add_block('simulink/Sources/Clock',[mdl '/Clk']);
    add_block('simulink/Sinks/To Workspace',[mdl '/TWT']); set_param([mdl '/TWT'],'VariableName','Tvec','SaveFormat','Array'); add_line(mdl,'Clk/1','TWT/1');
    add_block('simulink/Sinks/To Workspace',[mdl '/TWSv']); set_param([mdl '/TWSv'],'VariableName','Svx','SaveFormat','Array'); add_line(mdl,'Svx/1','TWSv/1');
    add_block('simulink/Sinks/To Workspace',[mdl '/TWSvy']); set_param([mdl '/TWSvy'],'VariableName','Svy','SaveFormat','Array'); add_line(mdl,'Svy/1','TWSvy/1');
    add_block('simulink/Sinks/To Workspace',[mdl '/TWFn']); set_param([mdl '/TWFn'],'VariableName','Fn','SaveFormat','Array'); add_line(mdl,'ClampFn/1','TWFn/1');

    set_param(mdl,'StopTime',sprintf('%.1f',sp.sim_time),'Solver','ode45','RelTol','1e-3','MaxStep','0.005','ZeroCross','off');
    fprintf('  Built fresh model for mu=%.4f\n', mu);
end
