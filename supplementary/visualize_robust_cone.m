%% visualize_robust_cone.m
% Refined animated visualization of the exact robust motion cone.
% Features:
%   • Left panel:  contact geometry with pusher, friction cone,
%                  push line, uncertainty box, and CoM velocity rays
%   • Right panel: polar motion cone plot with numerical angle
%                  annotations, stick/slide region labels, cone
%                  aperture arcs, and all candidate cones
%   • Bottom bar:  parameter sweep progress, live values,
%                  cone aperture comparison
%
% Parameter sweeps cycle:  δ (pusher offset) → α (normal tilt)
%                          → μ (friction) → c (push distance)
%
% Requires: MATLAB R2024b, core MATLAB only (no toolboxes needed).
% Output:   ../figures/robust_cone_visualization.mp4
% Run:      matlab -batch "run('supplementary/visualize_robust_cone.m')"

clear; close all; clc;
fprintf('=== Robust Motion Cone Visualization: Refined ===\n\n');

%% 1. Helper functions
    function v = perp(a); v = [-a(2); a(1)]; end
    function R = rot(theta)
        ct = cos(theta); st = sin(theta);
        R = [ct, -st; st, ct];
    end
    function a = angle_of(v); a = atan2(v(2), v(1)); end
    function a_w = wrap_to(a, ref)
        a_w = ref + mod(a - ref + pi, 2*pi) - pi;
    end
    function [fm, fp] = friction_edges(n, mu)
        g = atan(mu); fm = rot(-g)*n; fp = rot(+g)*n;
    end
    function M = push_map(p, c)
        s = 1.0/c^2; q = perp(p); M = eye(2) + s*(q*q');
    end
    function [um, up, ui] = cone_rays(p0, n0, mu, c, alpha, delta)
        t0 = perp(n0)/norm(perp(n0));
        p = p0 + delta*t0; n = rot(alpha)*n0;
        [fm, fp] = friction_edges(n, mu);
        M = push_map(p, c);
        um = M*fm; up = M*fp; ui = M*n;
    end
    function [lo, hi] = cone_angles(p0, n0, mu, c, alpha, delta, ref)
        [um, up, ~] = cone_rays(p0, n0, mu, c, alpha, delta);
        a1 = wrap_to(angle_of(um), ref); a2 = wrap_to(angle_of(up), ref);
        lo = min(a1, a2); hi = max(a1, a2);
    end

%% 2. Scenario parameters
p0 = [-0.30; 0.0]; n0 = [1.0; 0.0];
mu_lo = 0.25; mu_hi = 0.65;
c_lo = 0.045; c_hi = 0.075;
alpha_lo = deg2rad(-6); alpha_hi = deg2rad(6);
delta_lo = -0.015; delta_hi = 0.015;
mu0 = 0.45; c0 = 0.06; alpha0 = 0.0; delta0 = 0.0;

% Reference direction
[~, ~, ui_ref] = cone_rays(p0, n0, mu0, c0, alpha0, delta0);
ref = angle_of(ui_ref);

% Nominal cone
[lo_nom, hi_nom] = cone_angles(p0, n0, mu0, c0, 0, 0, ref);
aper_nom = rad2deg(hi_nom - lo_nom);

%% 3. Precompute robust cone and all 16 extreme candidates

% Robust cone (vertex grid + delta continuum)
lo_rob = -inf; hi_rob = inf;
for c_v = [c_lo, c_hi]
    for mu_v = [mu_lo, mu_hi]
        for a_v = [alpha_lo, alpha_hi]
            for d_v = linspace(delta_lo, delta_hi, 25)
                [l, h] = cone_angles(p0, n0, mu_v, c_v, a_v, d_v, ref);
                lo_rob = max(lo_rob, l); hi_rob = min(hi_rob, h);
            end
        end
    end
end
aper_rob = rad2deg(hi_rob - lo_rob);

% 16 box-vertex cones (illustration: the robust cone is contained in their
%   intersection; by Theorem 1 the mu/c/alpha optima are at these vertices,
%   while delta needs the continuum swept above). These are the box vertices,
%   NOT the paper's <=16 Algorithm-1 candidates (which use mu_lo only, both c,
%   per-edge alpha, and delta endpoints+roots).
extreme_deltas = [delta_lo, delta_hi];
extreme_alphas = [alpha_lo, alpha_hi];
extreme_mus    = [mu_lo, mu_hi];
extreme_cs     = [c_lo, c_hi];

candidate_lo = zeros(16,1); candidate_hi = zeros(16,1);
candidate_params = zeros(16,4);
idx = 1;
for d_e = extreme_deltas
    for a_e = extreme_alphas
        for m_e = extreme_mus
            for c_e = extreme_cs
                [l, h] = cone_angles(p0, n0, m_e, c_e, a_e, d_e, ref);
                candidate_lo(idx) = l; candidate_hi(idx) = h;
                candidate_params(idx,:) = [d_e, a_e, m_e, c_e];
                idx = idx + 1;
            end
        end
    end
end

fprintf('Scenario:\n');
fprintf('  Nominal cone: [%.2f, %.2f] deg  (aperture %.1f°)\n', ...
    rad2deg(lo_nom), rad2deg(hi_nom), aper_nom);
fprintf('  Robust  cone: [%.2f, %.2f] deg  (aperture %.1f°)\n', ...
    rad2deg(lo_rob), rad2deg(hi_rob), aper_rob);
fprintf('  16 candidates: lo ∈ [%.1f, %.1f]°, hi ∈ [%.1f, %.1f]°\n', ...
    rad2deg(min(candidate_lo)), rad2deg(max(candidate_lo)), ...
    rad2deg(min(candidate_hi)), rad2deg(max(candidate_hi)));
fprintf('  Reference:     0.00°\n\n');

%% 4. Video writer setup
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
video_file = fullfile(out_dir, 'robust_cone_visualization.mp4');
v = VideoWriter(video_file, 'MPEG-4');
v.FrameRate = 20; v.Quality = 90;
open(v);
fprintf('Writing video to %s ...\n', video_file);

%% 5. Figure layout
fig = figure('Position', [30, 30, 1480, 820], 'Color', [0.97,0.97,0.97]);

% Colour palette
nom_color    = [0.25, 0.45, 0.78];
rob_color    = [0.78, 0.08, 0.16];
rob_fill_c   = [0.78, 0.08, 0.16];
delta_color  = [0.56, 0.27, 0.68];  % purple
alpha_color  = [0.80, 0.25, 0.18];  % red
mu_color     = [0.12, 0.52, 0.26];  % green
c_color      = [0.90, 0.55, 0.10];  % orange
mode1_color  = [0.30, 0.60, 0.85];  % stick
mode2_color  = [0.85, 0.40, 0.20];  % slide
cand_color   = [0.50, 0.50, 0.50];  % grey for candidates

% Top title
annotation('textbox', [0.05, 0.94, 0.90, 0.05], ...
    'String', 'Exact Robust Motion Cone: Parameter Sweep', ...
    'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none', 'VerticalAlignment', 'middle');

% ========== LEFT PANEL: Contact Geometry ==========
ax1 = axes('Position', [0.04, 0.20, 0.42, 0.72]);
hold(ax1, 'on'); axis(ax1, 'equal');
xlim(ax1, [-0.48, 0.22]); ylim(ax1, [-0.22, 0.22]);
grid(ax1, 'on'); box(ax1, 'on');
xlabel(ax1, 'x (m)', 'FontSize', 10); ylabel(ax1, 'y (m)', 'FontSize', 10);
title(ax1, 'Contact Geometry', 'FontWeight', 'bold', 'FontSize', 12);

% Board halves (colored to show edge)
h_board_top  = patch(ax1, [-0.30, 0.30, 0.30, -0.30], [0, 0, 0.05, 0.05], ...
                     [0.85, 0.80, 0.65], 'EdgeColor', 'none');
h_board_bot  = patch(ax1, [-0.30, 0.30, 0.30, -0.30], [-0.05, -0.05, 0, 0], ...
                     [0.80, 0.75, 0.58], 'EdgeColor', 'none');
h_board_edge = plot(ax1, [-0.30, 0.30], [0, 0], 'Color', [0.35, 0.30, 0.15], ...
                    'LineWidth', 2.5);
text(ax1, 0.30, 0.055, 'board', 'FontSize', 8, 'HorizontalAlignment', 'right', ...
     'Color', [0.50, 0.45, 0.30]);

% CoM
h_com = plot(ax1, 0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(ax1, 0.01, 0.010, 'CoM', 'FontSize', 10, 'FontWeight', 'bold');

% Push line (dashed from contact to CoM)
h_push_line = plot(ax1, [0,0], [0,0], '--', 'Color', [0.30,0.30,0.30], 'LineWidth', 1.0);

% Pusher
h_pusher = plot(ax1, 0, 0, 'o', 'MarkerSize', 10, ...
                'MarkerFaceColor', [0.72,0,0.18], 'Color', [0.72,0,0.18]);
h_p_txt  = text(ax1, 0, 0, '', 'FontSize', 11, 'FontWeight', 'bold', ...
                'Color', [0.72,0,0.18]);

% Uncertainty box (delta+alpha range)
h_uncert = patch(ax1, [0,0,0,0], [0,0,0,0], [0.72,0,0.18], ...
                 'FaceAlpha', 0.06, 'EdgeColor', [0.72,0,0.18], ...
                 'LineWidth', 0.6, 'LineStyle', '--');

% Contact normal with quiver
h_normal_q = quiver(ax1, 0, 0, 0, 0, ...
                    'Color', [0.72,0,0.18], 'LineWidth', 1.8, 'MaxHeadSize', 0.45);
h_n_txt = text(ax1, 0, 0, '', 'Color', [0.72,0,0.18], ...
               'FontSize', 10, 'FontWeight', 'bold');

% Friction cone (two quivers + arc)
h_fric_p = quiver(ax1, 0, 0, 0, 0, ...
                  'Color', [0.16,0.44,0.75], 'LineWidth', 1.5, 'MaxHeadSize', 0.35);
h_fric_m = quiver(ax1, 0, 0, 0, 0, ...
                  'Color', [0.16,0.44,0.75], 'LineWidth', 1.5, 'MaxHeadSize', 0.35);
h_fric_arc = plot(ax1, 0, 0, 'Color', [0.16,0.44,0.75], 'LineWidth', 1.5);
h_fp_txt = text(ax1, 0, 0, '', 'Color', [0.16,0.44,0.75], 'FontSize', 9);
h_fm_txt = text(ax1, 0, 0, '', 'Color', [0.16,0.44,0.75], 'FontSize', 9);
h_gamma_txt = text(ax1, 0, 0, '', 'FontSize', 8, 'Color', [0.16,0.44,0.75]);

% CoM velocity rays (u_-, u_i, u_+)
h_ray_um = plot(ax1, [0,0], [0,0], 'Color', mode2_color, 'LineWidth', 2.0);
h_ray_ui = plot(ax1, [0,0], [0,0], 'Color', mode1_color, 'LineWidth', 2.0);
h_ray_up = plot(ax1, [0,0], [0,0], 'Color', mode2_color, 'LineWidth', 2.0);
% Arrow markers at tips
h_arrow_um = plot(ax1, 0, 0, '>', 'Color', mode2_color, ...
                  'MarkerSize', 7, 'MarkerFaceColor', mode2_color);
h_arrow_ui = plot(ax1, 0, 0, '>', 'Color', mode1_color, ...
                  'MarkerSize', 7, 'MarkerFaceColor', mode1_color);
h_arrow_up = plot(ax1, 0, 0, '>', 'Color', mode2_color, ...
                  'MarkerSize', 7, 'MarkerFaceColor', mode2_color);

% Ray labels
h_ray_label_um = text(ax1, 0, 0, '', 'FontSize', 8, 'FontWeight', 'bold', ...
                      'Color', mode2_color);
h_ray_label_ui = text(ax1, 0, 0, '', 'FontSize', 8, 'FontWeight', 'bold', ...
                      'Color', mode1_color);
h_ray_label_up = text(ax1, 0, 0, '', 'FontSize', 8, 'FontWeight', 'bold', ...
                      'Color', mode2_color);

% Parameter value box
h_param_box = text(ax1, -0.46, 0.20, '', 'FontSize', 8.5, ...
    'BackgroundColor', [1,1,1,0.88], 'EdgeColor', [0.50,0.50,0.50], ...
    'VerticalAlignment', 'top');

% Parameter description
h_param_desc = text(ax1, 0.20, -0.20, '', 'FontSize', 9, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'FontWeight', 'bold');

% Nominal reference
h_nom_ref_txt = text(ax1, -0.46, -0.20, '', 'FontSize', 7.5, ...
    'VerticalAlignment', 'bottom');

% ========== RIGHT PANEL: Motion Cone ==========
ax2 = axes('Position', [0.54, 0.20, 0.44, 0.72]);
hold(ax2, 'on'); axis(ax2, 'equal');
xlim(ax2, [-1.5, 1.5]); ylim(ax2, [-1.5, 1.5]);
box(ax2, 'on'); grid(ax2, 'off');
set(ax2, 'XTick', [], 'YTick', []);
title(ax2, 'Motion Cone', 'FontWeight', 'bold', 'FontSize', 12);

% Polar grid (concentric circles + radial lines)
for r = [0.25, 0.50, 0.75, 1.00]
    th = linspace(0, 2*pi, 100);
    plot(ax2, r*cos(th), r*sin(th), 'Color', [0.82,0.82,0.82], 'LineWidth', 0.4);
end
for th_d = -165:15:180
    th_r = deg2rad(th_d);
    plot(ax2, [0, 1.35*cos(th_r)], [0, 1.35*sin(th_r)], ...
         'Color', [0.87,0.87,0.87], 'LineWidth', 0.25);
end
text(ax2, 0, 1.10, '0°', 'HorizontalAlignment', 'center', ...
     'FontSize', 7, 'Color', [0.45,0.45,0.45]);
text(ax2, -1.10, 0, '±90°', 'HorizontalAlignment', 'center', ...
     'FontSize', 7, 'Color', [0.45,0.45,0.45]);
text(ax2, 1.12, 0, '0°', 'HorizontalAlignment', 'center', ...
     'FontSize', 7, 'Color', [0.45,0.45,0.45]);

% 16 candidate edge rays (thin grey lines, pre-created)
h_cand = gobjects(16,2);
for k = 1:16
    h_cand(k,1) = plot(ax2, [0,0], [0,0], 'Color', cand_color, ...
                       'LineWidth', 0.5, 'Visible', 'off');
    h_cand(k,2) = plot(ax2, [0,0], [0,0], 'Color', cand_color, ...
                       'LineWidth', 0.5, 'Visible', 'off');
end

% Nominal cone (fill + dashed edge)
h_nom_fill = fill(ax2, [0,0,0], [0,0,0], nom_color, ...
                  'FaceAlpha', 0.04, 'EdgeColor', 'none');
h_nom_edge = plot(ax2, [0,0], [0,0], 'Color', nom_color, ...
                  'LineWidth', 2.0, 'LineStyle', '--');

% Robust cone (fill + thick edge)
h_rob_fill = fill(ax2, [0,0,0], [0,0,0], rob_fill_c, ...
                  'FaceAlpha', 0.18, 'EdgeColor', 'none');
h_rob_edge = plot(ax2, [0,0], [0,0], 'Color', rob_color, ...
                  'LineWidth', 3.2);

% Current cone (fill + edge, color changes per parameter)
h_cur_fill = fill(ax2, [0,0,0], [0,0,0], [0,0,0], ...
                  'FaceAlpha', 0.12, 'EdgeColor', 'none');
h_cur_edge = plot(ax2, [0,0], [0,0], 'Color', [0,0,0], 'LineWidth', 2.2);

% Reference direction
h_ref_line = plot(ax2, [0, 0.30], [0, 0], 'k-', 'LineWidth', 1.0);

% Interior ray (separates two modes)
h_ui_polar = plot(ax2, [0,0], [0,0], ':', 'Color', [0.40,0.40,0.40], 'LineWidth', 1.0);

% Aperture arc for robust cone (arc with angle label)
h_aper_arc_rob = plot(ax2, 0, 0, 'Color', rob_color, 'LineWidth', 1.5);
h_aper_txt_rob = text(ax2, 0, 0, '', 'FontSize', 9, 'FontWeight', 'bold', ...
                      'Color', rob_color, 'HorizontalAlignment', 'center');

% Region labels
h_stick_label = text(ax2, 0, 0, '', 'FontSize', 9, 'FontWeight', 'bold', ...
                     'Color', mode1_color, 'HorizontalAlignment', 'center');
h_slide_label = text(ax2, 0, 0, '', 'FontSize', 9, 'FontWeight', 'bold', ...
                     'Color', mode2_color, 'HorizontalAlignment', 'center');

% Angle edge labels
h_lo_label = text(ax2, 0, 0, '', 'FontSize', 8, 'Color', rob_color, ...
                  'HorizontalAlignment', 'right');
h_hi_label = text(ax2, 0, 0, '', 'FontSize', 8, 'Color', rob_color, ...
                  'HorizontalAlignment', 'left');

% Legend box
h_legend_box = text(ax2, 1.30, 1.20, '', 'FontSize', 7.5, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1,1,1,0.85], 'EdgeColor', [0.50,0.50,0.50]);

% Cone description
h_cone_desc = text(ax2, 0, -1.42, '', 'FontSize', 9, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% ========== BOTTOM BAR ==========
% Sweep progress bar background
ax_bar = axes('Position', [0.04, 0.07, 0.52, 0.06]);
hold(ax_bar, 'on'); axis(ax_bar, 'off');
xlim(ax_bar, [0, 1]); ylim(ax_bar, [0, 1]);
% Progress bar fill
h_prog_fill = patch(ax_bar, [0,0,0,0], [0,1,1,0], [0.30,0.50,0.80], ...
                    'EdgeColor', 'none');
% Progress bar border
plot(ax_bar, [0,1,1,0,0], [0,0,1,1,0], 'k-', 'LineWidth', 1.0);
% Segment boundary markers
for s = 0:4
    xm = s / 4;
    plot(ax_bar, xm, 0, 'k|', 'MarkerSize', 10, 'LineWidth', 0.8);
end
% Segment labels
seg_labels = {'\delta', '\alpha', '\mu', 'c'};
for s = 1:4
    xm = (s-0.5) / 4;
    text(ax_bar, xm, -0.45, seg_labels{s}, 'FontSize', 8, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end
% Parameter labels under the bar
h_sweep_label_bar = text(ax_bar, 0.50, -0.60, '', ...
    'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Right side status panel
ax_status = axes('Position', [0.58, 0.07, 0.40, 0.06]);
hold(ax_status, 'on'); axis(ax_status, 'off');
xlim(ax_status, [0, 1]); ylim(ax_status, [0, 1]);
h_status_txt = text(ax_status, 0.50, 0.50, '', ...
    'FontSize', 8.5, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

%% 6. Animation loop
n_frames = 300;
fps = v.FrameRate;
n_param_frames = 72;

for frame = 1:n_frames
    
    % ---- Parameter cycling ----
    if frame <= n_param_frames
        frac = (frame-1) / max(1, n_param_frames-1);
        delta_v = delta_lo + frac*(delta_hi - delta_lo);
        alpha_v = 0; mu_v = mu0; c_v = c0;
        param_label = '\delta  offset';
        active_param = 'delta';
        sweep_desc = 'Pusher offset \delta varies: contact slides along the edge';
        param_symbol = sprintf('\\delta = %+.1f mm', delta_v*1000);
    elseif frame <= 2*n_param_frames
        frac = (frame - n_param_frames - 1) / max(1, n_param_frames-1);
        delta_v = 0; alpha_v = alpha_lo + frac*(alpha_hi - alpha_lo);
        mu_v = mu0; c_v = c0;
        param_label = '\alpha  tilt';
        active_param = 'alpha';
        sweep_desc = 'Normal tilt \alpha varies: contact face rotates';
        param_symbol = sprintf('\\alpha = %+.1f°', rad2deg(alpha_v));
    elseif frame <= 3*n_param_frames
        frac = (frame - 2*n_param_frames - 1) / max(1, n_param_frames-1);
        delta_v = 0; alpha_v = 0;
        mu_v = mu_lo + frac*(mu_hi - mu_lo);
        c_v = c0;
        param_label = '\mu  friction';
        active_param = 'mu';
        sweep_desc = 'Friction \mu varies: friction cone opens / closes';
        param_symbol = sprintf('\\mu = %.2f  (\\gamma = %.1f°)', mu_v, rad2deg(atan(mu_v)));
    else
        frac = (frame - 3*n_param_frames - 1) / max(1, n_frames - 3*n_param_frames - 1);
        delta_v = 0; alpha_v = 0; mu_v = mu0;
        c_v = c_lo + frac*(c_hi - c_lo);
        param_label = 'c  distance';
        active_param = 'c';
        sweep_desc = 'Push distance c varies: pressure distribution width';
        param_symbol = sprintf('c = %.1f mm', c_v*1000);
    end
    
    % Compute current geometry
    [um, up, ui] = cone_rays(p0, n0, mu_v, c_v, alpha_v, delta_v);
    [lo_cur, hi_cur] = cone_angles(p0, n0, mu_v, c_v, alpha_v, delta_v, ref);
    t0_vec = perp(n0)/norm(perp(n0));
    p_cur = p0 + delta_v * t0_vec;
    n_cur = rot(alpha_v) * n0;
    gamma_v = atan(mu_v);
    [fm_cur, fp_cur] = friction_edges(n_cur, mu_v);
    
    aper_cur  = rad2deg(hi_cur - lo_cur);
    
    % ======================== LEFT PANEL UPDATE ========================
    
    % Board (static, already drawn)
    
    % Push line: from contact to CoM
    set(h_push_line, 'XData', [p_cur(1), 0], 'YData', [p_cur(2), 0]);
    
    % Pusher at current position
    set(h_pusher, 'XData', p_cur(1), 'YData', p_cur(2));
    % Annotate with delta value
    set(h_p_txt, 'Position', [p_cur(1)-0.025, p_cur(2)-0.015], ...
        'String', sprintf('p (δ=%.0fmm)', delta_v*1000));
    
    % Uncertainty box (delta range × alpha-induced lateral shift range)
    % Compute the four corner positions due to (delta_lo/delta_hi, alpha_lo/alpha_hi)
    t0 = perp(n0)/norm(perp(n0));
    % For simplicity, show a box around the delta range with a small height
    % representing alpha uncertainty effect on lateral position (approximate)
    d_span = delta_hi - delta_lo;
    lat_span = 0.5 * sin(alpha_hi - alpha_lo) * norm(p0); % rough lateral shift
    bx = p0(1) + [delta_lo*t0(1)-lat_span, delta_hi*t0(1)+lat_span, ...
                  delta_hi*t0(1)+lat_span, delta_lo*t0(1)-lat_span];
    by = p0(2) + [delta_lo*t0(2)-lat_span, delta_hi*t0(2)-lat_span, ...
                  delta_hi*t0(2)+lat_span, delta_lo*t0(2)+lat_span];
    set(h_uncert, 'XData', bx, 'YData', by);
    
    % Contact normal
    n_len = 0.09;
    set(h_normal_q, 'XData', p_cur(1), 'YData', p_cur(2), ...
                    'UData', n_len*n_cur(1), 'VData', n_len*n_cur(2));
    set(h_n_txt, 'Position', [p_cur(1)+0.11*n_cur(1), p_cur(2)+0.11*n_cur(2)], ...
        'String', 'n');
    
    % Friction cone edges
    f_len = 0.11;
    f_p = rot(+gamma_v)*n_cur; f_m = rot(-gamma_v)*n_cur;
    set(h_fric_p, 'XData', p_cur(1), 'YData', p_cur(2), ...
                  'UData', f_len*f_p(1), 'VData', f_len*f_p(2));
    set(h_fric_m, 'XData', p_cur(1), 'YData', p_cur(2), ...
                  'UData', f_len*f_m(1), 'VData', f_len*f_m(2));
    
    % Friction cone arc
    theta_f = linspace(-gamma_v, gamma_v, 30);
    arc_x = p_cur(1) + 0.09*cos(theta_f + angle_of(n_cur));
    arc_y = p_cur(2) + 0.09*sin(theta_f + angle_of(n_cur));
    set(h_fric_arc, 'XData', arc_x, 'YData', arc_y);
    
    % Edge labels
    set(h_fp_txt, 'Position', [p_cur(1)+1.15*f_len*f_p(1), p_cur(2)+1.15*f_len*f_p(2)], ...
        'String', 'f_+');
    set(h_fm_txt, 'Position', [p_cur(1)+1.15*f_len*f_m(1), p_cur(2)+1.15*f_len*f_m(2)], ...
        'String', 'f_-');
    
    % γ annotation (friction angle)
    gamma_mid = angle_of(n_cur) - gamma_v/2;
    set(h_gamma_txt, 'Position', [p_cur(1)+0.13*cos(gamma_mid), ...
                                   p_cur(2)+0.13*sin(gamma_mid)], ...
        'String', sprintf('\\gamma = %.1f°', rad2deg(gamma_v)));
    
    % CoM velocity rays (scaled for visibility)
    ray_scale = 0.10;  % 10 cm on the board
    um_unit = um / norm(um); up_unit = up / norm(up); ui_unit = ui / norm(ui);
    set(h_ray_um, 'XData', [0, ray_scale*um_unit(1)], 'YData', [0, ray_scale*um_unit(2)]);
    set(h_ray_ui, 'XData', [0, ray_scale*ui_unit(1)], 'YData', [0, ray_scale*ui_unit(2)]);
    set(h_ray_up, 'XData', [0, ray_scale*up_unit(1)], 'YData', [0, ray_scale*up_unit(2)]);
    
    % Arrow markers
    set(h_arrow_um, 'XData', ray_scale*um_unit(1), 'YData', ray_scale*um_unit(2));
    set(h_arrow_ui, 'XData', ray_scale*ui_unit(1), 'YData', ray_scale*ui_unit(2));
    set(h_arrow_up, 'XData', ray_scale*up_unit(1), 'YData', ray_scale*up_unit(2));
    
    % Ray labels
    set(h_ray_label_um, 'Position', 1.05*ray_scale*um_unit', 'String', 'u_-');
    set(h_ray_label_ui, 'Position', 1.05*ray_scale*ui_unit', 'String', 'u_i');
    set(h_ray_label_up, 'Position', 1.05*ray_scale*up_unit', 'String', 'u_+');
    
    % Parameter info box
    param_str = sprintf(['\\mu = %.2f\n\\gamma = %.1f°\n' ...
                         'c = %.1f mm\n\\alpha = %+.1f°\n' ...
                         '\\delta = %+.1f mm'], ...
        mu_v, rad2deg(gamma_v), c_v*1000, rad2deg(alpha_v), delta_v*1000);
    aper_str = sprintf(['\n----- Apertures -----\n' ...
                        'Nominal:  %.0f°\nRobust:   %.0f°\nCurrent:  %.0f°'], ...
        aper_nom, aper_rob, aper_cur);
    set(h_param_box, 'String', [param_str, aper_str]);
    
    % Parameter description
    set(h_param_desc, 'String', sweep_desc);
    
    % Nominal reference
    set(h_nom_ref_txt, 'String', ...
        sprintf('Nominal: μ=%.2f, c=%.0fmm, α=0°, δ=0mm', mu0, c0*1000));
    
    % ======================== RIGHT PANEL UPDATE ========================
    
    % 16 extreme candidate edge rays
    for k = 1:16
        set(h_cand(k,1), 'XData', [0, cos(candidate_lo(k))], ...
                         'YData', [0, sin(candidate_lo(k))], 'Visible', 'on');
        set(h_cand(k,2), 'XData', [0, cos(candidate_hi(k))], ...
                         'YData', [0, sin(candidate_hi(k))], 'Visible', 'on');
    end
    
    % Nominal cone
    th_fill = linspace(lo_nom, hi_nom, 80);
    set(h_nom_fill, 'XData', [0, cos(th_fill), 0], 'YData', [0, sin(th_fill), 0]);
    th_edge = linspace(lo_nom, hi_nom, 80);
    set(h_nom_edge, 'XData', cos(th_edge), 'YData', sin(th_edge));
    
    % Robust cone
    if lo_rob < hi_rob
        th_fill = linspace(lo_rob, hi_rob, 80);
        set(h_rob_fill, 'XData', [0, cos(th_fill), 0], 'YData', [0, sin(th_fill), 0]);
        th_edge = linspace(lo_rob, hi_rob, 80);
        set(h_rob_edge, 'XData', cos(th_edge), 'YData', sin(th_edge));
    end
    
    % Current cone with colour
    switch active_param
        case 'delta'; cur_color = delta_color;
        case 'alpha'; cur_color = alpha_color;
        case 'mu';    cur_color = mu_color;
        case 'c';     cur_color = c_color;
    end
    th_fill = linspace(lo_cur, hi_cur, 60);
    set(h_cur_fill, 'XData', [0, cos(th_fill), 0], 'YData', [0, sin(th_fill), 0], ...
                    'FaceColor', cur_color);
    th_edge = linspace(lo_cur, hi_cur, 60);
    set(h_cur_edge, 'XData', cos(th_edge), 'YData', sin(th_edge), 'Color', cur_color);
    
    % Reference direction
    set(h_ref_line, 'XData', [0, 0.30*cos(ref)], 'YData', [0, 0.30*sin(ref)]);
    
    % Interior ray (u_i direction, separates the two motion modes)
    ui_angle = angle_of(ui);  % should be close to ref
    set(h_ui_polar, 'XData', [0, 0.85*cos(ui_angle)], 'YData', [0, 0.85*sin(ui_angle)]);
    
    % Robust cone aperture arc (at r=0.65)
    if lo_rob < hi_rob
        arc_r = 0.65;
        th_arc = linspace(lo_rob, hi_rob, 40);
        set(h_aper_arc_rob, 'XData', arc_r*cos(th_arc), 'YData', arc_r*sin(th_arc));
        % Aperture label at the midpoint
        mid_rob = (lo_rob + hi_rob) / 2;
        set(h_aper_txt_rob, 'Position', arc_r*1.12*[cos(mid_rob), sin(mid_rob)], ...
            'String', sprintf('%.0f°', aper_rob));
    end
    
    % Stick / Slide region labels.
    % The ENTIRE cone [lo_cur, hi_cur] is the sticking set (any command whose
    % direction lies inside it sticks); commands outside the cone slide. The
    % interior ray u_i is NOT a stick/slide boundary, so STICK is placed at the
    % cone midpoint and SLIDE just outside the upper edge.
    if lo_cur < hi_cur
        mid_cone  = (lo_cur + hi_cur) / 2;
        slide_ang = hi_cur + deg2rad(18);
        label_r = 0.55;
        set(h_stick_label, 'Position', label_r*[cos(mid_cone), sin(mid_cone)], ...
            'String', 'STICK');
        set(h_slide_label, 'Position', label_r*[cos(slide_ang), sin(slide_ang)], ...
            'String', 'SLIDE');
    end
    
    % Angle labels at robust cone edges
    if lo_rob < hi_rob
        label_r = 1.08;
        set(h_lo_label, 'Position', label_r*[cos(lo_rob), sin(lo_rob)], ...
            'String', sprintf('%.1f°', rad2deg(lo_rob)));
        set(h_hi_label, 'Position', label_r*[cos(hi_rob), sin(hi_rob)], ...
            'String', sprintf('+%.1f°', rad2deg(hi_rob)));
    end
    
    % Legend
    leg_str = sprintf(['  ---- Nominal cone (%3.0f°)\n' ...
                       '  ---- Robust cone  (%3.0f°)\n' ...
                       '  ---- Current cone (%3.0f°)\n' ...
                       '  . . .  Interior ray (u_i)\n' ...
                       '  · · ·  16 box-vertex cones'], ...
        aper_nom, aper_rob, aper_cur);
    set(h_legend_box, 'String', leg_str);
    
    % Cone description
    set(h_cone_desc, 'String', ...
        sprintf('Robust cone = intersection of all candidate cones (%.0f° aperture)', aper_rob));
    
    % ======================== BOTTOM BAR ========================
    seg_start = floor((frame-1) / n_param_frames);
    seg_frac = ((frame-1) - seg_start*n_param_frames) / max(1, n_param_frames-1);
    overall_frac = (frame-1) / (n_frames-1);
    
    set(h_prog_fill, 'XData', [0, overall_frac, overall_frac, 0], ...
                     'YData', [0,0,1,1]);
    
    % Progress bar color matches active sweep
    switch active_param
        case 'delta'; bar_color = delta_color;
        case 'alpha'; bar_color = alpha_color;
        case 'mu';    bar_color = mu_color;
        case 'c';     bar_color = c_color;
    end
    set(h_prog_fill, 'FaceColor', bar_color);
    
    set(h_sweep_label_bar, 'String', ...
        sprintf('Sweeping:  %s', param_label));
    
    % Status panel
    set(h_status_txt, 'String', ...
        sprintf(['Parameters  |  \\mu=%.2f  c=%.0fmm  \\alpha=%+.1f°  \\delta=%+.1fmm  |  ' ...
                 '|  Nom: %.0f°  Rob: %.0f°  Cur: %.0f°  |  %s'], ...
        mu_v, c_v*1000, rad2deg(alpha_v), delta_v*1000, ...
        aper_nom, aper_rob, aper_cur, param_symbol));
    
    % ======================== CAPTURE ========================
    drawnow;
    frame_img = getframe(fig);
    writeVideo(v, frame_img);
    
    if mod(frame, 50) == 0
        fprintf('  Frame %d/%d (%.0f%%)\n', frame, n_frames, 100*frame/n_frames);
    end
end

% Hold final frame for 1 second
for rep = 1:fps
    writeVideo(v, frame_img);
end

close(v);
fprintf('\nVideo written: %s\n', video_file);
fprintf('  %d frames at %d fps = %.1f seconds\n', n_frames+fps, fps, (n_frames+fps)/fps);
fprintf('Done.\n');
