%% MATLAB Version R2017b
clear all; clc; close all;

Time= 6;               % simulation time, sec
T   = 10e-3;            % step, sec
N   = floor(Time/T);    % steps
t   = 0;

% define the initial position
p1=[[-1.0; 1.0], ...        % position of agent-1
    [ 0.5; 0.5], ...
    [-0.5;-0.5]];           % position of obstacle-1
vp = [2, -2];
p2 = p1;

% agent and obstacle setting
na  = 1;              % number of mobile agents
no  = size(p1, 2)-1;  % number of obstacles
nao = na+no;
x1  = zeros(8, na);   % actuation state of agent1
x2  = zeros(8, na);   % actuation state of agent2
vp1 = zeros(2, na);   % ideal velocity
vp2 = zeros(2, na);   % ideal velocity
v_set1 = zeros(2, na);   % the input of the actuation
v_set2 = zeros(2, na);   % the input of the actuation
delta1 = zeros(1, na);   % the slack variable
delta2 = zeros(1, na);   % the slack variable

% quadraitc programming setting
Ap_i= getPolyAb(11, 1, 0);      %
ap_i= zeros(size(Ap_i, 1), na);
A1   = zeros(nao-1, 2, na);
b1   = zeros(nao-1, 1, na);
A2   = zeros(nao-1, 2, na);
b2   = zeros(nao-1, 1, na);

Ds  = 0.65;
% Ds  = sqrt(2)/2-10e-3;
Do  = 0.6;
delta_bar = 1e2;
alpha_v = 0.75;
vpd_bar = 5;
vp_bar  = 5;

% initialize graphic
H = myGraphic([], p1, x1, p2, x2, Ds, Ap_i, ap_i);

% simulation process
for k = 1:N
    % sample time
    t = k*T;
    
    % ideal velocity
    vp1(:, 1)= vp;
    vp2(:, 1)= vp;
    
    % collision avoidance
    [v_set1(:,1), delta1(1), ap_i(:,1), A1(:,:,1), b1(:,:,1)] = collision_avoidance_reshape(p1, vp1, Ds, Ap_i, delta_bar, alpha_v);
    [v_set2(:,1), delta2(1),            A2(:,:,1), b2(:,:,1)] = collision_avoidance(p2, vp2, Ds, delta_bar, alpha_v);
        
    % update model
    x1(:, 1) = rungekutta(@m_speed,    x1(:, 1), v_set1(  :, 1), T);     % ????-1
    x2(:, 1) = rungekutta(@m_speed,    x2(:, 1), v_set2(  :, 1), T);     % ????-2
    
    % Model 210918, res_06.mat
    C=[ 0.7761   -1.9815         0         0    2.1315    2.4144         0         0
             0         0   -2.2000   -2.8153         0         0   -1.5060   -0.9899];
    p1(:, 1) = rungekutta(@m_position, p1(:, 1),     C*x1, T);     % ????-1
    p2(:, 1) = rungekutta(@m_position, p2(:, 1),     C*x2, T);     % ????-2
    
    
    % update the graphic
    H = myGraphic(H, p1, x1, p2, x2, Ds, Ap_i, ap_i);
end
%% end of main.

%%
function [v_set, delta, A, b] = collision_avoidance(p, vp, Ds, delta_bar, alpha_v)
persistent opt;
if(isempty(opt))
%     opt = optimoptions('quadprog',  'Algorithm','interior-point-convex','Display','off', 'MaxIterations', 1e3);
    opt = optimoptions('quadprog',  'Display','off', 'MaxIterations', 1e3);
end

pm    = p(:, 1);
v_aim = vp(:, 1);
nao   = size(p, 2);

% generate conditions
A = zeros(nao-1, 2);
b = zeros(nao-1, 1);
count = 0;
for k = 2:nao
    count = count + 1;
    p_delta = pm - p(:, k);     % relative position
    dist    = norm(p_delta, 2); % distance
    h = dist - Ds;              % barrier function
    V = 1/(h+Ds);               % Lyapunov function
    % condition - k
    A(count, :) = -p_delta'/norm(p_delta);
    b(count, :) = -alpha_v*(V-1/Ds)/delta_bar;
%     b(count, :) = -alpha_v*(V-1/Ds);
%     b(count, :) = alpha_v*(dist^2-Ds^2);
end

% quadratic programming
f_aim = [v_aim; delta_bar];

% Standard Quadratic Programming with Slack Variables
v_set_det = quadprog(...
    diag([1, 1, 1]), -f_aim,...                             % cost function
    [A, -b], zeros(size(A, 1), 1), [], [],  ...                   % constraints
    [-inf; -inf; 0], [inf; inf; delta_bar], zeros(3, 1), opt);    % limits and setting

v_set  = v_set_det(1:end-1);
delta  = v_set_det(end);
end

function [v_set, delta, ap_i, A, b] = collision_avoidance_reshape(p, vp, Ds, Ap_i, delta_bar, alpha_v)
persistent opt;
if(isempty(opt))
%     opt = optimoptions('quadprog',  'Algorithm','interior-point-convex','Display','off', 'MaxIterations', 1e3);
    opt = optimoptions('quadprog',  'Display','off', 'MaxIterations', 1e3);
end

pm    =  p(:, 1);
v_aim = vp(:, 1);
nao   = size(p, 2);

% generate conditions
A = zeros(nao-1, 2);
b = zeros(nao-1, 1);
count = 0;

for k = 2:nao
    count = count + 1;
    p_delta = pm - p(:, k);     % relative position
    dist    = norm(p_delta, 2); % distance
    h = dist - Ds;              % barrier function
    V = 1/(h+Ds);               % Lyapunov function
    
    % condition - count
    A(count, :) = -p_delta'/norm(p_delta);
    b(count, :) = -alpha_v*(V-1/Ds)/delta_bar;
%     b(count, :) = alpha_v*(dist^2-Ds^2);
end

% quadratic programming
f_aim = [v_aim; delta_bar];
tau   = 2*sqrt(2)/delta_bar;
ap_i  = gen_a(Ap_i, A, -b, tau);  % Ap_i x + a delta \le 0
M     = [Ap_i, ap_i];

% Quadratic Programming with Slack Variables - Lipschitz
[v_set_det, ~, ~, ~, ~] = quadprog(...
    diag([1, 1, 1]), -f_aim,...                             % cost function
    M, zeros(size(M, 1), 1), [], [],  ...                   % constraints
    [-inf; -inf; 0], [inf; inf; delta_bar], zeros(3, 1), opt);    % limits and setting

v_set  = v_set_det(1:end-1);
delta  = v_set_det(end);
if(delta <=  1e-6)
    1;
end
end

% function H = myGraphic(H, p1, x1, Ds, Ap_i, ap_i, delta)
function H = myGraphic(H, p1, x1, p2, x2, Ds, Ap_i, ap_i)
persistent crashflag_2;
if(isempty(crashflag_2))
    crashflag_2 = false;
end

nao = size(p1, 2);
na  = 1;
no  = nao - na;
Do  = evalin('base', 'Do');
drawnow limitrate;
if(isempty(H))
    % initialize graphic
    figure(1); clf;
    axis equal;
    hold on;
    grid minor;
    set(gcf, 'position', [0,0, 800, 700 ], 'color', 'w');
    
    % obstacles
    for iO = (na+1):(na+no)
        H.shapeO(iO) = patch(p1(1,iO)+Ds*sin(linspace(0, 2*pi, 30)), p1(2,iO)+Ds*cos(linspace(0, 2*pi, 30)), 'c', 'facealpha', 0.1, 'HandleVisibility', 'off', 'linestyle', '--', 'EdgeColor', 'k');
        H.realO(iO)  = patch(p1(1,iO)+Do*sin(linspace(0, 2*pi, 30)), p1(2,iO)+Do*cos(linspace(0, 2*pi, 30)), 'c', 'facealpha', 1.0, 'HandleVisibility', 'off');
    end
    H.posO = plot(p1(1, na+1:end), p1(2, na+1:end), 'k.', 'markersize', 5, 'linestyle', 'none', 'HandleVisibility', 'off');
    
    % mobile agents
    HandleVisibility  = 'on';
    H.trajA(2) = animatedline('color', 'k', 'HandleVisibility', HandleVisibility, 'linewidth', 2, 'linestyle', ':', 'marker', 'none', 'Displayname', 'RP');
    H.trajA(3) = animatedline('color', 'k', 'HandleVisibility', HandleVisibility, 'linewidth', 2, 'linestyle', '--', 'marker', 'none', 'Displayname', 'RP after crash');
    H.trajA(1) = animatedline('color', 'b', 'HandleVisibility', HandleVisibility, 'linewidth', 2, 'linestyle', '-', 'marker', 'none', 'Displayname', 'RPRF');
    H.A = line([p1(1, 1), p2(1, 1)], [p1(2, 1), p2(2, 1)], 'Marker', 'o', 'linestyle', 'none', 'color','k', 'HandleVisibility', 'off');
    legend;
    text(-0.85, 1.2, 'Initial Position', 'fontsize', 18, 'HorizontalAlignment', 'center');
    text(-1.0, 1.02, '{\bf\times}', 'fontsize', 24, 'HorizontalAlignment', 'center');
    axis(1.3*[-1, 1, -1, 1]);
%     xlabel('$[p_1]_1$ (m)','Interpreter','latex'); ylabel('$[p_1]_2$ (m)','Interpreter','latex');
    xlabel('[p_1]_1 (m)'); ylabel('[p_1]_2 (m)');
%     set(gca, 'fontsize', 16, 'position', [0.06,0.13,0.85,0.85]);
    set(gca, 'fontsize', 19, 'box', 'on');
    fillwindows(gca);
    
    % feasible region
    figure(2); clf;
    axis equal;
    hold on;
    grid on;
    set(gcf, 'position', [0, 530, 800, 500 ], 'color', 'w');
    set(gca, 'fontsize', 18);
    H.feasibleRegion_1 = half_space_2D([], Ap_i, ap_i, 2*[-1, 1, -1, 1], 'b', 'domain', 'reshaped');
    H.feasibleRegion_2 = half_space_2D([], zeros(nao-1, 2), zeros(nao-1, 1), 2*[-1, 1, -1, 1], 'k', 'domain', 'original');
    
    %     H.v_p1_hdl   = quiver(0, 0, 0.0, 0.0, 1.0, 'linewidth', 2, 'color', 'b', 'MaxHeadSize',2, 'DisplayName', '$v^p_1$');
    H.v_set1_hdl = quiver(0, 0, 0.0, 0.0, 1.0, 'linewidth', 2, 'color', 'b', 'MaxHeadSize',2, 'linestyle', '-', 'DisplayName', '$v^*_1$');
    H.v_set2_hdl = quiver(0, 0, 0.0, 0.0, 1.0, 'linewidth', 2, 'color', 'k', 'MaxHeadSize',2, 'linestyle', '-', 'DisplayName', '$v^*_2$');
    lgd2_handle = legend('location', 'northwest');
    set(lgd2_handle,'Interpreter','latex', 'FontSize', 18);
    
    % trajectory
    figure(3); clf;
    set(gcf, 'position', [800, 530, 1000, 700 ], 'color', 'w');
    subplot(2, 1, 1);
    hold on;
    grid on;
    set(gca, 'fontsize', 16);
    H.Dist1= animatedline('color', 'b', 'HandleVisibility', HandleVisibility, 'linewidth', 2, 'linestyle', '-', 'marker', 'none', 'DisplayName', 'RPRF');
    H.Dist2= animatedline('color', 'k', 'HandleVisibility', HandleVisibility, 'linewidth', 2, 'linestyle', ':', 'marker', 'none', 'DisplayName', 'RP');
    H.D_safe= plot([0 evalin('base', 'Time')], Do*[1, 1], 'k-', 'linestyle', '--', 'HandleVisibility', 'off');
    lgd3_handle = legend('location', 'southwest','Orientation','horizontal');
    set(lgd3_handle,'Interpreter','latex', 'FontSize', 18)
    xlim([0 evalin('base', 'Time')]);
    xlabel('time (sec)');
    ylim([0 2]);
    ylabel({'minimal distance', '$\min\{|\tilde{p}_{12}|, |\tilde{p}_{12}|\}$  (m)'}, 'Interpreter' , 'latex');
    
    % trajectory
    subplot(2, 1, 2);
    hold on;
    grid on;
    set(gca, 'fontsize', 16);
    H.v_set1 = animatedline('color', 'b', 'HandleVisibility', HandleVisibility, 'linewidth', 2, 'linestyle', '-', 'marker', 'none', 'DisplayName', 'RPRF');
    H.v_set2 = animatedline('color', 'k', 'HandleVisibility', HandleVisibility, 'linewidth', 2, 'linestyle', ':', 'marker', 'none', 'DisplayName', 'RP');
%     H.D_safe= plot([0 evalin('base', 'Time')], Do*[1, 1], 'k-', 'linestyle', '--', 'HandleVisibility', 'off');
    lgd3_handle = legend('location', 'northwest','Orientation','horizontal');
    set(lgd3_handle,'Interpreter','latex', 'FontSize', 18)
    xlim([0 evalin('base', 'Time')]);
    xlabel('time (sec)');
    ylim([0 3]);
    ylabel({'velocity reference', '$|v^*_1|$ (m/s)'}, 'Interpreter' , 'latex');
    
    figure(4); clf;
    set(gcf, 'position', [800, 530, 800, 400 ], 'color', 'w');
    hold on;
    grid on;
    set(gca, 'fontsize', 18, 'box', 'on');
%     set(gca, 'fontsize', 16, 'position', [0.15,0.22,0.80,0.78]);
    H.v_set2 = animatedline('color', 'k', 'HandleVisibility', HandleVisibility, 'linewidth', 2, 'linestyle', ':', 'marker', 'none', 'DisplayName', 'RP');
    H.v_set2_crash = animatedline('color', 'k', 'HandleVisibility', HandleVisibility, 'linewidth', 2, 'linestyle', '--', 'marker', 'none', 'Displayname', 'RP after crash');
    H.v_set1 = animatedline('color', 'b', 'HandleVisibility', HandleVisibility, 'linewidth', 2, 'linestyle', '-', 'marker', 'none', 'DisplayName', 'RPRF');
    
%     H.D_safe= plot([0 evalin('base', 'Time')], Do*[1, 1], 'k-', 'linestyle', '--', 'HandleVisibility', 'off');
%     lgd3_handle = legend('location', 'northwest','Orientation','horizontal');
    lgd3_handle = legend('location', 'northwest');
    set(lgd3_handle, 'FontSize', 18)
    xlim([0 evalin('base', 'Time')]);
    xlabel('time (sec)');
    ylim([0 3.0]);
%     ylabel({'velocity command', '$|v^*_1|$ (m/s)'}, 'Interpreter' , 'latex');
    ylabel({'velocity reference', '|v^*_1| (m/s)'});
    fillwindows(gca);
else
    % get variables from workspace
    t     = evalin('base', 't');
    T     = evalin('base', 'T');
    
    delta1= evalin('base', 'delta1');   % reashaped algorighm
    v1_real= x1(3:4, :);
    A1    = evalin('base', 'A1');
    b1    = evalin('base', 'b1');
    v_set1 = evalin('base', 'v_set1');
    vp1   = evalin('base', 'vp1');
    
    delta2= evalin('base', 'delta2');   % original
    v2_real= x2(3:4, :);
    A2    = evalin('base', 'A2');
    b2    = evalin('base', 'b2');
    v_set2 = evalin('base', 'v_set2');
    vp2   = evalin('base', 'vp2');
    
    % update the graphic
    % mobile agents
    H.A.XData = [p1(1, 1), p2(1, 1)];
    H.A.YData = [p1(2, 1), p2(2, 1)];
    
    addpoints(H.trajA(1), p1(1, 1), p1(2, 1));    
    if(min(vecnorm(p2(:, 1) - p2(:, 2:end))) < Do || crashflag_2)
        addpoints(H.trajA(3), p2(1, 1), p2(2, 1));
        addpoints(H.v_set2_crash, t,  norm(v_set2) );
        crashflag_2 = true;
    else
        addpoints(H.trajA(2), p2(1, 1), p2(2, 1));
        addpoints(H.v_set2, t,  norm(v_set2) );
    end
    
    %   addpoints(H.V,         t, maxLog10_0(1/norm(p(:, 1)-p(:, 2))) );
    addpoints(H.Dist1, t,  min(vecnorm(p1(:, 1) - p1(:, 2:end))));
    addpoints(H.Dist2, t,  min(vecnorm(p2(:, 1) - p2(:, 2:end))));
    addpoints(H.v_set1, t,  norm(v_set1) );
    
    
    % update feasible region and arrow
    half_space_2D(H.feasibleRegion_1,      Ap_i,  -ap_i(:, 1)*delta1(1));
    half_space_2D(H.feasibleRegion_2, A2(:,:,1),    b2(:,:,1)*delta2(1));
    H.v_set1_hdl.UData = v_set1(1, 1); H.v_set1_hdl.VData = v_set1(2, 1);
    H.v_set2_hdl.UData = v_set2(1, 1); H.v_set2_hdl.VData = v_set2(2, 1);
end
end

function dx = m_speed(x, u)
% Model 210918, res_06.mat
A=[-1.5828    2.9188         0         0         0         0         0         0
   -2.9188   -1.5828         0         0         0         0         0         0
         0         0   -2.6833    7.1816         0         0         0         0
         0         0   -7.1816   -2.6833         0         0         0         0
         0         0         0         0   -2.5615    6.8558         0         0
         0         0         0         0   -6.8558   -2.5615         0         0
         0         0         0         0         0         0   -2.1391    3.7051
         0         0         0         0         0         0   -3.7051   -2.1391];
     
B=[ 1.6527         0
    0.6473         0
    1.4972         0
    0.9178         0
         0    1.5791
         0    0.8422
         0    1.5056
         0   -2.2906];
dx = A*x + B*u;
end

function dx = m_position(x, u)
A  = 0*[-1, 0; 0, -1];
B  = [1, 0; 0 1];
dx = A*x + B*u;
end

function x = rungekutta(fun, x0, u, h)
% FcnHandlesUsed  = isa(fun,'function_handle');
k1 = fun(x0       , u);
k2 = fun(x0+h/2*k1, u);
k3 = fun(x0+h/2*k2, u);
k4 = fun(x0+  h*k3, u);
x = x0 + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

function  res = limit(res, lower, upper)
res = min(max(res, lower), upper);
end

%%
function ap = gen_a(Ap, A, a, tau)
N  = size(Ap, 1);
ap = zeros(N, 1);
cos_tau = cos(2*pi/N);              % ????
A_unit     = unit_vector(A);        % R^{m x n}, a\in R{m x 1}
polyA_unit = unit_vector(Ap);       % R^{N x n}
cos_theta  = polyA_unit*A_unit';    % R^{N x m}
chi = cos_theta.*a' - limit(1*(cos_tau-cos_theta), 0, 1).*(tau+cos_theta.*a');
ap  = max(chi, [], 2); % ?????????????????
end

% ?????
function A = getPolyAb(N, Radius, theta_bias)
% 2???
theta    = linspace(theta_bias, theta_bias+2*pi, N+1)';
% position = Radius*[cos(theta(1:end-1)), sin(theta(1:end-1))];
% [A, ~]   = getlines(position);
A = [cos(theta(1:end-1)), sin(theta(1:end-1))];
% fprintf('??????????????? %f \n', min(matrix_minsvd(A)));
end

function [A, a] = getlines(p)
N = size(p, 1);
m = size(p, 2);
A = zeros(N, m);
a = -ones(N, 1);
for idx = 1:N
    p_index   = mod((idx:idx+(m-1))-1, N) + 1;
    p_tmp     = p(p_index, :);
    A(idx, :) = (p_tmp\ones(m, 1))';
    len = norm(A(idx, :));
    A(idx, :) = A(idx, :)/len;
    a(idx, :) = a(idx, :)/len;
end
end


function  M = unit_vector(M)

for idx = 1:size(M,1)
    M(idx, :) = M(idx, :)/ norm(M(idx, :));
end

end


function res = maxLog10_0(in)
res = max(log10(in), 0);
end

function dx = m_idealvelocity(x, u, vp_bar, vpd_bar)
T   = 1;
x_r = u/norm(u)*min(norm(u), vp_bar);
dx = -1/T*(x - x_r);
dx = dx/norm(dx)*min(norm(dx), vpd_bar);
end

function fillwindows(ax)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
end