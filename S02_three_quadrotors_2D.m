%% MATLAB Version R2017b
clear all; clc; close all;

Time= 20;               % simulation time, sec
T   = 50e-3;             % step, sec
N   = floor(Time/T);    % steps
t   = 0;

% agent and obstacle setting
na = 3;                 % number of mobile agents
x  = zeros(8, na);      % actuation state of agent
vp = zeros(2, na);      % ideal velocity
v_set = zeros(2, na);   % the input of the actuation
delta = zeros(1, na);   % the slack variable

% define the initial position
theta = linspace(0, 2*pi, 4);
p     = [cos(theta(1:3)+pi/6); sin(theta(1:3)+pi/6)];
p_aim = [cos(theta(1:3)+pi/1); sin(theta(1:3)+pi/1)];

% quadraitc programming setting
Ap_i= getPolyAb(11, 1, 0);      % 
ap_i= zeros(size(Ap_i, 1), na);
A   = zeros(na-1, 2, na);
b   = zeros(na-1, 1, na);
Ds  = 0.3;
Do  = 0.25;
k_delta = 100;
alpha_v = 0.1;
vpd_bar = 1;
vp_bar  = 1;

% initialize graphic
H = myGraphic([], p, x, Ds, Ap_i, ap_i);

% simulation process
for k = 1:N
    % sample time
    
    t = k*T;
    
    % update the graphic
    H = myGraphic(H, p, x, Ds, Ap_i, ap_i);
    
    x_tmp = x;
    p_tmp = p;

    for iA = 1:na
        % ideal velocity
        k_vp = 1.0;
        
        if(iA == 1)
            p_aim(:, iA) = -[cos(2*pi*t/10); cos(2*pi*t/10)];
        elseif(iA == 2)
            p_aim(:, iA) =  [cos(2*pi*t/10); cos(2*pi*t/10)] + [-1; 1]*Ds/sqrt(2);
        elseif(iA == 3)
            p_aim(:, iA) =  [cos(2*pi*t/10); cos(2*pi*t/10)] + [1; -1]*Ds/sqrt(2);
        end
            
        
        vp(:, iA)= rungekutta(@(x,u)m_idealvelocity(x, u, vp_bar, vpd_bar), vp(:, iA), k_vp*(p_aim(:, iA)-p(:, iA)), T);     % model of ideal velocity
        
        % collision avoidance
%         [v_set(:,iA), delta(iA), ap_i(:,iA), A(:,:,iA), b(:,:,iA)] = collision_avoidance_reshape(p_tmp, vp, Ds, Ap_i, k_delta, alpha_v, iA);
        [v_set(:,iA), delta(iA),             A(:,:,iA), b(:,:,iA)] = collision_avoidance(p_tmp, vp, Ds, k_delta, alpha_v, iA);
        %     v_set1(:,1) = vp1(:, 1);
        %     v_set2(:,1) = vp2(:, 1);
        
        % update model
        x(:, iA) = rungekutta(@m_speed,    x_tmp(:, iA), v_set(  :, iA), T);     % velocity model
        C=[ 0.7761   -1.9815         0         0    2.1315    2.4144         0         0
             0         0   -2.2000   -2.8153         0         0   -1.5060   -0.9899];
        p(:, iA) = rungekutta(@m_position, p_tmp(:, iA),  C*x(:, iA), T);     % position model
    end
end
%% end of main.

%% 
function [v_set, delta, A, b] = collision_avoidance(p, vp, Ds, delta_bar, alpha_v, iA)
persistent opt;
if(isempty(opt))
    opt = optimoptions('quadprog',  'Algorithm','interior-point-convex','Display','off', 'MaxIterations', 1e3);
end

pm    =  p(:, iA);
v_aim = vp(:, iA);
na    = size(p, 2);

% generate conditions
A = zeros(na-1, 2);
b = zeros(na-1, 1);
count = 0;

for k = 1:na
    if(k == iA)
        continue;
    end
    
    count = count + 1;
    p_delta = pm - p(:, k);     % relative position
    dist    = norm(p_delta, 2); % distance 
    h = dist - Ds;              % barrier function 
    V = 1/(h+Ds);               % Lyapunov function
    
    % condition - count
    A(count, :) = -p_delta'/norm(p_delta);
%     b(count, :) = alpha_v*h;
    b(count, :) = -alpha_v*(V-1/Ds)/delta_bar;
%     b(count, :) = alpha_v*(dist^2-Ds^2);
end

% Quadratic programming
f_aim = [v_aim; delta_bar];

% Standard Quadratic Programming with Slack Variables
v_set_det = quadprog(...
    diag([1, 1, 1]), -f_aim,...                             % cost function
    [A, -b], zeros(size(A, 1), 1), [], [],  ...             % constraints
    [-inf; -inf; 0], [inf; inf; inf], zeros(3, 1), opt);    % limits and setting

v_set  = v_set_det(1:end-1);
delta  = v_set_det(end);
end

function [v_set, delta, ap_i, A, b] = collision_avoidance_reshape(p, vp, Ds, Ap_i, delta_bar, alpha_v, iA)
persistent opt;
if(isempty(opt))
    opt = optimoptions('quadprog',  'Algorithm','interior-point-convex','Display','off', 'MaxIterations', 1e3);
end

pm    =  p(:, iA);
v_aim = vp(:, iA);
na    = size(p, 2);
xi    = 1.6;

% generate conditions
A = zeros(na-1, 2);
b = zeros(na-1, 1);
count = 0;

for k = 1:na
    if(k == iA)
        continue;
    end
    
    count = count + 1;
    p_delta = pm - p(:, k);     % relative position
    dist    = norm(p_delta, 2); % distance 
    h = dist - Ds;              % barrier function 
    V = 1/(h+Ds);               % Lyapunov function
    
    % condition - count
    A(count, :) = -p_delta'/norm(p_delta);
%     b(count, :) = alpha_v*h;
    b(count, :) = -alpha_v*(V-1/Ds)/delta_bar;
%     b(count, :) = alpha_v*(dist^2-Ds^2);
end

% quadratic programming
f_aim = [v_aim; delta_bar];
ap_i  = gen_a(Ap_i, A, -b, 1/delta_bar);  % Ap_i x + a delta \le 0
M     = [Ap_i, ap_i];

% Quadratic Programming with Slack Variables - Lipschitz
[v_set_det, ~, ~, ~, ~] = quadprog(...
    diag([1, 1, 1]), -f_aim,...                             % cost function
    M, zeros(size(M, 1), 1), [], [],  ...                   % constraints
    [-inf; -inf; 0], [inf; inf; inf], zeros(3, 1), opt);    % limits and setting

v_set  = v_set_det(1:end-1);
delta  = v_set_det(end);
end

% function H = myGraphic(H, p1, x1, Ds, Ap_i, ap_i, delta)
function H = myGraphic(H, p, x, Ds, Ap_i, ap_i)
na = size(p, 2);
Do  = evalin('base', 'Do');
drawnow;
if(isempty(H))
    % initialize graphic
    figure(1); clf; 
    axis equal;
    hold on; 
    grid on;
    set(gcf, 'position', [0,0, 800, 500 ], 'color', 'w');
    set(gca, 'fontsize', 16);
    axis([-2.5, 2.5, -2, 2]);
    
    % mobile agents
    HandleVisibility  = 'on';
    for iA = 1:na
        H.shapeDs = 0.5*Ds*[sin(linspace(0, 2*pi, 30)); cos(linspace(0, 2*pi, 30))];
        H.shapeDo = 0.5*Do*[sin(linspace(0, 2*pi, 30)); cos(linspace(0, 2*pi, 30))];
        H.shapeA(iA) = patch(p(1,iA)+H.shapeDs(1,:), p(2,iA)+H.shapeDs(2,:), 'c', 'facealpha', 0.2, 'HandleVisibility', 'off');
        H.realA(iA)  = patch(p(1,iA)+H.shapeDo(1,:), p(2,iA)+H.shapeDo(2,:), 'c', 'facealpha', 0.9, 'HandleVisibility', 'off');
        H.trajA(iA) = animatedline('color', 'r', 'HandleVisibility', HandleVisibility, 'linewidth', 1.5, 'linestyle', '-', 'marker', 'none', 'Displayname', 'reshaped'); 
    end
    H.A = line(p(1, 1:na), p(2, 1:na), 'Marker', 'o', 'linestyle', 'none', 'color','k', 'HandleVisibility', 'off');
    legend;
else
    % get variables from workspace
    t     = evalin('base', 't');
    
    % mobile agents
    for iA = 1:na
        H.shapeA(iA).XData = p(1,iA)+H.shapeDs(1,:);
        H.shapeA(iA).YData = p(2,iA)+H.shapeDs(2,:);
        H.realA(iA).XData = p(1,iA)+H.shapeDo(1,:);
        H.realA(iA).YData = p(2,iA)+H.shapeDo(2,:);
        addpoints(H.trajA(iA), p(1, iA), p(2, iA));
    end
    
    % position of mobile agents
    H.A.XData = p(1, 1:na);
    H.A.YData = p(2, 1:na);
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
cos_tau = cos(2*pi/N);              % 最低角度
A_unit     = unit_vector(A);        % R^{m x n}, a\in R{m x 1}
polyA_unit = unit_vector(Ap);       % R^{N x n}
cos_theta  = polyA_unit*A_unit';    % R^{N x m}
chi = cos_theta.*a' - limit(1*(cos_tau-cos_theta), 0, 1).*(tau+cos_theta.*a');
ap  = max(chi, [], 2); % 选取每一行的最小元素，组成列向量。
end

% 获取二维点
function A = getPolyAb(N, Radius, theta_bias)
% 2维情况
theta    = linspace(theta_bias, theta_bias+2*pi, N+1)';
position = Radius*[cos(theta(1:end-1)), sin(theta(1:end-1))];
[A, ~]   = getlines(position);
% fprintf('幂集中的有效集合的最小奇异值： %f \n', min(matrix_minsvd(A)));
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