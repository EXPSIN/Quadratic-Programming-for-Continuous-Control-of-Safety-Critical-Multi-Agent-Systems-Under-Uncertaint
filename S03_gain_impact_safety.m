clear all; close all; clc;

Time= 5;                 % simulation time, sec
T   = 10e-3;             % step, sec
N   = floor(Time/T);     % steps
t   = 0;
p   = 0.4;
x   = zeros(2, 1);
Ds  = 0.8;
alpha_v =[1,4,8,16];


figure(1);
set(gcf, 'position', [10, 10, 800, 400], 'color', 'w');
subplot(1, 5, [2, 5]); hold on;

px = [  0, Time, Time,  0,  0];
py = [-Ds,  -Ds,   Ds, Ds, -Ds];
patch(px, py+1, 'c', 'facealpha', 0.2, 'displayname', 'unsafe area');
patch(px, py-1, 'c', 'facealpha', 0.2, 'displayname', 'unsafe area', 'handlevisibility', 'off');
klinestyle = {'-', '--', ':', '-.'};
for i_Sim  = 1:length(alpha_v)
    H(i_Sim) = animatedline('color', 'k', 'displayname', sprintf('$p_1(t)~(k=%d)$', alpha_v(i_Sim)), 'linewidth', 2, 'linestyle', klinestyle{i_Sim}); 
end
plot([0, Time], [ 1,  1], 'r--', 'linewidth', 2, 'displayname', '$p_2$');
plot([0, Time], [-1, -1], 'r--', 'linewidth', 2, 'displayname', '$p_3$');
set(gca, 'position', [0.1, 0.2, 0.9, 0.7]);


% lH = legend('location', 'southwest','Orientation','horizontal');
lH = legend('location', 'northeastoutside');
set(lH, 'Interpreter','latex');
axis([0, Time, -1.5, 1.5]);
xlabel('time (sec)');
ylabel('position (m)');
set(gca, 'fontsize', 16);

% subplot(1, 5, 1); hold on;
% plot([0, 0], [1, -1], 'k*','linestyle', 'none');
% H2 = line('marker', 'o'); 
% set(gca, 'position', [0.1, 0.1, 0.1, 0.8]);
% axis equal;
% xticks([]);
% ylabel('position (meter)');
% set(gca, 'fontsize', 16);

% simulation process
for i_Sim  = 1:length(alpha_v)
    p   = 0.4;
    x   = zeros(2, 1);
    for k = 1:N
        % sample time
        t = k*T;

        % collision avoidance
        v_set = collision_avoidance(p, 0, Ds, alpha_v(i_Sim));

        % update model
        x = rungekutta(@m_speed,    x, v_set, T);     % velocity model
        C=[ 0, 1];
        p = rungekutta(@m_position, p,  C*x, T);     % position model

        addpoints(H(i_Sim), t, p);
    %     H2.XData = 0;
    %     H2.YData = p;
        drawnow limitrate;
    end
end



function dx = m_speed(x, u)
A=[-1, 2; -2 -1];
B=[-2.5; 0];
C=[0, 1];
dx = A*x + B*u;
end

function dx = m_position(x, u)
A  = 0;
B  = 1;
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


function v_set= collision_avoidance(p, vp, Ds, alpha_v)
persistent opt;
if(isempty(opt))
    opt = optimoptions('quadprog',  'Algorithm','interior-point-convex','Display','off', 'MaxIterations', 1e3);
end

% generate conditions
A = zeros(2, 1);
A(1, :) = -(p-1)'/norm(p-1);
A(2, :) = -(p+1)'/norm(p+1);
% a(1, :) = alpha_v*(1/norm(p-1)-1/Ds);
% a(2, :) = alpha_v*(1/norm(p+1)-1/Ds);
a(1, :) = -alpha_v*(norm(p-1)-Ds);
a(2, :) = -alpha_v*(norm(p+1)-Ds);
b = -a;

% Quadratic programming
% Standard Quadratic Programming with Slack Variables
v_set = quadprog(...
    1, -vp, ...                             % cost function
    A,   b, [], [],  ...             % constraints
    [], [], zeros(1, 1), opt);    % limits and setting

if(isempty(v_set))
    1;
end
end