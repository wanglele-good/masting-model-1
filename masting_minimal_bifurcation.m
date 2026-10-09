%% masting_minimal_bifurcation.m
% 最小连续 Resource-Reproduction 模型
%
% dR/dt = r*R*(1-R/K) - c*R*B/(1+h*R)
% dB/dt = e*c*R*B/(1+h*R) - d*B
%
% 目的：
% 1. 扫描 K，检查平衡点、trace(J)、det(J)
% 2. 找 Hopf 候选点
% 3. 看 Hopf 附近的时间序列和相图
% 4. 检查正平衡点是否出现折叠(saddle-node)
%
% 注意：
% 对这个“最小模型”，正平衡点只有一个显式分支，
% 因此改变 K 不会产生内部 saddle-node fold。
% K = R* 处正分支与 B=0 的平衡支相交，更接近
% transcritical，而不是 fold。
%
% 推荐先直接运行，不要改方程。

clear; clc; close all;

%% =========================
% 1. 基本参数
% ==========================
r = 1.0;
K = 3.0;
c = 1.0;
h = 0.5;
e = 0.5;
d = 0.2;

% 初始条件
R0 = 0.6;
B0 = 1.0;

%% =========================
% 2. 先计算正平衡点存在条件
% ==========================
den = e*c - d*h;

if den <= 0
    error('当前参数满足 e*c <= d*h，没有正的内部平衡点。');
end

Rstar = d / den;

if K <= Rstar
    fprintf('当前 K = %.4f <= R* = %.4f：正内部平衡点不存在。\n', K, Rstar);
else
    Bstar = (r/c)*(1 - Rstar/K)*(1 + h*Rstar);
    fprintf('正平衡点：R* = %.6f, B* = %.6f\n', Rstar, Bstar);
end

%% =========================
% 3. 当前参数下 Jacobian
% ==========================
if K > Rstar
    [J, trJ, detJ, eigJ] = jacobian_at_eq(Rstar,Bstar,r,K,c,h,e,d);

    fprintf('\n当前参数：\n');
    fprintf('trace(J) = %.8f\n', trJ);
    fprintf('det(J)   = %.8f\n', detJ);
    fprintf('eigenvalues:\n');
    disp(eigJ);

    if abs(trJ) < 1e-5 && detJ > 0
        fprintf('>>> 当前参数非常接近 Hopf 候选点。\n');
    elseif trJ < 0 && detJ > 0
        fprintf('>>> 正平衡点局部稳定（二维系统条件）。\n');
    elseif trJ > 0 && detJ > 0
        fprintf('>>> 正平衡点局部不稳定，可能进入极限环等动力学。\n');
    end
end

%% =========================
% 4. 扫描 K：寻找 Hopf
% ==========================
Kmin = max(Rstar*1.001, 0.55);
Kmax = 8.0;
NK = 2000;

Kvec = linspace(Kmin,Kmax,NK);

Rvec = nan(size(Kvec));
Bvec = nan(size(Kvec));
traceVec = nan(size(Kvec));
detVec = nan(size(Kvec));
lambda1 = nan(size(Kvec));
lambda2 = nan(size(Kvec));

for i = 1:length(Kvec)

    Ki = Kvec(i);

    Rsi = Rstar;

    if Ki > Rsi
        Bsi = (r/c)*(1 - Rsi/Ki)*(1 + h*Rsi);

        [~, tr_i, det_i, eig_i] = ...
            jacobian_at_eq(Rsi,Bsi,r,Ki,c,h,e,d);

        Rvec(i) = Rsi;
        Bvec(i) = Bsi;
        traceVec(i) = tr_i;
        detVec(i) = det_i;

        lambda1(i) = eig_i(1);
        lambda2(i) = eig_i(2);
    end
end

%% 理论 Hopf 点
% 对当前模型：
% K_H = 1/h + 2*R*
K_H = 1/h + 2*Rstar;

fprintf('\n=============================\n');
fprintf('理论 Hopf 候选参数：\n');
fprintf('K_H = %.8f\n', K_H);
fprintf('=============================\n');

%% =========================
% 5. 图：平衡点分支
% ==========================
figure;

plot(Kvec,Rvec,'LineWidth',1.8);
hold on;
plot(Kvec,Bvec,'LineWidth',1.8);
xline(K_H,'--','Hopf candidate');

xlabel('K');
ylabel('Equilibrium value');
legend('R^*','B^*','K_H','Location','best');
title('Equilibrium branches');
grid on;

%% =========================
% 6. 图：trace(J)
% ==========================
figure;

plot(Kvec,traceVec,'LineWidth',1.8);
hold on;
yline(0,'k--');
xline(K_H,'--','Hopf');

xlabel('K');
ylabel('trace(J)');
title('Trace of Jacobian');
grid on;

%% =========================
% 7. 图：det(J)
% ==========================
figure;

plot(Kvec,detVec,'LineWidth',1.8);
hold on;
yline(0,'k--');

xlabel('K');
ylabel('det(J)');
title('Determinant of Jacobian');
grid on;

%% =========================
% 8. 图：最大特征值实部
% ==========================
maxReal = max(real([lambda1;lambda2]),[],1);

figure;

plot(Kvec,maxReal,'LineWidth',1.8);
hold on;
yline(0,'k--');
xline(K_H,'--','Hopf');

xlabel('K');
ylabel('max Re(\lambda)');
title('Stability indicator');
grid on;

%% =========================
% 9. 在 Hopf 前后模拟
% ==========================
Ktest = [2.5, K_H, 3.5];

figure;

for j = 1:length(Ktest)

    Kj = Ktest(j);

    Rj = Rstar;
    Bj = (r/c)*(1 - Rj/Kj)*(1 + h*Rj);

    % 给一个很小扰动，避免恰好落在平衡点
    y0 = [Rj*1.01; Bj*0.99];

    tspan = [0 300];

    [t,y] = ode45(@(t,y) rhs(t,y,r,Kj,c,h,e,d), ...
                  tspan,y0);

    subplot(3,1,j);

    plot(t,y(:,1),'LineWidth',1.2);
    hold on;
    plot(t,y(:,2),'LineWidth',1.2);

    xlabel('Time');
    ylabel('State');

    title(sprintf('K = %.4f',Kj));

    legend('R','B');
    grid on;

end

%% =========================
% 10. 相平面
% ==========================
figure;

for j = 1:length(Ktest)

    Kj = Ktest(j);

    Rj = Rstar;
    Bj = (r/c)*(1 - Rj/Kj)*(1 + h*Rj);

    y0 = [Rj*1.01; Bj*0.99];

    [t,y] = ode45(@(t,y) rhs(t,y,r,Kj,c,h,e,d), ...
                  [0 300],y0);

    subplot(1,3,j);

    plot(y(:,1),y(:,2),'LineWidth',1.2);
    hold on;

    plot(Rj,Bj,'ko','MarkerFaceColor','k');

    xlabel('R');
    ylabel('B');

    title(sprintf('K = %.4f',Kj));

    grid on;

end

%% =========================
% 11. Nullclines at Hopf
% ==========================
K = K_H;

Rplot = linspace(0.001,K*1.1,1000);

% dR/dt = 0:
% B = r(1-R/K)(1+hR)/c
B_null_R = r*(1-Rplot/K).*(1+h*Rplot)/c;

% dB/dt = 0:
% B = 0 OR R = R*
B_null_B = nan(size(Rplot));
B_null_B(:) = 0;

figure;

plot(Rplot,B_null_R,'LineWidth',1.8);
hold on;

plot(Rplot,B_null_B,'LineWidth',1.8);

xline(Rstar,'--','R^*');

plot(Rstar, ...
     (r/c)*(1-Rstar/K)*(1+h*Rstar), ...
     'ko','MarkerFaceColor','k');

xlabel('R');
ylabel('B');
legend('dR/dt = 0','dB/dt = 0','R^*','Equilibrium');

title('Nullclines at Hopf candidate');
grid on;


%% =========================
% Local functions
% ==========================

function dydt = rhs(~,y,r,K,c,h,e,d)

    R = y(1);
    B = y(2);

    dR = r*R*(1-R/K) ...
         - c*R*B/(1+h*R);

    dB = e*c*R*B/(1+h*R) ...
         - d*B;

    dydt = [dR; dB];

end


function [J,trJ,detJ,eigJ] = ...
    jacobian_at_eq(R,B,r,K,c,h,e,d)

    J11 = r*(1-2*R/K) ...
          - c*B/(1+h*R)^2;

    J12 = -c*R/(1+h*R);

    J21 = e*c*B/(1+h*R)^2;

    J22 = e*c*R/(1+h*R) - d;

    J = [J11 J12;
         J21 J22];

    trJ = trace(J);
    detJ = det(J);

    eigJ = eig(J);

end
