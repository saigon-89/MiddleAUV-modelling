close all

%% ЗАГРУЗИТЬ ПАРАМЕТРЫ ИЗ ФАЙЛА
rov_param 

%% ТЕНЗОР ИНЕРЦИИ ПРЯМОУГОЛЬНОЙ ПРИЗМЫ
if ~exist("I0",'var') 
  I0  = diag(m*[(W^2+H^2) (L^2+H^2) (W^2+L^2)]/12).*10^-6;
end

%% ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
% преобразование в кососимметричную матрицу
S = @(x)([ 0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0 ]); % (1.5)

%% РАСЧЕТ ЯКОБИАНА
J_k_o = @(eta)[ 1  0           -sin(eta(5)); ...
                0  cos(eta(4))  cos(eta(5))*sin(eta(4)); ...
                0 -sin(eta(4))  cos(eta(5))*cos(eta(4)) ];
R_I_B = @(eta)[ cos(eta(6))*cos(eta(5))                                        sin(eta(6))*cos(eta(5))                                       -sin(eta(5)); ...
               -sin(eta(6))*cos(eta(4)) + cos(eta(6))*sin(eta(5))*sin(eta(4))  cos(eta(6))*cos(eta(4)) + sin(eta(6))*sin(eta(5))*sin(eta(4))  sin(eta(4))*cos(eta(5)); ...
                sin(eta(6))*sin(eta(4)) + cos(eta(6))*sin(eta(5))*cos(eta(4)) -cos(eta(6))*sin(eta(4)) + sin(eta(6))*sin(eta(5))*cos(eta(4))  cos(eta(4))*cos(eta(5)) ];
J = @(eta)[ R_I_B(eta) zeros(3); ...
            zeros(3)   J_k_o(eta) ];

%% РАСЧЕТ МАТРИЦЫ M
% расчет M_RB
M_RB = [ m*eye(3) -m*S(r_g_b); m*S(r_g_b) I0 ]; % (1.4)
% расчет M_A
M_A = rectangular_added_mass(L, H, W, rho, PF, PS, PT); % [6]
% расчет M
M = M_RB + M_A; % (1.3)

%% РАСЧЕТ C(v)
% расчет C_RB(v)
%C_RB = @(v)([ zeros(3) -m*S(v(1:3)); -m*S(v(1:3)) -S(diag(I0).*v(4:end)) ]); % (1.7)
C_RB = @(v)[ zeros(3)                           -m*S(v(1:3))-m*S(S(v(4:end))*r_g_b); ...
            -m*S(v(1:3))-m*S(S(v(4:end))*r_g_b)  m*S(S(v(1:3))*r_g_b)-S(I0*v(4:end)) ];
% расчет C_A(v)
diag_M_A = diag(M_A);
C_A = @(v)([ zeros(3) -S(diag_M_A(1:3).*v(1:3)); % (1.8)
    -S(diag_M_A(1:3).*v(1:3)) -S(diag_M_A(4:end).*v(4:end))]);
% расчет C(v)
C = @(v)(C_RB(v) + C_A(v)); % (1.3)
% M11 = M(1:3,1:3); M12 = M(1:3,4:6); M21 = M(4:6,1:3); M22 = M(4:6,4:6);
% C = @(v)[ zeros(3) -S(M11*v(1:3)+M12*v(4:6)); -S(M11*v(1:3)+M12*v(4:6)) -S(M21*v(1:3)+M22*v(4:6)) ];

%% РАСЧЕТ g(n)
x_g = r_g_b(1); y_g = r_g_b(2); z_g = r_g_b(3);
x_b = r_b_b(1); y_b = r_b_b(2); z_b = r_b_b(3);
g = @(eta)[ (m*9.81-B)*sin(eta(5)); 
           -(m*9.81-B)*cos(eta(5))*sin(eta(4)); 
           -(m*9.81-B)*cos(eta(5))*cos(eta(4));
           -(y_g*m*9.81-y_b*B)*cos(eta(5))*cos(eta(4)) + ...
              (z_g*m*9.81-z_b*B)*cos(eta(5))*sin(eta(4));
            (z_g*m*9.81-z_b*B)*sin(eta(5)) + ...
              (x_g*m*9.81-x_b*B)*cos(eta(5))*cos(eta(4));
           -(x_g*m*9.81-x_b*B)*cos(eta(5))*sin(eta(4)) - ...
              (y_g*m*9.81-y_b*B)*sin(eta(5)) ]; % (1.10)

%% РАСЧЕТ D(v)
[D_LIN, D_QUAD] = rectangular_damping(L, H, W, rho, PF, PS, PT, M_RB, M_A, B, r_g_b, r_b_b);
D = @(v)(D_LIN + D_QUAD.*abs(v)); % [6]

%% МОДЕЛИРОВАНИЕ ДВИЖЕНИЯ (FEEDBACK CONTROL)
eta0 = [0, 0, -2, 0, 0, 0]'; % начальная глубина 2 метра
v0 = [0, 0, 0, 0, 0, 0]';
set = [5 -5 0 0 pi/4 pi/6]'; % уставка для единичной обратной связи
t_end = 60; dt = 0.01;
[t,Y] = ode45(@(t,y)odefcn(t,y,M,C,D,g,J,set,'fb'), 0:dt:t_end, [eta0; v0]);

%% ПОСТРОЕНИЯ ГРАФИКОВ
figure
v = Y(:,7:end);
subplot(2,2,1), title('Скорости (линейные)'), hold on, grid on
plot(t, v(:,1:3)), xlabel('t, сек'), ylabel('Скорость, м/c'), xlim([0 t_end])
legend('u(t)', 'v(t)', 'w(t)', 'Location', 'Best')
subplot(2,2,2), title('Скорости (угловые)'), hold on, grid on
plot(t, v(:,4:end)), xlabel('t, сек'), ylabel('Скорость, рад/c'), xlim([0 t_end])
legend('p(t)', 'q(t)', 'r(t)', 'Location', 'Best')
eta = Y(:,1:6);
subplot(2,2,3), title('Положения (по осям)'), hold on, grid on
plot(t, eta(:,1:3)), xlabel('t, сек'), ylabel('Положения, м'), xlim([0 t_end])
legend('x(t)', 'y(t)', 'z(t)', 'Location', 'Best')
subplot(2,2,4), title('Положения (углы Эйлера)'), hold on, grid on
plot(t, eta(:,4:end)), xlabel('t, сек'), ylabel('Положения, рад'), xlim([0 t_end])
legend('\phi(t)', '\theta(t)', '\psi(t)', 'Location', 'Best')
sgtitle('Позиционирование при управлении обратной связью')

figure, plot3(eta(:,1), eta(:,2), eta(:,3)), title('Позиционирование AUV')
hold on
plot3(set(1), set(2), set(3), 'rO') 
legend('траектория AUV', 'желаемое положение', 'Location', 'Best')
xlabel('x(t)'), ylabel('y(t)'), zlabel('z(t)'), grid on
%set(gca, 'YDir', 'reverse');
view([-1,1,1])
axis equal
legend('off')
rectangular_plot(eta(1,:), L*10^-3, W*10^-3, H*10^-3, '--r')
rectangular_plot(eta(end,:), L*10^-3, W*10^-3, H*10^-3, 'r')

%% МОДЕЛИРОВАНИЕ ДВИЖЕНИЯ ПРИ ПРЕДОПРЕДЕЛЕННОМ СИЛОВОМ ВОЗДЕЙСТВИИ
eta0 = [0, 0, -2, 0, 0, pi]'; % начальная глубина 2 метра
v0 = [0, 0, 0, 0, 0, 0]';
tau = @(t)[0.05 0 0 0 0 60]';
t_end = 60; dt = 0.01;
[t,Y] = ode45(@(t,y)odefcn(t,y,M,C,D,g,J,tau,'simp'), 0:dt:t_end, [eta0; v0]);

%% ПОСТРОЕНИЯ ГРАФИКОВ
figure
v = Y(:,7:end);
subplot(2,2,1), title('Скорости (линейные)'), hold on, grid on
plot(t, v(:,1:3)), xlabel('t, сек'), ylabel('Скорость, м/c'), xlim([0 t_end])
legend('u(t)', 'v(t)', 'w(t)', 'Location', 'Best')
subplot(2,2,2), title('Скорости (угловые)'), hold on, grid on
plot(t, v(:,4:end)), xlabel('t, сек'), ylabel('Скорость, рад/c'), xlim([0 t_end])
legend('p(t)', 'q(t)', 'r(t)', 'Location', 'Best')
eta = Y(:,1:6);
subplot(2,2,3), title('Положения (по осям)'), hold on, grid on
plot(t, eta(:,1:3)), xlabel('t, сек'), ylabel('Положения, м'), xlim([0 t_end])
legend('x(t)', 'y(t)', 'z(t)', 'Location', 'Best')
subplot(2,2,4), title('Положения (углы Эйлера)'), hold on, grid on
plot(t, eta(:,4:end)), xlabel('t, сек'), ylabel('Положения, рад'), xlim([0 t_end])
legend('\phi(t)', '\theta(t)', '\psi(t)', 'Location', 'Best')
sgtitle('Позиционирование при силовом воздействии \tau')

figure, plot3(eta(:,1), eta(:,2), eta(:,3)), title('Позиционирование AUV')
hold on
plot3(eta(end,1), eta(end,2), eta(end,3), 'rO') 
legend('траектория AUV', 'желаемое положение', 'Location', 'Best')
xlabel('x(t)'), ylabel('y(t)'), zlabel('z(t)'), grid on
%set(gca, 'YDir', 'reverse');
view([-1,1,1])
axis equal
legend('off')
rectangular_plot(eta(1,:), L*10^-3, W*10^-3, H*10^-3, '--r')
rectangular_plot(eta(end,:), L*10^-3, W*10^-3, H*10^-3, 'r')

%auv_animation(eta,L,W,H,dt/1000)

figure
tau_val = zeros(numel(t),6);
for i=1:length(tau_val)
    tau_val(i,:) = tau(t(i));
end
subplot(1,2,1), title('Действующие силы'), hold on, grid on
plot(t, tau_val(:,1:3)), xlabel('t, сек'), ylabel('Силы, Н'), xlim([0 t_end])
legend('\tau_1(t)', '\tau_2(t)', '\tau_3(t)', 'Location', 'Best')
subplot(1,2,2), title('Прилагаемые моменты'), hold on, grid on
plot(t, tau_val(:,4:end)), xlabel('t, сек'), ylabel('Моменты, Н*м'), xlim([0 t_end])
legend('\tau_4(t)', '\tau_5(t)', '\tau_6(t)', 'Location', 'Best')
sgtitle('Позиционирование при силовом воздействии \tau(t)')

%% ARX MODELLING
u = tau_val;
y = eta;
z = iddata(y, u, dt, 'Name', 'AUV', 'TimeUnit', 's', ... 
    'OutputName', {'x','y','z','phi','theta','psi'}, ...
    'InputName', {'tau_1','tau_2','tau_3','tau_4','tau_5','tau_6'});
na = eye(6); nb = ones(6,6); nk = ones(6,6);
NN = [na, 2*nb, nk];
mw2 = nlarx(z, NN, idWaveletNetwork);
%getreg(mw2)

figure, compare(z,mw2)