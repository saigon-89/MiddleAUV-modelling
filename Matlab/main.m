close all
clear all

%% ��������� ��������� �� �����
auv_param 

%% ������ ������� ������������� ������
if ~exist("I0",'var') 
  l = L * 10^-3;
  w = W * 10^-3;
  h = H * 10^-3;
  I0  = (m/12).*diag([ (w^2 + h^2) (l^2 + h^2) (l^2 + w^2) ]);
end

%% ��������������� �������
% �������������� � ���������������� �������
S = @(x)([ 0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0 ]); 

%% ������ ��������
J_k_o = @(eta)[ 1  0           -sin(eta(5)); ...
                0  cos(eta(4))  cos(eta(5))*sin(eta(4)); ...
                0 -sin(eta(4))  cos(eta(5))*cos(eta(4)) ];
R_I_B = @(eta)[ cos(eta(6))*cos(eta(5))                                        sin(eta(6))*cos(eta(5))                                       -sin(eta(5)); ...
               -sin(eta(6))*cos(eta(4)) + cos(eta(6))*sin(eta(5))*sin(eta(4))  cos(eta(6))*cos(eta(4)) + sin(eta(6))*sin(eta(5))*sin(eta(4))  sin(eta(4))*cos(eta(5)); ...
                sin(eta(6))*sin(eta(4)) + cos(eta(6))*sin(eta(5))*cos(eta(4)) -cos(eta(6))*sin(eta(4)) + sin(eta(6))*sin(eta(5))*cos(eta(4))  cos(eta(4))*cos(eta(5)) ];
J = @(eta)[ R_I_B(eta) zeros(3); ...
            zeros(3)   J_k_o(eta) ];

%% ������ ������� M
% ������ M_RB
M_RB = [ m*eye(3) -m*S(r_g_b); m*S(r_g_b) I0 ]; 
% ������ M_A
M_A = rectangular_added_mass(L, H, W, rho, PF, PS, PT);
% ������ M
M = M_RB + M_A; 

%% ������ C(v)
% ������ C_RB(v)
%C_RB = @(v)([ zeros(3) -m*S(v(1:3)); -m*S(v(1:3)) -S(diag(I0).*v(4:end)) ]); 
C_RB = @(v)[ zeros(3)                           -m*S(v(1:3))-m*S(S(v(4:end))*r_g_b); ...
            -m*S(v(1:3))-m*S(S(v(4:end))*r_g_b)  m*S(S(v(1:3))*r_g_b)-S(I0*v(4:end)) ];
% ������ C_A(v)
diag_M_A = diag(M_A);
C_A = @(v)([ zeros(3) -S(diag_M_A(1:3).*v(1:3)); 
    -S(diag_M_A(1:3).*v(1:3)) -S(diag_M_A(4:end).*v(4:end))]);
% ������ C(v)
C = @(v)(C_RB(v) + C_A(v)); 
% M11 = M(1:3,1:3); M12 = M(1:3,4:6); M21 = M(4:6,1:3); M22 = M(4:6,4:6);
% C = @(v)[ zeros(3) -S(M11*v(1:3)+M12*v(4:6)); -S(M11*v(1:3)+M12*v(4:6)) -S(M21*v(1:3)+M22*v(4:6)) ];

%% ������ g(n)
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
              (y_g*m*9.81-y_b*B)*sin(eta(5)) ]; 

%% ������ D(v)
[D_LIN, D_QUAD] = rectangular_damping(L, H, W, rho, PF, PS, PT, M_RB, M_A, B, r_g_b, r_b_b);
D = @(v)(D_LIN + D_QUAD.*abs(v));

%% ������������� �������� ��� ���������������� ������� �����������
eta0 = [0, 0, 0, 0, 0, 0]'; 
v0 = [0, 0, 0, 0, 0, 0]';
tau = @(code)tau_c*[0.01;0.01;0.01;0.01];
t_end = 60; dt = 0.01;
[t,Y] = ode45(@(t,y)odefcn(t,y,M,C,D,g,J,tau,'simp'), 0:dt:t_end, [eta0; v0]);

%% ���������� ��������
figure
v = Y(:,7:end);
subplot(2,2,1), title('�������� (��������)'), hold on, grid on
plot(t, v(:,1:3)), xlabel('t, ���'), ylabel('��������, �/c'), xlim([0 t_end])
legend('u(t)', 'v(t)', 'w(t)', 'Location', 'Best')
subplot(2,2,2), title('�������� (�������)'), hold on, grid on
plot(t, v(:,4:end)), xlabel('t, ���'), ylabel('��������, ���/c'), xlim([0 t_end])
legend('p(t)', 'q(t)', 'r(t)', 'Location', 'Best')
eta = Y(:,1:6);
subplot(2,2,3), title('��������� (�� ����)'), hold on, grid on
plot(t, eta(:,1:3)), xlabel('t, ���'), ylabel('���������, �'), xlim([0 t_end])
legend('x(t)', 'y(t)', 'z(t)', 'Location', 'Best')
subplot(2,2,4), title('��������� (���� ������)'), hold on, grid on
plot(t, eta(:,4:end)), xlabel('t, ���'), ylabel('���������, ���'), xlim([0 t_end])
legend('\phi(t)', '\theta(t)', '\psi(t)', 'Location', 'Best')
sgtitle('���������������� ��� ����������� ����������')

figure, plot3(eta(:,1), eta(:,2), eta(:,3)), title('���������������� AUV')
hold on
plot3(eta(end,1), eta(end,2), eta(end,3), 'rO') 
legend('���������� AUV', '�������� ���������', 'Location', 'Best')
xlabel('x(t)'), ylabel('y(t)'), zlabel('z(t)'), grid on
set(gca, 'YDir', 'reverse');
set(gca, 'ZDir', 'reverse');
view([-1,1,1])
axis equal
legend('off')
rectangular_plot(eta(1,:), L*10^-3, W*10^-3, H*10^-3, '--r')
rectangular_plot(eta(end,:), L*10^-3, W*10^-3, H*10^-3, 'r')

% auv_animation(eta,L,W,H,dt/1000)
