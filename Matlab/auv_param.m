%% Middle AUV parameters
L = 213.799; % ����� AUV [mm]
H = 89.07; % ������ AUV [mm]
W = 193.866; % ������ AUV [mm]
rho = 1000; % ��������� �������� [kg/m^3]
PF = 10637.097; % ������� ��������� �������� [mm^2]
PS = 23902.237; % ������� ������� �������� [mm^2]
PT = 13724.628; % ������� ������� �������� [mm^2]
m = 1421.943 * 10^-3; % ����� AUV [kg]
r_b_b = [0; 0; -21]; % ����� ����������
r_g_b = [0; 0; 0]; % ����� ���� [mm]
r_g_b = r_g_b / 10^3; % ��������� [mm] � [m]
r_b_b = r_b_b / 10^3; % ��������� [mm] � [m]
B = m * 9.81;

I0 = [ 2.551E+06 -612.421   -1.730E+05;
      -612.421    3.174E+06  320.636;
      -1.730E+05  320.636    5.029E+06 ]; % [g mm^2]
I0 = I0 * 10^-9;

%% THRUSTER ALLOCATION
T = [ 9,   74,  -10; 
     -75,  74,   0; 
     -75, -74,   0;
      9,  -74,  -10 ] .* 10^-3;

F_max = 0.2 * 9.80665; % [���] � [�]

f1 = [ 0;     0; F_max ];
f2 = [ F_max; 0; 0     ];
f3 = [ F_max; 0; 0     ];
f4 = [ 0;     0; F_max ];

r1 = T(1,:)';
r2 = T(2,:)';
r3 = T(3,:)';
r4 = T(4,:)';

tau1 = [ f1; cross(r1, f1) ];
tau2 = [ f2; cross(r2, f2) ];
tau3 = [ f3; cross(r3, f3) ];
tau4 = [ f4; cross(r4, f4) ];

tau_c = [ tau1, tau2, tau3, tau4 ];