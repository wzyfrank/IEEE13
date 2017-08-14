function [Voltage_output, phasors] = ieee13_iter(loads, Zbus, Cbus)

tapARG60 = 1.0 + 0.00625 * 0;
tapBRG60 = 1.0 + 0.00625 * 0;
tapCRG60 = 1.0 + 0.00625 * 0;
alphaRG60 = [tapARG60; tapBRG60; tapCRG60];
alphaMRG60 = alphaRG60 * alphaRG60';

Z632645 = Zbus([30, 31],[19, 35]);
Z632633 = Zbus([29, 30, 31],[10, 11, 12]);
Z633634 = Zbus([10, 11, 12],[13, 14, 15]);
Z645646 = Zbus([19, 35],[20, 21]);
ZRG60632 = Zbus([7, 8, 9],[29, 30, 31]);
Z684652 = Zbus([37],[28]);
Z632671 = Zbus([29, 30, 31],[16, 17, 18]);
Z671684 = Zbus([16, 18],[37, 38]);
Z671680 = Zbus([16, 17, 18],[32, 33, 34]);
Z671692 = Zbus([16, 17, 18],[23, 36, 22]);
Z684611 = Zbus([38],[27]);
Z692675 = Zbus([23, 36, 22],[24, 25, 26]);
Z650RG60 = Zbus([4, 5, 6],[7, 8, 9]);
ZSOURCEBUS650 = Zbus([1, 2, 3],[4, 5, 6]);


% three phase voltage at slack bus
Vbase = 4160 / sqrt(3);
v0=1.05 * Vbase * [0,sqrt(3),0]';
% voltage upper and lower bounds
V_lb = 0.90 * Vbase;
V_ub = 1.10 * Vbase;
v_lb = V_lb * V_lb;
v_ub = V_ub * V_ub;

% sequential component parameters
a = -0.5 + 0.5 * 1i * sqrt(3);
A = 1/sqrt(3) * [1,1,1; 1, a*a, a; 1, a, a*a];
AH = 1/sqrt(3) * [1,1,1; 1, a, a*a; 1, a*a, a];



cvx_begin sdp quiet
% the solver: 
% cvx_solver SeDuMi;
cvx_solver Mosek;

% voltage square variables
variable v645(2,2) hermitian
variable v632_abc(3,3) hermitian
variable v633(3,3) hermitian
variable v634(3,3) hermitian
variable v646(2,2) hermitian
variable v632(3,3) hermitian
variable v652(1,1) hermitian
variable v671(3,3) hermitian
variable v684(2,2) hermitian
variable v671_abc(3,3) hermitian
variable v680(3,3) hermitian
variable v692(3,3) hermitian
variable v611(1,1) hermitian
variable v675(3,3) hermitian
variable vRG60(3,3) hermitian
variable v650(3,3) hermitian
variable vSOURCEBUS(3,3) hermitian

% complex power variables
variable S632645(2,2) complex
variable S632633(3,3) complex
variable S633634(3,3) complex
variable S645646(2,2) complex
variable SRG60632(3,3) complex
variable S684652(1,1) complex
variable S632671(3,3) complex
variable S671684(2,2) complex
variable S671680(3,3) complex
variable S671692(3,3) complex
variable S684611(1,1) complex
variable S692675(3,3) complex
variable S650RG60(3,3) complex
variable SSOURCEBUS650(3,3) complex

% current square variables
variable l632645(2,2) hermitian
variable l632633(3,3) hermitian
variable l633634(3,3) hermitian
variable l645646(2,2) hermitian
variable lRG60632(3,3) hermitian
variable l684652(1,1) hermitian
variable l632671(3,3) hermitian
variable l671684(2,2) hermitian
variable l671680(3,3) hermitian
variable l671692(3,3) hermitian
variable l684611(1,1) hermitian
variable l692675(3,3) hermitian
variable l650RG60(3,3) hermitian
variable lSOURCEBUS650(3,3) hermitian


minimize(trace(real(Z632645*l632645)) + trace(real(A * Z632633*l632633 * AH)) + trace(real(A * Z633634*l633634 * AH)) + trace(real(Z645646*l645646)) + trace(real(A * ZRG60632*lRG60632 * AH)) + trace(real(Z684652*l684652)) + trace(real(A * Z632671*l632671 * AH)) + trace(real(Z671684*l671684)) + trace(real(A * Z671680*l671680 * AH)) + trace(real(A * Z671692*l671692 * AH)) + trace(real(Z684611*l684611)) + trace(real(A * Z692675*l692675 * AH)) + trace(real(A * Z650RG60*l650RG60 * AH)) + trace(real(A * ZSOURCEBUS650*lSOURCEBUS650 * AH)) + 0)
subject to


% constraints: 
% (1): voltage lower and upper bounds 
v_lb <= diag(v645) <= v_ub;
v_lb <= diag(A * v633 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v634 * ctranspose(A)) <= v_ub;
v_lb <= diag(v646) <= v_ub;
v_lb <= diag(A * v632 * ctranspose(A)) <= v_ub;
v_lb <= diag(v652) <= v_ub;
v_lb <= diag(A * v671 * ctranspose(A)) <= v_ub;
v_lb <= diag(v684) <= v_ub;
v_lb <= diag(A * v680 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v692 * ctranspose(A)) <= v_ub;
v_lb <= diag(v611) <= v_ub;
v_lb <= diag(A * v675 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * vRG60 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v650 * ctranspose(A)) <= v_ub;
vSOURCEBUS == v0 * ctranspose(v0);

% (1): voltage across a line 
v645 == v632_abc([2, 3],[2, 3]) - S632645*ctranspose(Z632645) - Z632645*ctranspose(S632645) + Z632645*l632645*ctranspose(Z632645);
[v632_abc([2, 3],[2, 3]), S632645; ctranspose(S632645), l632645] >= 0;
diag(S632645-Z632645*l632645)- loads([19, 35]) + diag(v645 * Cbus([19, 35],[19, 35])) == diag(S645646) + 0;

v633 == v632 - S632633*ctranspose(Z632633) - Z632633*ctranspose(S632633) + Z632633*l632633*ctranspose(Z632633);
[v632, S632633; ctranspose(S632633), l632633] >= 0;
diag(A *(S632633-Z632633*l632633) * AH)- loads([10, 11, 12]) + diag(A * v633 * Cbus([10, 11, 12],[10, 11, 12]) * AH) == diag(A * S633634 * AH) + 0;

v634 == v633 - S633634*ctranspose(Z633634) - Z633634*ctranspose(S633634) + Z633634*l633634*ctranspose(Z633634);
[v633, S633634; ctranspose(S633634), l633634] >= 0;
diag(A *(S633634-Z633634*l633634) * AH)- loads([13, 14, 15]) + diag(A * v634 * Cbus([13, 14, 15],[13, 14, 15]) * AH) == 0;

v646 == v645 - S645646*ctranspose(Z645646) - Z645646*ctranspose(S645646) + Z645646*l645646*ctranspose(Z645646);
[v645, S645646; ctranspose(S645646), l645646] >= 0;
diag(S645646-Z645646*l645646)- loads([20, 21]) + diag(v646 * Cbus([20, 21],[20, 21])) == 0;

v632 == vRG60 - SRG60632*ctranspose(ZRG60632) - ZRG60632*ctranspose(SRG60632) + ZRG60632*lRG60632*ctranspose(ZRG60632);
[vRG60, SRG60632; ctranspose(SRG60632), lRG60632] >= 0;
diag(A *(SRG60632-ZRG60632*lRG60632) * AH)- loads([29, 30, 31]) + diag(A * v632 * Cbus([29, 30, 31],[29, 30, 31]) * AH) == [0; diag(S632645)] + diag(A * S632633 * AH) + diag(A * S632671 * AH) + 0;

v652 == v684([1],[1]) - S684652*ctranspose(Z684652) - Z684652*ctranspose(S684652) + Z684652*l684652*ctranspose(Z684652);
[v684([1],[1]), S684652; ctranspose(S684652), l684652] >= 0;
diag(S684652-Z684652*l684652)- loads([28]) + diag(v652 * Cbus([28],[28])) == 0;

v671 == v632 - S632671*ctranspose(Z632671) - Z632671*ctranspose(S632671) + Z632671*l632671*ctranspose(Z632671);
[v632, S632671; ctranspose(S632671), l632671] >= 0;
diag(A *(S632671-Z632671*l632671) * AH)- loads([16, 17, 18]) + diag(A * v671 * Cbus([16, 17, 18],[16, 17, 18]) * AH) == [S671684(1,1); 0; S671684(2,2) ] + diag(A * S671680 * AH) + diag(A * S671692 * AH) + 0;

v684 == v671_abc([1, 3],[1, 3]) - S671684*ctranspose(Z671684) - Z671684*ctranspose(S671684) + Z671684*l671684*ctranspose(Z671684);
[v671_abc([1, 3],[1, 3]), S671684; ctranspose(S671684), l671684] >= 0;
diag(S671684-Z671684*l671684)- loads([37, 38]) + diag(v684 * Cbus([37, 38],[37, 38])) == [diag(S684652); 0] + [0; diag(S684611)] + 0;

v680 == v671 - S671680*ctranspose(Z671680) - Z671680*ctranspose(S671680) + Z671680*l671680*ctranspose(Z671680);
[v671, S671680; ctranspose(S671680), l671680] >= 0;
diag(A *(S671680-Z671680*l671680) * AH)- loads([32, 33, 34]) + diag(A * v680 * Cbus([32, 33, 34],[32, 33, 34]) * AH) == 0;

v692 == v671 - S671692*ctranspose(Z671692) - Z671692*ctranspose(S671692) + Z671692*l671692*ctranspose(Z671692);
[v671, S671692; ctranspose(S671692), l671692] >= 0;
diag(A *(S671692-Z671692*l671692) * AH)- loads([23, 36, 22]) + diag(A * v692 * Cbus([23, 36, 22],[23, 36, 22]) * AH) == diag(A * S692675 * AH) + 0;

v611 == v684([2],[2]) - S684611*ctranspose(Z684611) - Z684611*ctranspose(S684611) + Z684611*l684611*ctranspose(Z684611);
[v684([2],[2]), S684611; ctranspose(S684611), l684611] >= 0;
diag(S684611-Z684611*l684611)- loads([27]) + diag(v611 * Cbus([27],[27])) == 0;

v675 == v692 - S692675*ctranspose(Z692675) - Z692675*ctranspose(S692675) + Z692675*l692675*ctranspose(Z692675);
[v692, S692675; ctranspose(S692675), l692675] >= 0;
diag(A *(S692675-Z692675*l692675) * AH)- loads([24, 25, 26]) + diag(A * v675 * Cbus([24, 25, 26],[24, 25, 26]) * AH) == 0;

A * vRG60 * AH == (A * v650([1, 2, 3],[1, 2, 3]) * AH) .* alphaMRG60;
[v650([1, 2, 3],[1, 2, 3]), S650RG60; ctranspose(S650RG60), l650RG60] >= 0;
diag(A *(S650RG60-Z650RG60*l650RG60) * AH)- loads([7, 8, 9]) + diag(A * vRG60 * Cbus([7, 8, 9],[7, 8, 9]) * AH) == diag(A * SRG60632 * AH) + 0;

v650 == vSOURCEBUS - SSOURCEBUS650*ctranspose(ZSOURCEBUS650) - ZSOURCEBUS650*ctranspose(SSOURCEBUS650) + ZSOURCEBUS650*lSOURCEBUS650*ctranspose(ZSOURCEBUS650);
[vSOURCEBUS, SSOURCEBUS650; ctranspose(SSOURCEBUS650), lSOURCEBUS650] >= 0;
diag(A *(SSOURCEBUS650-ZSOURCEBUS650*lSOURCEBUS650) * AH)- loads([4, 5, 6]) + diag(A * v650 * Cbus([4, 5, 6],[4, 5, 6]) * AH) == diag(A * S650RG60 * AH) + 0;

v671_abc == A * v671 * AH;
v632_abc == A * v632 * AH;


cvx_end


VSOURCEBUS = A * v0;
ISOURCEBUS650 = 1/trace(A * vSOURCEBUS([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * SSOURCEBUS650 * AH)*VSOURCEBUS([1, 2, 3]);
V650 = VSOURCEBUS([1, 2, 3]) - A * ZSOURCEBUS650* AH *ISOURCEBUS650;
I650RG60 = 1/trace(A * v650([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S650RG60 * AH)*V650([1, 2, 3]);
VRG60 = V650([1, 2, 3]) .* alphaRG60;
IRG60632 = 1/trace(A * vRG60([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * SRG60632 * AH)*VRG60([1, 2, 3]);
V632 = VRG60([1, 2, 3]) - A * ZRG60632* AH *IRG60632;
I632645 = 1/trace(v632_abc([2, 3],[2, 3]) ) * ctranspose(S632645)*V632([2, 3]);
V645 = V632([2, 3]) - Z632645*I632645;
I632633 = 1/trace(A * v632([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S632633 * AH)*V632([1, 2, 3]);
V633 = V632([1, 2, 3]) - A * Z632633* AH *I632633;
I632671 = 1/trace(A * v632([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S632671 * AH)*V632([1, 2, 3]);
V671 = V632([1, 2, 3]) - A * Z632671* AH *I632671;
I645646 = 1/trace(v645([1, 2],[1, 2]))*ctranspose(S645646)*V645([1, 2]);
V646 = V645([1, 2]) - Z645646*I645646;
I633634 = 1/trace(A * v633([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S633634 * AH)*V633([1, 2, 3]);
V634 = V633([1, 2, 3]) - A * Z633634* AH *I633634;
I671684 = 1/trace(v671_abc([1, 3],[1, 3]) ) * ctranspose(S671684)*V671([1, 3]);
V684 = V671([1, 3]) - Z671684*I671684;
I671680 = 1/trace(A * v671([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S671680 * AH)*V671([1, 2, 3]);
V680 = V671([1, 2, 3]) - A * Z671680* AH *I671680;
I671692 = 1/trace(A * v671([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S671692 * AH)*V671([1, 2, 3]);
V692 = V671([1, 2, 3]) - A * Z671692* AH *I671692;
I684652 = 1/trace(v684([1],[1]))*ctranspose(S684652)*V684([1]);
V652 = V684([1]) - Z684652*I684652;
I684611 = 1/trace(v684([2],[2]))*ctranspose(S684611)*V684([2]);
V611 = V684([2]) - Z684611*I684611;
I692675 = 1/trace(A * v692([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S692675 * AH)*V692([1, 2, 3]);
V675 = V692([1, 2, 3]) - A * Z692675* AH *I692675;


phasors=[];
phasors=[phasors;recover_voltage(VSOURCEBUS, 123)];
phasors=[phasors;recover_voltage(V650, 123)];
phasors=[phasors;recover_voltage(VRG60, 123)];
phasors=[phasors;recover_voltage(V632, 123)];
phasors=[phasors;recover_voltage(V645, 23)];
phasors=[phasors;recover_voltage(V633, 123)];
phasors=[phasors;recover_voltage(V671, 123)];
phasors=[phasors;recover_voltage(V646, 23)];
phasors=[phasors;recover_voltage(V634, 123)];
phasors=[phasors;recover_voltage(V684, 13)];
phasors=[phasors;recover_voltage(V680, 123)];
phasors=[phasors;recover_voltage(V692, 123)];
phasors=[phasors;recover_voltage(V652, 1)];
phasors=[phasors;recover_voltage(V611, 3)];
phasors=[phasors;recover_voltage(V675, 123)];

% change to per unit
phasors(:, 1) = phasors(:, 1) / Vbase;
phasors(:, 3) = phasors(:, 3) / Vbase;
phasors(:, 5) = phasors(:, 5) / Vbase;


Voltage_output=[];
Voltage_output = [Voltage_output; recover_voltage(V634, 123)];
Voltage_output = [Voltage_output; recover_voltage(V645, 23)];
Voltage_output = [Voltage_output; recover_voltage(V646, 23)];
Voltage_output = [Voltage_output; recover_voltage(V652, 1)];
Voltage_output = [Voltage_output; recover_voltage(V671, 123)];
Voltage_output = [Voltage_output; recover_voltage(V675, 123)];
Voltage_output = [Voltage_output; recover_voltage(V692, 123)];
Voltage_output = [Voltage_output; recover_voltage(V611, 3)];
Voltage_output = [Voltage_output; recover_voltage(V632, 123)];
Voltage_output = [Voltage_output; recover_voltage(V671, 123)];

% change to per unit
Voltage_output(:, 1) = Voltage_output(:, 1) / Vbase;
Voltage_output(:, 3) = Voltage_output(:, 3) / Vbase;
Voltage_output(:, 5) = Voltage_output(:, 5) / Vbase;

% output 
% disp(diag(S632645) / 1000);
% disp(diag(A * S632633 * AH) / 1000);
% disp(diag(A * S633634 * AH) / 1000);
% disp(diag(S645646) / 1000);
% disp(diag(A * SRG60632 * AH) / 1000);
% disp(diag(S684652) / 1000);
% disp(diag(A * S632671 * AH) / 1000);
% disp(diag(S671684) / 1000);
% disp(diag(A * S671680 * AH) / 1000);
% disp(diag(A * S671692 * AH) / 1000);
% disp(diag(S684611) / 1000);
% disp(diag(A * S692675 * AH) / 1000);
disp(diag(A * S650RG60 * AH) / 1000);
% disp(diag(A * SSOURCEBUS650 * AH) / 1000);

% disp(diag(A * S632671 * AH) / 1000);
% disp(loads([23, 36, 22]) / 1000);
% disp(loads([20, 21]) / 1000);

% disp(S632633 / 1000);
% disp(S633634 / 1000);
% disp(S632671 / 1000);
% disp(S671692 / 1000);
% 