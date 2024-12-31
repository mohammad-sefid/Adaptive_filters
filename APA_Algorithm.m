% Parameters
% ----------
% M = input('Enter The Number of Filter,s Tab : '); % Tap Filter
% D = input('Enter The Value of Delay : '); % Delay Value : C
% N_It = input('Enter The Number Process Itteration : '); % Num of Itteration
M = 35;
D = 15;
N = [1 2 3 7];
Md = M + 5; % Max of Delay Time
N_It = 500;
Time = 1000;
J = sqrt(-1);
H = [0.5 1.2 1.5 -1]; % Channel Response
L = length(H); % Length of Channel
Sig_S = 2; % Power of Signal
Sig_V = 0.01; % Power of Noise
r_h = conv(H,fliplr(H)); % AutoCorrelation of Channel
r_x = Sig_S * r_h; % AutoCorrelation of Signal
r_xp = [fliplr(r_x(1:L)),zeros(1,M - L)]; % AUX Parameter
Rx = toeplitz(r_xp); % AutoCorrelation Matrix of Signal
Rv = Sig_V * eye(M); % AutoCorrelation Matrix of Noise
Ru = Rx + Rv; % AutoCorrelation Matrix of Equ. Input
Sig_d = Sig_S; % Power of Desired Signal
P = zeros(M,Md); % AutoCorrelation of Input & Desired
J_min = zeros(1,Md); % Minimum Mean Square Error(MSE)
W_opt = zeros(M,Md); % Optimum Wiener Filter
for d=1:Md
% Generating P ( Correlation of Input & Desired )
for i=1:M
if d-i == 0 || d-i == 1 || d-i == 2 || d-i == 3
P(i,d) = 2*H(1,d-i+1);
end
end
W_opt(:,d) = inv(Ru) * P(:,d); % Optimim Wiener Filter
J_min(d) = Sig_d - P(:,d)'*W_opt(:,d); % MSE
end
% Finding Minimum of MSE and its Delay Time :
[Min, Index] = min(J_min);
disp(['Delay for Minimum MSE is :',num2str(Index)])
%% Part A : APA Algorithm.
disp('_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_')
disp('Part A :')
disp(' ')
Mu = 0.4; % Step Size Parameter
J_n_Tot = zeros(length(N),Time); % Total MSE
for i=1:length(N)
disp('----------------------------------------')
disp(['Part A : For N = ',num2str(N(i))])
disp(' ')
J_n = zeros(1,Time); % MSE
W_C = zeros(M,1); % Filter Coefficient
S_C = zeros(1,Time); % Main Signal
Y_C = zeros(1,Time); % Output of Filter
for j = 1:N_It % QAM Constellation
Re_C = 2*randi([0 1],1,Time) - 1; % Real Part
Im_C = 2*randi([0 1],1,Time) - 1; % Imaginary Part
S_C = Re_C + J*Im_C; % Main Signal
D_C = [zeros(1,D),S_C(1:end)]; % Desired Signal
X_C = conv(S_C,H); % Output of Channel
X_C = X_C(1:Time); % Get "Time" Elements
V_C = sqrt(Sig_V/2) * (randn(1,Time) + J*randn(1,Time)); % Noise
U_C = X_C + V_C; % Input of Filter
% Affine Projection Algorithm (APA)
W_old = zeros(M,1); % Initial Weights
U_Aux = zeros(M,1); % Seperating M Elements of U
D_n = zeros(N(i),1); % Aux Parameter
An_C = zeros(N(i),M); % Conditions
W_Aux = zeros(M,Time); % Aux Parameter
E_C = zeros(1,Time); % Error
for l=1:Time
% W'(n+1).u(n) = d(n)
% Dimensions --> W : M*1 , u : M*1 , d : 1*1
% A'(n) = [u(n) u(n-1) u(n-2) ... u(n-N+1)]
U_Aux = circshift(U_Aux,1); % Updating u(n)
U_Aux(1) = U_C(l);
An_C = circshift(An_C,1); % Updating A(n)
An_C(1,:) = U_Aux';
D_n = circshift(D_n,1); % Updating D(n)
D_n(1) = conj(D_C(l));
if l >= M
% W(n+1) = W(n) + Mu*A'(n)*inv[(A(n).A'(n)]*[d(n) - A(n).W(n)]
W_Aux(:,l) = W_old + Mu * An_C'*inv(An_C*An_C')*(D_n - An_C*W_old);
W_old = W_Aux(:,l); % W(n)
Y_Aux(1,l) = W_Aux(:,l-1)'*U_Aux; % y(n) = W'(n).u(n)
% e(n) = d(n) - y(n)
E_C(1,l) = D_C(1,l) - Y_Aux(1,l); % Error , e(n)
end
end
Y_C = Y_C + Y_Aux; % Final Output
J_n = J_n + abs(E_C).^2; % Final MSE for i'th Itteration
W_C = W_C + W_Aux; % Final Weights
end
J_n_Tot(i,:) = J_n/N_It; % Final MSE
% MissAdjustment = J_exc/J_min , J_exc = J(n) - J_min
J_exc_C(i) = J_n_Tot(i,end) - J_min(D); % J_excces
MisAdj_C(i) = J_exc_C(i)/J_min(D); % MissAdjustment
Exc_C = ['J_exc is Equal with : ',num2str(J_exc_C(i))];
Adj_C = ['MissAdjustment is Equal with : ',num2str(MisAdj_C(i))];
% Results in dB
% J_exc_dB = 10*log10(J_exc(i));
% MisAdj_dB = 10*log10(MisAdj(i));
% Exc = ['J_exc is Equal with(in dB) : ',num2str(J_exc_dB)];
% Adj = ['MissAdjustment is Equal with(in dB) : ',num2str(MisAdj_dB)];
disp(Exc_C)
disp(' ')
disp(Adj_C)
disp(' ')
end
% Plot Part
figure('name','Part A : Learning Curve')
for i=1:4
semilogy(J_n_Tot(i,:),'linewidth',1)
hold on
grid on
end
semilogy(J_min(D)*ones(1,Time),'linewidth',1)
title('Learning Curve','fontsize',13,'interpreter','latex','color','b')
xlabel('Discrete Time (n)','fontsize',13,'interpreter','latex')
ylabel('J(n)','fontsize',13,'interpreter','latex')
legend('N = 1','N = 2','N = 3','N = 7','$J_{min}$','interpreter','latex','fontsize',13)
%% Part B : APA Algorithm,Scattering Plot
disp('_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_')
disp('Part B :')
disp(' ')
N_D = [1 7]; % Condition Numbers
for j=1:length(N_D)
disp('----------------------------------------')
disp(['Part B : For N = ',num2str(N_D(j))])
disp(' ')
Re_D = 2*randi([0 1],1,Time) - 1; % Real Part of Signal
Im_D = 2*randi([0 1],1,Time) - 1; % Imaginary Part of Signal
S_D = Re_D + J*Im_D; % Main Input Signal
S_D_Tot(j,:) = S_D; % Store Signal
D_D = [zeros(1,D),S_D(1,1:end)]; % Desired Signal
X_D = conv(S_D,H); % Output of Channel
X_D = X_D(1:Time); % Get "Time" Elements
V_D = sqrt(Sig_V/2) * (randn(1,Time) + J*randn(1,Time)); % Noise
U_D = X_D + V_D; % Input of Filter
U_D_Tot(j,:) = U_D; % Store Input
% Affine Projection Algorithm (APA)
W_old_D = zeros(M,1); % Initial Weights
U_Aux_D = zeros(M,1); % Seperating M Elements of U
Dn_D = zeros(N_D(j),1); % Aux Parameter
An_D = zeros(N_D(j),M); % A(n)
W_Aux_D = zeros(M,Time); % Aux Parameter
for l=1:Time
% W'(n+1).u(n) = d(n)
% Dimensions --> W : M*1 , u : M*1 , d : 1*1
% A'(n) = [u(n) u(n-1) u(n-2) ... u(n-N+1)]
U_Aux_D = circshift(U_Aux_D,1); % Updating u(n)
U_Aux_D(1) = U_C(1,l);
An_D = circshift(An_D,1); % Updating A(n)
An_D(1,:) = U_Aux_D';
Dn_D = circshift(Dn_D,1); % Updating d(n)
Dn_D(1) = conj(D_C(l));
if l >= M
% W(n+1) = W(n) + Mu*A'(n)*inv[(A(n).A'(n)]*[d(n) - A(n).W(n)]
W_Aux_D(:,l) = W_old_D + Mu * An_D'*inv(An_D*An_D')*(Dn_D - An_D*W_old_D);
W_old_D = W_Aux_D(:,l); % W(n)
Y_D(j,l) = W_Aux_D(:,l-1)'*U_Aux_D; % y(n) = W'(n).u(n)
end
end
end
% Plot
figure('name','Part B : Scattering Diagram ')
subplot(2,2,1)
plot(real(U_D_Tot(1,:)),imag(U_D_Tot(1,:)),'.')
title('For N = 1 : Input Scattering','color','b')
xlabel('Real')
ylabel('Imaginary')
subplot (2,2,2)
plot(real(Y_D(1,:)),imag(Y_D(1,:)),'.r')
title('For N = 1 : Output Scattering','color','b')
xlabel('Real')
ylabel('Imaginary')
subplot (2,2,3)
plot(real(U_D_Tot(2,:)),imag(U_D_Tot(2,:)),'.')
title('For N = 7 : Input Scattering','color','b')
xlabel('Real')
ylabel('Imaginary')
subplot (2,2,4)
plot(real(Y_D(2,:)),imag(Y_D(2,:)),'.r')
title('For N = 7 : Output Scattering','color','b')
xlabel('Real')
ylabel('Imaginary')
%% Part C : APA Algorithm,Decision Directed
disp('_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_')
disp('Part C :')
disp(' ')
N_Tr = 200; % Training Set
N_E = N_Tr + 800; % Validation Set
W_E = zeros(M,1);
S_E_Tot = zeros(length(N),N_E); % Main Signal
Y_E = zeros(length(N),N_E); % Output of Filter
Jn_E_Tot = zeros(length(N),N_E); % MSE
for z=1:length(N)
disp('----------------------------------------')
disp(['Part C : For N = ',num2str(N(z))])
disp(' ')
J_n_E = zeros(1,N_E);
S_E = zeros(1,N_E);
for i = 1:N_It
Re_E = 2*randi([0 1],1,N_E) - 1; % Reai Part
Im_E = 2*randi([0 1],1,N_E) - 1; % Imaginary Part
S_E = Re_E + J*Im_E; % Main Signal
S_E_Tot(z,:) = S_E;
D_E = [zeros(1,D),S_E(1,1:end)]; % Desired Signal
X_E = conv(S_E,H); % Output of Channel
X_E = X_E(1:N_E); % Select N Samples
V_E = sqrt(Sig_V/2)*(randn(1,N_E) + J*randn(1,N_E)); % Noise
U_E = X_E + V_E; % Input of Filter
%%% Affine Projection Algorithm
W_old_E = zeros(M,1); % Initial Weights
U_Aux_E = zeros(M,1); % Input for Updating
Dn_E = zeros(N(z),1); % Desired Signal for Updating
An_E = zeros(N(z),M); % Conditions
W_Aux_E = zeros(M,N_E); % Aux Parameter
E_E = zeros(1,N_E); % Error
% QAM Constellation
S_e(1) = +1 + 1j;
S_e(2) = -1 + 1j;
S_e(3) = +1 - 1j;
S_e(4) = -1 - 1j;
for j=2:N_E
U_Aux_E = circshift(U_Aux_E,1);
U_Aux_E(1,1) = U_E(1,j);
An_E = circshift(An_E,1);
An_E(1,:) = U_Aux_E';
Dn_E = circshift(Dn_E,1);
if j<=200 % After Training Set
% e*(n) = d*(n) - W'(n).u(n)
Dn_E(1,1) = conj(D_E(1,j));
Y_Aux_E(1,j) = W_Aux_E(:,j-1)'*U_Aux_E; % y(n) = W'(n).u(n)
E_E(1,j) = D_E(1,j) - Y_Aux_E(1,j); % d*(n) - y(n)
else
Y_Aux_E(1,j) = W_Aux_E(:,j-1)'*U_Aux_E; % y(n) = W'(n).u(n)
for k=1:4
if k ==1
% Distance = |S(n) - y(n)|
Dis_E(k) = abs(S_e(k) - Y_Aux_E(1,j)) ;
Min_Dis_E = Dis_E(k); % Minimum Distance
S_new_E(j) = S_e(k); % Estimated Signal
else
Dis_E(k) = abs(S_e(k) - Y_Aux_E(1,j));
if Min_Dis_E > Dis_E(k)
Min_Dis_E = Dis_E(k);
S_new_E(j) = S_e(k); % Estimated Signal
end
end
Dn_E(1,1) = conj(S_new_E(j));
E_E(1,j) = D_E(1,j) - Y_Aux_E(1,j); % e(n) = d(n) - y(n)
end
end
if j >= M
% W(n+1) = W(n) + Mu*A'(n)*inv[(A(n).A'(n)]*[d(n) - A(n).W(n)]
W_Aux_E(:,j) = W_old_E + Mu * An_E'*inv(An_E*An_E')*(Dn_E - An_E*W_old_E);
W_old_E = W_Aux_E(:,j);
end
end
Y_E = Y_E + Y_Aux_E; % Final Output
J_n_E = J_n_E + abs(E_E).^2;
W_E = W_E + W_Aux_E; % Filnal Weights
end
Jn_E_Tot(z,:) = J_n_E/N_It; % Final MSE
% MissAdjustment = J_exc/J_min , J_exc = J(n) - J_min
J_exc_E(z) = Jn_E_Tot(z,end) - J_min(D); % J_excces
MisAdj_E(z) = J_exc_E(z)/J_min(D); % MissAdjustment
Exc_E = ['J_exc is Equal with : ',num2str(J_exc_E(z))];
Adj_E = ['MissAdjustment is Equal with : ',num2str(MisAdj_E(z))];
% Results in dB
% J_exc_dB = 10*log10(J_exc(i));
% MisAdj_dB = 10*log10(MisAdj(i));
% Exc = ['J_exc is Equal with(in dB) : ',num2str(J_exc_dB)];
% Adj = ['MissAdjustment is Equal with(in dB) : ',num2str(MisAdj_dB)];
disp(Exc_E)
disp(' ')
disp(Adj_E)
disp(' ')
end
%% Part C : APA Algorithm,Scattering Plot , Decision Directed
disp('_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_')
disp('Part C :')
disp(' ')
N_F = [1 7];
U_F_Tot = zeros(length(N_F),Time);
Y_F = zeros(length(N_F),Time);
%%
for i=1:length(N_F)
disp('----------------------------------------')
disp(['Part C : For N = ',num2str(N_F(i))])
disp(' ')
S_F = zeros(1,Time);
Re_F = 2*randi([0 1],1,Time) - 1; % Real Part
Im_F = 2*randi([0 1],1,Time) - 1; % Imaginary Part
S_F = Re_F + J*Im_F; % Main Signal
S_E_Tot(z,:) = S_F;
D_F = [zeros(1,D),S_F(1,1:end)]; % Desired Signal
X_F = conv(S_F,H);
X_F = X_F(1:Time); % Output of Channel
V_F = sqrt(Sig_V/2)*(randn(1,Time)+1j*randn(1,Time)); % Noise
U_F = X_F + V_F; % Input of Filter
U_F_Tot(i,:) = U_F;
% APA
W_old_F = zeros(M,1); % Initial Weights
U_Aux_F = zeros(M,1); % Aux Parameter, u(n)
Dn_F = zeros(N_F(i),1); % Aux Parameter, d(n)
An_F = zeros(N_F(i),M); % Conditions
W_Aux_F = zeros(M,Time); % Aux Parameter, W(n)
E_F = zeros(1,Time); % Error
% QAM Constellation
S_f(1) = +1 + 1j;
S_f(2) = -1 + 1j;
S_f(3) = +1 - 1j;
S_f(4) = -1 - 1j;
for j=2:Time
% Update & Shifting Process
U_Aux_F = circshift(U_Aux_F,1);
U_Aux_F(1) = U_F(1,j);
An_F = circshift(An_F,1);
An_F(1,:) = U_Aux_F';
Dn_F = circshift(Dn_F,1);
if j<=200
Dn_F(1) = conj(D_F(1,j));
Y_F(i,j) = W_Aux_F(:,j-1)'*U_Aux_F;
E_F(j) = D_F(1,j) - Y_F(i,j);
else
Y_F(i,j) = W_Aux_F(:,j-1)'*U_Aux_F;
for k=1:4
if k == 1
% Distance : |S(n) - y(n)|
Dis_F(k) = abs(S_f(k)-Y_F(i,j)) ;
Min_Dis_E = Dis_F(k); % Minimum Distance
S_new_F(j) = S_f(k); % Estimated Signal
else
Dis_F(k) = abs(S_f(k)-Y_F(i,j));
if Min_Dis_E > Dis_F(k)
Min_Dis_E = Dis_F(k); % Minimum Distance
S_new_F(j) = S_f(k); % Estimated Signal
end
end
Dn_F(1) = conj(S_new_F(j));
E_F(1,j) = D_F(1,j) - Y_F(i,j); % e(n) = d(n) - y(n)
end
end
if j >= M
% W(n+1) = W(n) + Mu*A'(n)*inv[(A(n).A'(n)]*[d(n) - A(n).W(n)]
W_Aux_F(:,j) = W_old_F + Mu * An_F'*inv(An_F*An_F')*(Dn_F - An_F*W_old_F);
W_old_F = W_Aux_F(:,j);
end
end
end
% Plot
figure('name','Part C : Scattering Diagram ')
subplot(2,2,1)
plot(real(U_F_Tot(1,:)),imag(U_F_Tot(1,:)),'.')
title('For N = 1 : Input Scattering','color','b')
xlabel('Real')
ylabel('Imaginary')
subplot (2,2,2)
plot(real(Y_F(1,:)),imag(Y_F(1,:)),'.r')
title('For N = 1 : Output Scattering','color','b')
xlabel('Real')
ylabel('Imaginary')
subplot (2,2,3)
plot(real(U_F_Tot(2,:)),imag(U_F_Tot(2,:)),'.')
title('For N = 7 : Input Scattering','color','b')
xlabel('Real')
ylabel('Imaginary')
subplot (2,2,4)
plot(real(Y_F(2,:)),imag(Y_F(2,:)),'.r')
title('For N = 7 : Output Scattering','color','b')
xlabel('Real')
ylabel('Imaginary')
%% Comparison
figure('name','Comparison')
Color = ['m','g','b','k'];
for i=1:length(N)
semilogy(J_n_Tot(i,:),'--','linewidth',1,'color',Color(i))
hold on
grid on
end
for i=1:4
semilogy(Jn_E_Tot(i,:),':','linewidth',1,'color',Color(i))
hold on
grid on
end
semilogy(J_min(D)*ones(1,N_E),'linewidth',2)
title('Learning Curve Comparison','fontsize',13,'interpreter','latex','color','b')
xlabel('Discrete Time (n)','fontsize',13,'interpreter','latex')
ylabel('J(n)','fontsize',13,'interpreter','latex')
Legend = {'Part A : N = 1','Part A : N = 2','Part A : N = 3','Part C : N = 7','Part C : N = 1',...
'Part C : N = 2','Part C : N = 3','Part C : N = 7','$J_{min}$'};
legend(Legend,'interpreter','latex','fontsize',13)
format loose