clear
close all
clc
format compact
disp(' ')
% Parameters
rng('default')
Var_Noise = input('Enter Variance of Noise : '); % Variance of Noise
Sig_d = input('Enter Power of Desired Signal : '); % Sigma^2_d
M = input('Enter Length of Filter : '); % Length of Filter
L = input('Enter Length of Channel : '); % Length of Channel
Delay = input('Enter The Value of Delay : '); % Value of Delay
N = input('Enter Time Length(n) : '); % Descrete Time Length
Rep = input('Enter The Number of Iterration : '); % Number of Iterration
B = input('Enter The Value of b(in Vector) : '); % Channel Coefficient
Beta = input('Enter The Value of Beta : '); % Value of Beta
Mu = input('Enter The Value of Mu : '); % Value of Mu
Lambda = input('Enter The Value of Lambda : '); % Value of Lambda
L_Ch = 1:L;
Sig_RLS = 0.004;
H_Tot = zeros(length(B),L); % Total Channel Response
J_min = zeros(1,length(B)); % Minimum MSE
P_Tot = zeros(M,length(B)); % Total Cross Correlation
J_Tot_LMS = zeros(N,length(B)); % MSE
J_Tot_DCT = zeros(N,length(B)); % MSE
J_Tot_RLS = zeros(N,length(B)); % MSE
J_LMS = zeros(N,1); % MSE
J_DCT = zeros(N,1); % MSE
J_RLS = zeros(N,1); % MSE
W_LMS = zeros(M,N); % Filter Weights
W_DCT = zeros(M,N); % Filter Weights
W_RLS = zeros(M,N); % Filter Weights
for z=1:length(B)
disp(['For b = ',num2str(B(z))])
disp(' ')
H = 0.5*(1 + cos((2*pi/B(z)).*(L_Ch - 2))); % Channel Response
H_Tot(z,:) = H;
r_h = conv(H,fliplr(H)); % Convlolution of Ch. Res.
r_x = [r_h(L:(2*L) - 1),zeros(1,M - L)];
R_x = toeplitz(r_x); % Correlation Matrix of Signal
R_v = Var_Noise*eye(M); % Correlation Matrix of Noise
% U(n) = X(n) + V(n) , X(n) = S(n)*H(n)
R_u = R_x + R_v; % Correlation Matrix of Input
EVD = eig(R_u); % EigenValue Decomposition
X_R = (max(EVD))/(min(EVD)); % EigenValue Spread
P = zeros(M,1); % Correlation of Input & Desired
for i = 1:M
if (Delay - i + 1 == 3) || (Delay - i + 1 == 2) || (Delay - i + 1 == 1)
P(i,1) = H(Delay - i + 1);
else
P(i,1) = 0;
end
end
P_Tot(:,z) = P;
W_opt = R_u\P; % Optimum Wiener Filter
J_min(z) = Sig_d - P.'*W_opt; % MSE
disp(['Jmin Equals With : ',num2str(J_min(z))])
disp(' ')
S = zeros(N,1); % Signal
T = zeros(M,M); % Transform Matrix
% Updating Transform Matrix
for k=1:M
for l=1:M
if k==1 % For k=0
T(k,l) = 1/sqrt(M);
else
T(k,l) = sqrt(2/M)*cos((pi*(k-1)*(2*(l-1)+1))/(2*M));
end
end
end
R_t = T*R_u*T'; % New Correlation Matrix
EVD_t = eig(R_t); % EigenValue Decomposition after Transform
X_R_T = (max(EVD_t))/(min(EVD_t)); % EigenValue Spread after Transform
disp(['EigenValue Spread Before Transform Equals With : ',num2str(X_R),...
' And After Transform Equals With : ',num2str(X_R_T)])
disp(' ')
for m=1:Rep
S = 2*randi([0 1],N,1)-1; % Signal
X = conv(S,H); % Output of Channel
X = X(1:N); % Set N samples
V = sqrt(Var_Noise)*randn(N,1); % Noise
U = X + V; % Input of Filter
D_D(Delay + 1:N,1) = S(1:N - Delay,1); % Desired Signal
% AUX Parameters
U_Aux_LMS = zeros(M,1); % LMS
W_AUX_LMS = zeros(M,N); % LMS
W_Old_LMS = zeros(M,1); % LMS
U_Aux_DCT = zeros(M,1); % DCT-LMS
W_AUX_DCT = zeros(M,N); % DCT-LMS
W_Old_DCT = zeros(M,1); % DCT-LMS
Sig_k = ones(1,M); % DCT-LMS
U_Aux_RLS = zeros(M,1); % RLS
W_Aux_RLS = zeros(M,N); % RLS
W_Old_RLS = zeros(M,1); % RLS
K_n = zeros(M,1); % RLS
Zeta = zeros(1,N); % RLS
P_RLS = (1/Sig_RLS)*eye(M); % RLS
% Calculating Error
for j=1:N
% ~> % For " LMS "
% e(n) = d(n) - W' * u(n)
% W(n+1) = W(n) + Mu(e*(n).u(n))
Err_LMS(j,1) = D_D(j,1) - W_Old_LMS'*U_Aux_LMS; % Err. Equivalent
W_AUX_LMS(:,j) = W_Old_LMS + Mu*Err_LMS(j,1)*U_Aux_LMS;
U_Aux_LMS = circshift(U_Aux_LMS,1); % Shifting
U_Aux_LMS(1,1) = U(j,1);
W_Old_LMS = W_AUX_LMS(:,j); % Save for Next Loop
% ~> % For " DCT - LMS "
% e(n) = d(n) - W' * x(n)
% W(n+1) = W(n) + Mu*D*(e*(n).x(n))
U_Aux_DCT = circshift(U_Aux_DCT,1); % Shifting
U_Aux_DCT(1,1) = U(j,1); % Updating
X_N = T*U_Aux_DCT; % Output of Transform Matrix
% Calculating Sigma^2_k
for q=1:M
Sig_k(q) = Beta*Sig_k(q) + (1 - Beta)*abs(X_N(q))^2;
end
D = diag(1./Sig_k); % Diagonal Matrix
Err_DCT(j,1) = D_D(j,1) - W_Old_DCT'*X_N; % Err. Equivalent
W_AUX_DCT(:,j) = W_Old_DCT + Mu*D*Err_DCT(j,1)*X_N; % Updating Coeeficient
W_Old_DCT = W_AUX_DCT(:,j); % Save for Next Loop
% ~> % For " RLS "
% W(n) = W(n-1) + K(n).Zeta*(n)
% K(n) = (inv(Lambda).P(n-1).u(n)) / (1 + inv(Lambda).u(n)'.P(n-1).u(n))
% Zeta(n) = d(n) - W'(n-1).u(n)
U_Aux_RLS = circshift(U_Aux_RLS,1); % Shifting
U_Aux_RLS(1,1) = U(j,1); % Updating
K_n = (Lambda^(-1)*P_RLS*U_Aux_RLS)/(1 + Lambda^(-1)*U_Aux_RLS'*P_RLS*U_Aux_RLS);
Zeta(1,j) = D_D(j,1) - W_Old_RLS'*U_Aux_RLS; % Previous Error
W_Aux_RLS(:,j) = W_Old_RLS + K_n*conj(Zeta(1,j)); % W(n+1)
Err_RLS(j,1) = D_D(j,1) - W_Old_RLS'*U_Aux_RLS; % Error
P_RLS = Lambda^(-1)*P_RLS - Lambda^(-1)*K_n*U_Aux_RLS'*P_RLS;
W_Old_RLS = W_Aux_RLS(:,j); % Save for Next Loop
end
W_LMS = W_LMS + W_AUX_LMS; % Updating, LMS
J_LMS = J_LMS + Err_LMS.^2; % Updating, LMS
W_DCT = W_DCT + W_AUX_DCT; % Updating, DCT-LMS
J_DCT = J_DCT + Err_DCT.^2; % Updating, DCT-LMS
W_RLS = W_RLS + W_Aux_RLS; % Updating, RLS
J_RLS = J_RLS + Err_RLS.^2; % Updating, RLS
end
% For LMS
W_LMS = W_LMS/Rep; % Final Weights , LMS
J_LMS = J_LMS/Rep; % Final MSE, LMS
J_Tot_LMS(:,z) = J_LMS; % Save MSE for Plot, LMS
% For DCT-LMS
W_DCT = W_DCT/Rep; % Final Weights, DCT-LMS
J_DCT = J_DCT/Rep; % Final MSE, DCT-LMS
J_Tot_DCT(:,z) = J_DCT; % Save MSE for Plot, DCT-LMS
% For RLS
W_RLS = W_RLS/Rep; % Final Weights, RLS
J_RLS = J_RLS/Rep; % Final MSE, RLS
J_Tot_RLS(:,z) = J_RLS; % Save MSE for Plot, RLS
Jexc_LMS(z) = J_LMS(end) - J_min(z); % Jexcess LMS
Jexc_DCT(z) = J_DCT(end) - J_min(z); % Jexcess DCT-LMS
Jexc_RLS(z) = J_RLS(end) - J_min(z); % Jexcess RLS
MisAdj_LMS(z) = Jexc_LMS(z)/J_min(z); % MissAdjustment LMS
MisAdj_DCT(z) = Jexc_DCT(z)/J_min(z); % MissAdjustment DCT-LMS
MisAdj_RLS(z) = Jexc_RLS(z)/J_min(z); % MissAdjustment RLS
disp(['MissAdjustment for LMS Equals With : ',num2str(MisAdj_LMS(z))])
disp(['MissAdjustment for DCT-LMS Equals With : ',num2str(MisAdj_DCT(z))])
disp(['MissAdjustment for RLS Equals With : ',num2str(MisAdj_RLS(z))])
disp(' ')
end
%% Plot
figure('name','LMS & DCT-LMS & RLS')
for z=1:length(B)
subplot(2,2,z)
semilogy(J_Tot_LMS(:,z),'linewidth',1.5) % LMS
hold on
semilogy(J_Tot_DCT(:,z),'linewidth',1.5) % DCT-LMS
hold on
semilogy(J_Tot_RLS(:,z),'linewidth',1.5) % RLS
hold on
semilogy(J_min(z)*ones(N,1),'linewidth',2); % Jmin
grid on
title(['Learning Curve for \mu = ',num2str(Mu),' and b = ', num2str(B(z))],'color','b','fontsize',13)
xlabel('$Time(n)$','interpreter','latex','fontsize',13);
ylabel('$J(n)$','interpreter','latex','fontsize',13);
legend('$J(n)-LMS$','$J(n) - DCT-LMS$','$J(n)-RLS$','$J_{min}$','interpreter','latex','fontsize',13)
end
figure('name','Learning Curves for RLS')
Color = ['r','b','g','m'];
for z=1:length(B)
semilogy(J_Tot_RLS(:,z),'linewidth',1.5,'color',Color(z)) % RLS
hold on
semilogy(J_min(z)*ones(N,1),'--','linewidth',2,'color',Color(z)); % Jmin
grid on
end
title('Learning Curve','color','b','fontsize',13)
xlabel('$Time(n)$','interpreter','latex','fontsize',13);
ylabel('$J(n)$','interpreter','latex','fontsize',13);
legend('b = 2.9','Jmin, b = 2.9','b = 3.1','Jmin, b = 3.1',...
'b = 3.3','Jmin, b = 3.3','b = 3.5','Jmin, b = 3.5')
format loose