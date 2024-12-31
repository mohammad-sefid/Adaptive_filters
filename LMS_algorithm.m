% Parameters
% ----------
rng('default')
Var_Noise = input('Enter Variance of Noise : '); % Variance of Noise
Sig_d = input('Enter Power of Desired Signal : '); % Sigma^2_d
M = input('Enter Length of Filter : '); % Length of Filter
L = input('Enter Length of Channel : '); % Length of Channel
Delay = input('Enter The Value of Delay : '); % Value of Delay
N = input('Enter Time Length(n) : '); % Descrete Time Length
Rep = input('Enter The Number of Iterration : '); % Number of Iterration
B = input('Enter The Value of b(in Vector) : '); % Channel Coefficient
L_Ch = 1:L;
% Main :
% ------
%% Part A : Correlation Matrix & EigenValue Spread.
disp('_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_')
disp('Parts A, B and C :')
disp(' ')
J_n_C = zeros(length(B),N); % MSE for Part C
H_Tot = zeros(length(B),L); % Total Channel Response
J_min = zeros(1,length(B)); % Minimum MSE
for z=1:length(B)
disp('--------------------------------------------------')
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
disp(['EigenValue Spread Equals With = ',num2str(X_R)])
disp(' ')
% Part B
P = zeros(M,1); % Correlation of Input & Desired
for i = 1:M
if (Delay - i + 1 == 3) || (Delay - i + 1 == 2) || (Delay - i + 1 == 1)
P(i,1) = H(Delay - i + 1);
else
P(i,1) = 0;
end
end
W_opt = inv(R_u)*P; % Optimum Wiener Filter
J_min(z) = Sig_d - P.'*W_opt; % MSE
disp(['Jmin Equals With : ',num2str(J_min)])
disp(' ')
% Part C
% Updating Process
Mu_C = 1/(max(EVD) + min(EVD)); % Optimum Value of Step-Size
W = zeros(M,N);
for k=2:N
W(:,k) = W(:,k-1) + Mu_C * (P - R_u * W(:,k-1));
end
Error_C = zeros(N,Rep);
for i=1:N
for j=1:Rep
S = 2*randi([0 1],N,1)-1; % Input Signal
X = conv(S,H); % Output of Channel
X = X(1:N); % Same Size of Signal & Noise
V = sqrt(Var_Noise)*randn(N,1); % Noise
U = X + V; % Input of Filter
D = [zeros(Delay,1);S(1:N-Delay)]; % Desired Signal
U_C = zeros(M,1); % AUX Parameter
for k=1:N
Error_C(k,j) = D(k,1) - W(:,i)'*U_C; % Error
U_C = circshift(U_C,1); % Shift for Convolution
U_C(1) = U(k); % Updating Shift
end
end
% J(n) = (1/N) Sum(|e(n)|^2)
J_n_C(z,i) = sum(sum(Error_C.^2,2)/Rep)/N; % MSE
end
W = W(:,end);
W_new = [W_opt,W];
disp(W_new)
end
%% Plot
figure('name','Parts A, B, and C')
for z=1:length(B)
subplot(2,2,z)
semilogy(J_n_C(z,:),'linewidth',1.5)
grid on
hold on
semilogy(J_min(z)*ones(N,1),'linewidth',2);
title(['Learning Curve for b = ', num2str(B(z))],'color','b','fontsize',13)
xlabel('$Time(n)$','interpreter','latex','fontsize',13);
ylabel('$J(n)$','interpreter','latex','fontsize',13);
legend('$J(n)$','$J_{min}$','interpreter','latex','fontsize',13)
end
figure('name','Channel Response H(n)')
for z=1:length(B)
subplot(2,2,z)
stem(H_Tot(z,:),'filled','linewidth',2)
grid on
title(['Channel Response for b = ', num2str(B(z))],'color','k','fontsize',13)
xlabel('$Time(n)$','interpreter','latex','fontsize',13,'color','r');
ylabel('$H(n)$','interpreter','latex','fontsize',13,'color','r');
xlim([0 L+1])
end
%% Part D
disp('_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_')
disp('Part D :')
disp(' ')
Mu_D = 0.075; % Step Size Value
S_D = zeros(N,1); % Signal
E_D = zeros(N,Rep); % Error
J_n_D = zeros(N,length(B)); % MSE
W_D = zeros(M,N); %Filter Weights
for y=1:length(B)
J_D = zeros(N,1);
for m=1:Rep
H_N = H_Tot(y,:); % Channel Response
S_D = 2*randi([0 1],N,1)-1; % Signal
X_D = conv(S_D,H_N); % Output of Channel
X_D = X_D(1:N); % Set N samples
V_D = sqrt(Var_Noise)*randn(N,1); % Noise
U_D = X_D + V_D; % Input of Filter
D_D(8:N,1) = S_D(1:N-7,1); % Desired Signal
% AUX Parameters
U_AD = zeros(M,1);
W_AUX = zeros(M,N);
W_old = zeros(M,1);
% Calculating Error
for j=1:N
% e(n) = d(n) - W' * u(n)
% W(n+1) = W(n) + Mu*(e(n).u(n))
Error_D(j,1) = D_D(j,1) - W_old'*U_AD; % Err. Equivalent
W_AUX(:,j) = W_old + Mu_D*Error_D(j,1)*U_AD;
U_AD = circshift(U_AD,1); % Shifting
U_AD(1,1) = U_D(j,1);
W_old = W_AUX(:,j); % Save for Next Loop
end
W_D = W_D + W_AUX; % Updating
J_D = J_D + Error_D.^2; % Updating
end
W_D = W_D/Rep; % Final Weights
J_D = J_D/Rep; % Final MSE
J_n_D(:,y) = J_D; % Save MSE for Plot
end
%% Plot
figure('name','Part D')
for z=1:length(B)
subplot(2,2,z)
semilogy(J_n_D(:,z),'linewidth',1.5)
grid on
hold on
semilogy(J_min(z)*ones(N,1),'linewidth',2);
title(['Learning Curve for \mu = 0.075 and b = ', num2str(B(z))],'color','b','fontsize',13)
xlabel('$Time(n)$','interpreter','latex','fontsize',13);
ylabel('$J(n)$','interpreter','latex','fontsize',13);
legend('$J(n)$','$J_{min}$','interpreter','latex','fontsize',13)
end
%% Part N
disp('_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_')
disp('Part N :')
disp(' ')
b = 3.1;
N_N = 3000; % New Time(n)
H_N = 0.5*(1 + cos((2*pi/b).*(L_Ch - 2))); % Channel Response
r_h_N = conv(H_N,fliplr(H_N));
r_x_N = [r_h_N(3:5),zeros(1,M-3)]; % Auto Correlation
R_x_N = toeplitz(r_x_N); % Auto Correlation Matrix
R_v_N = Var_Noise * eye(M); % Auto Correlation Matrix of Noise
R_u_N = R_x_N + R_v_N; % Auto Correlation Matrix of Input
% Generating P
for i = 1:M
if (Delay - i + 1 == 3) || (Delay - i + 1 == 2) || (Delay - i + 1 == 1)
P_N(i,1) = H_N(Delay - i + 1);
else
P_N(i,1) = 0;
end
end
W_opt_N = inv(R_u_N)*P_N; % Optimum Wiener Filter
J_min_N = Sig_d - P_N.'*W_opt_N; % Minimum MSE
disp(['Minimum MSE is Equals with :',num2str(J_min_N)])
Mo_N = [0.0075 0.025 0.075]; % Step Size
S_N = zeros(N_N,1); % Signal
E_N = zeros(N_N,Rep); % Error
J_n_N = zeros(N_N,length(Mo_N)); % MSE
W_N = zeros(M,N_N); % Filter Weights
for i=1:length(Mo_N)
J_N = zeros(N_N,1);
for j=1:Rep
S_N = 2*randi([0 1],N_N,1)-1; % Mail Signal
X_N = conv(S_N,H_N); % Output of Channel
X_N = X_N(1:N_N); % Set N samples
V_N = sqrt(Var_Noise)*randn(N_N,1); % Noise
U_N = X_N + V_N; % Input of Filter
D_N(8:N_N,1) = S_N(1:N_N-7,1); % Desired Signal
% Aux Parameters
In_N = zeros(M,1);
W_AUX__N = zeros(M,N_N);
W_old_N = zeros(M,1);
% Calculating Error
for k=1:N_N
% e(n) = d(n) - W' * u(n)
% W(n+1) = W(n) + Mu*(e(n).u(n))
Error_N(k,1) = D_N(k,1) - W_old_N'*In_N; % Err.
W_AUX__N(:,k) = W_old_N + Mo_N(i)*Error_N(k,1)*In_N; % Weights
In_N = circshift(In_N,1); % Shifting
In_N(1,1) = U_N(k,1);
W_old_N = W_AUX__N(:,k);
end
% Updating Weights and MSE
W_N = W_N + W_AUX__N;
J_N = J_N + Error_N.^2;
end
W_N = W_N/N_N; % Final Weights
J_N = J_N/Rep; % Final MSE
J_n_N(:,i) = J_N; % Save for Plot
end
%% Plot
figure('name','Part N')
for z=1:length(Mo_N)
subplot(2,2,z)
semilogy(J_n_N(:,z),'linewidth',1)
grid on
hold on
semilogy(J_min_N*ones(N_N,1),'linewidth',1.5);
title(['Learning Curve for \mu = ',num2str(Mo_N(z)),' and b = 3.1'],'color','b','fontsize',13)
xlabel('$Time(n)$','interpreter','latex','fontsize',13);
ylabel('$J(n)$','interpreter','latex','fontsize',13);
legend('$J(n)$','$J_{min}$','interpreter','latex','fontsize',13)
end
figure('name','Part N - Comparison')
for z=1:length(Mo_N)
semilogy(J_n_N(:,z),'linewidth',1)
grid on
hold on
end
semilogy(J_min_N*ones(N_N,1),'linewidth',1.5);
title(['Learning Curve for \mu = ',num2str(Mo_N(z)),' and b = 3.1'],'color','b','fontsize',13)
xlabel('$Time(n)$','interpreter','latex','fontsize',13);
ylabel('$J(n)$','interpreter','latex','fontsize',13);
legend('Jmin : \mu = 0.0075','Jmin : \mu = 0.025','Jmin : \mu = 0.075','Jmin','fontsize',12)
%% Comparison Curves
% Compare Part C and D
figure('name','Compare Parts C & D')
for z=1:length(B)
semilogy(J_n_C(z,:),'--','linewidth',1.5)
hold on
semilogy(J_n_D(:,z),'-.','linewidth',1.5)
hold on
semilogy(J_min(z)*ones(N,1),'linewidth',2);
end
grid on
legend('Jmin C: b=2.9','Jmin D: b=2.9','Jmin: b=3.1','Jmin C: b=3.1','Jmin D: b=2.9','Jmin: b=3.1',...
'Jmin C: b=3.3','Jmin D: b=3.3','Jmin: b=3.3','Jmin C: b=3.5','Jmin D: b=3.5','Jmin: b=3.5',...
'interpreter','latex','fontsize',12)
title('Learning Curve','color','b','fontsize',13)
xlabel('$Time(n)$','interpreter','latex','fontsize',13,'color','r');
ylabel('$J(n)$','interpreter','latex','fontsize',13,'color','r');
format loose