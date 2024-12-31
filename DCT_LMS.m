clear
close all
clc
format compact

% Parameters
rng('default')
Var_Noise = input('Enter Variance of Noise : ');        % Variance of Noise
Sig_d =  input('Enter Power of Desired Signal : ');     % Sigma^2_d
M     =  input('Enter Length of Filter : ');            % Length of Filter
L     =  input('Enter Length of Channel : ');           % Length of Channel
Delay =  input('Enter The Value of Delay : ');          % Value of Delay
N     =  input('Enter Time Length(n) : ');              % Descrete Time Length
Rep   =  input('Enter The Number of Iterration : ');    % Number of Iterration
B     =  input('Enter The Value of b(in Vector) : ');   % Channel Coefficient
Beta  =  input('Enter The Value of Beta : ');           % Value of Beta
Mu    =  input('Enter The Value of Mu : ');             % Value of Mu
L_Ch  =  1:L;
disp(' LMS & DCT-LMS : Numerical Results')
H_Tot = zeros(length(B),L);             % Total Channel Response
J_min = zeros(1,length(B));             % Minimum MSE
P_Tot = zeros(M,length(B));             % Total Cross Correlation
J_n1  = zeros(N,length(B));             % MSE
W_D   = zeros(M,N);                     %Filter Weights

for z=1:length(B)

    disp(['For b = ',num2str(B(z))])
    disp(' ')

    H   = 0.5*(1 + cos((2*pi/B(z)).*(L_Ch - 2)));   % Channel Response
    H_Tot(z,:) = H;
    r_h = conv(H,fliplr(H));                        % Convlolution of Ch. Res.
    r_x = [r_h(L:(2*L) - 1),zeros(1,M - L)];
    R_x = toeplitz(r_x);                    % Correlation Matrix of Signal
    R_v = Var_Noise*eye(M);                 % Correlation Matrix of Noise
    
    % U(n) = X(n) + V(n) , X(n) = S(n)*H(n)
    R_u = R_x + R_v;                        % Correlation Matrix of Input
    
    EVD = eig(R_u);                         % EigenValue Decomposition
    X_R = (max(EVD))/(min(EVD));            % EigenValue Spread
    
    
    P = zeros(M,1);                         % Correlation of Input & Desired
    
    for i = 1:M
        if (Delay - i + 1 == 3) || (Delay - i + 1 == 2) || (Delay - i + 1 == 1)
            P(i,1) = H(Delay - i + 1);
        else
            P(i,1) = 0;
        end
    end
    P_Tot(:,z) = P;
    
    W_opt = R_u\P;                          % Optimum Wiener Filter
    J_min(z) = Sig_d - P.'*W_opt;           % MSE 
    
    disp(['Jmin Equals With : ',num2str(J_min(z))])
    disp(' ')
    
       
   
    S    = zeros(N,1);               % Signal
    W_D1 = zeros(M,N);               % Filter Weights
    J_D1 = zeros(N,1);               % MSE
    T    = zeros(M,M);               % Transform Matrix
    J_D  = zeros(N,1);
    % Updating Transform Matrix 
    for k=1:M
        for l=1:M
            if k==1                     % For k=0
                T(k,l) = 1/sqrt(M);
            else
                T(k,l) = sqrt(2/M)*cos((pi*(k-1)*(2*(l-1)+1))/(2*M));
            end
        end
    end
    
    R_t   = T*R_u*T';
    EVD_t = eig(R_t);                         % EigenValue Decomposition after Transform
    X_R_T = (max(EVD_t))/(min(EVD_t));        % EigenValue Spread after Transform
    
    disp(['EigenValue Spread Before Transform Equals With : ',num2str(X_R),...
        ' And After Transform Equals With : ',num2str(X_R_T)])
    disp(' ')
    
    for m=1:Rep
        
        S = 2*randi([0 1],N,1)-1;         % Signal
        X = conv(S,H);                    % Output of Channel
        X = X(1:N);                       % Set N samples
        V = sqrt(Var_Noise)*randn(N,1);   % Noise
        U = X + V;                        % Input of Filter
        D_D(8:N,1) = S(1:N-7,1);          % Desired Signal
        
        % AUX Parameters
        U_AD1  = zeros(M,1);    % DCT-LMS
        W_AUX1 = zeros(M,N);    % DCT-LMS
        W_old1 = zeros(M,1);    % DCT-LMS
        Sig_k  = ones(1,M);     % DCT-LMS
        
        U_AD   = zeros(M,1);    % LMS
        W_AUX  = zeros(M,N);    % LMS
        W_old  = zeros(M,1);    % LMS
        
        
        % Calculating Error
        for j=1:N
            
            % For LMS
            % e(n) = d(n) - W' * u(n)
            % W(n+1) = W(n) + Mu*(e(n).u(n))
            
            Error_D(j,1) = D_D(j,1) - W_old'*U_AD;       % Err. Equivalent
            W_AUX(:,j)   = W_old + Mu*Error_D(j,1)*U_AD;
            U_AD         = circshift(U_AD,1);            % Shifting
            U_AD(1,1)    = U(j,1);
            W_old        = W_AUX(:,j);                   % Save for Next Loop
            
            
            % For DCT - LMS
            % e(n) = d(n) - W' * x(n)
            % W(n+1) = W(n) + Mu*D*(e*(n).x(n))
            
            U_AD1         = circshift(U_AD1,1);  % Shifting
            U_AD1(1,1)    = U(j,1);              % Updating
            X_N = T*U_AD1;                       % Output of Transform Matrix
            
            % Calculating Sigma^2_k
            for q=1:M
                Sig_k(q) = Beta*Sig_k(q) + (1 - Beta)*abs(X_N(q))^2;
            end
            
            D = diag(1./Sig_k);                              % Diagonal Matrix
            Error_D1(j,1) = D_D(j,1) - W_old1'*X_N;          % Err. Equivalent
            W_AUX1(:,j)   = W_old1 + Mu*D*Error_D1(j,1)*X_N; % Updating Coeeficient
            W_old1        = W_AUX1(:,j);                     % Save for Next Loop
        end 
        
        W_D1 = W_D1 + W_AUX1;               % Updating, DCT-LMS
        J_D1 = J_D1 + Error_D1.^2;          % Updating, DCT-LMS
        
        W_D = W_D + W_AUX;                  % Updating, LMS
        J_D = J_D + Error_D.^2;             % Updating, LMS
    end
    
    % For DCT-LMS
    W_D1 = W_D1/Rep;                        % Final Weights,     DCT-LMS
    J_D1 = J_D1/Rep;                        % Final MSE,         DCT-LMS
    J_n1(:,z) = J_D1;                       % Save MSE for Plot, DCT-LMS
    
    % For LMS
    W_D = W_D/Rep;                          % Final Weights ,    LMS
    J_D = J_D/Rep;                          % Final MSE,         LMS
    J_n_D(:,z) = J_D;                       % Save MSE for Plot, LMS
   
    Jexc_DCT(z) = J_D1(end) - J_min(z);     % Jexcess DCT-LMS
    Jexc_LMS(z) = J_D(end) - J_min(z);      % Jexcess LMS
    
    MisAdj_DCT(z) = Jexc_DCT(z)/J_min(z);   % MissAdjustment DCT-LMS
    MisAdj_LMS(z) = Jexc_LMS(z)/J_min(z);   % MissAdjustment LMS
    
    disp(['MissAdjustment for LMS Equals With : ',num2str(MisAdj_LMS(z)),...
        ' MissAdjustment for DCT-LMS Equals With : ',num2str(MisAdj_DCT(z))])
    disp(' ')
end

%% Plot

figure('name','LMS & DCT-LMS')
 
for z=1:length(B)
    subplot(2,2,z)
    semilogy(J_n_D(:,z),'linewidth',1.5)        % LMS
    hold on
    semilogy(J_n1(:,z),'linewidth',1.5)         % DCT-LMS
    hold on
    semilogy(J_min(z)*ones(N,1),'linewidth',2); % Jmin
    grid on
    title(['Learning Curve for \mu = ',num2str(Mu),' and b = ', num2str(B(z))],'color','b','fontsize',13)
    xlabel('$Time(n)$','interpreter','latex','fontsize',13);
    ylabel('$J(n)$','interpreter','latex','fontsize',13);
    legend('$J(n)-LMS$','$J(n) - DCT-LMS$','$J_{min}$','interpreter','latex','fontsize',13)
    
end

figure('name','Learning Curves')
Color = ['r','b','g','m'];

for z=1:length(B)
    
    semilogy(J_n1(:,z),'linewidth',1.5,'color',Color(z))              % DCT-LMS
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