clear;clc;close all;

% Loading in Data
A_26 = load('Aluminum_26V_250mA');
A_28 = load('Aluminum_28V_269mA');
B_26 = load('Brass_26V_245mA');
B_29 = load('Brass_29V_273mA');
S_21 = load('Steel_21V_192mA');

% Prelab Data
T = [17.6,21.61,25.13,29.22,34.92,38.10,45.21,47.01];
x_L = [1.375,1.875,2.375,2.875,3.375,3.875,4.375,4.875];

%% Steady State Distributions

A_28(:,2) = [];
B_26(:,2) = [];
B_29(:,2) = [];
S_21(:,2) = [];


% Steady State temperatures
T_a26 = A_26(end,2:end);
T_a28 = A_28(end,2:end);
T_b26 = B_26(end,2:end);
T_b29 = B_29(end,2:end);
T_s21 = S_21(end,2:end);


%Time vectors
time_a26 = A_26(:, 1);
time_a28 = A_28(:, 1);
time_b26 = B_26(:, 1);
time_b29 = B_29(:, 1);
time_s21 = S_21(:, 1);


% thermal conductivity values
ka = 130; %W/mK
kb = 115; %W/mk
ks = 16.2; %W/mK

A = pi*(0.5*0.0254)^2;

% Analytical temperature distribution slopes
H_an_a26 = (26*0.250*0.0254)/(ka*A); %C/in
H_an_a28 = (28*0.269*0.0254)/(ka*A); %C/in
H_an_b26 = (26*0.245*0.0254)/(kb*A); %C/in
H_an_b29 = (29*0.273*0.0254)/(kb*A); %C/in
H_an_s21 = (21*0.192*0.0254)/(ks*A); %C/in


H_an_a26_m = ((26*0.250*0.0254)/(ka*A))/0.0254; %C/m
H_an_a28_m = ((28*0.269*0.0254)/(ka*A))/0.0254; %C/m
H_an_b26_m = ((26*0.245*0.0254)/(kb*A))/0.0254; %C/m
H_an_b29_m = ((29*0.273*0.0254)/(kb*A))/0.0254; %C/m
H_an_s21_m = ((21*0.192*0.0254)/(ks*A))/0.0254; %C/m


% Linear fit of experimental data to determine Experimental temp distribution slope

p1 = polyfit(x_L,T_a26,1);
H_a26 = p1(1);
T0_a26 = p1(2);
v_a26 = H_a26*x_L + T0_a26;

p2 = polyfit(x_L,T_a28,1);
H_a28 = p2(1);
T0_a28 = p2(2);
v_a28 = H_a28*x_L + T0_a28;


p3 = polyfit(x_L,T_b26,1);
H_b26 = p3(1);
T0_b26 = p3(2);
v_b26 = H_b26*x_L + T0_b26;


p4 = polyfit(x_L,T_b29,1);
H_b29 = p4(1);
T0_b29 = p4(2);
v_b29 = H_b29*x_L + T0_b29;



p5 = polyfit(x_L,T_s21,1);
H_s21 = p5(1);
T0_s21 = p5(2);
v_s21 = H_s21*x_L + T0_s21;

H_a26_m = H_a26/0.0254; %C/m
H_a28_m = H_a28/0.0254; %C/m
H_b26_m = H_b26/0.0254; %C/m
H_b29_m = H_b29/0.0254; %C/m
H_s21_m = H_s21/0.0254; %C/m


%% Steady State Distributions


% analytical temperature distributions
u_an_a26 = H_an_a26*x_L + T0_a26;
u_an_a28 = H_an_a28*x_L + T0_a28;
u_an_b26 = H_an_b26*x_L + T0_b26;
u_an_b29 = H_an_b29*x_L + T0_b29;
u_an_s21 = H_an_s21*x_L + T0_s21;


%Steady State Distributions
figure()
plot(x_L,u_an_a26)
hold on
plot(x_L,v_a26)
plot(x_L,T_a26)
xlabel('Linear Distance (in)')
ylabel('Temperature Reading (deg C)')
legend('Analytical Model','Experimental Model','Experimental Data','Location','southeast')
title('Steady State Distribution for A26')
hold off

figure()
plot(x_L,u_an_a28)
hold on
plot(x_L,v_a28)
plot(x_L,T_a28)
xlabel('Linear Distance (in)')
ylabel('Temperature Reading (deg C)')
legend('Analytical Model','Experimental Model','Experimental Data','Location','southeast')
title('Steady State Distribution for A28')
hold off

figure()
plot(x_L,u_an_b26)
hold on
plot(x_L,v_b26)
plot(x_L,T_b26)
xlabel('Linear Distance (in)')
ylabel('Temperature Reading (deg C)')
legend('Analytical Model','Experimental Model','Experimental Data','Location','southeast')
title('Steady State Distribution for B26')
hold off

figure()
plot(x_L,u_an_b29)
hold on
plot(x_L,v_b29)
plot(x_L,T_b29)
xlabel('Linear Distance (in)')
ylabel('Temperature Reading (deg C)')
legend('Analytical Model','Experimental Model','Experimental Data','Location','southeast')
title('Steady State Distribution for B29')
hold off

figure()
plot(x_L,u_an_s21)
hold on
plot(x_L,v_s21)
plot(x_L,T_s21)
xlabel('Linear Distance (in)')
ylabel('Temperature Reading (deg C)')
legend('Analytical Model','Experimental Model','Experimental Data','Location','southeast')
title('Steady State Distribution for S21')
hold off



%% Time-Dependent Temperature Profiles


% Model 1A


alpha_A26 = 4.82*10^(-5); %m^2/s
alpha_A28 = alpha_A26;
alpha_B26 = 3.56*10^(-5);
alpha_B29 = alpha_B26;
alpha_S21 = 4.05*10^(-6);


%Model 1A

time_dependent_temp(H_an_a26,T0_a26, A_26, time_a26, alpha_A26, 100, 'Aluminum 26V')
time_dependent_temp(H_an_a28,T0_a28, A_28, time_a28, alpha_A28, 200, 'Aluminum 28V')
time_dependent_temp(H_an_b26,T0_b26, B_26, time_b26, alpha_B26, 300, 'Brass 26V')
time_dependent_temp(H_an_b29,T0_b29, B_29, time_b29, alpha_B29, 400, 'Brass 29V')
time_dependent_temp(H_an_s21,T0_s21, S_21, time_s21, alpha_S21, 500, 'Steele 21V')



%Model 1B

time_dependent_temp(H_a26,T0_a26, A_26, time_a26, alpha_A26, 600, 'Aluminum 26V')
time_dependent_temp(H_a28,T0_a28, A_28, time_a28, alpha_A28, 700, 'Aluminum 28V')
time_dependent_temp(H_b26,T0_b26, B_26, time_b26, alpha_B26, 800, 'Brass 26V')
time_dependent_temp(H_b29,T0_b29, B_29, time_b29, alpha_B29, 900, 'Brass 29V')
time_dependent_temp(H_s21,T0_s21, S_21, time_s21, alpha_S21, 1000, 'Steele 21V')

%% Initial State Distributions (Model 2)

% Initial temperatures
T_in_a26 = A_26(1,2:end);
T_in_a28 = A_28(1,2:end);
T_in_b26 = B_26(1,2:end);
T_in_b29 = B_29(1,2:end);
T_in_s21 = S_21(1,2:end);


% Linear fit of experimental data to determine Initial state slope

p1_in = polyfit(x_L,T_in_a26,1);
M_exp_a26 = p1_in(1);
T0_exp_a26 = p1_in(2);
v_exp_a26 = M_exp_a26*x_L + T0_exp_a26;

p2_in = polyfit(x_L,T_in_a28,1);
M_exp_a28 = p2_in(1);
T0_exp_a28 = p2_in(2);
v_exp_a28 = M_exp_a28*x_L + T0_exp_a28;

p3_in = polyfit(x_L,T_in_b26,1);
M_exp_b26 = p3_in(1);
T0_exp_b26 = p3_in(2);
v_exp_b26 = M_exp_b26*x_L + T0_exp_b26;

p4_in = polyfit(x_L,T_in_b29,1);
M_exp_b29 = p4_in(1);
T0_exp_b29 = p4_in(2);
v_exp_b29 = M_exp_b29*x_L + T0_exp_b29;

p5_in = polyfit(x_L,T_in_s21,1);
M_exp_s21 = p5_in(1);
T0_exp_s21 = p5_in(2);
v_exp_s21 = M_exp_s21*x_L + T0_exp_s21;


M_exp_a26_m = M_exp_a26/0.0254;
M_exp_a28_m = M_exp_a28/0.0254;
M_exp_b26_m = M_exp_b26/0.0254;
M_exp_b29_m = M_exp_b29/0.0254;
M_exp_s21_m = M_exp_s21/0.0254;


%Steady State Distributions
figure()
hold on
plot(x_L,v_exp_a26)
plot(x_L,T_in_a26)
xlabel('Linear Distance (in)')
ylabel('Temperature Reading (deg C)')
legend('Linear Fit Model','Experimental Data','Location','southeast')
title('Initial State Distribution A26')
hold off

figure()
hold on
plot(x_L,v_exp_a28)
plot(x_L,T_in_a28)
xlabel('Linear Distance (in)')
ylabel('Temperature Reading (deg C)')
legend('Linear Fit Model','Experimental Data','Location','southeast')
title('Initial State Distribution A28')
hold off

figure()
hold on
plot(x_L,v_exp_b26)
plot(x_L,T_in_b26)
xlabel('Linear Distance (in)')
ylabel('Temperature Reading (deg C)')
legend('Linear Fit Model','Experimental Data','Location','southeast')
title('Initial State Distribution B26')
hold off

figure()
hold on
plot(x_L,v_exp_b29)
plot(x_L,T_in_b29)
xlabel('Linear Distance (in)')
ylabel('Temperature Reading (deg C)')
legend('Linear Fit Model','Experimental Data','Location','southeast')
title('Initial State Distribution B29')
hold off

figure()
hold on
plot(x_L,v_exp_s21)
plot(x_L,T_in_s21)
xlabel('Linear Distance (in)')
ylabel('Temperature Reading (deg C)')
legend('Linear Fit Model','Experimental Data','Location','southeast')
title('Initial State Distribution S21')
hold off


%Model 2 Plots


model_2(H_a26, M_exp_a26, T0_exp_a26, A_26, time_a26, alpha_A26, 1100, 'Model 2 for Aluminum 26V')
model_2(H_a28, M_exp_a28, T0_exp_a28, A_28, time_a28, alpha_A28, 1200, 'Model 2 for Aluminum 28V')
model_2(H_b26, M_exp_b26, T0_exp_b26, B_26, time_b26, alpha_B26, 1300, 'Model 2 for Brass 26V')
model_2(H_b29, M_exp_b29, T0_exp_b29, B_29, time_b29, alpha_B29, 1400, 'Model 2 for Brass 29V')
model_2(H_s21, M_exp_s21, T0_exp_s21, S_21, time_s21, alpha_S21, 1500, 'Model 2 for Steele 21V')


%% Variane in Thermal Diffusivities (Model 3)

alpha_A26_adj = 4.82*10^(-5) - 2.5*10^(-5); %m^2/s
alpha_A28_adj = alpha_A26;
alpha_B26_adj = 3.56*10^(-5) - 2*10^(-5);
alpha_B29_adj = alpha_B26;
alpha_S21_adj = 4.05*10^(-6) + 0.5*10^(-6);

time_dependent_temp(H_a26, T0_a26, A_26, time_a26, alpha_A26_adj, 1600, 'Model 3 for Aluminum 26V')
time_dependent_temp(H_a28, T0_a28, A_28, time_a28, alpha_A28_adj, 1700, 'Model 3 for Aluminum 28V')
time_dependent_temp(H_b26, T0_b26, B_26, time_b26, alpha_B26_adj, 1800, 'Model 3 for Brass 26V')
time_dependent_temp(H_b29, T0_b29, B_29, time_b29, alpha_B29_adj, 1900, 'Model 3 for Brass 29V')
time_dependent_temp(H_s21, T0_s21, S_21, time_s21, alpha_S21_adj, 2000, 'Model 3 for Steele 21V')

%% Time to Steady State

tss_A26 = 1000;
tss_A28 = 1500;
tss_B26 = 2500;
tss_B29 = 2500;
tss_S21 = 8000;

L = 5.875*0.0254; %m

Fo_A26 = (alpha_A26_adj*tss_A26)/L^2;
Fo_A28 = (alpha_A28_adj*tss_A28)/L^2;
Fo_B26 = (alpha_B26_adj*tss_B26)/L^2;
Fo_B29 = (alpha_B29_adj*tss_B29)/L^2;
Fo_S21 = (alpha_S21_adj*tss_S21)/L^2;


%% Error of Models 2 and 3

error_bars_analysis(H_a26, M_exp_a26, T0_a26, A_26, time_a26,  alpha_A26, alpha_A26_adj, 3000, 'A26 Comparison of Model 2 and 3')
error_bars_analysis(H_a28, M_exp_a28, T0_a28, A_28, time_a28,  alpha_A28, alpha_A28_adj, 3100, 'A28 Comparison of Model 2 and 3')
error_bars_analysis(H_b26, M_exp_b26, T0_b26, B_26, time_b26,  alpha_B26, alpha_B26_adj, 3200, 'B26 Comparison of Model 2 and 3')
error_bars_analysis(H_b29, M_exp_b29, T0_b29, B_29, time_b29,  alpha_B29, alpha_B29_adj, 3300, 'B29 Comparison of Model 2 and 3')
error_bars_analysis(H_s21, M_exp_s21, T0_s21, S_21, time_s21,  alpha_S21, alpha_S21_adj, 3400, 'S21 Comparison of Model 2 and 3')







function time_dependent_temp(H, T0, experimental_data, time, alpha, fig_numb, fig_title)

x = [1.375,1.875,2.375,2.875,3.375,3.875,4.375,4.875]*0.0254; %m

H_m = H/0.0254; %C/m
L = (4 + (7/8) + 1)*0.0254; %m

n_axis = [1:1:10];

figure(fig_numb)
plot(experimental_data(:,1),experimental_data(:,2),'-b', 'LineWidth', 2)
hold on
plot(experimental_data(:,1),experimental_data(:,3),'-b', 'LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,4),'-b', 'LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,5),'-b', 'LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,6),'-b','LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,7),'-b','LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,8),'-b','LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,9),'-b','Linewidth', 2)
ylabel('Temeprature (deg C)')
xlabel('Time (s)')
title(fig_title)

linestyles = {'-r', '-r', '-r', '-r', '-r', '-r', '-r', '-r'};

figure(fig_numb)
hold on
for j = 1:length(x)

        for i=1:length(time)
             
            for n = 1:10
               
                %n = 1;
                bn = (8*H_m*L*((-1)^n))/((pi^2)*((2*n - 1)^2));
                lambda = ((2*n - 1)*pi)/(2*L);
        
                sums1(n) = bn*(sin(lambda*x(j)))*exp(-(lambda^2)*alpha*time(i));
            end
            
            u1(j,i) = (H_m*x(j) + T0) + sum(sums1(1:n));
             
        end
        plot(time, u1(j,:), linestyles{j}, 'LineWidth', 2)

    end

figure(fig_numb)
legend('T1','T2','T3','T4','T5','T6','T7','T8','T1_{an}','T2_{an}','T3_{an}','T4_{an}','T5_{an}','T6_{an}','T7_{an}','T8_{an}','Location','northwest')


end




function model_2(H, M, T0, experimental_data, time, alpha, fig_numb, fig_title)

x = [1.375,1.875,2.375,2.875,3.375,3.875,4.375,4.875]*0.0254; %m

H_m = H/0.0254; %C/m
M_m = M/0.0254;
L = (4 + (7/8) + 1)*0.0254; %m

n_axis = [1:1:10];

figure(fig_numb)
plot(experimental_data(:,1),experimental_data(:,2),'-b', 'LineWidth', 2)
hold on
plot(experimental_data(:,1),experimental_data(:,3),'-b', 'LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,4),'-b', 'LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,5),'-b', 'LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,6),'-b','LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,7),'-b', 'LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,8),'-b', 'LineWidth', 2)
plot(experimental_data(:,1),experimental_data(:,9),'-b', 'LineWidth', 2)
ylabel('Temeprature (deg C)')
xlabel('Time (s)')
title(fig_title)

linestyles = {'-r', '-r', '-r', '-r', '-r', '-r', '-r', '-r'};



figure(fig_numb)
hold on
for j = 1:length(x)

        for i=1:length(time)
             
            for n = 1:10
               
                %n = 1;
                bn = -(8*(M_m - H_m)*L*((-1)^n))/((pi^2)*((2*n - 1)^2));
                lambda = ((2*n - 1)*pi)/(2*L);
        
                sums1(n) = bn*(sin(lambda*x(j)))*exp(-(lambda^2)*alpha*time(i));
            end
            
            w = H_m*x(j) + T0;
            u1(j,i) = (w) + sum(sums1(1:n)) + M_m*x(j);
             
        end
        plot(time, u1(j,:), linestyles{j}, 'LineWidth', 2)

    end

figure(fig_numb)
legend('T1','T2','T3','T4','T5','T6','T7','T8','T1_{an}','T2_{an}','T3_{an}','T4_{an}','T5_{an}','T6_{an}','T7_{an}','T8_{an}','Location','northwest')



end




function error_bars_analysis(H, M, T0, experimental_data, time, alpha2, alpha3, fig_numb, fig_title)

x = [1.375,1.875,2.375,2.875,3.375,3.875,4.375,4.875]*0.0254; %m

H_m = H/0.0254; %C/m
M_m = M/0.0254;
L = (4 + (7/8) + 1)*0.0254; %m

n_axis = [1:1:10];

err = 2*zeros(1, length(experimental_data));
for z = 1:(length(err)/20)
    err(z*20) = 2;
end

figure(fig_numb)
hold on
errorbar(experimental_data(:,1),experimental_data(:,9), err,'-b', 'LineWidth', 2)
ylabel('Temeprature (deg C)')
xlabel('Time (s)')
title(fig_title)

linestyles = {'-r', '-r', '-r', '-r', '-r', '-r', '-r', '-r'};

figure(fig_numb)
hold on
for j = 1:length(x)

        for i=1:length(time)
             
            for n = 1:10
               
                %n = 1;
                bn = (8*H_m*L*((-1)^n))/((pi^2)*((2*n - 1)^2));
                lambda = ((2*n - 1)*pi)/(2*L);
        
                sums1(n) = bn*(sin(lambda*x(j)))*exp(-(lambda^2)*alpha3*time(i));
            end
            
            u1(j,i) = (H_m*x(j) + T0) + sum(sums1(1:n));
             
        end

end

figure(fig_numb)
plot(time, u1(8,:), linestyles{j}, 'LineWidth', 2)

for j = 1:length(x)

        for i=1:length(time)
             
            for n = 1:10
               
                %n = 1;
                bn = -(8*(M_m - H_m)*L*((-1)^n))/((pi^2)*((2*n - 1)^2));
                lambda = ((2*n - 1)*pi)/(2*L);
        
                sums1(n) = bn*(sin(lambda*x(j)))*exp(-(lambda^2)*alpha2*time(i));
            end
            
            w = H_m*x(j) + T0;
            u1(j,i) = (w) + sum(sums1(1:n)) + M_m*x(j);
             
        end
    end

figure(fig_numb)
plot(time, u1(8,:), '-g', 'LineWidth', 2)
legend('Th8 Experimental Data','Th8 Model 3','Th8 Model 2','Location','northwest')


end

