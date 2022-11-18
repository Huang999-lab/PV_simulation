%% 迭代计算 k1--Stunden
SOC_sk1= zeros(365*24*12,1);
SOC_sk1(1) = Batch_SOC(1,1)/100*SOC_max;

initial = 1;
iterations = 500;
stepsize = 0.001;
recordf = zeros(iterations,365*24); % MSE2.0k1  每小时一个k_1
recordk1 = []; % k1_2.0

for k = 1:iterations
    k1 = initial - (k-1)*stepsize;
    recordk1 = [recordk1,k1];
    for i = 2:365*24*12
        SOC_sk1(i) = SOC_sk1(i-1) + k1*Batch_EIBE(i-1) - Batch_EABB(i-1); 
%         SOC_m2(i) = SOC_m2(i-1) + k1*EIBE_sim(i-1) - EABB_sim(i-1);      %Batch_EIBE-->EIBE_sim,   Batch_EABB-->EABB_sim,   with simlation EI and EB
    end
    SOC_sk1 = SOC_sk1*100/SOC_max;
    for s = 1:365*24
        recordf(k,s) = sum((SOC_sk1((12*s-11):12*s) - Batch_SOC((12*s-11):12*s)).^2)/12;
    end
end

optimumf = zeros(365*24,1);
optimumsk1 = zeros(365*24,1);
for s = 1:365*24
    optimumf(s) = min(recordf(:,s));
    optimumsk1(s) = recordk1(find(recordf(:,s) == optimumf(s)));
end


%% 绘制 三个模型Day.k-Stunde.s的比较
s=1200; % Studen
SOC_sk1= zeros(365*24*12,1);
SOC_sk1(1) = Batch_SOC(1,1)/100*SOC_max;

for i = 2:365*24*12
    SOC_sk1(i) = SOC_sk1(i-1) + optimumsk1(s)*Batch_EIBE(i-1) - Batch_EABB(i-1); %   
end
SOC_sk1 = SOC_sk1*100/SOC_max;

C = reshape(SOC,[365*24*12,1]);
figure
plot(C(s*12:(s+1)*12))
hold on
grid on
plot(SOC_sk1(s*12:(s+1)*12))
hold off
xlabel('t / 5min')
ylabel('SOC / %')
legend('data','sim')%'sim1.0',
title(['vergleich Studen.',num2str(s)]) 

%%
figure
plot(optimumsk1)
grid on
xlabel('$ Stunden $','interpreter','latex', 'FontSize', 20)
ylabel('$ k_1  $','interpreter','latex', 'FontSize', 20)

figure
plot(optimumk1)
grid on
xlabel('$ Tag $','interpreter','latex', 'FontSize', 20)
ylabel('$ k_1  $','interpreter','latex', 'FontSize', 20)

%%
figure
plot(recordk1,recordf(:,5000))
grid on
xlabel('$ k_1 $','interpreter','latex', 'FontSize', 20)
ylabel('$ MSE \%  $','interpreter','latex', 'FontSize', 20)
