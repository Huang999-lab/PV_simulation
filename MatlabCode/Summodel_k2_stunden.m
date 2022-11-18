%% 迭代计算 Summodel_2.0k1k2
SOC_k1k2 = zeros(288*365,1);
SOC_k1k2(1) = Batch_SOC(1,1)/100*SOC_max;
SOC_s = reshape(SOC,12,8760);

initial = 1;
iterations = 600;
stepsize = 0.0001;
recordk2 = []; % k2_2.0
recordf = zeros(iterations,1); % MSE2.0k1k2 一天24个k_2
optimumsfk1k2 = zeros(8760,1);% Min-Costen from everyday 
optimumsk2 = zeros(8760,1);% k2 from everyday

tic
for s = 1:8760
    recordk2 = []; % recordk2清零
    
    for k = 1:iterations
        k2 = initial + (k-1)*stepsize;
        recordk2 = [recordk2,k2];
        SOC_k1k2 = zeros(288*365,1);% SOC清零
        SOC_k1k2(1) = SOC(1,1)/100*SOC_max;
        for i = 2:288*365
            SOC_k1k2(i) = SOC_k1k2(i-1) + optimumsk1(s)*Batch_EIBE(i-1) - k2*Batch_EABB(i-1);% absolute wert
        end
        SOC_k1k2 = SOC_k1k2*100/SOC_max; % radio wert
        SOC_k1k2 = reshape(SOC_k1k2,[12,8760]);
        recordf(k) = sum((SOC_k1k2(:,s) - SOC_s(:,s)).^2)/12;
    end
    
    optimumsfk1k2(s) = min(recordf);
    index_k2 = recordk2(recordf == optimumsfk1k2(s));
    optimumsk2(s) = index_k2(1);

end
toc  

%%
figure
plot(optimumsk2)
grid on
xlabel('$ Stunden $','interpreter','latex', 'FontSize', 20)
ylabel('$ k_2  $','interpreter','latex', 'FontSize', 20)

figure
plot(optimumk2)
grid on
xlabel('$ Tag $','interpreter','latex', 'FontSize', 20)
ylabel('$ k_2  $','interpreter','latex', 'FontSize', 20)
%% 绘制 三个模型Day.j的比较
j = 141; % Tag
SOC_sk1k2= zeros(365*24*12,1);
SOC_sk1k2(1) = Batch_SOC(1,1)/100*SOC_max;
SOC_simtag = [];

for s = (j-1)*24+1:j*24
    for i = 2:365*24*12
        SOC_sk1k2(i) = SOC_sk1k2(i-1) + optimumsk1(s)*Batch_EIBE(i-1) - optimumsk2(s)*Batch_EABB(i-1); %
    end
    SOC_simtag = [SOC_simtag; SOC_sk1k2(((s-1)*12+1):((s-1)*12+12))];
end
SOC_simtag = SOC_simtag*100/SOC_max;
SOC = reshape(SOC,[288,365]);

figure
plot(SOC(:,j),'--','linewidth',2)
hold on
grid on
plot(SOC_simtag,'-','linewidth',2)
hold off
xlabel('$ t/5min$','interpreter','latex', 'FontSize', 18)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 18)
legend('$ data $','$ sim $','$ diff $','interpreter','latex', 'FontSize', 20)
title(['Tag. ',num2str(j)],'interpreter','latex', 'FontSize', 14)
%% MSE plot
tic
recordf = zeros(8760,1); 
for s = 1:8760
    SOC_k1k2 = zeros(288*365,1);% SOC清零
    SOC_k1k2(1) = SOC(1,1)/100*SOC_max;
    for i = 2:288*365
        SOC_k1k2(i) = SOC_k1k2(i-1) + optimumsk1(s)*Batch_EIBE(i-1) - optimumsk2(s)*Batch_EABB(i-1);% absolute wert
    end
    SOC_k1k2 = SOC_k1k2*100/SOC_max; % radio wert
    SOC_k1k2 = reshape(SOC_k1k2,[12,8760]);
    recordf(s) = sqrt(sum((SOC_k1k2(:,s) - SOC_s(:,s)).^2)/12);
end
toc 

recordf = sum(reshape(recordf,24,365))/24;
figure
plot(recordf,'--','linewidth',2)
grid on
xlabel('$ Tag $','interpreter','latex', 'FontSize', 18)
ylabel('$ MSE \% $','interpreter','latex', 'FontSize', 18)