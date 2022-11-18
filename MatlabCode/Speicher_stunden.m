%% minimun costen_S, optimun Autarkiegrad_S,Eigenverbrauchsgrad_S:  Mit größerem Speicher bei verschiedenen Laufzeiten 
% Kosten des Speichers 1000€/kWh 
% Einspeisevergütung 0,078€/kWh
% Stromkosten für bezogene Energie 0,3€/kWh
% Laufzeit 10 oder 15 oder 20 Jahre
% Verbrauch 7000 KWh pro jahr.

SOC = reshape(SOC,[12,8760]);
V = reshape(Verbrauch,[12,8760]);
PV = reshape(PV,[12,8760]);
EVNB = reshape(EVNB,[12,8760]);
EINE = reshape(EINE,[12,8760]);
EABB = reshape(EABB,[12,8760]);
EIBE = reshape(EIBE,[12,8760]);

Faktor_V = 7e+06/sum(sum(V));
% V = V*Faktor_V;
% V = V_EAuto1;
% V = V_EAuto2;

SOC_lowerlimit = 5;% SOC_lowerlimit*SOC_max/100 = 545.0166 Wh
Laufzeit = linspace(1,20,20); % Jahre
% iterations = length(Laufzeit);
% Cos_S = zeros(iterations,1); % Kosten
Kap = linspace(-100,100,201);
profitS = zeros(8760,length(Kap));

SOC_radioS = zeros(length(Kap),12,8760);
SOC_absS = zeros(length(Kap),12,8760);
EINE_Ssim = zeros(length(Kap),12,8760);
EVNB_Ssim = zeros(length(Kap),12,8760);
EIBE_Ssim = zeros(length(Kap),12,8760);
EABB_Ssim = zeros(length(Kap),12,8760);
DV_sim = zeros(length(Kap),12,8760);

tic
for m = 1:length(Kap)
    SOC_Smax = SOC_max + Kap(m)*100  ; % mehr 2.5 kWh Akku 此模型只适用于增加电池容量的模拟   -+ 2500
    SOC_absS(m,1,1) = SOC(1,1)*SOC_max/100;
    EABB_Ssim(m,1,1) = EABB(1,1);
    EIBE_Ssim(m,1,1) = EIBE(1,1);
    %%%
    for s = 1:8760
        if s == 1 % initial die erste Stunden
            SOC_absS(m,1,1) = SOC(1,1)*SOC_max/100;
            EABB_Ssim(m,1,1) = EABB(1,1);
            EIBE_Ssim(m,1,1) = EIBE(1,1);
            EINE_Ssim(m,1,1) = EINE(1,1);
            EVNB_Ssim(m,1,1) = EVNB(1,1);
            for i = 2:12 % SOC berechnen ersten Stunden
                SOC_absS(m,i,s) = SOC_absS(m,i-1,s) + optimumsk1(s)*EIBE_Ssim(m,i-1,s) - optimumsk2(s)*EABB_Ssim(m,i-1,s);
                if PV(i,s) > V(i,s)
                    DV_sim(m,i,s) = V(i,s);
                    EVNB_Ssim(m,i,s) = 0;
                    EABB_Ssim(m,i,s) = 0;
                    if SOC_absS(m,i,s) < SOC_Smax
                        EIBE_Ssim(m,i,s) = PV(i,s) - V(i,s);
                        EINE_Ssim(m,i,s) = 0;
                    else
                        SOC_absS(m,i,s) = SOC_Smax;
                        EIBE_Ssim(m,i,s) = 0;
                        EINE_Ssim(m,i,s) = PV(i,s) - V(i,s);
                    end
                else
                    DV_sim(m,i,s) = PV(i,s);
                    EIBE_Ssim(m,i,s) = 0;
                    EINE_Ssim(m,i,s) = 0;
                    if (SOC_absS(m,i,s) > SOC_lowerlimit*SOC_max/100) && (SOC_absS(m,i,s)< SOC_Smax)  % >=
                        EVNB_Ssim(m,i,s) = 0;
                        EABB_Ssim(m,i,s) = V(i,s) - PV(i,s);
                    elseif SOC_absS(m,i,s) >= SOC_Smax
                        SOC_absS(m,i,s) = SOC_Smax;
                        EVNB_Ssim(m,i,s) = 0;
                        EABB_Ssim(m,i,s) = V(i,s) - PV(i,s);
                    else
                        SOC_absS(m,i,s) = SOC_lowerlimit*SOC_max/100;
                        EABB_Ssim(m,i,s) = 0;
                        EVNB_Ssim(m,i,s) = V(i,s) - PV(i,s);
                    end
                end
            end
        else
            SOC_absS(m,1,s) = SOC_absS(m,end,s-1) + optimumsk1(s)*EIBE_Ssim(m,end,s-1) - optimumsk2(s)*EABB_Ssim(m,end,s-1); % initial die erste Punkt 会出现低于SOC_lowerlimit
            if PV(1,s) > V(1,s)
                EABB_Ssim(m,1,s) = 0;
                if SOC_absS(m,1,s) < SOC_Smax
                    EIBE_Ssim(m,1,s) = PV(1,s) - V(1,s);
                    EINE_Ssim(m,1,s) = 0;
                else
                    SOC_absS(m,1,s) = SOC_Smax;
                    EIBE_Ssim(m,1,s) = 0;
                    EINE_Ssim(m,1,s) = PV(1,s) - V(1,s);
                end
            else
                EIBE_Ssim(m,1,s) = 0;
                if (SOC_absS(m,1,s) > SOC_lowerlimit*SOC_max/100) && (SOC_absS(m,1,s) < SOC_Smax) % >=
                    EVNB_Ssim(m,1,s) = 0;
                    EABB_Ssim(m,1,s) = V(1,s) - PV(1,s);
                elseif SOC_absS(m,1,s) >= SOC_Smax
                    SOC_absS(m,1,s) = SOC_Smax;
                    EVNB_Ssim(m,1,s) = 0;
                    EABB_Ssim(m,1,s) = V(1,s) - PV(1,s);
                else
                    SOC_absS(m,1,s) = SOC_lowerlimit*SOC_max/100;
                    EABB_Ssim(m,1,s) = 0;
                    EVNB_Ssim(m,1,s) = V(1,s) - PV(1,s);
                end
            end
            for i = 2:12
                SOC_absS(m,i,s) = SOC_absS(m,i-1,s) + optimumsk1(s)*EIBE_Ssim(m,i-1,s) - optimumsk2(s)*EABB_Ssim(m,i-1,s);
                if PV(i,s) > V(i,s)
                    DV_sim(m,i,s) = V(i,s);
                    EVNB_Ssim(m,i,s) =0;
                    EABB_Ssim(m,i,s) = 0;
                    if SOC_absS(m,i,s) < SOC_Smax
                        EIBE_Ssim(m,i,s) = PV(i,s) - V(i,s);
                        EINE_Ssim(m,i,s) = 0;
                    else
                        SOC_absS(m,i,s) = SOC_Smax;
                        EIBE_Ssim(m,i,s) = 0;
                        EINE_Ssim(m,i,s) = PV(i,s) - V(i,s);
                    end
                else
                    DV_sim(m,i,s) = PV(i,s);
                    EIBE_Ssim(m,i,s) = 0;
                    EINE_Ssim(m,i,s) = 0;
                    if (SOC_absS(m,i,s) > SOC_lowerlimit*SOC_max/100) && (SOC_absS(m,i,s) < SOC_Smax) % >=
                        EVNB_Ssim(m,i,s) = 0;
                        EABB_Ssim(m,i,s) = V(i,s) - PV(i,s);
                    elseif SOC_absS(m,i,s) >= SOC_Smax
                        SOC_absS(m,i,s) = SOC_Smax;
                        EVNB_Ssim(m,i,s) = 0;
                        EABB_Ssim(m,i,s) = V(i,s) - PV(i,s);
                    else
                        SOC_absS(m,i,s) = SOC_lowerlimit*SOC_max/100;
                        EABB_Ssim(m,i,s) = 0;
                        EVNB_Ssim(m,i,s) = V(i,s) - PV(i,s);
                    end
                end
            end
        end
%         仅仅对EVNB_Ssim 修正, 低容量阶段SOC模型误差较大,SOC模型要先于真实模型到最低点
%         idx1 = [];
%         for i = 1:12
%             if SOC_absS(m,i,s) == SOC_lowerlimit*SOC_max/100
%                 idx1 = [idx1,i];
%             end
%         end
%         for i = 1:length(idx1)
%             if (EVNB_Ssim(m,idx1(i),s) ~= 0) && (EVNB(idx1(i),s) == 0)
%                 EVNB_Ssim(m,idx1(i),s) = EVNB(idx1(i),s);
%             end
%         end
        profitS(s,m) = -(sum(EVNB_Ssim(m,:,s))*0.3 - sum(EINE_Ssim(m,:,s))*0.078)/1000; % 1-Jahre-Profit
    end
    %%%
    SOC_radioS(m,:,:) = SOC_absS(m,:,:)*100/SOC_Smax; % absluto --> relativ SOC_Smax
    %         Cos_S(k) = 2500 - Laufzeit(k)*sum(profitS); % Laufzeit-Jahre-Kosten
end
toc


profit_S = sum(profitS); % ohne Kosten des Speichers
figure
plot((SOC_max + Kap*100)/1000,profit_S,'o','linewidth',1)
grid on
xlabel('$ kW/h$','interpreter','latex', 'FontSize', 18)
ylabel('$ Euro $','interpreter','latex', 'FontSize', 18)
title('Profit ohne Kosten des Speichers')
profit_S = sum(profitS) - 100*Kap; % mit Kosten des Speichers
figure
plot((SOC_max + Kap*100)/1000,profit_S,'o','linewidth',1)
grid on
xlabel('$ kW/h$','interpreter','latex', 'FontSize', 18)
ylabel('$ Euro $','interpreter','latex', 'FontSize', 18)
title('Profit mit Kosten des Speichers')

%% 检查天 EABB EIBE
EABB_Ssim = reshape(EABB_Ssim,3,288,365);
EIBE_Ssim = reshape(EIBE_Ssim,3,288,365);

Kap = 2;
k = 177;
figure
plot(EABB_Ssim(Kap,:,k),'--','linewidth',2)
hold on
plot(EIBE_Ssim(Kap,:,k),'-','linewidth',2)
grid on
xlabel('$ t/5min$','interpreter','latex', 'FontSize', 18)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 18)
legend('$ EABB_{sim} $','$ EIBE_{sim} $','interpreter','latex', 'FontSize', 20)
title(['Tag. ',num2str(k)],'interpreter','latex', 'FontSize', 18)

%% EVNB EINE
EVNB_85 = sum(EVNB_Ssim(1,:))
EVNB_11 = sum(EVNB_Ssim(2,:))
EVNB_135 = sum(EVNB_Ssim(3,:))

EINE_85 = sum(EINE_Ssim(1,:))
EINE_11 = sum(EINE_Ssim(2,:))
EINE_135 = sum(EINE_Ssim(3,:))

%% Profit
EVNB = reshape(Energie_vom_Netz_bezogen,[288,365]);
EINE = reshape(Energie_ins_Netz_eingespeist,[288,365]);
profit = zeros(365,1);
for j =1:365
    profit(j) = -(sum(EVNB(:,j))*0.3 - sum(EINE(:,j))*0.078)/1000; % 1-Jahre-Profit
end
sum(profit)
profit_S = sum(profitS)

A = reshape(profitS,24,365,3);
b = sum(A);
profitS_11_11 = b(1,:,2);
profitS_11_85 = b(1,:,1);
profitS_11_135 = b(1,:,3);

figure
plot(profitS_11_11,'o','linewidth',1)
hold on
plot(profit,'*','linewidth',1)
grid on
xlabel('$ Tag $','interpreter','latex', 'FontSize', 18)
ylabel(' Gewinn / €', 'interpreter','latex', 'FontSize', 18)
legend('$ 11 \; kWh \;sim $','$ 11 \; kWh \; data $','interpreter','latex', 'FontSize', 20)


%%% 盈利去除8个异常天
delete = [47 79 115 146 177 208 239 270];
profit(delete) = []; 
profitS_11_11(delete) = [];
sum(profit)
sum(profitS_11_11)

%% Profit 13,5 与 11 kWh 模拟值比较
a = profitS_11_135-profitS_11_11;
profit_blue = [];
profit_red = [];
for j =1:365
    if a(j) >= 0
        profit_blue = [profit_blue,j];
    else
        profit_red = [profit_red,j];
    end
end
figure
plot(profit_blue,a(profit_blue),'o')
hold on
plot(profit_red,a(profit_red),'o')
grid on
xlabel('$ Tag $','interpreter','latex', 'FontSize', 18)
ylabel('\Delta Gewinn_5 / €','FontSize', 18)

%% Profit-- 绘图月份柱状图
monate = [31 28 31 30 31 30 31 31 30 31 30 31];
% PV-Monate-Verlauf
profitS_M = zeros(3,length(monate));
profit_M = zeros(1,length(monate));
profitS = reshape(profitS,24,365,3);
A = sum(profitS);

mm = zeros(3,length(monate));
for m = 1:3
    for i = 1:length(monate)
        mm(i) = sum(monate(1:i));
        if i == 1
            profitS_M(m,i) = sum(sum(A (1,1:mm(i),m)));
            profit_M(i) = sum(profit(1:mm(i)));
        else
            profitS_M(m,i) = sum(sum(A (1,1+mm(i-1):mm(i),m)));
            profit_M(i) = sum(profit(1+mm(i-1):mm(i)));
        end
    end
end

a = [profitS_M;profit_M];
figure
bar(a.')
grid on
xlabel('$ Monat $','interpreter','latex', 'FontSize', 18)
ylabel('Gewinn / €','FontSize', 18)
legend('$ 8,5 \; kWh \; sim $','$11 \; kWh \; sim$','$13,5 \; kWh \; sim$','$11 \; kWh \; data$','interpreter','latex', 'FontSize', 16)

%% Cos
Cos_S = zeros(1,length(Laufzeit));
for k = 1:length(Laufzeit)
    Cos_S(k) = 2500 - Laufzeit(k)*sum(profitS_11_135); % Laufzeit-Jahre-Kosten
end
figure
plot(Cos_S,'o')
grid on
xlabel('$ Jahr $','interpreter','latex', 'FontSize', 18)
ylabel('Kosten / €','FontSize', 18)
% ylabel('$ Kosten $','interpreter','latex', 'FontSize', 18)

%% EVNB EINE-- 绘图月份柱状图
monate = [31 28 31 30 31 30 31 31 30 31 30 31]*24;
EVNB_SimM = zeros(3,length(monate));
EINE_SimM = zeros(3,length(monate));
mm = zeros(3,length(monate));
for m = 1:3
    for i = 1:length(monate)
        mm(i) = sum(monate(1:i));
        if i == 1
            EVNB_SimM(m,i) = sum(sum(EVNB_Ssim(m,:,1:mm(i))));
            EINE_SimM(m,i) = sum(sum(EINE_Ssim(m,:,1:mm(i))));
        else
            EVNB_SimM(m,i) = sum(sum(EVNB_Ssim(m,:,1+mm(i-1):mm(i))));
            EINE_SimM(m,i) = sum(sum(EINE_Ssim(m,:,1+mm(i-1):mm(i))));
        end
    end
end

figure
bar(EVNB_SimM.'/1000)
grid on
xlabel('$ Monat $','interpreter','latex', 'FontSize', 18)
ylabel('$ kWh $','interpreter','latex', 'FontSize', 18)
legend('8,5 kWh','11 kWh','13,5 kWh','interpreter','latex', 'FontSize', 18)
title('$ EVNB $','interpreter','latex', 'FontSize', 18)


figure
bar(EINE_SimM.'/1000)
grid on
xlabel('$ Monat $','interpreter','latex', 'FontSize', 18)
ylabel('$ kWh $','interpreter','latex', 'FontSize', 18)
legend('8,5 kWh','11 kWh','13,5 kWh','interpreter','latex', 'FontSize', 18)
title('$ EINE $','interpreter','latex', 'FontSize', 18)

sum(EVNB_SimM.'/1000)
sum(EINE_SimM.'/1000)
%% 检查第91与第21,22天
k = 146;
Kap = 2;
figure
plot(reshape(SOC(:,(k-1)*24+1:k*24),288,1),'--','linewidth',3)
hold on
plot(reshape(SOC_radioS(Kap,:,(k-1)*24+1:k*24),288,1),'-','linewidth',3)
grid on
xlabel('$ t/5min$','interpreter','latex', 'FontSize', 14)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 14)
legend('$ 11 kWh $','$ 13,5 kWh $','interpreter','latex', 'FontSize', 22)
title('$ Tag.146 $','interpreter','latex', 'FontSize', 14) 
% title(['Tag. ',num2str(k),' und Tag,',num2str(k+1)],'interpreter','latex', 'FontSize', 14)

figure
plot(reshape(EIBE(:,(k-1)*24+1:k*24),288,1),'--','linewidth',3)
hold on
plot(reshape(EIBE_Ssim(Kap,:,(k-1)*24+1:k*24),288,1),'-','linewidth',3)
grid on
xlabel('$ t / 5min $','interpreter','latex', 'FontSize', 14)
ylabel('$ wH \%  $','interpreter','latex', 'FontSize', 14)
legend('$ EIBE $','$ EIBE sim $','interpreter','latex', 'FontSize', 22)
title('$ Tag.146 $','interpreter','latex', 'FontSize', 14) 
% title(['Tag. ',num2str(k),' und Tag,',num2str(k+1)],'interpreter','latex', 'FontSize', 14)

%%
k = 277;
Kap = 2;
figure
plot(reshape(SOC(:,(k-1)*24+1:(k+1)*24),576,1),'--','linewidth',3)
hold on
plot(reshape(SOC_radioS(Kap,:,(k-1)*24+1:(k+1)*24),576,1),'-','linewidth',3)
grid on
xlabel('$ t/5min$','interpreter','latex', 'FontSize', 14)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 14)
legend('$ data $','$ sim $','interpreter','latex', 'FontSize', 22)
title('$ Tag.21,22 $','interpreter','latex', 'FontSize', 14) 
title(['Tag. ',num2str(k),' und Tag,',num2str(k+1)],'interpreter','latex', 'FontSize', 14)


%% 检查所有天
A = reshape(SOC_radioS(2,:,:),288,365);
SOC = reshape(SOC,288,365);

SOC_radioS_mean = sum(A)/288;
SOC_mean = sum(SOC)/288;
MSE_S = sqrt(sum((A-SOC).^2)/288);

figure
plot(SOC_mean,'--','linewidth',2)
hold on
plot(SOC_radioS_mean,'-','linewidth',2)
hold on
plot(abs(SOC_mean-SOC_radioS_mean),':','linewidth',2)
% hold on
% plot(MSE_S,':','linewidth',2)
grid on
xlabel('$ t/5min$','interpreter','latex', 'FontSize', 14)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 14)
legend('$ data $','$ sim $','$ diff $','interpreter','latex', 'FontSize', 20)

sum(abs(SOC_mean-SOC_radioS_mean))/365
%% 检查异常天的k1k2
k = 207;

figure
plot(optimumsk1(24*k-23:24*k+24),'--','linewidth',2)
grid on
xlabel('$ Stunden $','interpreter','latex', 'FontSize', 18)
ylabel('$ k_1 $','interpreter','latex', 'FontSize', 18)
legend('$ k_1 $','interpreter','latex', 'FontSize', 20)
title(['Tag. ',num2str(k),' und Tag,',num2str(k+1)],'interpreter','latex', 'FontSize', 14)

figure
plot(optimumsk2(24*k-23:24*k+24),'--','linewidth',2)
grid on
xlabel('$ Stunden$','interpreter','latex', 'FontSize', 18)
ylabel('$ k_2 $','interpreter','latex', 'FontSize', 18)
legend('$ k_2 $','interpreter','latex', 'FontSize', 20)
title(['Tag. ',num2str(k),' und Tag,',num2str(k+1)],'interpreter','latex', 'FontSize', 14)

%% 绘图SOC--12:00 und 00:00 Verlauf
% 12:00 Verlauf
SOC_radioS_year = reshape(SOC_radioS(3,:,:),288*365,1);
A = reshape(SOC,288*365,1);
C1 = zeros(365,1);
C2 = zeros(365,1);
for i = 1:365
    C1(i) = A(144+(i-1)*288);
    C2(i) = SOC_radioS_year(144+(i-1)*288);
end
figure
plot(C1,'-','linewidth',1.5)
hold on
plot(C2,'--','linewidth',1.5)
grid on
xlabel('$ Tage $','interpreter','latex', 'FontSize', 20)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 20)
legend('$ 11 kWh $','$ 13,5 kWh $','interpreter','latex', 'FontSize', 16)
% title('SOC 1 Jahrverlauf mit 2.5 KWh mehr AKK')


% 00:00 Verlauf
C1 = zeros(365,1);
C2 = zeros(365,1);
for i = 1:365
    C1(i) = A(1+(i-1)*288);
    C2(i) = SOC_radioS_year(1+(i-1)*288);
end
figure
plot(C1,'-','linewidth',1.5)
hold on
plot(C2,'--','linewidth',1.5)
grid on
xlabel('$ Tage $','interpreter','latex', 'FontSize', 20)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 20)
legend('$ 11 kWh $','$ 13,5 kWh $','interpreter','latex', 'FontSize', 16)
% title('SOC 1 Jahrverlauf mit 2.5 KWh mehr AKK')

%% 绘图SOC--Sommer und Winter Tage-Verlauf
% Sommer: 6,7,8. Monaten
SOC = reshape(SOC,288,365);
SOC_radioS = reshape(SOC_radioS,3,288,365);
k = 151:243; 
SOC_sommer = zeros(92,1);
SOC_Ssommer = zeros(92,1);
for i = 1:length(k)
    SOC_sommer(i) = sum(SOC(:,k(i)))/288; %mean-wert
    SOC_Ssommer(i) = sum(SOC_radioS(3,:,k(i)))/288; %mean-wert
end

figure
plot(k,SOC_sommer,'--','linewidth',2)
hold on
plot(k,SOC_Ssommer,'-','linewidth',2)
grid on
xlabel('$ Tage $','interpreter','latex', 'FontSize', 20)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 20)
legend('$ 11 kWh $','$ 13,5 kWh $','interpreter','latex', 'FontSize', 16)

% Winter: 12,1,2. Monaten
k1 = 1:62;
k2 = 334:365;
k3 = [k2,k1];
SOC_winter = zeros(92,1);
SOC_Swinter = zeros(92,1);
for i = 1:length(k3)
    SOC_winter(i) = sum(SOC(:,k3(i)))/288; %mean-wert
    SOC_Swinter(i) = sum(SOC_radioS(3,:,k3(i)))/288; %mean-wert
end

a = -31:62;
figure
plot(a,SOC_winter,'--','linewidth',2)
hold on
plot(a,SOC_Swinter,'-','linewidth',2)
grid on
xlabel('$ Tage $','interpreter','latex', 'FontSize', 20)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 20)
legend('$ 11 kWh $','$ 13,5 kWh $','interpreter','latex', 'FontSize', 16)

%% 绘图SOC--Sommer und Winter Fragment-Verlauf
% Sommer: 7.1-7.7
k = 181:187; 
A = reshape(SOC,288*365,1);
B = reshape(SOC_radioS,3,288*365);
SOC_sommer = A(k(1)*288:k(end)*288);
SOC_Ssommer = B(3,k(1)*288:k(end)*288);

figure
plot(SOC_sommer,'--','linewidth',2)
hold on
plot(SOC_Ssommer,'-','linewidth',2)
grid on
xlabel('$ t/5min $','interpreter','latex', 'FontSize', 18)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 18)
legend('$ 11 kWh $','$ 13,5 kWh $','interpreter','latex', 'FontSize', 20)

% Winter: 1.1-1.7
k = 1:7; 
SOC_winter = A(k(1)*288:k(end)*288);
SOC_Swinter = B(3,k(1)*288:k(end)*288);

figure
plot(SOC_winter,'--','linewidth',2)
hold on
plot(SOC_Swinter,'-','linewidth',2)
grid on
xlabel('$ t/5min $','interpreter','latex', 'FontSize', 18)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 18)
legend('$ 11 kWh $','$ 13,5 kWh $','interpreter','latex', 'FontSize', 20)
