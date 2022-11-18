%%  V -> V_EAuto
% V 改变: WochenTag von 18:00 bis 0:00 Uhr  mehr Stromverbrauch
% durchschnittlichen Verbrauch von circa 15 kWh/100 km. (https://www.verivox.de/elektromobilitaet/themen/verbrauch-elektroauto/)
% ich nehme das E-Auto jeden Tag 10km fahren, damit 1.5 kWh oder 2.5kWh zu laden
% jeden Stunden 250 Wh oder 412 Wh laden, insgesamt 6 Stunden, jeder 5 minuten 20.8333 Wh oder 34.7222 mehr Verbrauch 
% 1.1.2020 ist Mittwoche

V = reshape(Verbrauch,[288,365]);
V_EAuto1 = V;
V_EAuto2 = V;
% set Montag,Dienstag,Mittwoche
for j = 1:3
    V_EAuto1(217:end,j:7:end) = V_EAuto1(217:end,j:7:end) + 20.8333;
    V_EAuto2(217:end,j:7:end) = V_EAuto2(217:end,j:7:end) + 34.7222;
end
% set Donerstag,Freitag
for j = 1:2
    V_EAuto1(217:end,(j+5):7:end) = V_EAuto1(217:end,(j+5):7:end) + 20.8333;
    V_EAuto2(217:end,(j+5):7:end) = V_EAuto2(217:end,(j+5):7:end) + 34.7222;
end

% Verbrauch 7000 KWh pro jahr.
Faktor_V = 7e+06/sum(sum(V));
V_EAuto3 = V*Faktor_V;
%% 最小costen_E, 最优Autarkiegrad_E,Eigenverbrauchsgrad_E:  Mit höherem Verbrauch und Simulierten Ladezeiten Auto (不考虑电网波动)
% Kosten des Speichers 1000€/kWh pro 1% Vergrößerung: (l-1)*100*1000
% Einspeisevergütung 0,078€/kWh
% Stromkosten für bezogene Energie 0,3€/kWh
% Laufzeit 10 oder 15 oder 20 Jahre

SOC_lowerlimit = 5;% SOC_lowerlimit*SOC_max/100 = 545.0166 Wh
Laufzeit = linspace(1,20,20); % Jahre
iterations = length(Laufzeit);
SOC_radioE = zeros(3,288,365);
Cos_E = zeros(3,iterations); % Kosten
recordk = []; % record years
profitE = zeros(3,365);
refer = [V_EAuto1,V_EAuto2,V_EAuto3];

tic
for n = 1:3
    V = refer(:,365*(n-1)+1:365*n);
    for k = 1:iterations
        recordk = [recordk,k];
        SOC_Smax = SOC_max + 2500  ; % mehr 2.5 kWh Akku 此模型只适用于增加电池容量的模拟+ 2500 
        SOC_absS = zeros(size(SOC));
        EINE_Ssim = zeros(size(SOC));
        EVNB_Ssim = zeros(size(SOC));
        EIBE_Ssim = zeros(size(SOC));
        EABB_Ssim = zeros(size(SOC));
        SOC_absS(1,1) = SOC(1,1)*SOC_max/100;
        EABB_Ssim(1,1) = EABB(1,1);
        EIBE_Ssim(1,1) = EIBE(1,1);
        DV_sim = zeros(size(DV));
        %%%
        for j = 1:365
            if j == 1 % initial die erste Punkt
                SOC_absS(1,1) = SOC(1,1)*SOC_max/100;
                EABB_Ssim(1,1) = EABB(1,1);
                EIBE_Ssim(1,1) = EIBE(1,1);
                EINE_Ssim(1,1) = EINE(1,1);
                EVNB_Ssim(1,1) = EVNB(1,1);
                for i = 2:288 % SOC berechnen ersten Tag
                    SOC_absS(i,j) = SOC_absS(i-1,j) + optimumk1(j)*EIBE_Ssim(i-1,j) - optimumk2(j)*EABB_Ssim(i-1,j);
                    if PV(i,j) > V(i,j)
                        DV_sim(i,j) = V(i,j);
                        EVNB_Ssim(i,j) = 0;
                        EABB_Ssim(i,j) = 0;
                        if SOC_absS(i,j) < SOC_Smax
                            EIBE_Ssim(i,j) = PV(i,j) - V(i,j);
                            EINE_Ssim(i,j) = 0;
                        else
                            SOC_absS(i,j) = SOC_Smax;
                            EIBE_Ssim(i,j) = 0;
                            EINE_Ssim(i,j) = PV(i,j) - V(i,j);
                        end
                    else
                        DV_sim(i,j) = PV(i,j);
                        EIBE_Ssim(i,j) = 0;
                        EINE_Ssim(i,j) = 0;
                        if (SOC_absS(i,j) > SOC_lowerlimit*SOC_max/100) && (SOC_absS(i,j) < SOC_Smax)  % >=
                            EVNB_Ssim(i,j) = 0;
                            EABB_Ssim(i,j) = V(i,j) - PV(i,j);
                        elseif SOC_absS(i,j) >= SOC_Smax
                            SOC_absS(i,j) = SOC_Smax;
                            EVNB_Ssim(i,j) = 0;
                            EABB_Ssim(i,j) = V(i,j) - PV(i,j);
                        else
                            SOC_absS(i,j) = SOC_lowerlimit*SOC_max/100;
                            EABB_Ssim(i,j) = 0;
                            EVNB_Ssim(i,j) = V(i,j) - PV(i,j);
                        end
                    end
                end
            else
                SOC_absS(1,j) = SOC_absS(end,j-1) + optimumk1(j)*EIBE_Ssim(end,j-1) - optimumk2(j)*EABB_Ssim(end,j-1); % initial die erste Punkt 会出现低于SOC_lowerlimit
                if PV(1,j) > V(1,j)
                    EABB_Ssim(1,j) = 0;
                    if SOC_absS(1,j) < SOC_Smax
                        EIBE_Ssim(1,j) = PV(1,j) - V(1,j);
                        EINE_Ssim(1,j) = 0;
                    else
                        SOC_absS(1,j) = SOC_Smax;
                        EIBE_Ssim(1,j) = 0;
                        EINE_Ssim(1,j) = PV(1,j) - V(1,j);
                    end
                else
                    EIBE_Ssim(1,j) = 0;
                    if (SOC_absS(1,j) > SOC_lowerlimit*SOC_max/100) && (SOC_absS(1,j) < SOC_Smax) % >=
                        EVNB_Ssim(1,j) = 0;
                        EABB_Ssim(1,j) = V(1,j) - PV(1,j);
                    elseif SOC_absS(1,j) >= SOC_Smax
                        SOC_absS(1,j) = SOC_Smax;
                        EVNB_Ssim(1,j) = 0;
                        EABB_Ssim(1,j) = V(1,j) - PV(1,j);
                    else
                        SOC_absS(1,j) = SOC_lowerlimit*SOC_max/100;
                        EABB_Ssim(1,j) = 0;
                        EVNB_Ssim(1,j) = V(1,j) - PV(1,j);
                    end
                end
                for i = 2:288
                    SOC_absS(i,j) = SOC_absS(i-1,j) + optimumk1(j)*EIBE_Ssim(i-1,j) - optimumk2(j)*EABB_Ssim(i-1,j);
                    if PV(i,j) > V(i,j)
                        DV_sim(i,j) = V(i,j);
                        EVNB_Ssim(i,j) =0;
                        EABB_Ssim(i,j) = 0;
                        if SOC_absS(i,j) < SOC_Smax
                            EIBE_Ssim(i,j) = PV(i,j) - V(i,j);
                            EINE_Ssim(i,j) = 0;
                        else
                            SOC_absS(i,j) = SOC_Smax;
                            EIBE_Ssim(i,j) = 0;
                            EINE_Ssim(i,j) = PV(i,j) - V(i,j);
                        end
                    else
                        DV_sim(i,j) = PV(i,j);
                        EIBE_Ssim(i,j) = 0;
                        EINE_Ssim(i,j) = 0;
                        if (SOC_absS(i,j) > SOC_lowerlimit*SOC_max/100) && (SOC_absS(i,j) < SOC_Smax) % >=
                            EVNB_Ssim(i,j) = 0;
                            EABB_Ssim(i,j) = V(i,j) - PV(i,j);
                        elseif SOC_absS(i,j) >= SOC_Smax
                            SOC_absS(i,j) = SOC_Smax;
                            EVNB_Ssim(i,j) = 0;
                            EABB_Ssim(i,j) = V(i,j) - PV(i,j);
                        else
                            SOC_absS(i,j) = SOC_lowerlimit*SOC_max/100;
                            EABB_Ssim(i,j) = 0;
                            EVNB_Ssim(i,j) = V(i,j) - PV(i,j);
                        end
                    end
                end
            end
            % 仅仅对EVNB_Ssim 修正, 低容量阶段SOC模型误差较大,SOC模型要先于真实模型到最低点
            idx1 = [];
            for i = 1:288
                if SOC_absS(i,j) == SOC_lowerlimit*SOC_max/100
                    idx1 = [idx1,i];
                end
            end
            for i = 1:length(idx1)
                if (EVNB_Ssim(idx1(i),j) ~= 0) && (EVNB(idx1(i),j) == 0)
                    EVNB_Ssim(idx1(i),j) = EVNB(idx1(i),j);
                end
            end
            profitE(n,j) = -(sum(EVNB_Ssim(:,j))*0.3 - sum(EINE_Ssim(:,j))*0.078)/1000; % 1-Jahre-Profit
        end
        %%%
        SOC_radioE(n,:,:) = SOC_absS*100/SOC_Smax; % absluto --> relativ SOC_Smax
        Cos_E(n,k) = 2500 - Laufzeit(k)*sum(profitE(n,:)); % Laufzeit-Jahre-Kosten
    end
end
toc

%% 绘图E_V_vergleich
% k.day
V = reshape(Verbrauch,[288,365]);
k = 1;
figure
plot(V(:,k),'-','linewidth',2)
hold on
plot(V_EAuto1(:,k),'--','linewidth',2)
hold on
plot(V_EAuto2(:,k),':','linewidth',2)
hold on
plot(V_EAuto3(:,k),'-.','linewidth',2)
grid on
legend('$ Rohrdaten $','$ Laden \; E-Auto \; mit \; 250 W/h  $','$ Laden \; E-Auto \; mit \; 412 W/h $','$ Haushaltsverbrauch $','interpreter','latex', 'FontSize', 10) %  
xlabel('$ t/5min $','interpreter','latex', 'FontSize', 20)
ylabel('$ W/h $','interpreter','latex', 'FontSize', 20)
title(['Tag. ',num2str(k)],'interpreter','latex', 'FontSize', 20)



%% profitE-- 绘图月份柱状图
monate = [31 28 31 30 31 30 31 31 30 31 30 31];
% PV-Monate-Verlauf
profitE_M = zeros(3,length(monate));
mm = zeros(3,length(monate));
for m = 1:3
    for i = 1:length(monate)
        mm(i) = sum(monate(1:i));
        if i == 1
            profitE_M(m,i) = sum(sum(profitE(m,1:mm(i))));
        else
            profitE_M(m,i) = sum(sum(profitE(m,1+mm(i-1):mm(i))));
        end
    end
end

figure
bar(profitE_M.')
grid on
xlabel('$ Monat $','interpreter','latex', 'FontSize', 18)
ylabel('Gewinn / €','FontSize', 18)
legend('$ Laden \; E-Auto \; mit \; 250 W/h $','$ Laden \; E-Auto \; mit \; 412 W/h $','$ Haushaltsverbrauch $','interpreter','latex', 'FontSize', 14)


%% profitE-- 绘图
figure
plot(profitE(1,:),'o','linewidth',1)
hold on
plot(profitE(2,:),'*','linewidth',1)
hold on
plot(profitE(3,:),'d','linewidth',1)
grid on
legend('$ Laden \; E-Auto \; mit \; 250 W/h $','$ Laden \; E-Auto \; mit \; 412 W/h $','$ Haushaltsverbrauch $','interpreter','latex', 'FontSize', 13)
xlabel('$ Tag $','interpreter','latex', 'FontSize', 18)
ylabel('Gewinn / €','FontSize', 18)

figure
plot(Cos_E(1,:),'-','linewidth',2)
hold on
plot(Cos_E(2,:),'--','linewidth',2)
hold on
plot(Cos_E(3,:),':','linewidth',2)
grid on
legend('$ Laden \; E-Auto \; mit \; 250 W/h $','$ Laden \; E-Auto \; mit \; 412 W/h $','$ Haushaltsverbrauch $','interpreter','latex', 'FontSize', 14)
xlabel('$ Jahr $','interpreter','latex', 'FontSize', 18)
ylabel('Kosten / €','FontSize', 18)

A1 = sum(profitE(1,:));
A2 = sum(profitE(2,:));
A3 = sum(profitE(3,:));

%% 绘图EINE-EVNB-SOC--Tagesverlauf
k = 91;% Tag.91
figure
plot(EINE_Ssim(:,k))
hold on
plot(EINE(:,k))
hold on
plot(EVNB_Ssim(:,k),'*')
hold on
plot(EVNB(:,k))
hold on
plot(SOC_radioE(:,k))
hold on
plot(SOC(:,k))
grid on
xlabel('$ t/5min$','interpreter','latex', 'FontSize', 14)
ylabel('$ Wh  $','interpreter','latex', 'FontSize', 14)
legend('EINE sim ','EINE ','EVNB sim ','EVNB','SOC sim','SOC','interpreter','latex', 'FontSize', 10)%
title(['Tag. ',num2str(k)],'interpreter','latex', 'FontSize', 14)

%% Tag 21,22
k = 21;
figure
plot([SOC_radioE(1,:,k),SOC_radioE(1,:,k+1)],'-','linewidth',2)
hold on
plot([SOC_radioE(2,:,k),SOC_radioE(2,:,k+1)],'--','linewidth',2)
hold on
plot([SOC_radioE(3,:,k),SOC_radioE(3,:,k+1)],':','linewidth',2)
grid on
xlabel('$ t/5min$','interpreter','latex', 'FontSize', 18)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 18)
legend('$ Laden \; E-Auto \; mit \; 250 W/h $','$ Laden \; E-Auto \; mit \; 412 W/h $','$ Haushaltsverbrauch $','interpreter','latex', 'FontSize', 12)
title(['Tag. ',num2str(k),' und ',num2str(k+1)],'interpreter','latex', 'FontSize', 18)
%% Tag 91
k = 91;
figure
plot(SOC_radioE(1,:,k),'-','linewidth',2)
hold on
plot(SOC_radioE(2,:,k),'--','linewidth',2)
hold on
plot(SOC_radioE(3,:,k),':','linewidth',2)
grid on
xlabel('$ t/5min$','interpreter','latex', 'FontSize', 18)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 18)
legend('$ Laden \; E-Auto \; mit \; 250 W/h $','$ Laden \; E-Auto \; mit \; 412 W/h $','$ Haushaltsverbrauch $','interpreter','latex', 'FontSize', 12)
title(['Tag. ',num2str(k)],'interpreter','latex', 'FontSize', 18)


%% Tag 21,22 und 91的 V与V1的SOC比较
k = 21;
figure
plot([SOC(:,k);SOC(:,k+1)],'-','linewidth',2)
hold on
plot([SOC_radioE(1,:,k),SOC_radioE(1,:,k+1)],'--','linewidth',2)
hold on
plot([V_EAuto1(:,k) - V(:,k);V_EAuto1(:,k+1) - V(:,k+1)],':','linewidth',2)
grid on
xlabel('$ t/5min$','interpreter','latex', 'FontSize', 18)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 18)
legend('$ Rohrdaten $','$ Laden \; E-Auto \; mit \; 250 W/h $','$ Lastzusatz $','interpreter','latex', 'FontSize', 12)
title(['Tag. ',num2str(k),' und ',num2str(k+1)],'interpreter','latex', 'FontSize', 18)

k = 91;
figure
plot(SOC(:,k),'-','linewidth',2)
hold on
plot(SOC_radioE(1,:,k),'--','linewidth',2)
hold on
plot(V_EAuto1(:,k) - V(:,k),':','linewidth',2)
grid on
xlabel('$ t/5min$','interpreter','latex', 'FontSize', 18)
ylabel('$ SOC \%  $','interpreter','latex', 'FontSize', 18)
legend('$ Rohrdaten $','$ Laden \; E-Auto \; mit \; 250 W/h $','$ Lastzusatz $','interpreter','latex', 'FontSize', 12)
title(['Tag. ',num2str(k)],'interpreter','latex', 'FontSize', 18)