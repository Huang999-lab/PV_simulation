%% Dataset
data1 = xlsread('PV-Anlage_Zack_01012020-31122020.xlsx');
% [C,ia,ic] = unique(data1,'rows');
% ia = ic :没有重复行
% data1: 2020年有366天
% 每5分钟测一次数据, 一天有288个样本，共有105630个样本(366*24*12=105408)? 数据多余222？
% data2 = xlsread('PV-Anlage_Zack_01012021-14112021.xlsx');
% data3 = xlsread('PV-Anlage_Zack_01092019-31122019.xlsx');

% error_idx = find(isnan(SOC)); % error-index
% data1(error_idx,:) = [];

SOC = data1(:,47);
Verbrauch = data1(:,46);
PV_Produktion = data1(:,45);
Energie_vom_Netz_bezogen = data1(:,44);
Energie_ins_Netz_eingespeist = data1(:,43);
Energie_in_Batterie_gespeichert = data1(:,42);
Energie_aus_Batterie_bezogen = data1(:,41);
Direkt_verbraucht = data1(:,40);

%%删除NaN元素,Dataset(data,day)
error_idx = find(isnan(SOC));
SOC(error_idx) = [];

error_idx = find(isnan(Verbrauch));
Verbrauch(error_idx) = [];

error_idx = find(isnan(PV_Produktion));
PV_Produktion(error_idx) = [];

error_idx = find(isnan(Energie_vom_Netz_bezogen));
Energie_vom_Netz_bezogen(error_idx) = [];

error_idx = find(isnan(Energie_ins_Netz_eingespeist));
Energie_ins_Netz_eingespeist(error_idx) = [];

error_idx = find(isnan(Energie_in_Batterie_gespeichert));
Energie_in_Batterie_gespeichert(error_idx) = [];

error_idx = find(isnan(Energie_aus_Batterie_bezogen));
Energie_aus_Batterie_bezogen(error_idx) = [];

error_idx = find(isnan(Direkt_verbraucht));
Direkt_verbraucht(error_idx) = [];

%%以天为单位取365天数据
year = 1:1:105120;
SOC = SOC(year);
Verbrauch = Verbrauch(year);
PV_Produktion = PV_Produktion(year);
Energie_vom_Netz_bezogen = Energie_vom_Netz_bezogen(year);
Energie_ins_Netz_eingespeist = Energie_ins_Netz_eingespeist(year);
Energie_in_Batterie_gespeichert = Energie_in_Batterie_gespeichert(year);
Energie_aus_Batterie_bezogen = Energie_aus_Batterie_bezogen(year);
Direkt_verbraucht = Direkt_verbraucht(year);

SOC = reshape(SOC,[288,365]);
V = reshape(Verbrauch,[288,365]);
PV = reshape(PV_Produktion,[288,365]);
EVNB = reshape(Energie_vom_Netz_bezogen,[288,365]);
EINE = reshape(Energie_ins_Netz_eingespeist,[288,365]);
EIBE = reshape(Energie_in_Batterie_gespeichert,[288,365]);
EABB = reshape(Energie_aus_Batterie_bezogen,[288,365]);
DV = reshape(Direkt_verbraucht,[288,365]);
%% 计算28次完整充放电的wirkungsgard
t = 1:1:288*365;%全年数据
Batch_SOC = SOC(t).';
Batch_V = Verbrauch(t);
Batch_PV = PV_Produktion(t);
Batch_EVNB = Energie_vom_Netz_bezogen(t);
Batch_EINE = Energie_ins_Netz_eingespeist(t);
Batch_EIBE = Energie_in_Batterie_gespeichert(t);
Batch_EABB = Energie_aus_Batterie_bezogen(t);
Batch_DV = Direkt_verbraucht(t);
%%% 简单求和求功率因素
wirkungsgard1 = sum(Batch_EABB)/sum(Batch_EIBE);

%%% 分charge与discharge阶段
idx_leer = find(Batch_SOC <= 5);
idx_voll = find(Batch_SOC == 100);
%%% charge model
idx_volleli_charge = [idx_voll(1); idx_voll(find(diff(idx_voll) ~= 1) + 1)];
idx_leeleli_charge = idx_leer(find(diff(idx_leer)~= 1));

%寻找每次Vollladung之前的首次leer状态的时间点 1434个vollladung对应1434个leer起点 很多无效起点
[row,~] = size(idx_volleli_charge);
idx_leercharge = zeros(row,1);
for i = 1:row
    A = find(idx_leeleli_charge < idx_volleli_charge(i));
    A = idx_leeleli_charge(A(end));
    idx_leercharge(i) = A;
end

%消除无效起点
A = find(diff(idx_leercharge) == 0);
idx_leercharge(A) = [];
idx_vollcharge = zeros(length(idx_leercharge),1);
for i = 1:length(idx_leercharge)
    B = find(idx_volleli_charge > idx_leercharge(i));
    idx_vollcharge(i) = B(1);
end

%一共28次Vollladung,每次存入能量储存在Energie_charge
[row,~] = size(idx_vollcharge);
Energie_Fullycharge = zeros(row,1);
for i = 1:row
    Energie_Fullycharge(i) = sum(Batch_EIBE(idx_leercharge(i):idx_volleli_charge(idx_vollcharge(i))) - Batch_EABB(idx_leercharge(i):idx_volleli_charge(idx_vollcharge(i))));
end


%%% discharge model
idx_volleli_discharge = idx_voll(find(diff(idx_voll) ~= 1));
idx_leeleli_discharge = idx_leer(find(diff(idx_leer)~= 1) + 1);

%寻找每次Vollladung之后的首次leer状态的时间点
[row,~] = size(idx_volleli_discharge);
idx_leerdischarge = zeros(row,1);
for i = 1:row
    A = find(idx_leeleli_discharge > idx_volleli_discharge(i));
    idx_leerdischarge(i) = A(1);
end

A = find(diff(idx_leerdischarge) == 0);
idx_leerdischarge(A) = [];
idx_volldischarge = idx_volleli_discharge;
idx_volldischarge(A) = [];

%一共28次Vollentladung,每次能量消耗储存在Energie_discharge
[row,~] = size(idx_volldischarge);
Energie_Fullydischarge = zeros(row,1);
for i = 1:row
    Energie_discharge(i) = sum(Batch_EABB(idx_volldischarge(i):idx_leeleli_discharge(idx_leerdischarge(i))) - Batch_EIBE(idx_volldischarge(i):idx_leeleli_discharge(idx_leerdischarge(i))));
end

wirkungsgard = zeros(28,1);
for i = 1:28
    wirkungsgard(i) = Energie_discharge(i)/Energie_Fullycharge(i);
end
wirkungsgard2 = sum(wirkungsgard)/length(wirkungsgard);

%%% SOC_max definition
SOC_max = sum(Energie_Fullycharge)/length(Energie_Fullycharge);

figure
plot(Energie_Fullycharge,'*')
grid on
hold on 
plot(Energie_discharge,'o')
hold off
xlabel('Nummer.')
ylabel('Energie / Wh')
legend('fullcharge','fulldischarge')
title('Energie in and aus SOC ')

%%  28次完整充放电绘图
figure
subplot(4,2,1);
plot(idx_leercharge(1):idx_volleli_charge(idx_vollcharge(1)),Batch_SOC(idx_leercharge(1):idx_volleli_charge(idx_vollcharge(1))))
grid on
title('1.fullcharge ')
subplot(4,2,2);
plot(idx_volldischarge(1):idx_leeleli_discharge(idx_leerdischarge(1)),Batch_SOC(idx_volldischarge(1):idx_leeleli_discharge(idx_leerdischarge(1))))
grid on
title('1.fulldischarge ')

subplot(4,2,3);
plot(idx_leercharge(19):idx_volleli_charge(idx_vollcharge(19)),Batch_SOC(idx_leercharge(19):idx_volleli_charge(idx_vollcharge(19))))
grid on
title('19.fullcharge ')
subplot(4,2,4);
plot(idx_volldischarge(19):idx_leeleli_discharge(idx_leerdischarge(19)),Batch_SOC(idx_volldischarge(19):idx_leeleli_discharge(idx_leerdischarge(19))))
grid on
title('19.fulldischarge ')

subplot(4,2,5);
plot(idx_leercharge(20):idx_volleli_charge(idx_vollcharge(20)),Batch_SOC(idx_leercharge(20):idx_volleli_charge(idx_vollcharge(20))))
grid on
title('20.fullcharge ')
subplot(4,2,6);
plot(idx_volldischarge(20):idx_leeleli_discharge(idx_leerdischarge(20)),Batch_SOC(idx_volldischarge(20):idx_leeleli_discharge(idx_leerdischarge(20))))
grid on
title('20.fulldischarge ')

subplot(4,2,7);
plot(idx_leercharge(22):idx_volleli_charge(idx_vollcharge(22)),Batch_SOC(idx_leercharge(22):idx_volleli_charge(idx_vollcharge(22))))
grid on
title('22.fullcharge ')
subplot(4,2,8);
plot(idx_volldischarge(22):idx_leeleli_discharge(idx_leerdischarge(22)),Batch_SOC(idx_volldischarge(22):idx_leeleli_discharge(idx_leerdischarge(22))))
grid on
title('22.fulldischarge ')

%% 每天检查实际与模型数据 EIBE,EABB:完全正确
for i = 1:365
    for j = 1:288
        EIBE_sim(j,i) = PV(j,i) - DV(j,i) - EINE(j,i);
        EABB_sim(j,i) = V(j,i) - DV(j,i) - EVNB(j,i);
    end
end

%%%
MSE_EIBE = zeros(365,1);
MSE_EABB = zeros(365,1);
for i = 1:365
    MSE_EIBE(i) = sum((EIBE(:,i) - EIBE_sim(:,i)).^2)/288;
    MSE_EABB(i) = sum((EABB(:,i) - EABB_sim(:,i)).^2)/288;
end

figure
plot(MSE_EIBE)
hold on
plot(MSE_EABB)
legend('MSE EIBE','MSE EABB')
%% 计算每天模型数据 DV 并与实际模型比较
DV_sim = zeros(size(DV));
for j = 1:365
    for i = 1:288
        if PV(i,j) > V(i,j)
            DV_sim(i,j) = V(i,j);
        else
            DV_sim(i,j) = PV(i,j);
        end
    end
end

MSE_DV = zeros(365,1);
for j = 1:365
    MSE_DV(j) = sum((DV_sim(:,j) - DV(:,j)).^2)/288;
end

figure
plot(MSE_DV)
grid on
legend('MSE DV')

%查看异常Day.k = 47,114,145,176,267,312
k = [47,114,145,176,267,312];
figure
for i = 1:length(k)
    subplot(2,3,i)
    plot(DV_sim(:,k(i)))
    hold on
    plot(DV(:,k(i)))
    grid on
    xlabel('t / 5min')
    ylabel('DV W/h')
    legend('sim','data')%,'PV',,,'EINE'
    title(['Tag.',num2str(k(i)),' MSE:',num2str(MSE_DV(k(i)))])
end

%% 绘图EINE-EVNB-SOC--Tagesverlauf
SOC = reshape(SOC,[288,365]);
V = reshape(Verbrauch,[288,365]);
PV = reshape(PV_Produktion,[288,365]);
EVNB = reshape(Energie_vom_Netz_bezogen,[288,365]);
EINE = reshape(Energie_ins_Netz_eingespeist,[288,365]);
EIBE = reshape(Energie_in_Batterie_gespeichert,[288,365]);
EABB = reshape(Energie_aus_Batterie_bezogen,[288,365]);
DV = reshape(Direkt_verbraucht,[288,365]);

k = 177;% Tag.47
figure
% plot(PV(:,k),'linewidth',2)
% hold on
% plot(V(:,k),'linewidth',2)
% hold on
% plot(EINE(:,k),'o','linewidth',2)
% hold on
% plot(EVNB(:,k),'linewidth',2)
% hold on
% plot(SOC(:,k),'linewidth',2)
% hold on
plot(EABB(:,k),'--','linewidth',2)
hold on
plot(EIBE(:,k),'-','linewidth',2)
grid on
xlabel('$ t/5min$','interpreter','latex', 'FontSize', 18)
ylabel('$ W/h  $','interpreter','latex', 'FontSize', 18)
legend('$ EABB_{data} $','$ EIBE_{data} $','interpreter','latex', 'FontSize', 20) %,'EINE','EVNB','SOC','EABB','EIBE'
title(['Tag. ',num2str(k)],'interpreter','latex', 'FontSize', 18)


%% 绘图月份柱状图
monate = [31 28 31 30 31 30 31 31 30 31 30 31];
% PV-Monate-Verlauf
PV_M = zeros(length(monate),1);
mm = zeros(length(monate),1);
for i = 1:length(monate)
    mm(i) = sum(monate(1:i));
    if i == 1
        PV_M(i) = sum(sum(PV(:,1:mm(i))));
    else
        PV_M(i) = sum(sum(PV(:,1+mm(i-1):mm(i))));
    end
end

PV_M = PV_M/1000;
figure
bar(PV_M)
xlabel('$ Monate $','interpreter','latex', 'FontSize', 20)
ylabel('$ PV/kWh $','interpreter','latex', 'FontSize', 20)


% V-Monate-Verlauf
V_M = zeros(length(monate),1);
for i = 1:length(monate)
    mm(i) = sum(monate(1:i));
    if i == 1
        V_M(i) = sum(sum(V(:,1:mm(i))));
    else
        V_M(i) = sum(sum(V(:,1+mm(i-1):mm(i))));
    end
end

V_M = V_M/1000;
figure
bar(V_M)
xlabel('$ Monate $','interpreter','latex', 'FontSize', 20)
ylabel('$ V/kWh $','interpreter','latex', 'FontSize', 20)

% EVNB-Monate-Verlauf
EVNB_M = zeros(length(monate),1);
for i = 1:length(monate)
    mm(i) = sum(monate(1:i));
    if i == 1
        EVNB_M(i) = sum(sum(EVNB(:,1:mm(i))));
    else
        EVNB_M(i) = sum(sum(EVNB(:,1+mm(i-1):mm(i))));
    end
end

EVNB_M = EVNB_M/1000;
figure
bar(EVNB_M)
xlabel('$ Monate $','interpreter','latex', 'FontSize', 20)
ylabel('$ EVNB/kWh $','interpreter','latex', 'FontSize', 20)

% EINE-Monate-Verlauf
EINE_M = zeros(length(monate),1);
for i = 1:length(monate)
    mm(i) = sum(monate(1:i));
    if i == 1
        EINE_M(i) = sum(sum(EINE(:,1:mm(i))));
    else
        EINE_M(i) = sum(sum(EINE(:,1+mm(i-1):mm(i))));
    end
end

EINE_M = EINE_M/1000;
figure
bar(EINE_M)
xlabel('$ Monate $','interpreter','latex', 'FontSize', 20)
ylabel('$ EINE/kWh $','interpreter','latex', 'FontSize', 20)


%% 绘图6月Tag柱状图: juni:Tag von 151:180
tag = 151:180;
Juni_PV = zeros(30,1);
Juni_V = zeros(30,1);
Juni_EINE = zeros(30,1);
Juni_EVNB = zeros(30,1);
Juni_EIBE = zeros(30,1);
Juni_EABB = zeros(30,1);

for i = 1:length(tag)
    Juni_PV(i) = sum(PV(:,tag(i)));
    Juni_V(i) = sum(V(:,tag(i)));
    Juni_EINE(i) = sum(EINE(:,tag(i)));
    Juni_EVNB(i) = sum(EVNB(:,tag(i)));
    Juni_EIBE(i) = sum(EIBE(:,tag(i)));
    Juni_EABB(i) = sum(EABB(:,tag(i)));
end

figure
bar([Juni_PV,Juni_V])
xlabel('$ Tag $','interpreter','latex', 'FontSize', 20)
ylabel('$ kWh $','interpreter','latex', 'FontSize', 20)
legend('PV','V','interpreter','latex', 'FontSize', 20)

figure
bar([Juni_EINE,Juni_EVNB])
xlabel('$ Tag $','interpreter','latex', 'FontSize', 20)
ylabel('$ kWh $','interpreter','latex', 'FontSize', 20)
legend('EINE','EVNB','interpreter','latex', 'FontSize', 20)

figure
bar([Juni_EIBE,Juni_EABB])
xlabel('$ Tag $','interpreter','latex', 'FontSize', 20)
ylabel('$ kWh $','interpreter','latex', 'FontSize', 20)
legend('EIBE','EABB','interpreter','latex', 'FontSize', 20)


% Juni-SOC-Verlauf in Stunden
Juni_SOC = SOC(:,tag);
Juni_SOC_tag = zeros(30,1);
for i = 1:30
    Juni_SOC_tag(i) = sum(Juni_SOC(:,i))/288;
end

figure
plot(Juni_SOC_tag,'-','linewidth',2)
grid on
xlabel('$ Tag $','interpreter','latex', 'FontSize', 20)
ylabel('$ SOC/ \% $','interpreter','latex', 'FontSize', 20)



%% 原始Profit
profit = zeros(365,1);
profit_blue = [];
profit_red = [];
for j =1:365
    profit(j) = -(sum(EVNB(:,j))*0.3 - sum(EINE(:,j))*0.078)/1000; % 1-Jahre-Profit
    if profit(j) >= 0
        profit_blue = [profit_blue,j];
    else
        profit_red = [profit_red,j];
    end
end
% profit in Juni
figure
plot(profit_blue,profit(profit_blue),'o')
hold on
plot(profit_red,profit(profit_red),'o')
grid on
xlabel('$ Tag $','interpreter','latex', 'FontSize', 20)
ylabel('Gewinn / €','FontSize', 20)


figure
plot(profit,'o')
hold on 
plot(profitS,'*')
grid on
xlabel('$ Tag $','interpreter','latex', 'FontSize', 20)
ylabel('$ Gewinn $','interpreter','latex', 'FontSize', 20)
legend('$ ohne \; Zusatzspeicher $','$ mit \; Zusatzspeicher \; 13,5 kWh $','interpreter','latex', 'FontSize', 14)







