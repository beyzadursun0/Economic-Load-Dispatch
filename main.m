clc;
clear;
global a b c B B0 B00 PD Pmin Pmax penalty_factor

PD = 1263;

Pmin = [100 50 80 50 50 50];
Pmax = [500 200 300 150 200 120];

a = [240 200 220 200 220 190];
b = [0.7 10 8.5 11 10.5 12];
c = [0.007 0.0095 0.009 0.009 0.008 0.0075];

B = [0.0017  0.0012  0.0007 -0.0001 -0.0005 -0.0002;
     0.0012  0.0011  0.0009  0.0001 -0.0006 -0.0001;
     0.0007  0.0009  0.0031  0.0000 -0.0010 -0.0006;
    -0.0001  0.0001  0.0000  0.0024 -0.0006 -0.0008;
    -0.0005 -0.0006 -0.0010 -0.0006  0.0129 -0.0002;
    -0.0002 -0.0001 -0.0006 -0.0008 -0.0002  0.0150];

B0 = [-0.3908 -0.1297 0.7047 0.0591 0.2161 -0.6635]*1e-3;
B00 = 0.056;

nVar = 6;            
VarSize = [1 nVar];  
MaxIt = 1000;          
nPop = 500;           
penalty_factor = 1e2;  

sol = repmat(struct('Position', [], 'Cost', []), nPop, 1);
for i=1:nPop
    start_pos = Pmin + (Pmax - Pmin).*rand(1,nVar);
    lower_bound = max(Pmin, start_pos*0.70);
    upper_bound = min(Pmax, start_pos*1.30);
    sol(i).Position = unifrnd(lower_bound, upper_bound);
    sol(i).Cost = CostFunction(sol(i).Position);
end

[~, best_idx] = min([sol.Cost]);
BestSol = sol(best_idx);

BestCost = zeros(MaxIt,1);
BestPositions = zeros(MaxIt, nVar);


for it=1:MaxIt
    for i=1:nPop
        j = randi([1 nPop]);
        while j == i
            j = randi([1 nPop]);
        end
        BF1 = randi([1 2]);
        MutualVector = (sol(i).Position + sol(j).Position)/2;
        NewSol1 = sol(i).Position + rand(1,nVar).*(BestSol.Position - BF1*MutualVector);
        NewSol2 = sol(j).Position + rand(1,nVar).*(BestSol.Position - BF1*MutualVector);
        NewSol1 = ApplyBounds(NewSol1, Pmin, Pmax);
        NewSol2 = ApplyBounds(NewSol2, Pmin, Pmax);
        NewCost1 = CostFunction(NewSol1);
        NewCost2 = CostFunction(NewSol2);
        if NewCost1 < sol(i).Cost
            sol(i).Position = NewSol1;
            sol(i).Cost = NewCost1;
        end
        if NewCost2 < sol(j).Cost
            sol(j).Position = NewSol2;
            sol(j).Cost = NewCost2;
        end
        j = randi([1 nPop]);
        while j == i
            j = randi([1 nPop]);
        end
        NewSol = sol(i).Position + (rand(1,nVar)*2 - 1).*(BestSol.Position - sol(j).Position);
        NewSol = ApplyBounds(NewSol, Pmin, Pmax);
        NewCost = CostFunction(NewSol);
        if NewCost < sol(i).Cost
            sol(i).Position = NewSol;
            sol(i).Cost = NewCost;
        end
        parasite = sol(i).Position;
        rand_index = randi(nVar);
        parasite(rand_index) = unifrnd(Pmin(rand_index), Pmax(rand_index));
        parasite = ApplyBounds(parasite, Pmin, Pmax);
        parasite_cost = CostFunction(parasite);
        j = randi([1 nPop]);
        while j == i
            j = randi([1 nPop]);
        end
        if parasite_cost < sol(j).Cost
            sol(j).Position = parasite;
            sol(j).Cost = parasite_cost;
        end
    end
    
    [~, best_idx] = min([sol.Cost]);
    BestSol = sol(best_idx);
    BestPositions(it,:) = BestSol.Position;
    BestCost(it) = BestSol.Cost;
    
    fprintf('Iterasyon: %d | En İyi Maliyet: %.4f\n', it, BestCost(it));
end

disp('Optimal Güç Dağılımı (MW):');
disp(BestSol.Position);

disp(['Toplam Maliyet (Yakıt): ' num2str(BestSol.Cost)]);

P_opt = BestSol.Position(:);
B0_col = B0(:);
PL_opt = P_opt' * B/100 * P_opt + B0_col' * P_opt + B00*100;
P_total = sum(P_opt);

fprintf('Toplam Güç Kaybı (MW): %.4f\n', PL_opt);
fprintf('Toplam Üretim Gücü (MW) (Üretim - Kayıp): %.4f\n', P_total - PL_opt);
fprintf('Talep Gücü (MW): %.4f\n', PD);

figure;
bar(BestSol.Position);
title('Jeneratörlerin Ürettiği Güç (MW)');
xlabel('Jeneratör No');
ylabel('Güç (MW)');
grid on;
ylim([0 max(Pmax)*1.1]);

figure;
plot(1:MaxIt, BestCost, 'LineWidth', 2);
title('İterasyonlara Göre Toplam Maliyet');
xlabel('İterasyon');
ylabel('Maliyet ($)');
grid on;

MMBTU = 3.412; 
TotalPower = zeros(MaxIt,1);
PowerLoss = zeros(MaxIt,1);
for k=1:MaxIt
    Pk = BestPositions(k,:)';
    PowerLoss(k) = Pk' * B/100 * Pk + B0_col' * Pk + B00*100;
    TotalPower(k) = sum(Pk);
end
BTU_per_iter = TotalPower * MMBTU;
CostPerBTU = BestCost ./ BTU_per_iter;

figure;
plot(1:MaxIt, CostPerBTU, 'm-', 'LineWidth', 2);
title('İterasyonlara Göre BTU Başına Maliyet ($/MMBTU)');
xlabel('İterasyon');
ylabel('Maliyet ($/MMBTU)');
grid on;

figure;
plot(1:MaxIt, TotalPower, 'b-', 'LineWidth', 2);
hold on;
plot(1:MaxIt, PowerLoss, 'r--', 'LineWidth', 2);
yline(PD,'k-.','Talep Gücü (PD)');
title('İterasyonlara Göre Üretim Gücü ve Güç Kayıpları (MW ve Maliyet Birimleri)');
xlabel('İterasyon');
ylabel('Güç (MW) / Maliyet ($, BTU)');
legend('Toplam Üretim Gücü (MW)','Toplam Güç Kaybı (MW)','Talep Gücü (MW)','Location','best');
grid on;
hold off;

Efficiency = (PD / P_total) * 100;
Regulation = ((P_total - PD) / PD) * 100;

fprintf('Sistem Verimi (%%): %.2f\n', Efficiency);
fprintf('Regülasyon Oranı (%%): %.2f\n', Regulation);

GenLossContribution = zeros(nVar,1);
for i = 1:nVar
    Pi = P_opt(i);
    Bi = B(i,:)/100;
    GenLossContribution(i) = Pi * sum(Bi .* P_opt');
end

disp('Her Jeneratörün Yaklaşık Güç Kaybı Katkısı (MW):');
for i = 1:nVar
    fprintf('Jeneratör %d: %.4f MW\n', i, GenLossContribution(i));
end

function cost = CostFunction(P)
    global a b c B B0 B00 PD penalty_factor
    P_col = P(:);
    B0_col = B0(:);
    PL = P_col' * B/100 * P_col + B0_col' * P_col + B00*100;
    penalty = penalty_factor * abs(sum(P) - (PD + PL));
    cost = sum(a + b.*P + c.*(P.^2)) + penalty;
end

function P = ApplyBounds(P, Pmin, Pmax)
    P = max(P, Pmin);
    P = min(P, Pmax);
end
