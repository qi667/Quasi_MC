function [S_cB] = Process_sensitivity_B_c(c_t)
% process B sensitivity index
N = 300;
nx = 11;
MB(1) = 0.5;
MB(2) = 0.5;
MA(1) = 0.5;
MA(2) = 0.5;
for i = 1:21
     c_total(i,:) = reshape(c_t(i,:,:,:,:),[4*N^2,1]);
     var_t_c(i) = nanmean(c_total(i,:).^2) - (nanmean(c_total(i,:)))^2;
end
clearvars c_total;

for i = 1:21
    for j = 1:2
        for k = 1:N
            for l = 1:2
                E_ta_c(i,l,j,k) = nanmean(c_t(i,l,:,j,k));
            end
            E_a_c(i,j,k) = MA(1)*E_ta_c(i,l,j,k) + MA(2)*E_ta_c(i,l,j,k);
        end
        E_tb_c(i,j) = nanmean(E_a_c(i,j,:));
        E_tb_c2(i,j) = nanmean(E_a_c(i,j,:).^2);
    end
    E_b_c(i) = MB(1)*E_tb_c(i,1) + MB(2)*E_tb_c(i,2);
    E_b_c2(i) = MB(1)*E_tb_c2(i,1) + MB(2)*E_tb_c2(i,2);
    Var_cB(i) = E_b_c2(i) - (E_b_c(i))^2;
    S_cB(i) = Var_cB(i)/var_t_c(i);
end

