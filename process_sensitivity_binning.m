% calculated process sensitivity through binning algorithm of case
% 
function [S_cA_bin] = Process_sensitivity_B_c(c_t)
for iii = 1:100
    clearvars a_index c_total var_t_c E_tb_c E_b_c; 
    N = 100*iii;
    MB(1) = 0.5;
    MB(2) = 0.5;
    MA(1) = 0.5;
    MA(2) = 0.5;
    bin = 10;
    for i = 1:21
        c_total(i,:) = reshape(c_t(i,:,1:N,:),[4*N,1]);
        var_t_c(i) = nanvar(c_total(i,:));
    end
    clearvars c_total;
    
    for i = 1:2
        if i == 1
            pd=makedist('Normal',3.35,1);
            y = linspace(0.01,0.99,11);
            abin = icdf(pd,y);
        else
            abin = linspace(min(a(i,:)),max(a(i,:)),bin +1 );
        end
        for m = 2:bin+1
            mm(m-1) = 0;
            for j = 1:N
                if a(i,j)>=abin(m-1) && a(i,j)<=abin(m)
                    mm(m-1) = mm(m-1) + 1;
                    a_index{i,m-1}(mm(m-1)) = j;
                end
            end
        end
    end
    
    for i = 1:21
        for j = 1:2
            for k = 1:bin
                for l = 1:2
                    E_tb_c(i,j,k,l) = nanmean(c_t(i,j,a_index{j,k},l));
                end
                E_b_c(i,j,k) = MB(1)*E_tb_c(i,j,k,1) + MB(2)*E_tb_c(i,j,k,2);
            end
            E_ta_c(i,j) = nanmean(E_b_c(i,j,:));
            E_ta_c2(i,j) = nanmean(E_b_c(i,j,:).^2);
        end
        E_a_c(i) = MA(1)*E_ta_c(i,1) + MA(2)*E_ta_c(i,2);
        E_a_c2(i) = MA(1)*E_ta_c2(i,1) + MA(2)*E_ta_c2(i,2);
        Var_cA(i) = E_a_c2(i) - (E_a_c(i))^2;
        S_cA_bin(iii,i) = Var_cA(i)/var_t_c(i);
    end
end
