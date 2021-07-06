function [S_cB] = Process_sensitivity_B_c(c_t)

MB(1) = 0.5;
MB(2) = 0.5;
MA(1) = 0.5;
MA(2) = 0.5;
for i = 1:21   
   c_total(i,:) = reshape(c_t(i,:,:,:),[8*N,1]);
   var_t_c(i) = nanvar(c_total(i,:));
   E_t_c(i) = nanmean(c_total(i,:));
end

for jj = 1:21
    for i = 1:2
        for j = 1:2
            for r = 1:N
                if j == 1
                    Fijj_c(jj,i,j,r) = c_t(jj,i,r,j).*c_t(jj,i,r,3);
                else
                    Fijj_c(jj,i,j,r) = c_t(jj,i,r,j).*c_t(jj,i,r,4);
                end
            end
              MFijj_c(jj,i,j) = nanmean(Fijj_c(jj,i,j,:));
        end
    end
    SFijjj_c(jj) = MFijj_c(jj,1,1)+MFijj_c(jj,1,2)+MFijj_c(jj,2,1)+MFijj_c(jj,2,2); % eqq
end
V1_c = SFijjj_c.*(0.5^3);
for jj = 1:21
    for i = 1:2
        for j = 1:2
            for r = 1:N
                if j == 1
                    Fijl_c(jj,i,j,r) = c_t(jj,i,r,j).*c_t(jj,i,r,2);
                else
                    Fijl_c(jj,i,j,r) = c_t(jj,i,r,j).*c_t(jj,i,r,1);
                end
            end
              MFijl_c(jj,i,j) = nanmean(Fijl_c(jj,i,j,:));

        end
            SFijl_c(jj,i) = MFijl_c(jj,i,1)+MFijl_c(jj,i,2);
    end
    SFijjl_c(jj) = SFijl_c(jj,1)+SFijl_c(jj,2);  
end
V2_c = SFijjl_c.*0.5^3;
E2c = E_t_c.^2;
Var_cB = V1_c+V2_c-E2c;
S_cB = Var_cB./var_t_c;
