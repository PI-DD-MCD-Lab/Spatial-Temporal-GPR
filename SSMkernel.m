function Cov = SSMkernel(F,H,G,Sigma_0,m,T)
Cov = zeros(m,m);
        for i = 1:m
            for j = 1:m
                temp = 0;
                for k = 0:1:min(i,j)-1
                    temp = temp + H*F^(i-1-k)*(G*G')*(F^(j-1-k))'*H';
                end
                Cov(i,j) = H*(F^i)*Sigma_0*(F^j)'*H' + temp;
            end
        end