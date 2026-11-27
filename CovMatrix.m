function Cov = CovMatrix(x,y,kernel)

if sum(double(size(x)==[1,2]))==2
    m = 1;
else
    m = length(x);
end
if sum(double(size(y)==[1,2]))==2
    n = 1;
else
    n = length(y);
end

kernelfunc = kernelFunctionNEW(kernel);
Cov = zeros(m,n);
for i = 1:m
    for j = 1:n
        Cov(i,j) = kernelfunc(x(i,:),y(j,:));
    end
end
% CovD = zeros(m,n);
% for i = 1:m
%         CovD(i,i) = exp(kernel.alt*x(i,4));
% end
% Cov = diag(diag(Cov).*diag(CovD));
% Cov = Cov + Cov' - diag(diag(Cov));