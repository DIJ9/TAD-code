V = 7.69;
alpha = 12.2*pi/180;
% mu = V*cos(alpha)/(omega*R);
mu = 0.20227;
Ttotal = Th/cos(alpha);
Ct = Ttotal/(rho*A*Vtip^2);
lambda = sqrt(Ct/2);
lambda_old = 10*lambda;
iter = 0;
while abs(lambda_old-lambda) > 10^-3 && lambda ~= 0
    lambda_old = lambda;
    lambda = lambda - (lambda-mu*tan(alpha)-Ct/2/sqrt(mu^2+lambda^2))...
             /(1+Ct/2*lambda/((mu^2+lambda^2)^(3/2)));
    iter = iter + 1;
end

disp(lambda)
