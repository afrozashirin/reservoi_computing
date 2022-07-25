function Ldot = Lorenz(t,L)

rho=28;
sigma=10;
beta=8/3;

Ldot = zeros(3,1);% you should set to zero

Ldot(1) = sigma*(L(2) - L(1));
Ldot(2) = L(1)*(rho - L(3)) - L(2);
Ldot(3) = L(1)*L(2) - beta*L(3);

%Ldot = Ldot*0.02;
end

