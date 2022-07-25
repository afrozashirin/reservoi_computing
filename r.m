function rdot = r(t,r,flag,p1,p2,p3,A,W,meanX,stdX,V,n)


rho=28;
sigma=10;
beta=8/3;
z = 5;

rdot = zeros(n+3,1); % you did set some of the columns/values to zero. I dont know why. you should set all to zero

rdot(1) = (sigma*(r(2) - r(1)));
rdot(2) = (r(1)*(rho - r(3)) - r(2));
rdot(3) = (r(1)*r(2) - beta*r(3));

X = (r(1)-meanX)/stdX; % there was an absolute value sign here which was very wrong.

drive_sig = X;
gamma = p2;    

%linear
%rdot(4:n+3) = gamma*(-r(4:n+3) + (A*r(4:n+3) + drive_sig.*W));

%%Nonlinear
rdot(4:n+3) = gamma*(-r(4:n+3) + tanh(A*r(4:n+3) + drive_sig.*W));




end