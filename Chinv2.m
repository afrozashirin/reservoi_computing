function Ominv = Chinv2(Om,beta)
n = length(Om(1,:));
M = Om'*Om+beta*eye(n);
%M = network_out'*network_out;
L = chol(M,'lower');
U = L\eye(n);
Minv = L'\U;
Ominv = Minv*Om';
end