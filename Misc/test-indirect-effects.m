%# Compute marginal effects for SAR probit
randn('seed',0)

d = 10;
n = 10;
I_n = speye(n);

%intercept = ones(n, 1);
%x1 = normrnd(0, 1, n, 1);
%x2 = normrnd(0, 1, n, 1);
%X = [intercept x1 x2 ];
X = [intercept (-2:.2:-0.1)' (3:-0.5:-1.5)']


W = [   0   1   0   0   0   0   0   0   0   0;
      0.5   0 0.5   0   0   0   0   0   0   0;
        0 0.5   0 0.5   0   0   0   0   0   0;
        0   0 0.5   0 0.5   0   0   0   0   0;
        0   0   0 0.5   0 0.5   0   0   0   0;
        0   0   0   0 0.5   0 0.5   0   0   0;
        0   0   0   0   0 0.5   0 0.5   0   0;
        0   0   0   0   0   0 0.5   0 0.5   0;
        0   0   0   0   0   0   0 0.5   0 0.5;
        0   0   0   0   0   0   0   0   1   0;
    ];

beta = [ 0;  1;  -1 ];  % (3 x 1)
rho = 0.75;

%% compute traces for W^i
iiter=10000;
o=100;

diag_ests=zeros(n,o);  % SW: (n x o)
for iii=1:iiter

u=randn(n,1);

umat=u(:,ones(1,o));
wumat=zeros(n,o);
wu=u;
wumat(:,1)=wu;
for ii=2:o
    wu=W*wu;
    wumat(:,ii)=wu;
end
   
diag_estimates_iii=(umat.*wumat);

diag_ests=diag_ests+diag_estimates_iii;
end

estimated_diags=diag_ests/iiter;

estimated_diags(1:5,1:5)

%%
S = I_n - rho.*W;
mu = S\X*beta;
s = S\speye(n);

%% compute marginal effects
beff=[1;  -1 ];  % (2 x 1)
avg_total    = zeros(2,1);
avg_direct   = zeros(2,1);
avg_indirect = zeros(2,1);

rhovec=(rho.^(0:o-1))';
pdfz = normpdf(mu,0,1);
for kk=1:length(beff);
   avg_direct(kk,1) = (pdfz'*estimated_diags*rhovec)*beff(kk,1)/n; %av direct effect
   dd = spdiags(pdfz, 0, n, n);
   avg_total(kk,1) = mean(sum(dd*s*beff(kk,1),2));
   %dd = spdiags(pdfzs,0,n,n);
   %totals = sum(dd*ss*beff(kk,1),2); 
end;
avg_indirect = avg_total' - avg_direct';

avg_direct
avg_indirect
avg_total

%avg_direct =
%    0.0898
%   -0.0898
%avg_indirect
%avg_total
%avg_indirect =
%    0.1237   -0.1237
%avg_total =
%    0.2136
%   -0.2136



