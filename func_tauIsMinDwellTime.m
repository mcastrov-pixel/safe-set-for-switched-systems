





function flag_out = func_tauIsMinDwellTime(cell_A,tauTest,verbose)
% solves the LMI in J. C. Geromel and P. Colaneri, “Stability and stabilization of discrete
% time switched systems,” International Journal of Control, 2006 to see if
% tauTest is suitable as a min dwell time.

%feasibility problem is solved using Yalmip

if nargin<2 ||isempty(verbose)
	verbose = 0
end
% ---------- SETUP ----------
M = numel(cell_A);
n = size(cell_A{1},1);

for it = 2:M
	assert(all(size(cell_A{it}) == [n n]), 'All A{i} must be n-by-n.');
end

% Precompute A_i^D
Atau_cell = cell(M,1);
for it = 1:M
	Atau_cell{it} = cell_A{it}^tauTest;   % mpower
end


%form the LMI

% Decision variables P_i
P = cell(M,1);
Constraints = [];
eps_margin = 1e-6;     % strictness margin (since solvers handle non-strict)
for it = 1:M
	P{it} = sdpvar(n,n,'symmetric');
	Constraints = [Constraints, P{it} >= eps_margin*eye(n)];
	Constraints = [Constraints, cell_A{it}'*P{it}*cell_A{it} - P{it} <= -eps_margin*eye(n)];
end

% Build the (i,j) set
	[I,J] = ndgrid(1:M, 1:M);
	pairs = [I(:) J(:)];
% Cross constraints: (A_i^D)' P_j (A_i^D) - P_i < 0
for k = 1:size(pairs,1)
	it = pairs(k,1);
	j = pairs(k,2);
	Constraints = [Constraints, Atau_cell{it}'*P{j}*Atau_cell{it} - P{it} <= -eps_margin*eye(n)];
end

% ---------- SOLVE (feasibility) ----------
opts = sdpsettings('solver','mosek','verbose',0);   % or 'sdpt3' / 'sedumi'
sol = optimize(Constraints, [], opts);

% ---------- RESULTS ----------
if sol.problem ~= 0
	if verbose
		disp(sol.info);
		sprintf('Infeasible or solver error (code=%d).', sol.problem);
	end
	flag_out = 0;
end
flag_out = 1;
Psol = cellfun(@value, P, 'UniformOutput', false);
if verbose
	disp('Feasible. Smallest eigenvalue of each P_i:');
	for it = 1:M
		fprintf('i=%d: min eig(P_i) = %.3e\n', it, min(eig(Psol{it})));
	end
end
end