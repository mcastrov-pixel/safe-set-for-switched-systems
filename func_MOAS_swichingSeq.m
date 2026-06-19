function [A_infty,b_infty,it]= func_MOAS_swichingSeq(H,h,ACell,SwitchSequence,epsVal)
% computes moas for a given sequence. Assumes the sequence to be composed
% of a finite number of switches between different modes before  becoming
% constant (given by the last element of the sequence)
% the usual tewrmination criteria is used once we reach the end of the
% sequence.

% !! No tightening is  performed, the A matrices need to be Schur!!!

% constraints are given as H x<= h
%origin should be in the interior (h>0)

if nargin<5
	epsVal = 1e-4;
end
nx = length(ACell{1});
nSwitches = length(SwitchSequence);
nMax = 2e3+nSwitches;
nElimRedundant = 30;
%initiate the MOAS
A_infty = H;
b_infty = h;%*(1-epsVal);
%dynamics propagation given the sequence
ASeq = eye(nx);
flg_termination = 0;

for it =1:nMax
	if it<= nSwitches
		curMod = SwitchSequence(it);
	end
	
	%propagate one step forward
	ASeq = ACell{curMod}*ASeq;

	%new constraints
	ACstr_nw = H*ASeq; 
	bCstr_nw = h;

	if it>nSwitches
	flg_termination = isIncluded(ACstr_nw,bCstr_nw,A_infty,b_infty,epsVal);
	end

	if flg_termination
		break;
	end
	A_infty = [A_infty;ACstr_nw];
	b_infty = [b_infty;bCstr_nw];

	%eliminate almost redundant constraints
	if mod(it,nElimRedundant) == 0
		[A_infty,b_infty,~]=elimRedundant(A_infty,b_infty,epsVal);
	end

	if it == nMax
		sprintf("reached the max number of iterations without fulfilling the stopping criterion \n !!!!!!!!!!!!!!!!!!!!!!!! output might not be the MOAS!!!!!!!!!!!!!!!!!!!!!!!!")
	end
end
[A_infty,b_infty,~]=elimRedundant(A_infty,b_infty,epsVal);

end



function [A1,b1,IDX]=elimRedundant(A,b,epsVal)
%eliminate redundant constraints from a polytope in H form

if nargin==2
	epsVal=1e-6;
end

A1=A;
b1=b;
sizeA=size(A);

n=sizeA(1);
linProg_opts = optimoptions('linprog','Display','off');
IDX = [2:n];
n_del = 0;
for it=1:1:n
	[~,fVal]=linprog(-A(it,:),A(IDX,:),b(IDX),[],[],[],[],linProg_opts);
	IDX(it-n_del) =  it;
	if -fVal<b(it)+epsVal
		IDX(it-n_del) = [];
		n_del = n_del+1;
	end
	if n_del+1 == n
		IDX = n;
		break; disp('Warning: only one active constraint remains');
	end
end

A1 = A(IDX,:);
b1 = b(IDX);

end



function issub=isIncluded(A,b,A1,b1,epsVal)
% Checks whether the set of inequalities defined by A1,b1 is included in
% A,b issub takes values 'y' for yes or 'n' for no
% The parameter epsVal specifies the tolerance in determining set
% inclusion; 

if nargin<5
	epsVal=1e-8;
end


szA1=size(A1);
m=szA1(1); % number of inequalitites which define A1

issub=1; %
linProg_opts = optimoptions('linprog','Display','off');

for it=1:1:m
	f=-A1(it,:);   %select the ith constraint, i=1,...,m, of A1, b1
	[~,fVal] = linprog(f,A,b,[],[],[],[],linProg_opts); %solve an lp problem to determine if redundant
	hv = -fVal-b1(it);
	if hv>=epsVal
		issub=0; break,
	end % if not redundant stop

end

end

