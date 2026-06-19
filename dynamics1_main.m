clear all



% parameters
nx = 2;
% dynamics
cell_A{1} = [0.8 0 ;5 0.8];
cell_A{2} = [0.5 0;0 0.56];
% Constraint set 
H = kron(eye(2),[-1;1]);
h =[0.2;0.2;1;1];


%% Compute min dwell time
% look for a good minimum dwell time using the LMI in 
% in J. C. Geromel and P. Colaneri, “Stability and stabilization of discrete
% time switched systems,” International Journal of Control, 2006 
%
% simple for loop this could be changed to a bisection type of problem for
% faster convergence
tau_val = [];
for tau_test = 1:20
	flag_out = func_tauIsMinDwellTime(cell_A,tau_test,0);
	if flag_out
		tau_val = tau_test;
		sprintf("%d is a valid minimum dwell time", tau_val)
		break
	end
end

if isempty(tau_val)
	sprintf("valid min dwell time not found, assigning the value of %d", tau_test)
	tau_val = tau_test;
end


%% check that the usual stopping criterion might not work: (Figure 1 in paper)
% we consider the sequence {2,1,1,1,1,1,1...} for example
% we note the following O_{2} = O_{2,1} = O_{2,2} so we might think we are
% done but it turns out that O_{2,1,1} is a subset of O_{2} showing that
% the moas of the sequence under consideration is not reached

AMat = {};bMat = {};
AMat{1} = [H;H*cell_A{2}];
bMat{1} = [h;h];
AMat{2} = [AMat{1};H*cell_A{2}*cell_A{2}];bMat{2} = [bMat{1};h];
AMat{3} = [AMat{1};H*cell_A{1}*cell_A{2}];bMat{3} = [bMat{1};h];
AMat{4} = [AMat{1};H*cell_A{1}^2*cell_A{2}];bMat{4} = [bMat{1};h];

figure(1)
temPopo = {};
markerStyles={'o','+','*'};
markerStyles={'o','+','*'};
lineWidths = {8,5,2,.75};
for it =1:length(AMat)
	
		temPopo{it} = polytope(AMat{it},bMat{it});
		if it == length(AMat)
			temPopo{it}.plot([1,2],'LineWidth',1.5,'LineStyle','-.');%,'Color',[0 0 0.05]); hold on;
		else
			temPopo{it}.plot([1,2],'LineWidth',lineWidths{it}); hold on;
		end
end
xlim([-1.1*h(1),1.1*h(1)]);ylim([-1.1*h(3),1.1*h(3)]);
	hold off;
legend({'$\mathcal O_{\{2\}}$','$\mathcal O_{\{2,2\}}$','$\mathcal O_{\{2,1\}}$','$\mathcal O_{\{2,1,1\}}$',},'Interpreter','latex','FontSize',16,'Location','best');
xlabel('$ x_1$','Interpreter','latex','FontSize',14);
ylabel('$ x_2$','Interpreter','latex','FontSize',14);



%% compute the MOAS for different matrices (just for fun)
cell_seqs = {};
cell_MOAS_A = {};
cell_MOAS_b ={};
cell_MOAS_it = {};
cell_legend ={};

% generate desired sequences

cell_seqs{1} = 1;
cell_seqs{2} = 2;
cell_legend = {'$\mathcal O_{\{1\}^\infty}$','$\mathcal O_{\{2\}^\infty}$'};

% switching sequence: {1,2,1,2,1,2,2,2,...}
nSwitches = 6;
switchSeq = [];
for it =1:floor(nSwitches)/tau_val
	if mod(it,2) == 0
	switchSeq = [switchSeq ones(1,tau_val)*2];
	else
		switchSeq = [switchSeq ones(1,tau_val)];
	end
end
cell_seqs{3} = switchSeq;
cell_legend{end+1} = sprintf('$\\mathcal O_{\\{2,1\\}^%d\\cap\\{2\\}^\\infty}$',nSwitches);

% switching sequence: {2,1,2,1,2,1,1,1,...}
switchSeq = [];
for it =1:floor(nSwitches/tau_val)
	if mod(it,2) == 0
	switchSeq = [switchSeq ones(1,tau_val)];
	else
		switchSeq = [switchSeq ones(1,tau_val)*2];
	end
end
cell_seqs{4} = switchSeq;
cell_legend{end+1} = sprintf('$\\mathcal O_{\\{1,2\\}^%d\\cap\\{1\\}^\\infty}$',nSwitches);

%a random switching sequence
switchSeq = [];
for it =1:4
	nOnes = floor(4*rand(1));
	switchSeq = [switchSeq 2 ones(1,nOnes)];
end
cell_seqs{5} = switchSeq;
% get the name for this sequence
seqStr = '$\mathcal O_{';
compactSeqDescription = diff(find(diff([0 switchSeq 0]) ~=0));
for it =1:length(compactSeqDescription)-1;
	if mod(it,2)==1
		seqStr = strcat(seqStr,sprintf('\\{2\\}^%d\\cap',compactSeqDescription(it)));
	else
		seqStr = strcat(seqStr,sprintf('\\{1\\}^%d\\cap',compactSeqDescription(it)));
	end
end
if mod(it+1,2)==1
	seqStr = strcat(seqStr,'\{2\}^\infty}$');
else
	seqStr = strcat(seqStr,'\{1\}^\infty}$');
end
cell_legend{end+1} = seqStr; 


% compute MOASes

for it =1:length(cell_seqs)
	[cell_MOAS_A{it},cell_MOAS_b{it},cell_MOAS_it{it}] = func_MOAS_swichingSeq(H,h,cell_A,cell_seqs{it});
end

%plot all MOAS
figure(101)

for it =1:length(cell_MOAS_A)
	temPopo = polytope(cell_MOAS_A{it},cell_MOAS_b{it});
	temPopo.plot([1,2],'LineWidth',1.5);hold on
	% legendVals{it} = strcat("$O_{seq ",num2str(it),"}$");
end
legend(cell_legend,'interpreter','latex','fontsize',18)
hold off







