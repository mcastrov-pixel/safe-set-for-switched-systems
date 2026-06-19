clear all
% parameters
nx = 2;
% dynamics
% these are such that the MOAS increases as we include more swithces!
cell_A{1} = [0.9 .65;0 0.9 ];
cell_A{2} = [0.9 -.65;0 .9 ];
% Constraint set
H = kron(eye(2),[-1;1]);
h =[1;1;1;1];


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

%% compute the MOAS for different matrices (Figure 2 in paper)
cell_seqs = {};
cell_MOAS_A = {};
cell_MOAS_b ={};
cell_MOAS_it = {};
cell_legend = {};

% generate desired sequences

cell_seqs{1} = 1;
cell_seqs{2} = 2;

cell_legend = {'$\mathcal O_{\{1\}^\infty}$','$\mathcal O_{\{2\}^\infty}$'};

% switching sequence: {1,2,1,2,1,2,2,2,...}
for it =[2 4 6]
	nSwitches = it;
	cell_seqs{end+1} = kron(ones(1,nSwitches),kron([1 2],ones(1,tau_val)));
	cell_seqs{end+1} = kron(ones(1,nSwitches),kron([2 1],ones(1,tau_val)));
	cell_legend{end+1} = sprintf('$\\mathcal O_{\\{1,2\\}^%d\\cap\\{2\\}^\\infty}$',it);
	cell_legend{end+1} = sprintf('$\\mathcal O_{\\{2,1\\}^%d\\cap\\{2\\}^\\infty}$',it);
end
cell_legend{end+1} = '$\mathcal O_{\mathcal S}$';

% compute MOASes
for it =1:length(cell_seqs)
	[cell_MOAS_A{it},cell_MOAS_b{it},cell_MOAS_it{it}] = func_MOAS_swichingSeq(H,h,cell_A,cell_seqs{it});
end

% plot properties
[lineStyles, colors, lineColors, lineWidths] = plotCparamsAndColors();
orderToPlot = ([1 2 3 4 5 6 7 8]);

%plot all MOAS
for zt = 1:length(cell_MOAS_A)/2
	figure(20+zt);clf
	temPopo = polytope(H,h);
	temPopo.plot([1,2],'LineWidth',7.2,'Color',[0, 0, 128]/255,'LineStyle','-');
	hold on;
	for jt =1:2
		it = orderToPlot(jt+(zt-1)*2);
		it
		temPopo = polytope(cell_MOAS_A{it},cell_MOAS_b{it});
		temPopo.plot([1,2],'LineWidth',lineWidths(4-jt),'Color',lineColors(jt,:),'LineStyle',lineStyles(7),'FaceColor',[246, 38, 129]/255,'FaceAlpha',.18);
	end
	ylim([-1.05 1.05]);yticks([-1 -0.5 0 0.5 1]);
	xlim([-1.05 1.05]);xticks([-1 -0.5 0 0.5 1]);

	set(gca, 'FontSize', 30)

	xlabel('$ x_1$','Interpreter','latex','FontSize',34);
	ylabel('$ x_2$','Interpreter','latex','FontSize',34);
	hold off
end


hold off

%% control example considering only stable modes (Figure 4 in paper)

x0Vals = [-1 1 1 -1 1 1 -1 -1; -1 1 -1 1 -0.53 0.53 -0.53 0.53];
NSim = 100;
% Here we assume \mathcal S = "the last two sequences".
cell_OS_A = cell_MOAS_A([7 8]);
cell_OS_b = cell_MOAS_b([7 8]);


curFig = figure(4);clf;
for subfig = 1:2
	subplot(1,2,subfig)
	%plot the O_{infty} (we know it is equal to the constraint set from previous section)
	temPopo = polytope(H,h);
	temPopo.plot([1,2],'LineStyle','none','LineWidth',2,'FaceColor',[250 218 221]/250);hold on;
	% plot the {1}^\infty MOAS (because one of the cost functions tries to
	% apply only the first mode)
	temPopo = polytope(cell_MOAS_A{1},cell_MOAS_b{1});
	temPopo.plot([1,2],'LineStyle','-','Color',[1 1 1]*0,'LineWidth',2,'FaceColor',[0.2 0.2 0.2]*3)

end

% run simulations with different objective functions in the OCP defined in
% (13)


cell_sigma = cell(2,size(x0Vals,2));
cell_xHist = cell(2,size(x0Vals,2));

for mt = 1:size(x0Vals,2);
	x0 = x0Vals(:,mt);


	% run simulation considering both costs in parallel
	xHist1 = zeros(size(x0,1),NSim+1);
	xHist2 = zeros(size(x0,1),NSim+1);
	sigmaHist1 = zeros(1,NSim);
	sigmaHist2 = zeros(1,NSim);

	xHist1(:,1) = x0;
	xHist2(:,1) = x0;


	flg_ok1 = 0;
	for it = 1:NSim

		[sigmaCur1,flg_ok1] = solveOCP_cost1(xHist1(:,it),cell_A,cell_OS_A,cell_OS_b);
		[sigmaOut2,flg_ok2] = solveOCP_cost2(xHist2(:,it),cell_A,cell_OS_A,cell_OS_b);
		
		if flg_ok1
			sigmaHist1(it) = sigmaCur1;
			xHist1(:,it+1) = cell_A{sigmaHist1(it)}*xHist1(:,it);
		end

		if flg_ok2
			sigmaHist2(it) = sigmaOut2;
			xHist2(:,it+1) = cell_A{sigmaHist2(it)}*xHist2(:,it);
		end
		if ~flg_ok1 || ~flg_ok2
			keyboard
		end

	end
	

	subplot(1,2,1)
	plot(xHist1(1,[sigmaHist1 ==1 false]),xHist1(2,[sigmaHist1 ==1 false]),'r<','LineWidth',1)
	plot(xHist1(1,[sigmaHist1 ==2 false]),xHist1(2,[sigmaHist1 ==2 false]),'k<','LineWidth',1)
	plot(x0(1,1),x0(2,1),'bo','LineWidth',5,'MarkerSize',6)
	plot(xHist1(1,:),xHist1(2,:),'Color',lineColors(6,:),'LineWidth',1.5);

	subplot(1,2,2)
	if mt ==1
		plot(xHist2(1,:),xHist2(2,:),'Color',lineColors(6,:),'LineWidth',1.5);
		plot(xHist2(1,[sigmaHist2 ==1 false]),xHist2(2,[sigmaHist2 ==1 false]),'r<','LineWidth',1)
		plot(xHist2(1,[sigmaHist2 ==2 false]),xHist2(2,[sigmaHist2 ==2 false]),'k<','LineWidth',1)
		plot(xHist2(1,1),xHist2(2,1),'bo','LineWidth',2.5,'MarkerSize',3)
	end
	plot(xHist2(1,[sigmaHist2 ==1 false]),xHist2(2,[sigmaHist2 ==1 false]),'r>','LineWidth',1)
	plot(xHist2(1,[sigmaHist2 ==2 false]),xHist2(2,[sigmaHist2 ==2 false]),'k>','LineWidth',1)
	plot(x0(1,1),x0(2,1),'bo','LineWidth',5,'MarkerSize',6)
	plot(xHist2(1,:),xHist2(2,:),'Color',lineColors(6,:),'LineWidth',1.5)

	cell_sigma{1,mt} = sigmaHist1;
	cell_sigma{2,mt} = sigmaHist2;
	cell_xHist{1,mt} = xHist1;
	cell_xHist{2,mt} = xHist2;
end
curAx = findobj(gcf, 'Type', 'axes');
babies = curAx(1).Children;

legend(curAx(1),flip(babies(end-5:end)),{'$\mathcal O_{\mathcal S}$','$\mathcal O_{\{1\}^\infty}$',...'$\mathcal O_{\{2\}^\infty}$', ...
	'traj','$\sigma_k = 1$','$\sigma_k = 2$','$I.C.$'},'Interpreter','latex','FontSize',13.5,'Location','best','NumColumns',2);hold off;
xlabel(curAx(1),'$ x_1$','Interpreter','latex','FontSize',14);
ylabel(curAx(1),'$ x_2$','Interpreter','latex','FontSize',14);
xlabel(curAx(2),'$ x_1$','Interpreter','latex','FontSize',14);
ylabel(curAx(2),'$ x_2$','Interpreter','latex','FontSize',14);

xlim([-1,1])
ylim([-1,1])
title(curAx(1),'$J(\sigma,k) = \|A_{\sigma}x_k\|^2$','Interpreter','latex','FontSize',16);
title(curAx(2),'$J(\sigma,k) = \sigma$','Interpreter','latex','FontSize',16);


%% control example considering also unstable modes (Figure 5 in paper)

x0Vals = [-1 1 1 -1 1 1 -1 -1; -1 1 -1 1 -0.53 0.53 -0.53 0.53];
NSim = 100;
% Here we assume \mathcal S = "the last two sequences".
cell_OS_A = cell_MOAS_A([7 8]);
cell_OS_b = cell_MOAS_b([7 8]);

cell_A_aug = cell_A;
cell_A_aug{end+1} = [    1.4000    0.1700;
    0.060    0.8100]; 
% cell_A_aug = flip(cell_A_aug);

curFig = figure(5);clf;

% run simulations with different objective functions in the OCP defined in
% (13)

cell_sigma = cell(2,size(x0Vals,2));
cell_xHist = cell(2,size(x0Vals,2));

for mt = 1:size(x0Vals,2);
	x0 = x0Vals(:,mt);


	% run simulation considering both costs in parallel
	xHist1 = zeros(size(x0,1),NSim+1);
	xHist2 = zeros(size(x0,1),NSim+1);
	sigmaHist1 = zeros(1,NSim);
	sigmaHist2 = zeros(1,NSim);

	xHist1(:,1) = x0;
	xHist2(:,1) = x0;


	flg_ok1 = 0;
	for it = 1:NSim

		[sigmaCur1,flg_ok1] = solveOCP_cost2(xHist1(:,it),cell_A,cell_OS_A,cell_OS_b);
		[sigmaOut2,flg_ok2] = solveOCP_cost2(xHist2(:,it),cell_A_aug,cell_OS_A,cell_OS_b);
		
		if flg_ok1
			sigmaHist1(it) = sigmaCur1;
			xHist1(:,it+1) = cell_A{sigmaHist1(it)}*xHist1(:,it);
		end

		if flg_ok2
			sigmaHist2(it) = sigmaOut2;
			xHist2(:,it+1) = cell_A_aug{sigmaHist2(it)}*xHist2(:,it);
		end
		if ~flg_ok1 || ~flg_ok2
			keyboard
		end

	end
	
plot(vecnorm(xHist1,2,1),'c','LineWidth',1.5);hold on;plot(vecnorm(xHist2,2,1),'m-.','LineWidth',1.5);

	cell_sigma{1,mt} = sigmaHist1;
	cell_sigma{2,mt} = sigmaHist2;
	cell_xHist{1,mt} = xHist1;
	cell_xHist{2,mt} = xHist2;

end
figure(34);
legend({'$\Sigma = \{1,2\}$','$\Sigma = \{1,2,3\}$'},'Interpreter','latex','FontSize',14,'Location','best');hold off;
ylabel('$\|x_k\|_2$','Interpreter','latex','FontSize',14);
xlabel('$k$','Interpreter','latex','FontSize',14);




%% functions



function [sigmaOut,flg_ok] = solveOCP_cost1(xCur,cell_A,cell_OS_A,cell_OS_b)
% cost considered here is J(sigma,x) = sigma --> try to prioritize keeping
% mode 1
sigmaOut = [];
flg_ok = 0;
for sigmaCand =1:length(cell_A)
	ACand = cell_A{sigmaCand};
	for jt = 1:length(cell_OS_A)
		if max(cell_OS_A{jt}*ACand*xCur-cell_OS_b{jt})<=0
			flg_ok = 1;
			break
		end
	end
	if flg_ok
		break
	end
end
sigmaOut = sigmaCand;
end


function [sigmaOut,flg_ok] = solveOCP_cost2(xCur,cell_A,cell_OS_A,cell_OS_b)
% cost considered here is the norm of 1 step forward J(sigma,x) = ||A_sigma xCur||^2 
flg_ok = 0;
sigmaOut=[];
closest_d = Inf;
for sigmaCur =1:length(cell_A)
	ACand = cell_A{sigmaCur};
	cur_d = norm(ACand*xCur);
	for jt = 1:length(cell_OS_A)
		if max(cell_OS_A{jt}*ACand*xCur-cell_OS_b{jt})<=0
			if cur_d<=closest_d
				closest_d = cur_d;
				sigmaOut = sigmaCur;
			end
			flg_ok =1;
			break
		end
	end
end

end


function [lineStyles, colors, lineColors, lineWidths] = plotCparamsAndColors()


lineStyles = {':','-.',':',':','-.','-.','-'};
colors = [0.2 0.2 0.2;
	0.2 0.2 0.2;
	0.5 0.5 0.5;
	0.7 0.7 0.7;
	0.5 0.5 0.5;
	0.7 0.7 0.7;
	0.8 0.7 0.3
	0.9 1 0];

lineColors = [
	[255, 225, 13]/255;%gold
	0 1 1; %electric cyan
	0.0660, 0.4430, 0.7450 %Orange: [
	0.8660, 0.3290, 0.0000;%Yellow
	0.9290, 0.6940, 0.1250;%Purple:
	0.5210, 0.0860, 0.8190;%Green:
	0.2310, 0.6660, 0.1960;%Cyan:
	0.1840, 0.7450, 0.9370;% Pink:
	0.8190, 0.0150, 0.5450;
	[250 218 221]/250];


lineWidths=[2;1.9;3.4;8;4;8];

end