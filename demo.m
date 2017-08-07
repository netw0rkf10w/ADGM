%%% Demo for the Alternating Direction Graph Matching algorithm
%%% I (D. Khue Le-Huu) adapted this code from Duchenne's
%%% If you find this code useful, please cite:
%%% D. Khuê Lê-Huu and Nikos Paragios. Alternating Direction Graph Matching. arXiv preprint arXiv:1611.07583 (2016).

addpath('./ann_mwrapper');
addpath('./ADGM');
addpath('./TM/mex');
addpath('./TM');

close all;

%% Create a third-order graph matching problem

%number of points
nP1=30; 
nP2=40;

nPmax = max(nP1, nP2);

%randomly generate them
P1=randn(2,nPmax);

%generate modified version of points 1
scale= 0.5+rand();
theta = 0.5*rand();
Mrot = [cos(theta) -sin(theta) ; sin(theta) cos(theta) ];
P2=Mrot*P1*scale+0.05*randn(size(P1));

P1 = P1(:,1:nP1);
P2 = P2(:,1:nP2);

%number of used triangles (results can be bad if too low)
nT=nP1*50;
t1=floor(rand(3,nT)*nP1);
while 1
  probFound=false;
  for i=1:3
    ind=(t1(i,:)==t1(1+mod(i,3),:));
    if(nnz(ind)~=0)
      t1(i,ind)=floor(rand(1,nnz(ind))*nP1);
      probFound=true;
    end
  end
  if(~probFound)
    break;
  end
end

%generate features
t1=int32(t1);
[feat1,feat2] = mexComputeFeature(P1,P2,int32(t1),'simple');

%number of nearest neighbors used for each triangle (results can be bad if
%too low)
nNN=300;

%find the nearest neighbors
[inds, dists] = annquery(feat2, feat1, nNN, 'eps', 10);

%build the tensor
[i j k]=ind2sub([nP2,nP2,nP2],inds);
tmp=repmat(1:nT,nNN,1);
indH = double(t1(:,tmp(:)))'*nP2 + [k(:)-1 j(:)-1 i(:)-1];
valH = exp(-dists(:)/mean(dists(:)));


%% Initiatial solution
X0 = 1/nP2*ones(nP2,nP1);



%%
%% Solve with Tensor Matching
%%
[X_TM, ~]=tensorMatching(X0,[],[],[],[],indH,valH);

% Discretize and compute the objective score
X_TM = asgHun(X_TM);
objective_TM = getMatchingScore(X_TM, indH, valH);
fprintf('Objective score TM = %f\n', objective_TM);

%draw
if 1
  figure;
  imagesc(X_TM);
  title('TM');
  figure;
  hold on;
  plot(P1(1,:),P1(2,:),'r x');
  plot(P2(1,:),P2(2,:),'b o');
  [~, match] = max(X_TM);
  for p=1:nP1
    plot([P1(1,p),P2(1,match(p))],[P1(2,p),P2(2,match(p))],'k- ');
  end
  title('TM');
end


%%
%% Solve with ADGM
%%

% The following are typical values, but one may tune them to get better performance:
% + Increasing rho or eta usually results in faster convergence but lower 
% objective (and vice-versa: decreasing them usually offer higher objective values) 
% + Decreasing iter1 or iter2 usually results in faster convergence but 
% lower objective (and vice-versa: increasing them usually offer higher objective values)
rho = nP1*nP2/1000;
MAX_ITER = 5000;
eta = 2.0; 
iter1 = 200;
iter2 = 50;

%% ADGM1

[X_ADGM1] = ADGM1(X0, [], [], [], [], indH, valH, rho, MAX_ITER, false, eta, iter1, iter2);
%[X_ADGM1] = ADGM1_SYMMETRIC(X0, [], [], [], [], indH, valH, rho, MAX_ITER, true, eta, iter1, iter2);
X_ADGM1 = asgHun(X_ADGM1);
objective_ADGM1 = getMatchingScore(X_ADGM1, indH, valH);
fprintf('Objective score ADGM1 = %f\n', objective_ADGM1);           

%draw
if 1
  figure;
  imagesc(X_ADGM1);
  title('ADGM1');
  figure;
  hold on;
  plot(P1(1,:),P1(2,:),'r x');
  plot(P2(1,:),P2(2,:),'b o');
  [~, match] = max(X_ADGM1);
  for p=1:nP1
    plot([P1(1,p),P2(1,match(p))],[P1(2,p),P2(2,match(p))],'k- ');
  end
  title('ADGM1');
end

%% ADGM2
[X_ADGM2] = ADGM2(X0, [], [], [], [], indH, valH, rho, MAX_ITER, false, eta, iter1, iter2);
%[X_ADGM2] = ADGM2_SYMMETRIC(X0, [], [], [], [], indH, valH, rho, MAX_ITER, true, eta, iter1, iter2);
X_ADGM2 = asgHun(X_ADGM2);
objective_ADGM2 = getMatchingScore(X_ADGM2, indH, valH);
fprintf('Objective score ADGM2 = %f\n', objective_ADGM2);           

%draw
if 1
  figure;
  imagesc(X_ADGM2);
  title('ADGM2');
  figure;
  hold on;
  plot(P1(1,:),P1(2,:),'r x');
  plot(P2(1,:),P2(2,:),'b o');
  [~, match] = max(X_ADGM2);
  for p=1:nP1
    plot([P1(1,p),P2(1,match(p))],[P1(2,p),P2(2,match(p))],'k- ');
  end
  title('ADGM2');
end


