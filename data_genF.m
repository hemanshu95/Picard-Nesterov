clear, clc;
% This is an example for running the function overlapping
%  Problem:
%  min  1/2 || x - v||^2 + z_1 \|x\|_1 + z_2 * sum_i w_i ||x_{G_i}||
%   Minimizes 1/2 ||Xw-y||^2 + gamma ||Bw||_GL

%%
load('vant.mat')
m=size(X,1);  n=size(X,2);       % The data matrix is of size m x n
s1=floor(m*0.65)
s2=m-s1;
rander_train=randperm(m,s1);
rander_train=sort(rander_train);

k=1;
rander_test=[];
for i=1:m
    if(sum(i==rander_train)==0)
        rander_test(k)=i;
        k=k+1;
    end
end
X_train=X(rander_train,:);
X_test=X(rander_test,:);

Y_train=Y(rander_train,2);
Y_test=Y(rander_test,2);


eps_acc = 10e-8;

gamma = 0.5;
w0 = zeros(n,1);
iters_fp = 1000;
eps_fp = 10e-4;
iters_acc = 1000;
eps_acc = 10e-6;
svdX = 1;
svdB = 1;
kappa=0.5;



fileID=fopen('pathways.txt','r');
f1='%d %d';
file1=fscanf(fileID,f1,[2,15126])';
groups=[];
gn=max(file1(:,1));
disp('ma');
disp(gn);
for i=1:gn
    groups{i}=[];
end
for i=1:15126
    
    groups{file1(i,1)}= [groups{file1(i,1)} find(entrez==file1(i,2))];
end

disp('Created Groups' );
% Create group indexing matrix
B = [];
d = size(X_train,2);

for i=1:length(groups)
    temp = zeros(length(groups{i}),d);
    temp(:, groups{i} ) = eye(length(groups{i}));
	B = [B;temp];
end
disp('done Initialisation' );
[w,costs,fp_diff,vs, B, t1,theta,result] = sq_gl_B_fp_acc2(X_test,Y_test,X_train, Y_train,B, gamma, groups, svdX, svdB,kappa, w0, iters_fp, eps_fp, iters_acc, eps_acc);
