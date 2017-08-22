function [w, costs, fp_diff, vs, fun, B,theta,lambdas,results,precalls,recalls,fscores] = sq_gl_B_fp_acc2(X_test,Y_test,X, y,B, gamma, groups, svdX, svdB, kappa, ...
    w0, iters_fp, eps_fp, iters_acc, eps_acc)

% Minimizes 1/2 ||Xw-y||^2 + gamma ||Bw||_GL
% where ||.||_GL is the group Lasso norm
% using fixed point prox & accelerated
% groups : cell array of group indicators
% svdX : 1 to compute L using SVD of X , 0 to use Frobenius estimate
% svdB : 1 to compute lambda using SVD of B , 0 to use Frobenius estimate

d = size(X,2);
disp('size');
disp(d);
G = length(groups);
ind = [];
cnt = 0;
for i=1:G
  ind{i} = cnt+(1:length(groups{i}));
	cnt = cnt+length(groups{i});
end

f = @(z) (norm((X*z-y),2)^2/2);

    function res = phi(z)
        % Group Lasso norm with varying size groups
        res = 0;
        for i=1:G		
            res = res + norm(z(ind{i}));
        end
        res = gamma* res;
    end

gradf = @(z) (X'*(X*z-y));

    function res = prox(z,lam)
        % Group Lasso prox map
        res = [];
        for i=1:G
            nrm = norm(z(ind{i}));
            coeff = max( [ nrm-gamma/lam , 0 ] ,[],2 ) / (nrm+eps);
            res = [res; coeff*z(ind{i})];
		
        end
    end

fun = @phi;
disp('start');

L=norm(X)^2;;
disp(size(w0));
disp('Lambda and L');
disp(L);

lambdas=[20000 50000 100000 200000 500000 100000 ];
results=(1:6)*0;
precalls=(1:6)*0;
recalls=(1:6)*0;
fscores=(1:6)*0;
for j=1:6
    disp('j=');
    disp(j);
    
    [w, costs, fp_diff, vs,theta] = fp_acc(f, @phi, B, gradf, L, lambdas(j), kappa, @prox, w0, ...
				 iters_fp, eps_fp, iters_acc, eps_acc);
    Y3=X_test*w;
    Y5=Y3>=-0.2;
    Y5=(Y5*2)-1;
    results(j)=(sum(Y5==Y_test)*100)/size(Y_test,1);
    Z=[0 0;0 0]
    for i=1:95
        if(Y_test(i)==1)
           if(Y5(i)==1)
               Z(1,1)=Z(1,1)+1;
           else
               Z(1,2)=Z(1,2)+1;
           end
       else
           if(Y5(i)==1)
               Z(2,1)=Z(2,1)+1;
           else
               Z(2,2)=Z(2,2)+1;
           end
       end
   end
   disp(Z);
   Z=(Z*100)/size(Y_test,1);
   disp(Z);
   precalls(j) = Z(1,1) / ( Z(1,1)+ Z(2,1) );
   precalls(j)=precalls(j)*100;
   recalls(j)=Z(1,1)/(Z(1,1)+Z(1,2));
   recalls(j)=recalls(j)*100;
   fscores(j)=(2*precalls(j)*recalls(j))/(precalls(j)+recalls(j));
   disp('for lambda =');
   disp(lambdas(j));
   disp(results(j));

end

    
   
end


            