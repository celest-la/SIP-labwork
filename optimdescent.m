function [xh,result,xval] = optimdescent(critfun,params,options,x0)
% critfun : name of fonction .m evaluating
% the objective function, its gradient and hessien
% params : parameter for objective function evaluation
% options : necessary options for the algorithms
% options.method : direction type, ’gradient’, ’Newton’, ’Quasi-Newton’, ...
% options.step : ’fixed’, ’variable’
% options.const : Armijo constant
% options.beta : backtracking rate
% options.TolX, options.TolF, options.TolG, options.MaxIter : stopping thresholds
% x0 : initial point
% xh : final point
% result : structure contenant les r´esultats
% result.iter : nombre d’iterations r´ealis´ees
% result.crit : valeur du crit`ere par iteration
% result.grad : norm du gradient par iteration
% result.temp : temps de calcul
% result.stop : condition d’arret (’TolX’, ’TolG’, ’TolF’, ’Maxiter’)
% xval : valeurs des it´er´ees

%options.step

tic;
xa = x0;
xval=x0;
xb=0;
d=0;
Lambda=0;
c=options.const;
alpha=options.alpha;
alpha0=options.alpha;

       [f,g,H] = feval(critfun,xa,params);
       result.crit(1)=f;
       result.grad(1)=norm(g);

        k=1;
       condition=true;

invH=eye(params.n);
    
        while(condition)
            xb=xa;
            
            fb=f;
 
        switch options.method
             case 'gradient'
            
            [f,g] = feval(critfun,xa,params);    
            d=-g/sqrt(g'*g);
            if g'*d>0
                d=-d;
            end
         
         case 'Fletcher'
              
            gb=g;
            db=d;
            [f,g] = feval(critfun,xa,params);          
                  
            bf=(norm(g)^2)/(norm(gb)^2);
            d=-g+bf*db;
            if g'*d>0
               d=-d;

           end
    case 'Newton'
        
        [f,g,H] = feval(critfun,xa,params);
        

        d=-inv(H)*g;
       if g'*d>0
           d=-d;
          
       end
     case 'Quasi-Newton'
       [f,g] = feval(critfun,xa,params);
       d=-invH*g;
       if g'*d>0
           d=-d;
        end
     case 'Gauss-Newton'
     [f,g,H,J]=feval(critfun,xa,params);
     d=inv(J'*J+Lambda*(eye(params.n)))*g;
      if g'*d>0
           d=-d;
      end
    case 'Levenberg'
      [f,g,H,J]=feval(critfun,xa,params);
      Lambda=0.5*norm(g);
      d=inv(J'*J+Lambda*diag(J'*J))*g;
      if g'*d>0
           d=-d;
      end
      end
        switch options.step
            case 'fixed'
                alpha=options.alpha;
            case 'variable'
               
              stop = 0;
              while stop==0
       
               [fj,gj] = feval(critfun,xa+alpha*d,params);
                   if (fj - f > f + alpha*d'*g)
                   %if (fj - f > c*alpha*g'*d)
                       alpha = alpha*options.beta;
                       else
                           stop = 1;
                   end
              end
            case 'logbarrier'
                 
                Ik=find(d<=0);
                if isempty(Ik)>0
                else
                alpha0=min(-xa(Ik)./d(Ik));
                end
                stop = 0;
                alpha=alpha0;
              while stop==0
       
               [fj,gj] = feval(critfun,xa+alpha*d,params);
                   if (fj - f > f + alpha*d'*g)
                   %if (fj - f > c*alpha*g'*d)
                       alpha = alpha*options.beta;
                       else
                           stop = 1;
                   end
              end
        end
        
        xa=xa+alpha*d;
       switch options.method
         case 'Quasi-Newton'
         yQN=gj-g;
          invH=(eye(params.n)-(yQN*(alpha*d)')/((alpha*d)'*yQN))'*invH*(eye(params.n)-(yQN*(alpha*d)')/((alpha*d)'*yQN))+(alpha*d)*(alpha*d)'/((alpha*d)'*yQN);
        end
        
            xval=[xval xa];
            result.crit(k)=f;
            result.grad (k)= norm(g);
k=k+1;

  %%%%%%%%%%%%%test the tolerances%%%%%%%%%%
        if k>2
            if (norm(g)<params.TolG) 
                result.stop='TolG';
                condition=false;
                
            elseif (k>=params.kmax)
                condition=false;
                result.stop='kmax';
                
            elseif ((abs(f-fb)/fb)<params.TolF)
                condition=false;
                result.stop='TolF'; 

            elseif (((norm(xa-xb))/norm(xb))<params.TolX)         
               condition=false;
               result.stop='TolX';
            end
        end
        end
  
xh=xa;
result.iter=k;
result.time=toc;

end
