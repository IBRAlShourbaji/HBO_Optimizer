% ---------------------------------------------------------------------------------------------------------------------- % 
%  Author Programmer: Dr. Mohammed Jameel  
%  E-Mail: moh.jameel@su.edu.ye, mohjameel555@gmail.com
% ---------------------------------------------------------------------------------------------------------------------- %
% ---------------------------------------------------------------------------------------------------------------------- %
% Please Refer the Following:
% Heavy Ball Optimizer: A Momentum-Driven Metaheuristic with Superior 
% Scalability and Real-World Applications
% ---------------------------------------------------------------------------------------------------------------------- %

function [f_best, Xbest, Convergence_curve] = HBO(N, max_iter, lb, ub, dim, fobj) 
  % Input arguments:
    % N     - Number of particles in the population
    % max_iter - Maximum number of iterations
    % lb    - Lower bound of the search space
    % ub    - Upper bound of the search space
    % dim   - Dimensionality of the search space
    % fobj  - Objective function to minimize/maximize
    %% Initialization phase
    X = initialization(N, dim, lb, ub);  % Current population
    Fitness = zeros(N, 1);  % Fitness of current population
    Convergence_curve = zeros(1, max_iter);  % To track the best cost
    GSR = zeros(N,dim);
    GSRP = zeros(N,dim);
    Xnew = zeros(N,dim);

    % Evaluate initial fitness for X
    for i = 1:N
               Fitness(i) = fobj(X(i,:));
    end  

    % Determine initial best solution for X
    [~, Ind] = sort(Fitness);
    f_best = Fitness(Ind(1));  % Best fitness
    Xbest = X(Ind(1), :);      % Best solution

    % Independent initialization for Xprev (t-1)
    Xprev = initialization(N, dim, lb, ub);  % Initial previous position
    FitnessP = zeros(N, 1);  % Fitness for Xprev

      for i = 1:N
            FitnessP(i) = fobj(Xprev(i,:));
      end  

    % Determine best  solution for Xprev
    [~, Ind] = sort(FitnessP);
    XprevBest = Xprev(Ind(1), :);
    %% Main Optimization Loop
  iter=0;
  while iter < max_iter    
       X_M = mean(X); % mean postion of the current solutions
       Xprev_M = mean(Xprev); % mean postion of the previous solutions
       F_XM = fobj(X_M);
       F_XprevM = fobj(Xprev_M);

       for i = 1:N

%-------------------------Heavy Ball Search Rule (HBSR)-------------------%
         if iter==0
            GSRP(i,:) = initialization(1, dim, lb, ub);
         else
            GSRP(i,:) = GSR(i,:);
         end
          %%% Updating the Previous Solution (Xprev)
          k = Xprev(i, :);            
          y1 = randi(dim);
          k(y1) = rand * k(y1);
          Xprev(i, :) = k;      
          XR = initialization(1, dim, lb, ub);
          GSR(i,:) = gradientOperator(X(i,:),Xbest,Xprev(i,:),XR);      
          gamma = 1 + rand;
          h = (GSR(i,:) -rand.*  GSRP(i,:));
          DeltXX = (X(i, :) -rand.* Xprev(i, :)); % Delta X for search rule
          L = gamma*(norm(h) / (norm(DeltXX) + eps));
          mu = norm(h - L .* DeltXX) / (norm(DeltXX)+ eps);
          alpha = (rand *(sqrt(L)-sqrt(mu)).^2) ;
          beta = (rand *(sqrt(L)-sqrt(mu)) / (sqrt(L)+sqrt(mu))+eps)^2;

         % Adaptive Momentum Mechanism  
         HBM = beta.* (rand.* X(i, :) - Xprev(i, :)); 

         X_HBSR = X(i, :)- alpha .* (2.* GSR(i,:)-rand.*GSRP(i,:))+ HBM;  % Eq.(5)         

%-------------------Random Heavy Ball Strategy (RHBS) -----------------------%       
         omega1 = -1 + 2 * rand(); omega2 = -0.5 + rand(); omega3= 0.5 * rand;          
         %% Random indices
        a = randi(N);
        b = randi(N);

        while a==i || b==i || a==b
            a = randi(N);
            b = randi(N);
        end
        RB = Brownian(dim);
         if F_XprevM < F_XM % Eq(21)
             X_p = XprevBest+omega1.*(X_HBSR- rand.* Xprev_M)+rand.*(X(a,:)-X(b,:)) ; 
             X_EBSR = X_p+omega1.*(XprevBest-rand.*Xprev(i,:))+omega2.* RB.*(Xprev_M-X(a,:));        
         else
             if Fitness(i) < F_XM
                 X_c = Xbest+omega1.*(X_HBSR- rand.* X_M)+rand.*(X(a,:)-X(b,:)) ; 
                 X_EBSR = X_c+omega1.*(Xbest-X(i,:))+omega2.*RB.*(X_M-X(a,:));
             else
                 X_EBSR = Xbest-rand*X_M+omega3.*(X(i, :)- X(a, :))+omega2.*(X(i, :)- X(b, :));
             end
         end

%-------------------Accelerated Convergence Mechanism (ACM)----------------------%
          if rand < rand
             phi1 = rand;
            elseif rand < rand
               phi1 = RB(1);
          else
                phi1 = 0.05*levy(1,1,1.5);   %Levy random number vector
          end

          if rand < rand
              phi2 = rand;
          else
              phi2 = 1;
          end   
     rho1 = 0.5 - 0.5 * rand(); 
     if  f_best  < fobj(X_EBSR)    
           S3 = phi2.*X_HBSR-rand.*X_M; S2=phi1.* X_EBSR; S1= Xbest;
          X_ACM = S1- rand .* (S2- rho1.*S1).^2 ./( rand.* S3-2 .* rand(1,dim).*S2+rand(1,dim).*S1+eps);
     else
          S3 = (phi2.*X_HBSR-rand.* X_M); S2= phi1.* Xbest;  S1= X_EBSR;
           X_ACM = S1- rand .* (S2- rho1.*S1).^2 ./(rand.*S3-2.*rand(1,dim).*S2+rand(1,dim).*S1+eps);
     end 

       if rand < 0.5    % Eq(27)
           Xnew(i,:) = X_ACM ;
       else
           if i==N
              Xnew(i,:) = (lb+rand*(ub-lb));
           else
               if mod(iter,2) == 0  
                    RR = 0.05*levy(1,dim,1.5);
               else
                    RR = zeros(1,dim);
               end
             omega22 = -1 + 2 * rand();
               a1 = randperm(N,1);
                 b1 = randperm(N,1);
                   c1 = randperm(N,1);  
                   XX = (X(a1,:)+X(b1,:)+X(c1,:) )./3;
                 Xnew(i,:) = X(i,:)+rand.*(Xbest - XX)+omega22.*(Xbest- X_M)+RR;   
            end
      end
              % Enforce boundary constraints
            Xnew(i,:) = min(max(Xnew(i,:) , lb), ub);

            Xprev(i,:) =  X(i,:);       % Current X becomes previous Xprev

            XprevBest = Xbest;  % Current best becomes previous best

            Xnew_Cost = fobj(Xnew(i,:));

             % Update best positions
            if Xnew_Cost < Fitness(i)
                X(i, :) = Xnew(i,:);
                Fitness(i) = Xnew_Cost;

                % Update global best
                if Fitness(i) < f_best
                    f_best = Fitness(i);
                    Xbest = X(i, :);
                end
            end
       end
        %%% Record convergence
          iter = iter+1; 
          Convergence_curve(iter) = f_best ;
          disp(['Iteration ', num2str(iter), ': Best Cost = ', num2str(f_best)]);   
  end
end


function GSR = gradientOperator(X,Xbest,Xprev,XR)
    dim = size(X,2);
    S0 = rand(1,dim).*(X-Xbest)./(((Xprev-rand.*X).*(Xprev-rand.*Xbest))+eps);
    S1 = rand(1,dim).*(2.*X-Xprev-Xbest)./(((X-rand.*Xprev).*(X-rand.*Xbest))+eps);
    S2 = rand(1,dim).*(X-Xprev)./(((Xbest-rand.*Xprev).*(Xbest-rand.*X))+eps);
    Y0 = rand(1,dim).* (Xbest-rand.*Xprev);
    Y1 = rand(1,dim).* (Xbest-rand.*X);
    Y2 = rand(1,dim) .* (Xbest-rand.*XR);
   GSR = randn.*( Y0.*S0+ Y1.*S1+ Y2.* S2);
end


function o = Brownian(dim)
   T = 1;
   r = T/dim;
   dw = sqrt(r)*randn(1,dim);
   o = cumsum(dw);
end






 