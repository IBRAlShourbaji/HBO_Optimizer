% ---------------------------------------------------------------------------------------------------------------------- % 
%  Developed in MATLAB R2024a
%  programmer: Mohammed Jameel  
%  E-Mail: moh.jameel@su.edu.ye, mohjameel555@gmail.com
%  Paper:
%  Mohammed Jameel, Ahmed R. El-Saeed, Anis Elgabli, Ibrahim Al-Shourbaji, Pramod H. Kachare, Mohamed Abouhawwash 
%  Heavy Ball Optimizer: A Momentum-Driven Metaheuristic with Superior 
% Scalability and Real-World Applications          
% ---------------------------------------------------------------------------------------------------------------------- %

% This function initialize the first population 
function X=initialization(nP,dim,ub,lb)
Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb

if Boundary_no==1
    X=rand(nP,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        X(:,i)=rand(nP,1).*(ub(i)-lb(i))+lb(i);
    end
end