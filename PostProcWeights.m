% Calculate post-processor weights
clear all;
close all
clear

%load IMEX-EIS+ method
load('3s4pType2ImExEIS+.mat')
p=P;

% calculate the number of intervals needed to stretch polynomial basis.
if rank(tau(:,[p-1,end-1]),1.e-15)==1
    rank_tau=1;
else
    rank_tau=2;
end
 m = 2;                                
 while m*s-(p+rank_tau)<0 
     m=m+1;
 end
 %
 % form the postprocessor
 gp = c-(m-1:-1:0);                    % Identify Gridpoints
 S = (vander(gp(:)));                  % Build Polynomial Basis
if rank_tau==1
      TAU = repmat(tau(:,p-1),m,1);  % Build Extended Truncation Vector;
      S(:,1) = TAU;                       % Introduce Truncation Vector into Basis
      DD = diag([0,ones(1,length(S)-1)]); % Zero Out Truncation Error
      
 else
      TAU = repmat(tau(:,[p-1,end-1]),m,1);  % Build Extended Truncation Vector;
      S(:,1:2) = TAU;                       % Introduce Truncation Vector into Basis
      DD = diag([0,0,ones(1,length(S)-2)]); % Zero Out Truncation Error
end
 Phi = S*DD*inv(S); 
 weights=Phi(end-s+1,:)

 
 
 