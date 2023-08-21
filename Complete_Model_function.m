

function [Rates,time,DT,Stim] = Complete_Model_function(trialtype,VA,VB,varargin)

% Simulates the full model for a given trial type
%
% INPUTS:
% trialtype='Ab', 'aB', 'Ba', or 'bA'
% Va and VB : object values (values between 0.5 and 1)
% The function plots the dynamics of some neurons
% if you don't need the plots set Complete_Model_function(trialtype,VA,VB,'off')
%
% OUTPUTS:
% Rates : activity of the different neurons
%         Rates(:,1) : V1
%         Rates(:,2) : M1
%         Rates(:,3) : V2
%         Rates(:,4) : M2
%         Rates(:,5) : A1st
%         Rates(:,6) : B1st
%         Rates(:,7) : Switch
%         Rates(:,8) : C1
%         Rates(:,9) : C2
%         Rates(:,10:13) : four combination neurons
%         Rates(:,14) : CA
%         Rates(:,15) : CB
%
% time : time in seconds
% DT : resolution time
% Stim : indicates when stimuli were on
%
% From:
% Grabenhorst et al. (2023)
% "A view-based decision mechanism for rewards in the primate amygdala"
%--------------------------------------------------------------------------

if strncmp(trialtype,'Ab',2) == 1
    ttype = 1;
    Afirst = 1;
    if VA < VB
        warning('VA should be larger than VB for ''Ab'' trials')
        return
    end
elseif strncmp(trialtype,'aB',2) == 1
    ttype = 2;
    Afirst = 1;
    if VB < VA
        warning('VA should be lower than VB for ''aB'' trials')
        return
    end    
elseif strncmp(trialtype,'Ba',2) == 1
    ttype = 3;
    Afirst = 0;
    if VB < VA
        warning('VA should be lower than VB for ''Ba'' trials')
        return
    end        
elseif strncmp(trialtype,'bA',2) == 1
    ttype = 4;
    Afirst = 0;
    if VA < VB
        warning('VA should be larger than VB for ''bA'' trials')
        return
    end    
end

if nargin<4
   flag_plot_all = 1;
elseif nargin == 4
   if strncmp(varargin,'off',3) == 1 
   flag_plot_all = 0;
   else
   flag_plot_all = 1;    
   end
end


% Time/signals parameters:
tau = .01; % in [s]
tau1 = 2*tau;
dt=0.01*tau;
tmax = 4; % duration in [s]
tspan=0:dt:tmax;
L = length(tspan);

t1 = find(tspan==1); % 1st stim on 1s
t2 = find(tspan==1.5); % 1st stim off 1.5s
t3 = find(tspan==2); % 2nd stim on  2s
t4 = find(tspan==2.5); % 2nd stim off 2.5s

ds = 23; % downsampling
Tds = length(0:ds*dt:tmax)-1;
time = 0:ds*dt:tmax-ds*dt;
time = time-1;
DT = ds*dt;


% noise:
sigma = 0.03; 

% sequence of inputs:
It = zeros(6,L);
if Afirst == 1
It(3,t1:t2) = VA;
It(3,t3:t4) = VB;
It(1,t1:t2) = 1-VA;
It(1,t3:t4) = 1-VB;
It(5,t1:t2) = 1; % input to A neuron
It(6,t3:t4) = 1; % input to B meuron
else
It(3,t1:t2) = VB;
It(3,t3:t4) = VA;
It(1,t1:t2) = 1-VB;
It(1,t3:t4) = 1-VA;    
It(6,t1:t2) = 1; % input to A neuron
It(5,t3:t4) = 1; % input to B meuron
end

% Parameters
% Sequence neurons:
%--------------------------------------------------------------------------
% Connectivity:
k = 5;
wCC = 0;
wMC = .6;
wCM = .6;
wMM = k;
wI  = 0;
W = [wCC -wMC 0 0; wCM wMM 0 -wI; 0 0 wCC -wMC; -wI 0 wCM wMM];

I0 = 0;
%activation function:
RC = @(x) (x>-1).*tanh(x);
RM = @(x) (x>I0*wCM).*x/k;
F_seq = @(x) [feval(RC,x(1));feval(RM,x(2));feval(RC,x(3));feval(RM,x(4))];
r0 = .8; % baseline rate of V1 and V2 neurons

% A-first and B-first neurons
%--------------------------------------------------------------------------

% syn. depression:
k1 = 0.015;
k2 = 0.22;
tauq = 0.5;

% Connectivity:
wxA = 1;
wxB = 1;

%activation function:
F_first = @(x) (x>0).*x;

% E/I switch
%--------------------------------------------------------------------------

% E/I Connectivity:
wII=7; 
wIE=10; 
wEI=9; 
wEE=16; 

Wb = [wEE -wEI; wIE -wII];
Is = [-3*ones(1,L);-1*ones(1,L)];

% Transfer functions
aE = 1;
aI = 3;
FE = @(x) 1./(1 + exp(-x*aE) );
FI = @(x) 1./(1 + exp(-x*aI) );
Fb = @(x) [feval(FE,x(1));feval(FI,x(2))];

% Order-based decision:
%--------------------------------------------------------------------------    
% Connectivity:
wrec = 2.5; 
w11=wrec; 
w22=w11;
w12=2;
w21=2;

WdO = [w11 -w12; -w21 w22];

% Transfer function
F = @(x) 1./(1 + exp(-x/1) );


% Initialize:
%---------------------------------------------------------------------------
rV = zeros(4,1);  
RV = zeros(Tds,4);
II = zeros(1,Tds);
U = zeros(Tds,4);

a = 0;
b = 0;
xa = 0;
xb = 0;
qb = 1;
qa = 1;
qS = 1;

A1st = zeros(Tds,1);
B1st = zeros(Tds,1);
rS = zeros(2,1);
RS = zeros(Tds,2);
rDO = zeros(2,1);
RDO = zeros(Tds,2);
% rr = zeros(4,1);
% RR = zeros(Tds,4);
rrI = zeros(4,1);
RRI = zeros(Tds,4);
rDOb = zeros(2,1);
RDOb = zeros(Tds,2);


% Simulate:
%--------------------------------------------------------------------------
tt = 0;
                
for t = 1:L   
                   
    % Sequence Neurons ----------------------------------------------------
    rVth = rV.*(rV>0);
    u = W*rVth + It(1:4,t);   
    rV = rV + dt*(-rV + feval(F_seq,u))./tau + sqrt(dt)*sigma*randn(4,1);  
    
    % Afirst and Bfirst neurons -------------------------------------------
    ua = It(5,t);
    a = a + dt*( -a + feval(F_first,ua) )/tau1 + sqrt(dt)*sigma*randn;  

    ub = It(6,t);
    b = b + dt*( -b + feval(F_first,ub) )/tau1 + sqrt(dt)*sigma*randn;  

    ux = wxA*a + qb*wxB*b;
    xb = xb + dt*( -xb + feval(F_first,ux) )/tau1 + sqrt(dt)*sigma*randn;  

            %syn depression:
            qb = qb + dt*( k1*(1-qb) - k2*xb*qb )/tauq;
            
    ux = qa*wxA*a + wxB*b;
    xa = xa + dt*( -xa + feval(F_first,ux) )/tau1 + sqrt(dt)*sigma*randn;  

            %syn depression:
            qa = qa + dt*( k1*(1-qa) - k2*xa*qa )/tauq;
                     
    % Switch --------------------------------------------------------------        

    ss = xa + xb;
    u = Wb*rS + 0.1*qS*[ss; 0] + Is(:,t) ; 
    K = feval(Fb,u);            
    rS = rS + dt*(-rS + K)/(tau/3) + sqrt(dt)*sigma*randn(2,1);

            %syn facilitation:
            qS = qS + dt*( k1*(1-qS) + 5*k2*rS(1)*qS )/tauq;
    
    % order-based WTA decision --------------------------------------------
    w = 0.1;
    I = w*([rV(1);rV(3)]) + rS(1)  - 1.6;
    u = WdO*rDO + I;            
    K = feval(F,u);            
    rDO = rDO + dt*(-rDO + K)/tau + sqrt(dt)*sigma*randn(2,1);
    
    % Intermediate mix neurons --------------------------------------------
    
    u1 = rDO(1) + xa; %  1st wins and A was the 1st => A
    u2 = rDO(2) + xb; %  2nd wins and B was the 1st => A
    u3 = rDO(2) + xa; %  2nd wins and A was the 1st => B
    u4 = rDO(1) + xb; %  1st wins and B was the 1st => B
    u = [u1;u2;u3;u4];
    
    %rr = rr + dt*( -rr + u )/tau1; % + sqrt(dt)*sigma*randn(4,1);
    
    % Intermediate threshold neurons
    rrI = rrI + dt*( -rrI + feval(F_first,u-1.35) )/tau1;
    
    % object-based WTA decision --------------------------------------------
    w1 = .5;
    I = w1*[rrI(1)+rrI(2);rrI(3)+rrI(4)] + rS(1)  - 1.6;
    uu = WdO*rDOb + I;            
    K = feval(F,uu);            
    rDOb = rDOb + dt*(-rDOb + K)/tau1 + sqrt(dt)*sigma*randn(2,1);

    % Store variables -----------------------------------------------------                    
    if mod(t,ds)==0
       tt=tt+1;
       RV(tt,:)=rV;      
       B1st(tt) = xb;
       A1st(tt) = xa;             
       RS(tt,:) = rS;     
       RDO(tt,:) = rDO;    
       %RR(tt,:) = rr;
       RRI(tt,:) = rrI;   
       RDOb(tt,:) = rDOb;           
       II(tt)=It(1,t)+It(3,t);    
       U(tt,:) = u;    
    end            
    
end
        
RV(:,[1,3]) = RV(:,[1,3]) + r0;
        
scale = 20; %Hz
scale2 = 40; %Hz
        
% group all time-courses:

Rates = [scale*RV scale*A1st scale*B1st scale*RS(:,1) scale2*RDO scale*RRI scale2*RDOb];
% nb of var:  4       1         1           1               2         4         2

Stim = II;

% FIGURES
%--------------------------------------------------------------------------
if flag_plot_all == 1

        figure   
        area(time,scale*(II>0)*2.5,'facecolor',[.9 .9 .9],'edgecolor','none')
        hold on    
        h1=plot(time,scale*RV(:,1),'b-','linewidth',1);
        h2=plot(time,scale*RV(:,3),'r-','linewidth',1);
        box off
        xlabel('time [s]')
        ylabel('Firing rate [Hz]')
        set(gca,'xlim',[-.5 2.5])      
        legend([h1,h2],'V1','V2')
        
        if ttype == 1
         text( .25,.8,'A','fontsize',14,'units','normalized')
         text( .56,.8,'b','fontsize',14,'units','normalized')
        elseif ttype == 2
         text( .25,.8,'a','fontsize',14,'units','normalized')
         text( .56,.8,'B','fontsize',14,'units','normalized')
        elseif ttype == 4
         text( .25,.8,'b','fontsize',14,'units','normalized')
         text( .56,.8,'A','fontsize',14,'units','normalized')
        elseif ttype == 3
         text( .25,.8,'B','fontsize',14,'units','normalized')
         text( .56,.8,'a','fontsize',14,'units','normalized')   
        end
        
        
        figure   
        area(time,scale*(II>0)*2.5,'facecolor',[.9 .9 .9],'edgecolor','none')
        hold on    
        h1=plot(time,scale2*RDO(:,1),'linewidth',2);
        h2=plot(time,scale2*RDO(:,2),'linewidth',2);
        box off
        xlabel('time [s]')
        ylabel('Firing rate [Hz]')
        set(gca,'xlim',[-.5 2.5])      
        legend([h1,h2],'C1 1rst wins','C2 2nd wins')
        
        if ttype == 1
         text( .25,.8,'A','fontsize',14,'units','normalized')
         text( .56,.8,'b','fontsize',14,'units','normalized')
        elseif ttype == 2
         text( .25,.8,'a','fontsize',14,'units','normalized')
         text( .56,.8,'B','fontsize',14,'units','normalized')
        elseif ttype == 4
         text( .25,.8,'b','fontsize',14,'units','normalized')
         text( .56,.8,'A','fontsize',14,'units','normalized')
        elseif ttype == 3
         text( .25,.8,'B','fontsize',14,'units','normalized')
         text( .56,.8,'a','fontsize',14,'units','normalized')   
        end
        
        
        figure   
        area(time,scale*(II>0)*2.5,'facecolor',[.9 .9 .9],'edgecolor','none')
        hold on    
        h1=plot(time,scale*A1st,'k-','linewidth',1);
        h2=plot(time,scale*B1st,'color',[.5 .5 .5],'linewidth',1);
        box off
        xlabel('time [s]')
        ylabel('Firing rate [Hz]')
        set(gca,'xlim',[-.5 2.5])             
        legend([h1,h2],'A1st','B1st')

        if ttype == 1
         text( .25,.8,'A','fontsize',14,'units','normalized')
         text( .56,.8,'b','fontsize',14,'units','normalized')
        elseif ttype == 2
         text( .25,.8,'a','fontsize',14,'units','normalized')
         text( .56,.8,'B','fontsize',14,'units','normalized')
        elseif ttype == 4
         text( .25,.8,'b','fontsize',14,'units','normalized')
         text( .56,.8,'A','fontsize',14,'units','normalized')
        elseif ttype == 3
         text( .25,.8,'B','fontsize',14,'units','normalized')
         text( .56,.8,'a','fontsize',14,'units','normalized')   
        end
        
        

        figure   
        area(time,scale*(II>0)*2.5,'facecolor',[.9 .9 .9],'edgecolor','none')
        hold on    
        h1=plot(time,scale2*RDOb(:,1),'linewidth',2);
        h2=plot(time,scale2*RDOb(:,2),'linewidth',2);
        box off
        xlabel('time [s]')
        ylabel('Firing rate [Hz]')
        set(gca,'xlim',[-.5 2.5])             
        legend([h1,h2],'Choose A','Choose B')

        if ttype == 1
         text( .25,.8,'A','fontsize',14,'units','normalized')
         text( .56,.8,'b','fontsize',14,'units','normalized')
        elseif ttype == 2
         text( .25,.8,'a','fontsize',14,'units','normalized')
         text( .56,.8,'B','fontsize',14,'units','normalized')
        elseif ttype == 4
         text( .25,.8,'b','fontsize',14,'units','normalized')
         text( .56,.8,'A','fontsize',14,'units','normalized')
        elseif ttype == 3
         text( .25,.8,'B','fontsize',14,'units','normalized')
         text( .56,.8,'a','fontsize',14,'units','normalized')   
        end
        
        
        
end
