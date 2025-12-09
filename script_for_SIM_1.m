
clc
clear all
close all 

%% Add path 
%%
% To run this code, the SEED-G toolbox is required 

% Anzolin A, Toppi J, Petti M, Cincotti F and Astolfi L 2021 
% Seed-g: Simulated eeg data generator for testing connectivity 
% algorithms Sensors 21

% change the path 
addpath(genpath('C:\Users\marti\OneDrive - Politecnico di Milano\Progetto Dottorato\SEED-G-toolbox-main'));
 
%% Load of real sources from the toolbox for AR components 
% change the path
path='C:\Users\marti\OneDrive - Politecnico di Milano\Progetto Dottorato\SEED-G-toolbox-main\real data\'; 
name='sLOR_cortical_sources.mat'; 
load([path name]);
data=EEG.samp; 
fs=EEG.fSamp; 

%% Parameters 
% Parameters can be change according to the desired dataset properties 

Nodes=3; % number of nodes 
p=5; % Order of the model 
DataLength=250; % Length of the simulated data
Trials =100; % Number of trials 
SNR =10; % Signal to noise ratio 


%% Model definition 
% connection positions and values 
Model=zeros(Nodes,Nodes,2); 

Model(2,3,1)=0; 
Model(2,1,1)=0; 
Model(3,1,1)=0.9; 

Model(2,3,2)=0; 
Model(2,1,2)=0.9;
Model(1,3,2)=0; 

% Delay matrix --> number of samples lag 
DelayMatrix=zeros(Nodes,Nodes); 

DelayMatrix(2,3)=0;
DelayMatrix(3,1)=3;
DelayMatrix(2,1)=3;

%%
trans_null=0; % this parameter defines a step transition (see the SEED-G toolbox for further details)

%%
%%% Time-Varying Implementation Parameters (see SEED-G TOOLBOX for more
%%% details)

% These lines are taken from the SEED-G toolbox and adapted to our simulated data 

CampIntPerc=[0.5 0.5]; 
CampInt=round(CampIntPerc*DataLength);
NumIntCost=length(CampInt);
Del_Range = 1:p;

%%
%%%Input for MVGC-Toolbox function 'var_to_tsdata'

mtrunc=     0;       %default
decayfac=   100;     %default
Singtr=     1;

%%% Threshold for the signal amplitude
SigLim = 80;

%sig_num =       size(Model,1);
Sw =            eye(Nodes);                  % residual covariance matrix - White Noise (input Barnett Toolbox)
AR_num =      1; %round(sig_num*AR_perc);        % number of Auto-Regressive components
ARpos=1; 


max_att =       Trials*100;                    
% maximum number of attempts to generate the chosen num. of trials with the same model

%%% Inputs for function "arfit_v3"
selector=      'AIK';
no_const=      'zero';

%%% AR parameters evaluation from real data
for r=1:AR_num
    for tt=1:size(data,3)
        EEG_ch = data(:,r,tt);
        orgEEG{tt}=EEG_ch;
    end %cycle on real trial
    clear tt
    ar_coef(:,r) = arfit_v3(orgEEG, p, p, selector, no_const);
end

flag=1;
cont=0;

trans_samp=zeros(1, NumIntCost-1);  
while flag==1
    cont=cont+1;
    flag=0;
    
    %%% 1.TV model generation
    %   1.1 Model generation in the constant intervals
    for i=1:NumIntCost
        startInt(i)=sum(CampInt(1:i-1))+sum(trans_samp(1:i-1))+1;
        endInt(i)=startInt(i)-1+CampInt(i);
        aux=reshape(repmat(Model(:,:,i),1,CampInt(i)),[Nodes Nodes CampInt(i)]);
        TVmod(:,:,startInt(i):endInt(i))=aux;
    end
    clear i
    
    
    %% Time-varying model 
 
    for mo=1:size(TVmod,3)
        ModelDel(:,:,:,mo)=rearrangeModel(Del_Range,TVmod(:,:,mo),DelayMatrix);
        
        %%% Add Auto-Regressive coefficients on the main diagonal
        for ii=1:AR_num
            ModelDel(ARpos(ii),ARpos(ii),:,mo)=ar_coef(1:p,ii)';
        end
        
        clear ii
    end %cycle on the number of models
    
    clear mo
    
    Tr=1;
    num_iter=0;
    
    while Tr~=Trials+1
        num_iter=num_iter+1;
        
        if num_iter>max_att %Max number of attempts
            Tr=1;
            flag=1;
            num_iter=0;
            break
        end
        
        [Y,E]=var_to_tsdata_complete(ModelDel,Sw,DataLength,Singtr,mtrunc,decayfac);
        Ctr=find(abs(Y)>SigLim);
        
        if ~isempty(Ctr)
            continue
        else
            for ch=1:Nodes
                nc=randn(1,size(Y,2));
                Ynorm=norm(Y(ch,:));
                ncnorm=norm(nc);
                noise=(Ynorm/(ncnorm*sqrt(SNR)))*nc;  
                Ynoise(ch,:)=Y(ch,:)+noise;
            end
            Y_gen(:,:,Tr)=Ynoise';
            E_gen(:,:,Tr)=E';
            Tr=Tr+1;
        end
    end
    
    clear Tr Ctr num_iter
    flag_out=1;
    
    if cont>10
        flag_out=0;
        Y_gen=[];
        E_gen=[];
        break
    end
end

%% Plot of the simulated signals and noise

Nreal=size(data,1); 
t=(0:DataLength-1)/fs; 
t1=(0:Nreal-1)/fs; 
figure 
plot(t,Y_gen(:,:,1)); %sto plottando un trial 
title('Simulated signal');
figure
plot(t1,data(:,:,1)); 
title('Real signal'); 
figure
plot(t,E_gen(:,:,1)); 
title('Noise'); 


%% Spectral properties 
[p, fx]=pwelch(Y_gen(:,:,1),[],[],[],fs);
figure
plot(fx,p)
title('Spectrum')

%%
[p, fx] = pwelch(Y_gen(:,:,1), [], [], [], fs);

figure;
plot(fx, 10*log(p), 'LineWidth', 2); 
grid on;
legend('y1','y2','y3'); 
xlim([0 125]); 


%% save the results 

save('DelayMatrix','DelayMatrix');
save('E_gen','E_gen'); 
save('Model','Model');
save('ModelDel','ModelDel');
save('TVmod','TVmod');
save('Y_gen','Y_gen');

