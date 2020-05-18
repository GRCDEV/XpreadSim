function Sim_EpidemicXn_ROMA()
% This is a sample script for using the Xpread epidemic simulation
% It uses the Roma trace and the simulator 
% It takes some time to complete all the simulations!!


% Obtains the delivery time and delivery nodes
function [DT, average_DT, sum_delivered ] = getDeliveryTimeRatio(DN) 
    DT = ones(1,Nmax)*-1;
    sum_delivered = 0;
    sum_DT = 0;
    for ii=NF+1:size(DN,1)
        if min(DN(ii,:)) > 0
            sum_delivered = sum_delivered + 1;
            DT(ii) = max(DN(ii,:));
            sum_DT = sum_DT + DT(ii);
        end
    end
    average_DT = sum_DT/sum_delivered;
end

% Set up initial values for simulation
function [X0, BO] = SetupSimulation(N,Q)      
    X0 = zeros(N,Q);
    % Set sender nodes with message a 1.
    for ii = 1:NF
        X0(ii,:) = ones(1,Q);
    end
    BO = ones(N,Q);
    for ii = 1:N
        if mod(ii,2) == 1
            BO(ii,:) = 1:Q;
        else
            BO(ii,:) = Q:-1:1;
        end
    end
end

% Simulate with the values indicated, 
% Returns:
%    NR: Real Number of nodes involved,
%    DM: Vector with the messages delivered
%    RDM: Ratio of delivered messages
%    ADT: Average delay time
%    SMT: Messages succesusfully transmitted
%    NDM: Non delivel messages
%    TMB: Total MB transmitted
function [NR, DM, RDM, ADT, SMT, NDM, TMB] = Simulate()
    
    DM = zeros(1,length(Tam_ITEMs));
    RDM = zeros(1,length(Tam_ITEMs));
    ADT = zeros(1,length(Tam_ITEMs));
    SMT = zeros(1,length(Tam_ITEMs));
    NDM = zeros(1,length(Tam_ITEMs));
    TMB = zeros(1,length(Tam_ITEMs));

    fprintf('  Q   TamItem  (MB)     TxTime(s)   N   Deliv.  Ratio Deliv_time   TOT_MSG   *TAM (GB)\n');
    fprintf('--------------------------------------------------------------------------------------\n');
    for i = 1:length(Tam_ITEMs)

        Q = ceil(Tam_ITEMs(i)/MSJ_SIZE); 
        TX_TIME = Tam_ITEMs(i)*8/TX_SPEED/Q;
        T_items = TX_TIME*ones(1,Q);
        fprintf('%3d %7.3f -> %7.3f   %7.3f', Q, Tam_ITEMs(i)/MB, Tam_ITEMs(i)/MB/Q, TX_TIME );
        [X0, BO] = SetupSimulation(Nmax,Q);
        [DN, NDM(i), DMT, WTT ] = Epidemic_Simulation_ONETrace(trace,c_radio, Q, NF, X0, T_setup, T_items, BO);

        [DT, ADT(i), DM(i) ] = getDeliveryTimeRatio(DN); 

        SMT(i) = sum(DMT); %  ONLY messages succesfully transmitted.
        TMB(i) = WTT*TX_SPEED/8;
        % TMB(i) = (SMT(i) + 0.5*NDM(i))*Tam_ITEMs(i)/Q/MB;  % We assume that for non delivered message half message is transmitted.
        NR = size(DN,1);
        RDM(i)=DM(i)/(NR-NF);


        fprintf('  %3d     %4d  %5.3f   %8.2f     %5d    %8.2f\n', NR,  DM(i), RDM(i), ADT(i), SMT(i)+NDM(i), TMB(i)/(1024^3));

    end

end

% Configuration of the simulation
StartHour = 7;
EndHour = 12;
MB = 1024*1024;
c_radio = 50;
TX_SPEED = 10e6;
T_setup = 5;
Tam_ITEMs = [1 2.5:2.5:100]*MB;  

% Here you load a trace with a format compatible with The One Simulator,
% that is: time, node, pos_x, pos_y
% As an example, we use the Roma Trace
load('Roma_trace_2014_02_04.mat','Roma_trace_2014_02_04');

% This function extract from the trace a given time intervalñ
trace = ONETrace_extract_from_interval(Roma_trace_2014_02_04,3600*StartHour,3600*EndHour);
trace(:,1) = trace(:,1)-trace(1,1);

% Fixed Points defined for this trace
% We add these points as nodes 1 and 2
trace(:,2) = trace(:,2)+2;
trace = [0 1 60820 51230; 0 2 59093 501742; trace]; 
NF = 2;

Nmax = max(trace(:,2));

fprintf('  OF %s FixedNodes = %d  Interval time = [%d,%d] \n', 'ROME' , NF, StartHour, EndHour);

    
MSJ_SIZE = 1*MB;    % Fix the message size
[NR1, DM1, RDM1, ADT1, STM1, NDM1, TMB1] = Simulate();

MSJ_SIZE = 2*MB;    % Fix the message size
[NR2, DM2, RDM2, ADT2, STM2, NDM2, TMB2] = Simulate();

MSJ_SIZE = 5*MB;    % Fix the message size
[NR5, DM5, RDM5, ADT5, STM5, NDM5, TMB5] = Simulate();

MSJ_SIZE = 10*MB;    % Fix the message size
[NR10, DM10, RDM10, ADT10, STM10, NDM10, TMB10] = Simulate();

MSJ_SIZE = 20*MB;    % Fix the message size
[NR20, DM20, RDM20, ADT20, STM20, NDM20, TMB20] = Simulate();

%MSJ_SIZE = 50*MB;    % Fix the message size
%[NR50, DM50, RDM50, ADT50, STM50, NDM50, TMB50] = Simulate();

MSJ_SIZE = 100*MB;    % Fix the message size
[NR100, DM100, RDM100, ADT100, STM100, NDM100, TMB100] = Simulate();

MSJ_SIZE = 500*MB;    % Fix the message size
[NR0, DM0, RDM0, ADT0, STM0, NDM0, TMB0] = Simulate();    

SUBPLOTS = true;

if SUBPLOTS
    subplot(2,3,1);
else
    figure;
end

plot(Tam_ITEMs/MB,RDM0,'k');
hold on;
plot(Tam_ITEMs/MB,RDM20,'--');
hold on;
plot(Tam_ITEMs/MB,RDM10,'-.');
hold on;
plot(Tam_ITEMs/MB,RDM5,'--');
hold on;
plot(Tam_ITEMs/MB,RDM2,'-.');
hold on;
plot(Tam_ITEMs/MB,RDM1,'--');
hold on;
set(gca,'FontSize',22);
ylabel('Diffusion ratio');
xlabel('Message size (MB)')
legend('Epi.','20MB','10MB','5MB','2MB','1MB');
xlim([0 Tam_ITEMs(end)/MB]);
ylim([0.5 1]);

if SUBPLOTS
    subplot(2,3,2);
else
    figure;
end

plot(Tam_ITEMs/MB,RDM20./RDM0,'--');
hold on;
plot(Tam_ITEMs/MB,RDM10./RDM0,'-.');
hold on;
plot(Tam_ITEMs/MB,RDM5./RDM0,'--');
hold on;
plot(Tam_ITEMs/MB,RDM2./RDM0,'-.');
hold on;
plot(Tam_ITEMs/MB,RDM1./RDM0,'--');
hold on;
set(gca,'FontSize',22);
ylabel('Improvement factor');
xlabel('Message size (MB)')
legend('20MB','10MB','5MB','2MB','1MB');
xlim([0 Tam_ITEMs(end)/MB]);
ylim([1 2]);

if SUBPLOTS
    subplot(2,3,3);
else
    figure;
end


plot(Tam_ITEMs/MB,ADT0,'k');
hold on;
plot(Tam_ITEMs/MB,ADT20,'--');
hold on;
plot(Tam_ITEMs/MB,ADT10,'-.');
hold on;
plot(Tam_ITEMs/MB,ADT5,'--');
hold on;
plot(Tam_ITEMs/MB,ADT2,'-.');
hold on;
plot(Tam_ITEMs/MB,ADT1,'--');
hold on;
set(gca,'FontSize',22);
ylabel('Average delivery time (s)');
xlabel('Message size (MB)')
legend('Epi.','20MB','10MB','5MB','2MB','1MB');
xlim([0 Tam_ITEMs(end)/MB]);

if SUBPLOTS
    subplot(2,3,4);
else
    figure;
end

plot(Tam_ITEMs/MB,STM0,'k');
hold on;
plot(Tam_ITEMs/MB,STM20,'--');
hold on;
plot(Tam_ITEMs/MB,STM10,'-.');
hold on;
plot(Tam_ITEMs/MB,STM5,'--');
hold on;
plot(Tam_ITEMs/MB,STM2,'-.');
hold on;
plot(Tam_ITEMs/MB,STM1,'--');
hold on;
set(gca,'FontSize',22);
ylabel('Transmitted Messages');
xlabel('Message size (MB)')
legend('Epi.','20MB','10MB','5MB','2MB','1MB');
xlim([0 Tam_ITEMs(end)/MB]);

if SUBPLOTS
    subplot(2,3,5);
else
    figure;
end


plot(Tam_ITEMs/MB,NDM0,'k');
hold on;
plot(Tam_ITEMs/MB,NDM20,'--');
hold on;
plot(Tam_ITEMs/MB,NDM10,'-.');
hold on;
plot(Tam_ITEMs/MB,NDM5,'--');
hold on;
plot(Tam_ITEMs/MB,NDM2,'-.');
hold on;
plot(Tam_ITEMs/MB,NDM1,'--');
hold on;
set(gca,'FontSize',22);
ylabel('No transmitted messages');
xlabel('Message size (MB)')
legend('Epi.','20MB','10MB','5MB','2MB','1MB');
xlim([0 Tam_ITEMs(end)/MB]);

if SUBPLOTS
    subplot(2,3,6);
else
    figure;
end


plot(Tam_ITEMs/MB,TMB0/1024^3,'k');
hold on;
plot(Tam_ITEMs/MB,TMB20/1024^3,'--');
hold on;
plot(Tam_ITEMs/MB,TMB10/1024^3,'-.');
hold on;
plot(Tam_ITEMs/MB,TMB5/1024^3,'--');
hold on;
plot(Tam_ITEMs/MB,TMB2/1024^3,'-.');
hold on;
plot(Tam_ITEMs/MB,TMB1/1024^3,'--');
hold on;
set(gca,'FontSize',22);
ylabel('Bytes transmitted (GB)');
xlabel('Message size (MB)')
legend('Epi.','20MB','10MB','5MB','2MB','1MB');
xlim([0 Tam_ITEMs(end)/MB]);

figure;
plot(Tam_ITEMs/MB,max((TMB0./Tam_ITEMs)./(RDM0*NR0),1),'k');
hold on;
plot(Tam_ITEMs/MB,(TMB20./Tam_ITEMs)./(RDM20*NR20),'--');
hold on;
plot(Tam_ITEMs/MB,(TMB10./Tam_ITEMs)./(RDM10*NR10),'-.');
hold on;
plot(Tam_ITEMs/MB,(TMB5./Tam_ITEMs)./(RDM5*NR5),'--');
hold on;
plot(Tam_ITEMs/MB,(TMB2./Tam_ITEMs)./(RDM2*NR2),'-.');
hold on;
plot(Tam_ITEMs/MB,(TMB1./Tam_ITEMs)./(RDM1*NR1),'--');
hold on;
set(gca,'FontSize',22);
ylabel('Overhead');
xlabel('Message size (MB)')
legend('Epi.','20MB','10MB','5MB','2MB','1MB');
xlim([0 Tam_ITEMs(end)/MB]);
ylim([1 2.5]);


fprintf('Mean values\n');
aRDM(1) = mean(RDM0);
aRDM(2) = mean(RDM100);
aRDM(3) = mean(RDM20);
aRDM(4) = mean(RDM10);
aRDM(5) = mean(RDM5);
aRDM(6) = mean(RDM2);
aRDM(7) = mean(RDM1);

mRDM(1) = min(RDM0);
mRDM(2) = min(RDM100);
mRDM(3) = min(RDM20);
mRDM(4) = min(RDM10);
mRDM(5) = min(RDM5);
mRDM(6) = min(RDM2);
mRDM(7) = min(RDM1);

aADT(1) = mean(ADT0);
aADT(2) = mean(ADT100);
aADT(3) = mean(ADT20);
aADT(4) = mean(ADT10);
aADT(5) = mean(ADT5);
aADT(6) = mean(ADT2);
aADT(7) = mean(ADT1);

aTMB(1) = mean(TMB0);
aTMB(2) = mean(TMB100);
aTMB(3) = mean(TMB20);
aTMB(4) = mean(TMB10);
aTMB(5) = mean(TMB5);
aTMB(6) = mean(TMB2);
aTMB(7) = mean(TMB1);


M = [500,100,20,10,5,2,1];
fprintf('  M  aRDM  mRDM     aADT    aTMB\n');
fprintf('--------------------------------\n');
for i=1:length(M)
    fprintf('%3d %5.3f %5.3f %5.2f %5.1f\n', M(i), aRDM(i),mRDM(i),aADT(i),aTMB(i)/1024^3); 
end

end


