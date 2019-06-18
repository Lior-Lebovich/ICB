%% Fig. 5: The recurrent spiking network model

addpath('analyses','data');
figure;


%% Fig 5b - traces:

% load data: (col 1 = time;  cols 2-5 = voltage of 4 different neurons)
filName_voltage_E = 'v1e-densesparse-test100-k=400-gee=0.3-gei=1.5-gie=2.0-gii=2-ie=0.06-ii=0.04-ne=32000-ni=8000-t=100s-epsilon=3.4-traces-spont-net01-r01';
filName_voltage_I = 'v1i-densesparse-test100-k=400-gee=0.3-gei=1.5-gie=2.0-gii=2-ie=0.06-ii=0.04-ne=32000-ni=8000-t=100s-epsilon=3.4-traces-spont-net01-r01';
voltage_E = load( filName_voltage_E );
voltage_I = load( filName_voltage_I );
% plot:
subplot(8,8,[33 34]);
plot(voltage_E(:,1)-0.65-3,voltage_E(:,2),'r','lineWidth',0.5);
xlim([0,3]); ylim([-85,25]); box off; set(gca, 'FontName', 'Times', 'FontSize', 5);
yticks([-60,0]); ylabel('Voltage [mV]'); xticks(0:1:3);
subplot(8,8,[41 42]);
plot(voltage_I(:,1)-0.65-3,voltage_I(:,5),'b','lineWidth',0.5); 
xlim([0,3]); ylim([-85,25]); box off; set(gca, 'FontName', 'Times', 'FontSize', 5);
yticks([-60,0]); xlabel('Time [s]'); xticks(0:1:3);


%% Fig 5c - firing rate distribution:

% load data:
fileName_rateE = 'fire-densesparse-test100-k=400-gee=0.3-gei=1.5-gie=2.0-gii=2-ie=0.06-ii=0.04-ne=32000-ni=8000-t=100s-epsilon=3.0-FIRESTAT-SPONT-NET01-R01';
fileName_rateI = 'firi-densesparse-test100-k=400-gee=0.3-gei=1.5-gie=2.0-gii=2-ie=0.06-ii=0.04-ne=32000-ni=8000-t=100s-epsilon=3.0-FIRESTAT-SPONT-NET01-R01';
rateE_all = load( fileName_rateE );
rateI_all = load( fileName_rateI );
rateE = log10( rateE_all(:,2) ); 
rateI = log10( rateI_all(:,2) );

[h_rate, edges_FR] = histcounts(rateE,110,'Normalization', 'pdf', 'binEdges', linspace(-1,1,26));
rate_x = edges_FR(1:(end-1)) + 0.5 * (edges_FR(2)-edges_FR(1)); 
[h_rate_I, edges_FR_I] = histcounts(rateI,110,'Normalization', 'pdf', 'binEdges', linspace(-1,1,26));
rate_x_I = edges_FR_I(1:(end-1)) + 0.5 * (edges_FR_I(2)-edges_FR_I(1)); 
subplot(8,8,[49 50 57 58]);
area(rate_x, h_rate, 'FaceColor','red','EdgeColor','none' ); hold on;
area(rate_x_I, h_rate_I, 'FaceColor','blue','EdgeColor','none' );
box off; xlim([-1,1]); ylim([0,2]); 
ggg = gca;
xticks(-1:1:1); yticks(0:1:2);
ggg.XMinorTick = 'on';
set(gca, 'FontName', 'Times', 'FontSize', 5);
alpha 0.5;
xlabel('log_{10}(FR) [Hz]'); ylabel('PDF');


%% Fig 5d-e - example impossible trials:

random_indices_E = (1:1000);
random_indices_I = (1:1000);
col_U_E = 0.01 * [50 0 0];
col_D_E = 0.01 * [85 0 0];
col_U_I = 0.01 * [0 0 50];
col_D_I = 0.01 * [0 0 85];


%% Fig 5d - example impossible trials:

% read data:
fileName_E_raster_1 = 'spike-choice-test100-k=400-gee=0.3-gei=1.5-gie=2.0-gii=2-ie=0.06-ii=0.04-ne=32000-ni=8000-epsilon=3.0-DECISION-1TRIAL-NET01-R01';
fileName_I_raster_1 = 'spiki-choice-test100-k=400-gee=0.3-gei=1.5-gie=2.0-gii=2-ie=0.06-ii=0.04-ne=32000-ni=8000-epsilon=3.0-DECISION-1TRIAL-NET01-R01';
fileName_activity_1 = 'ux-choice-test100-k=400-gee=0.3-gei=1.5-gie=2.0-gii=2-ie=0.06-ii=0.04-ne=32000-ni=8000-epsilon=3.0-DECISION-1TRIAL-NET01-R01';
E_raster_1 = load( fileName_E_raster_1 );
I_raster_1 = load( fileName_I_raster_1 );
activity_1 = load( fileName_activity_1 );
decisionTime_1 = activity_1(min(find(((activity_1(:,3)-activity_1(:,2))./(activity_1(:,3)+activity_1(:,2)))>=0.4)),1)-0.5;

%plot:
% E:
subplot(8,8,3:5);
stem(E_raster_1(ismember(E_raster_1(:,2),random_indices_E),1)-2,...
E_raster_1(ismember(E_raster_1(:,2),random_indices_E),2),'marker','.','lineStyle','none','color',col_U_E,'MarkerSize',1)
xlim([-1,1]); ylim([0,max(random_indices_E)]); 
set(gca,'ytick',[]); xticks(-1:0.5:1); xticklabels([]); set(gca, 'FontName', 'Times', 'FontSize', 5);
ylabel('U_E');
subplot(8,8,11:13);
stem(E_raster_1(ismember(E_raster_1(:,2),random_indices_E+16000),1)-2,...
-16000+E_raster_1(ismember(E_raster_1(:,2),random_indices_E+16000),2),'marker','.','lineStyle','none','color',col_D_E,'MarkerSize',1)
xlim([-1,1]); ylim([0,max(random_indices_E)]); 
xticks(-1:0.5:1); xticklabels([]); set(gca,'ytick',[]); set(gca, 'FontName', 'Times', 'FontSize', 5);
ylabel('D_E');
% acticity:
subplot(8,8,[19:21, 27:29]);
fill([0,decisionTime_1-1.5,decisionTime_1-1.5,0,0],[0,0,9,9,0],[0.85 0.85 0.85],'lineStyle','none'); hold on;
plot(activity_1(:,1)-2,activity_1(:,2),'Color',col_U_E,'lineWidth',1); hold on;
plot(activity_1(:,1)-2,activity_1(:,3),'Color',col_D_E,'lineWidth',1)
xlim([-1,1]); ylim([0,6]); box off;  set(gca, 'FontName', 'Times', 'FontSize', 5);
ylabel('Activity [Hz]');
% I:
subplot(8,8,35:37);
stem(I_raster_1(ismember(I_raster_1(:,2),random_indices_I),1)-2,...
I_raster_1(ismember(I_raster_1(:,2),random_indices_I),2),'marker','.','lineStyle','none','color',col_U_I,'MarkerSize',1)
xlim([-1,1]); ylim([0, max(random_indices_I)]); 
xticks(-1:0.5:1); xticklabels([]); set(gca,'ytick',[]); set(gca, 'FontName', 'Times', 'FontSize', 5);
ylabel('U_I');
subplot(8,8,43:45);
stem(I_raster_1(ismember(I_raster_1(:,2),random_indices_I+4000),1)-2,...
-4000+I_raster_1(ismember(I_raster_1(:,2),random_indices_I+4000),2),'marker','.','lineStyle','none','color',col_D_I,'MarkerSize',1)
xlim([-1,1]); ylim([0, max(random_indices_I)]); 
xticks(-1:0.5:1); xticklabels([]); set(gca,'ytick',[]); set(gca, 'FontName', 'Times', 'FontSize', 5);
ylabel('D_I');
% acticity:
subplot(8,8,[51:53, 59:61]);
fill([0,decisionTime_1-1.5,decisionTime_1-1.5,0,0],[0,0,9,9,0],[0.85 0.85 0.85],'lineStyle','none'); hold on;
plot(activity_1(:,1)-2,activity_1(:,4),'Color',col_U_I,'lineWidth',1); hold on;
plot(activity_1(:,1)-2,activity_1(:,5),'Color',col_D_I,'lineWidth',1)
xlim([-1,1]); ylim([0,9]); box off; yticks(0:3:9);
ylabel('Activity [Hz]'); xlabel('Time [s]'); set(gca, 'FontName', 'Times', 'FontSize', 5);


%% Fig 5e - example impossible trial:

% read data:
fileName_E_raster_2 = 'spike-choice-test100-k=400-gee=0.3-gei=1.5-gie=2.0-gii=2-ie=0.06-ii=0.04-ne=32000-ni=8000-epsilon=3.0-decision-another1trial-net01-r01';
fileName_I_raster_2 = 'spiki-choice-test100-k=400-gee=0.3-gei=1.5-gie=2.0-gii=2-ie=0.06-ii=0.04-ne=32000-ni=8000-epsilon=3.0-decision-another1trial-net01-r01';
fileName_activity_2 = 'ux-choice-test100-k=400-gee=0.3-gei=1.5-gie=2.0-gii=2-ie=0.06-ii=0.04-ne=32000-ni=8000-epsilon=3.0-decision-another1trial-net01-r01';
E_raster_2 = load( fileName_E_raster_2 );
I_raster_2 = load( fileName_I_raster_2 );
activity_2 = load( fileName_activity_2 );
decisionTime_2 = activity_2(min(find(((activity_2(:,2)-activity_2(:,3))./(activity_2(:,2)+activity_2(:,3)))>=0.4)),1)-0.5;

%plot:
% E:
subplot(8,8,6:8);
stem(E_raster_2(ismember(E_raster_2(:,2),random_indices_E),1)-2,...
E_raster_2(ismember(E_raster_2(:,2),random_indices_E),2),'marker','.','lineStyle','none','color',col_U_E,'MarkerSize',1)
xlim([-1,1]); ylim([0,max(random_indices_E)]); 
xticks(-1:0.5:1); xticklabels([]); set(gca,'ytick',[]); set(gca, 'FontName', 'Times', 'FontSize', 5);
ylabel('U_E');
subplot(8,8,14:16);
stem(E_raster_2(ismember(E_raster_2(:,2),random_indices_E+16000),1)-2,...
-16000+E_raster_2(ismember(E_raster_2(:,2),random_indices_E+16000),2),'marker','.','lineStyle','none','color',col_D_E,'MarkerSize',1)
xlim([-1,1]); ylim([0,max(random_indices_E)]); 
xticks(-1:0.5:1); xticklabels([]); set(gca,'ytick',[]);
ylabel('D_E');set(gca, 'FontName', 'Times', 'FontSize', 5); set(gca, 'FontName', 'Times', 'FontSize', 5);
% acticity:
subplot(8,8,[22:24, 30:32]);
fill([0,decisionTime_2-1.5,decisionTime_2-1.5,0,0],[0,0,9,9,0],[0.85 0.85 0.85],'lineStyle','none'); hold on;
plot(activity_2(:,1)-2,activity_2(:,2),'Color',col_U_E,'lineWidth',1); hold on;
plot(activity_2(:,1)-2,activity_2(:,3),'Color',col_D_E,'lineWidth',1)
xlim([-1,1]); ylim([0,6]); box off; set(gca, 'FontName', 'Times', 'FontSize', 5);
ylabel('Activity [Hz]');
% I:
subplot(8,8,38:40);
stem(I_raster_2(ismember(I_raster_2(:,2),random_indices_I),1)-2,...
I_raster_2(ismember(I_raster_2(:,2),random_indices_I),2),'marker','.','lineStyle','none','color',col_U_I,'MarkerSize',1)
xlim([-1,1]); ylim([0, max(random_indices_I)]); 
xticks(-1:0.5:1); xticklabels([]); set(gca,'ytick',[]); set(gca, 'FontName', 'Times', 'FontSize', 5);
ylabel('U_I');
subplot(8,8,46:48);
stem(I_raster_2(ismember(I_raster_2(:,2),random_indices_I+4000),1)-2,...
-4000+I_raster_2(ismember(I_raster_2(:,2),random_indices_I+4000),2),'marker','.','lineStyle','none','color',col_D_I,'MarkerSize',1)
xlim([-1,1]); ylim([0, max(random_indices_I)]); 
xticks(-1:0.5:1); xticklabels([]); set(gca,'ytick',[]); set(gca, 'FontName', 'Times', 'FontSize', 5);
ylabel('D_I');
% acticity:
subplot(8,8,[54:56, 62:64]);
fill([0,decisionTime_2-1.5,decisionTime_2-1.5,0,0],[0,0,9,9,0],[0.85 0.85 0.85],'lineStyle','none'); hold on;
plot(activity_2(:,1)-2,activity_2(:,4),'Color',col_U_I,'lineWidth',1); hold on;
plot(activity_2(:,1)-2,activity_2(:,5),'Color',col_D_I,'lineWidth',1)
xlim([-1,1]); ylim([0,9]); box off; yticks(0:3:9);
ylabel('Activity [Hz]'); xlabel('Time [s]'); set(gca, 'FontName', 'Times', 'FontSize', 5);
