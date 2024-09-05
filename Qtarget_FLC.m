%% Statements
% ------------------------------------------------------------------------------------------------------------
% This .m file is the main function to train a fuzzy logic controller in the PFL-RBC,
% by tune its membership funtions (MFs) of inputs and outputs.
% The input variables are set as the total rainfall depth (PP1, PP2, …, PP9) and the storage capacitys (small, medium, large),
% and the output variable is set as the target flow (QT1, QT2, …, QT12) .
% The training data are retrieved by analysing the history rainfall-runoff processes.
% ------------------------------------------------------------------------------------------------------------

%% Reading training data
path = 'D:\PFL-RBC\Data\';
Filename = 'Target flow training data.xlsx';
xls_data = readmatrix(strcat(path,Filename), 'Sheet', 1);
data = xls_data(1:end,[1,2,3]);
N = length(data);

%% GA-optimization
for k = 1 : 18
    params{k} = char(strcat("x", string(k)));
end

lb(1:8) = [0.25,0.5,0.75,1.1,1.5,1.8,2.2,2.8];
ub(1:8) = [0.3,0.6,1,1.2,1.6,2,2.5,3];

lb(9:18) = [0.1,0.2,0.3,0.5,0.6,0.75,1.5,2,4,6];
ub(9:18) = [0.2,0.3,0.4,0.6,0.75,1.5,2,4,6,8];

nvars = length(params);
popSize = 100;
crossRate = 0.8;
maxGen = 50;
gaOpts = optimoptions('ga', 'PopulationSize', popSize, 'CrossoverFraction', crossRate, 'MaxGenerations', maxGen, 'PlotFcns', @gaplotbestf);
functn = @(opt)ObjFunc(opt, data, N);

disp('Start optimising! ...');
[opt, J] = ga(functn, nvars, [], [], [], [], lb, ub, [], gaOpts);

%-----------------------------------------------------------------------------------------
% Output the tuned MF parameters
%-----------------------------------------------------------------------------------------
disp(['PP: ', num2str(opt(1:8))]);
disp(['QT: ', num2str([0,opt(9:18)])]);
disp(['Fitness value: ', num2str(J)]);

%-----------------------------------------------------------------------------------------
% The objective function
%-----------------------------------------------------------------------------------------
function J = ObjFunc(opt, data, N)
for i = 1 : 8
    x(i) = opt(i) * 100;
    y(i+1) = opt(i+8) * 10;
end
y(1) = 0;
y(10) = opt(17) * 10;
y(11) = opt(18) * 10;

a = newfis('fuzzy target flow');
a = addvar(a,'input','PP',[0,x(8)]);
a = addmf(a,'input',1,'PP1','zmf',[0,x(1)]);
a = addmf(a,'input',1,'PP2','trimf',[x(1),x(1),x(2)]);
a = addmf(a,'input',1,'PP3','trimf',[x(1),x(2),x(3)]);
a = addmf(a,'input',1,'PP4','trimf',[x(2),x(3),x(4)]);
a = addmf(a,'input',1,'PP5','trimf',[x(3),x(4),x(5)]);
a = addmf(a,'input',1,'PP6','trimf',[x(4),x(5),x(6)]);
a = addmf(a,'input',1,'PP7','trimf',[x(5),x(6),x(7)]);
a = addmf(a,'input',1,'PP8','trimf',[x(6),x(7),x(8)]);
a = addmf(a,'input',1,'PP9','smf',[x(7),x(8)]);

a = addvar(a,'input','V',[0,70]);
a = addmf(a,'input',2,'V1','zmf',[0,17.5]);
a = addmf(a,'input',2,'V2','trimf',[0,17.5,35]);
a = addmf(a,'input',2,'V3','trimf',[17.5,35,70]);
a = addmf(a,'input',2,'V4','smf',[35,70]);


a = addvar(a,'output','QT',[0,y(11)]);
a = addmf(a,'output',1,'QT1','trimf',[0,0,y(1)]);
a = addmf(a,'output',1,'QT2','zmf',[0,y(2)]);
a = addmf(a,'output',1,'QT3','trimf',[y(1),y(2),y(3)]);
a = addmf(a,'output',1,'QT4','trimf',[y(2),y(3),y(4)]);
a = addmf(a,'output',1,'QT5','trimf',[y(3),y(4),y(5)]);
a = addmf(a,'output',1,'QT6','trimf',[y(4),y(5),y(6)]);
a = addmf(a,'output',1,'QT7','trimf',[y(5),y(6),y(7)]);
a = addmf(a,'output',1,'QT8','trimf',[y(6),y(7),y(8)]);
a = addmf(a,'output',1,'QT9','trimf',[y(7),y(8),y(9)]);
a = addmf(a,'output',1,'QT10','trimf',[y(8),y(9),y(10)]);
a = addmf(a,'output',1,'QT11','trimf',[y(9),y(10),y(11)]);
a = addmf(a,'output',1,'QT12','smf',[y(10),y(11)]);

% fuzzy rules set up
rulelist=[1 4 1 1 1;
    2 4 1 1 1;
    3 4 2 1 1;
    4 4 4 1 1;
    5 4 4 1 1;
    6 4 5 1 1;
    7 4 7 1 1;
    8 4 9 1 1;
    9 4 10 1 1;
    1 3 1 1 1;
    2 3 2 1 1;
    3 3 3 1 1;
    4 3 5 1 1;
    5 3 6 1 1;
    6 3 8 1 1;
    7 3 9 1 1;
    8 3 10 1 1;
    9 3 11 1 1;
    1 2 1 1 1;
    2 2 3 1 1;
    3 2 4 1 1;
    4 2 6 1 1;
    5 2 7 1 1;
    6 2 9 1 1;
    7 2 10 1 1;
    8 2 11 1 1;
    9 2 12 1 1;
    1 1 2 1 1;
    2 1 5 1 1;
    3 1 7 1 1;
    4 1 8 1 1;
    5 1 9 1 1;
    6 1 10 1 1;
    7 1 10 1 1;
    8 1 11 1 1;
    9 1 12 1 1];
a = addRule(a,rulelist);

% Setting up the defuzzification algorithm
a1 = setfis(a,'DefuzzMethod','mom');
writeFIS(a1,'target flow');
a2 = readfis('target flow');

%figure(1);
%plotfis(a2);
%figure(2);
%plotmf(a2,'input',1);
%figure(3);
%plotmf(a2,'output',1);
%figure(4);
%plotmf(a2,'input',2);

%showrule(a2);
%ruleview('target flow');


PP_train = data(:,1);
V_train = data(:,2);
QT_train = data(:,3);

for i = 1 : N
    QT_flc(i, 1) = evalfis([PP_train(i),V_train(i)],a2);
    if QT_train(i, 1) == 0
        if QT_flc(i, 1) <= 0.1
            error(i, 1) = 0;
        else
            error(i, 1) = 1;
        end
    else
        error(i, 1) = abs(QT_train(i, 1) - QT_flc(i, 1))/QT_train(i, 1);
    end
end
J = sum(error) / N;
end
