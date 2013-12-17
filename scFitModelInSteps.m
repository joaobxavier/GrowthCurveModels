% scFitModelInSteps - fit model to data
% 1. Fit muMax, ii0 and Yi using backbone data (series 1)
% 2. Fit Yc and kd using carbon data series (2 through 9)
% 3. Fit Yn and ni0 using nitrogen data series (10 through 17)


% preliminary parameters (already pretty well optimized
p.muMax = 0.358057;
p.ni0 = 0.818074;
p.Yc = 0.507823;
p.Yn = 2.446186;
p.ii0 = 516.96034;
p.Yi = 0.04959;
p.kd = 0.011320;
p.qR = 4710.879865;
p.qRb = 31289.796830;
p.kg = 0.225869;
p.h1 = 2.878829;
p.qG = 0.000126;
p.kdGFP = 0.024259;
p.x0 = 0.000280;
p.gi0 = 73200.000000;
p.C0Value = 3; %default
p.N0Value = 0.5; %default
p.I0Value = 1; %default

% create object to store all the time series and models
models = ModelsOfManyTimeSeries;

% load backbone and add it to 'models'
load('backbone');
p.C0Value = 3;
p.N0Value = 0.5;
p.I0Value = 1; % initial value of external iron
% y0 = [p.x0  p.C0Value  p.N0Value  p.I0Value  p.ni0  p.ii0  p.gi0];
models = models.addTimeSeries(timehoursBackBone', odBackBone,...
    p, gfpBackBone);
modeltest = ModelWithGfp(p);
%%

% add all other time series
% basedir = '/Users/xavierj/Joao/lab_results/kerry/mathematical_model_paper/aligned_mat_files/manual_alignment_finetuned/';
basedir = '../Iron Model V1_112713/';
matfile{1} = [basedir '05092013CarbonTitration1.mat'];
%matfile{1} = [basedir '05242013CarbonTitration2.mat'];
matfile{2} = [basedir '11042013CarbonTitrationIRON.mat'];
matfile{3} = [basedir '05092013NitrogenTitration1.mat'];
%matfile{3} = [basedir '05162013NitrogenTitration2.mat'];
matfile{4} = [basedir '110713NitrogenWithFe.mat'];

C0Array = [0.5 0.25 0.125 0.0625;...
    0.5 0.25 0.125 0.0625;...
    3 3 3 3;...
    3 3 3 3]; %gC
N0Array = [0.5 0.5  0.5  0.5;...
    0.5 0.5  0.5  0.5;...
    0.0625 0.02836 0.0156 0.0078;...
    0.0625 0.02836 0.0156 0.0078]; %gN
I0Array = [p.I0Value p.I0Value p.I0Value p.I0Value;...
    100 100 100 100;...
    p.I0Value p.I0Value p.I0Value p.I0Value;...
    100 100 100 100];

for j = 1:length(matfile)
    load(matfile{j});
    try 
        clear odData2 gfpData2
    catch
    end
    for i = 1:4
        odData2(:, i) = median(oddata(:, (1:7) + (i-1)*7), 2);
        gfpData2(:, i) = median(gfpdata(:, (1:7) + (i-1)*7), 2);
    end
    % background correction
    for i = 1:84
        backIndex = floor((i-1)/7) + 1;
        backIndex = mod(backIndex, 4) + 1;
        gfpdata(:,i) = gfpdata(:, i) - gfpData2(:,backIndex);
    end
    % add titration to the object 'models'
    for column = 1:4
        p.C0Value = C0Array(j, column);
        p.N0Value = N0Array(j, column);
        p.I0Value = I0Array(j, column); % initial value of external iron
%         y0 = [p.x0, p.C0Value, p.N0Value,...
%             p.I0Value, p.ni0, p.ii0, p.gi0];
        models = models.addTimeSeries(timehours',...
            oddata(:, (1:7) + 56 + (column-1)*7),...
            p, gfpdata(:, (1:7) + 56 + (column-1)*7));
    end
end

models.solveModel;

%%
figure(1)
subplot(2, 2, 1);
models.plotAllModels(1:5);
title('Carbon titration, no iron');
subplot(2, 2, 2);
models.plotAllModels([1 6:9])
title('Carbon titration + iron');
subplot(2, 2, 3);
models.plotAllModels([1 10:13]);
title('Nitrogen titration, no iron');
subplot(2, 2, 4);
models.plotAllModels([1 14:17])
title('Nitrogen titration + iron');

%% Fit muMax, ii0 and Yi using backbone data (series 1)
models = models.fitParameters({'ii0', 'Yi'}, 1);
% models = models.fitParameters({'muMax', 'ii0', 'Yi'}, 1);
%%
figure(2)
models.plotAllModels(1)
title('Fit of muMax, ii0 and Yi using backbone data')

%% Fit Yc and kd using also the carbon data series (1 + 2 through 9)
models = models.fitParameters({'muMax', 'Yc', 'kd'}, 1:9);


figure(3)
subplot(2, 1, 1);
models.plotAllModels(1:5);
title('Carbon titration, no iron');
subplot(2, 1, 2);
models.plotAllModels([1 6:9])
title('Carbon titration + iron');


%% Fit Yn and ni0 using all data series (1 through 17)
models = models.fitParameters({'Yn', 'ni0'}, [1:17]);

figure(4)
subplot(2, 1, 1);
models.plotAllModels([1 10:13]);
title('Nitrogen titration, no iron');
subplot(2, 1, 2);
models.plotAllModels([1 14:17])
title('Nitrogen titration + iron');


%% fit gfp parameters (qRb qR kg h1 kdGFP) using all data
models = models.optimizeParametersForGfpRate(...
    {'qRb','qR', 'kg', 'h1', 'kdGFP'}, [1:17]);

%%
sim = 2;
figure(5)

subplot(3, 1, 1);
models.plotAllModels(sim, 'dgidt');
set(gca, 'YLim', [-10000 20000]);
set(gca, 'YScale', 'lin');
set(gca, 'Color', [0.2 0.2 0.2]);
subplot(3, 1, 2);
models.plotAllModels(sim, 'dlogxdt');
set(gca, 'YLim', [-0.2 0.8]);
set(gca, 'YScale', 'lin');
set(gca, 'Color', [0.2 0.2 0.2]);
subplot(3, 1, 3);

models.plotAllModels(sim, 'gfp');
set(gca, 'Color', [0.2 0.2 0.2]);
set(gca, 'YLim', [5 8e4]);
set(gca, 'YScale', 'lin');

%%
figure(6);
models.arrayOfModels(18).model.plotDataAndModel([0 1 1], 'gi')
set(gca, 'Color', [0.2 0.2 0.2]);
set(gca, 'YLim', [5 8e4]);
set(gca, 'YScale', 'lin');

%%
models = models.changeParameterValue('C0Value',[.49,.49], [2,6]);
models = models.fitParametersSpecModels('C0Value',[2,6]);
