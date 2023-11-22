clc
clear all
close all
warning off
load '0.9/tumour0.9(50)ts.mat'
imp=predictorImportance(model);
figure;
bar(imp);
title('Predictor Importance Estimates');
ylabel('Estimates');
xlabel('Predictors');
grid on;
h=gca;
h.XTickLabel=model.PredictorNames;
Tree10 = model.Trained{10};
view(Tree10,'Mode','graph');
