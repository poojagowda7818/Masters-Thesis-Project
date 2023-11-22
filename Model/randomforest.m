clc;
clear all;
close all;
warning off
data=readtable('0.1/tumour0.1(50)ts.csv');
data=removevars(data,["Unnamed_0","x_index","parent","gval","alpha","flux","TumourVolume","TumourFraction","XExtent","YExtent","ZExtent","Unnamed_0_1"]);
% data=removevars(data,"Unnamed_0"); 
k=["High","Low"];
l=[1,0];
% here we encode the group high, low to 1,0 respectively
g=data.TumourPercent;
number=zeros(length(g),1);
for i=1:length(k)
    rs=ismember(g,k(i));
    number(rs)=l(i);
end
data.category_encoded=number;
data.TumourPercent=[];
cv = cvpartition(size(data,1),'k',10);
for i=1:10
    idx = cv.test(i);
    label = data(~idx,end);
    dataTrain = data(~idx,1:end);
    dataTest=data(idx,1:end);
end
testing=dataTest(1:end,1:end-1);
model = fitensemble(dataTrain,'category_encoded','Bag',100,'Tree','Type','classification');
prediction=predict(model,testing);
ms=(sum(prediction==table2array(dataTest(:,end)))/size(dataTest,1))*100;
disp(ms);
% Getting the confusion matrix and classification performance
cm = confusionchart(table2array(dataTest(:,end)),prediction);
cp = classperf(table2array(dataTest(:,end)));
cp2 = classperf(cp,prediction);
cm.Title = 'Tumour Classification Using Random forest';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
% saveas(model,'tumour0.8os.mat');

