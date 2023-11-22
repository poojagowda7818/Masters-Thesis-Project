clc;
clear all;
close all;
warning off
load 'tumour 20/tumour(20)180.mat'  
% Read the CSV file data
data2 = readtable('0.9/Test.csv');
% data2=removevars(data2,["Var1","TumourVolume","TumourFraction","XExtent","YExtent","ZExtent"]);
data2=removevars(data2,["Unnamed_0","x_index","parent","gval","alpha","flux","TumourVolume","TumourFraction","XExtent","YExtent","ZExtent"]);
k=["High","Low"];
l=[1,0];
g= data2.TumourPercent;
number=zeros(length(g),1);
    for j=1:length(k)
        rs=ismember(g,k(j));
        number(rs)=l(j);
    end
data2.category_encoded=number;
data2.TumourPercent=[];
testing2=data2(1:end,1:end-1);
prediction=predict(model,testing2);
ms=(sum(prediction==table2array(data2(:,end)))/size(data2,1))*100; 
cm = confusionchart(table2array(data2(:,end)),prediction);
cm.FontSize = 12;
cm2 = cfmatrix2(table2array(data2(:,end)),prediction,[0,1],0,1);
cm.Title = 'Tumour Classification Using Random forest';
