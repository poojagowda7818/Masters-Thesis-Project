clc;
clear all;
close all;
warning off
load '0.1/tumour0.1(50)ts.mat'
inputFolder = '0.5/50';
% Set the folder path for the output CSV files
outputFolder = '0.5/result50';
% Get a list of all the CSV files in the input folder
fileList = dir(fullfile(inputFolder, '*.csv'));
% Loop through each CSV file and read its data
for i = 1:length(fileList)
    % Construct the input file path
    inputFilePath = fullfile(inputFolder, fileList(i).name);
    
    % Read the CSV file data
    data2 = readtable(inputFilePath);
%     data2=removevars(data2,["TumourVolume","TumourFraction","XExtent","YExtent","ZExtent"]);
    data2=removevars(data2,["Var1","x_index","parent","gval","alpha","flux","TumourVolume","TumourFraction","XExtent","YExtent","ZExtent"]);
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
    % testing2 = normalize(testing2);
    prediction=predict(model,testing2);
    ms=(sum(prediction==table2array(data2(:,end)))/size(data2,1))*100; 
    cm = confusionchart(table2array(data2(:,end)),prediction);
    cm2 = cfmatrix2(table2array(data2(:,end)),prediction,[0,1],0,1);
    cm.Title = 'Tumour Classification Using Random forest';
    % Construct the output file path
    outputFilePath = fullfile(outputFolder, fileList(i).name);
    % Write the CSV file data to the output file
    writematrix(prediction, outputFilePath);
%     writematrix(prediction,'files/total/results 90%/028-1.csv')
     d=readtable(inputFilePath);
     d2=readtable(outputFilePath);
     writetable(horzcat(d(:,2:3),d2),outputFilePath); 
end

