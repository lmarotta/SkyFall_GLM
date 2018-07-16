clear 
close all
figure(1);
%set the number of cycles for the test

n_cycles=500; % for rapid inspections, 10-20 cycles are enough to have an idea

for nindex=1:n_cycles
load('RealFalls.mat'); %% FallsH
load('DeeptableNoFallsFinal.mat');
load('DeepTableNoSide5sec.mat');
load('index.mat');

%reshuffle ADLs for each cycle
l=randperm(size(DeepTableNoFalls,1),size(DeepTableNoFalls,1));
DeepTableNoFalls=DeepTableNoFalls(l,1:end);

%reshuffle simulated falls for each cycle
l=randperm(size(DeepTableNoSide5Sec,1),size(DeepTableNoSide5Sec,1));
DeepTableNoSide5Sec=DeepTableNoSide5Sec(l,1:end);


% create Training Set

TestTableNoFalls=DeepTableNoFalls(1:8000,:); %(1:size(DeepTableNoSide5Sec,1),:);

% create feature evaluation set

Ta1=table2array(TestTableNoFalls(:,2:end));
Sub1=table2array(TestTableNoFalls(:,1));

% randomly select as many falls as no falls we have
DeepTableNoSide5Sec2=DeepTableNoSide5Sec(1:600,:);
Test1=table2array(DeepTableNoSide5Sec2(:,5:end));
SubjectARRAY1=table2array(DeepTableNoSide5Sec2(:,2));
Sub_Arr=[Sub1;SubjectARRAY1];
F_0=zeros(size(Ta1,1),1);
F_1=ones(size(Test1,1),1);
F=[F_0;F_1];
TestData= [Ta1;Test1];
for i=1:size(TestData,2)
    media(:,i)=mean(TestData(:,i));
    sd(:,i)=std(TestData(:,i));
end


TableTest=array2table(TestData);
a1={'Variance','Angle','Kurtosis_a','Skewness_a','SD_a','Mean_a','Median_a','IQR_a','Kurtosis_g','Skewness_g','SD_g','Mean_g','Median_g','IQR_g','RMS','EnergyX','EnergyY','EnergyZ','EnergyG','Acc_steepness_afterpeak','Gyro_steepness_afterpeak','Acc_steepness_afterpeak_X','Acc_steepness_afterpeak_Y','Acc_steepness_afterpeak_Z','Kurtosis_x_a','Skewness_x_a','IQR_x_a','Kurtosis_y_a','Skewness_y_a','IQR_y_a','Kurtosis_z_a','Skewness_z_a','IQR_z_a','Kurtosis_x_g','Skewness_x_g','IQR_x_g','Kurtosis_y_g','Skewness_y_g','IQR_y_g','Kurtosis_z_g','Skewness_z_g','IQR_z_g','Max_f','Periodogram_maxf','Skewness_fft','Kurtosis_fft','S_Entropy','NF1','NF2','NF3','NF4','NF5','NF6','NF7','NF8','NF9','NF10','NF11','NF12','NF13','NF14','NF15','NF16','NF17','NF18','maxOrienx','varOrienx','maxOrieny','varOrieny','maxOrienz','varOrienz'};
TableTest.Properties.VariableNames = a1; 

a2={'Fall_Outcome','Subject'};

Table21=table(F);
a2={'Fall_Outcome'};
Table21.Properties.VariableNames = a2; 

EvalTable1=[Table21 TableTest];

% load('in_1.mat');
index1=[1 index];
EvalTable1=EvalTable1(:,index1);


RandomForrest = TreeBagger(200 ,EvalTable1,'Fall_Outcome','Cost',[0 0.85 ;0.15 0]);

%validation table
Variables=DeepTableNoFalls.Properties.VariableNames(index);
TestADLTable=DeepTableNoFalls(8001:11001,index);
TT=table2array(TestADLTable);
Sub1=table2array(DeepTableNoFalls(8001:11001,1));
% index=index+3 % only for laboratory falls!!
DeepTableNoSide5Sec2=DeepTableNoSide5Sec(601:680,index+3);
SubjectARRAY1=table2array(DeepTableNoSide5Sec2(:,1));
Sub_Arr=[Sub1;SubjectARRAY1];
F_0=zeros(size(TT,1),1);
F_1=ones(size(DeepTableNoSide5Sec2,1),1);
F=[F_0;F_1];
DeepTableNoSide5Sec2=table2array(DeepTableNoSide5Sec2);
TestData= [TT;DeepTableNoSide5Sec2];

X2=TestData;
yval2=F;

[prediction2, classifScore2]=RandomForrest.predict(X2);



[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.04; 0.91 0]);
SP1(nindex)=1-OPTROCPT(1);
SE1(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'ro')

T1=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr1(nindex)=T1;

[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.09; 0.90 0]);
SP2(nindex)=1-OPTROCPT(1);
SE2(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'bo')
T2=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr2(nindex)=T2;


[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.12; 0.88 0]);
SP3(nindex)=1-OPTROCPT(1);
SE3(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'co')
T3=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr3(nindex)=T3;

[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.14; 0.85 0]);
SP4(nindex)=1-OPTROCPT(1);
SE4(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'ko')
T4=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr4(nindex)=T4;


[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.15; 0.84 0]);
SP5(nindex)=1-OPTROCPT(1);
SE5(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'yo')
T5=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr5(nindex)=T5;


[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.16; 0.86 0]);
SP6(nindex)=1-OPTROCPT(1);
SE6(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'wo')
T6=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr6(nindex)=T6;


[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.18; 0.84 0]);
SP7(nindex)=1-OPTROCPT(1);
SE7(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'wo')
T7=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr7(nindex)=T7;

[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.20; 0.80 0]);
SP8(nindex)=1-OPTROCPT(1);
SE8(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'wo')
T8=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr8(nindex)=T8;

[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.25; 0.75 0]);
SP9(nindex)=1-OPTROCPT(1);
SE9(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'ro')

T9=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr9(nindex)=T9;

[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.32; 0.68 0]);
SP10(nindex)=1-OPTROCPT(1);
SE10(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'bo')
T10=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr10(nindex)=T10;


[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.50; 0.50 0]);
SP11(nindex)=1-OPTROCPT(1);
SE11(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'co')
T11=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr11(nindex)=T11;

[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.60; 0.40 0]);
SP12(nindex)=1-OPTROCPT(1);
SE12(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'ko')
T12=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr12(nindex)=T12;


[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.65; 0.35 0]);
SP13(nindex)=1-OPTROCPT(1);
SE13(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'yo')
T13=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr13(nindex)=T13;


[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.75; 0.25 0]);
SP14(nindex)=1-OPTROCPT(1);
SE14(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'wo')
T14=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr14(nindex)=T14;


[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.85; 0.15 0]);
SP15(nindex)=1-OPTROCPT(1);
SE15(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'wo')
T15=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr15(nindex)=T15;

[XA,YA,TA,AUC,OPTROCPT,suby] = perfcurve(yval2,classifScore2(:,2),1,'Cost',[0 0.90; 0.10 0]);
SP16(nindex)=1-OPTROCPT(1);
SE16(nindex)=OPTROCPT(2);
AUC;
hold on
plot(XA,YA,'m')
hold on;
plot(OPTROCPT(1),OPTROCPT(2),'wo')
T16=TA((XA==OPTROCPT(1))&(YA==OPTROCPT(2)));
Tr16(nindex)=T16;

%test data

Variables=DeepTableNoFalls.Properties.VariableNames(index);
TestADLTable=DeepTableNoFalls(end-4830:end,index);
TT=table2array(TestADLTable);
Sub1=table2array(DeepTableNoFalls(end-4830:end,1));
Test1=table2array(RealFalls(:,index));
SubjectARRAY1=table2array(RealFalls(:,1));
Sub_Arr=[Sub1;SubjectARRAY1];
F_0=zeros(size(TT,1),1);
F_1=ones(size(Test1,1),1);
F=[F_0;F_1];
TestData= [TT;Test1];

TableTest=array2table(TestData);
TableTest.Properties.VariableNames = Variables;

Table2=table(F,Sub_Arr);
a2={'Fall_Outcome','Subject'};
Table2.Properties.VariableNames = a2; 

TestingTable=[Table2 TableTest];
X2=table2array(TestingTable(:,3:end));
yval2=table2array(TestingTable(:,1));

[prediction2, classifScore2]=RandomForrest.predict(X2);

fp = sum((classifScore2(:,2) >= T1) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T1) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T1)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T1) & (yval2 == 0));
% precision1=tp/(tp+fp);
% recall1=tp/(tp+fn);
% myF1Tree(nindex)=2*precision1*recall1/(precision1+recall1);
myse1(nindex)=tp/(tp+fn);
mysp1(nindex)=tn/(tn+fp);
 
fp = sum((classifScore2(:,2) >= T2) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T2) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T2)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T2) & (yval2 == 0));
myse2(nindex)=tp/(tp+fn);
mysp2(nindex)=tn/(tn+fp);

fp = sum((classifScore2(:,2) >= T3) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T3) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T3)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T3) & (yval2 == 0));
myse3(nindex)=tp/(tp+fn);
mysp3(nindex)=tn/(tn+fp);
 
fp = sum((classifScore2(:,2) >= T4) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T4) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T4)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T4) & (yval2 == 0));
myse4(nindex)=tp/(tp+fn);
mysp4(nindex)=tn/(tn+fp);

fp = sum((classifScore2(:,2) >= T5) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T5) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T5)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T5) & (yval2 == 0));
myse5(nindex)=tp/(tp+fn);
mysp5(nindex)=tn/(tn+fp);
 
fp = sum((classifScore2(:,2) >= T6) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T6) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T6)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T6) & (yval2 == 0));
myse6(nindex)=tp/(tp+fn);
mysp6(nindex)=tn/(tn+fp);

fp = sum((classifScore2(:,2) >= T7) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T7) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T7)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T7) & (yval2 == 0));
myse7(nindex)=tp/(tp+fn);
mysp7(nindex)=tn/(tn+fp);

fp = sum((classifScore2(:,2) >= T8) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T8) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T8)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T8) & (yval2 == 0));
myse8(nindex)=tp/(tp+fn);
mysp8(nindex)=tn/(tn+fp);

fp = sum((classifScore2(:,2) >= T9) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T9) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T9)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T9) & (yval2 == 0));
% precision1=tp/(tp+fp);
% recall1=tp/(tp+fn);
% myF1Tree(nindex)=2*precision1*recall1/(precision1+recall1);
myse9(nindex)=tp/(tp+fn);
mysp9(nindex)=tn/(tn+fp);
 
fp = sum((classifScore2(:,2) >= T10) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T10) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T10)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T10) & (yval2 == 0));
myse10(nindex)=tp/(tp+fn);
mysp10(nindex)=tn/(tn+fp);

fp = sum((classifScore2(:,2) >= T11) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T11) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T11)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T11) & (yval2 == 0));
myse11(nindex)=tp/(tp+fn);
mysp11(nindex)=tn/(tn+fp);
 
fp = sum((classifScore2(:,2) >= T12) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T12) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T12)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T12) & (yval2 == 0));
myse12(nindex)=tp/(tp+fn);
mysp12(nindex)=tn/(tn+fp);

fp = sum((classifScore2(:,2) >= T13) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T13) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T13)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T13) & (yval2 == 0));
myse13(nindex)=tp/(tp+fn);
mysp13(nindex)=tn/(tn+fp);
 
fp = sum((classifScore2(:,2) >= T14) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T14) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T14)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T14) & (yval2 == 0));
myse14(nindex)=tp/(tp+fn);
mysp14(nindex)=tn/(tn+fp);

fp = sum((classifScore2(:,2) >= T15) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T15) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T15)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T15) & (yval2 == 0));
myse15(nindex)=tp/(tp+fn);
mysp15(nindex)=tn/(tn+fp);

fp = sum((classifScore2(:,2) >= T16) & (yval2 == 0));
tp = sum((classifScore2(:,2) >= T16) & (yval2 == 1));
fn = sum((classifScore2(:,2) <= T16)& (yval2 == 1));
tn = sum((classifScore2(:,2) <= T16) & (yval2 == 0));
myse16(nindex)=tp/(tp+fn);
mysp16(nindex)=tn/(tn+fp);

myAUCTREE(nindex)=AUC;

clearvars -except myAUCTREE myF1Tree Tr1 Tr2 Tr3 Tr4 Tr5 Tr6 Tr7 Tr8 myse1 mysp1 myse2 mysp2 myse3 mysp3 myse4 mysp4 myse5 mysp5 myse6 mysp6 myse7 mysp7 myse8 mysp8 Tr9 Tr10 Tr11 Tr12 Tr13 Tr15 Tr14 Tr16 myse9 mysp9 myse10 mysp10 myse11 mysp11 myse12 mysp12 myse13 mysp13 myse14 mysp14 myse15 mysp15 myse16 mysp16 SP1 SE1 SP2 SE2 SP3 SE3 SP4 SE4 SP5 SE5 SP6 SE6 SP7 SE7 SP8 SE8 SP9 SE9 SP10 SE10 SP11 SE11 SP12 SE12 SP13 SE13 SP14 SE14 SP15 SE15 SP16 SE16

end
title('ROC curves for each classifier');
xlabel('True positive rate');
ylabel('True negative rate');
legend('Tree','optimalpointTree');

meanAUC=mean(myAUCTREE);

SPTest = [ mysp1;  mysp2; mysp3; mysp4;  mysp5 ; mysp6;  mysp7; mysp8 ;mysp9 ;mysp10; mysp11 ; mysp12; mysp13; mysp14 ; mysp15 ; mysp16];
SETest=  [ myse1;  myse2; myse3; myse4 ; myse5  ;myse6 ; myse7; myse8; myse9; myse10; myse11 ; myse12; myse13 ;myse14 ; myse15 ; myse16];

SPval = [SP1; SP2; SP3; SP4; SP5; SP6; SP7; SP8; SP9 ;SP10 ;SP11; SP12 ;SP13; SP14; SP15; SP16];
SEval = [SE1; SE2; SE3 ;SE4; SE5; SE6; SE7; SE8 ;SE9; SE10; SE11 ;SE12; SE13; SE14 ;SE15; SE16];

for i=1:16
    SD_SEVal(:,i)=std(SEval(i,:));
    SD_SPVal(:,i)=std(SPval(i,:));
    Mean_SEVal(:,i)=mean(SEval(i,:));
    Mean_SPVal(:,i)=mean(SPval(i,:));
end

for i=1:16
    SD_SETest(:,i)=std(SETest(i,:));
    SD_SPTest(:,i)=std(SPTest(i,:));
    Mean_SETest(:,i)=mean(SETest(i,:));
    Mean_SPTest(:,i)=mean(SPTest(i,:));
end

% meanF1score=mean(myF1Tree);

meanSE1=mean(myse1);

meanSP1=mean(mysp1);

meanSE2=mean(myse2);

meanSP2=mean(mysp2);

meanSE3=mean(myse3);

meanSP3=mean(mysp3);

meanSE4=mean(myse4);

meanSP4=mean(mysp4);

meanSE5=mean(myse5);

meanSP5=mean(mysp5);

meanSE6=mean(myse6);

meanSP6=mean(mysp6);

meanSE7=mean(myse7)

meanSP7=mean(mysp7)

meanSE8=mean(myse8)

meanSP8=mean(mysp8)

%standard dev of sensitivity 

sdSE1=std(myse1);

sdSE2=std(myse2);

sdSE3=std(myse3);

sdSE4=std(myse4);

sdSE5=std(myse5);

sdSE6=std(myse6);

sdSE7=std(myse7);

sdSE8=std(myse8);

meanSE9=mean(myse9);

meanSP9=mean(mysp9);

meanSE10=mean(myse10);

meanSP10=mean(mysp10);

meanSE11=mean(myse11);

meanSP11=mean(mysp11);

meanSE12=mean(myse12);

meanSP12=mean(mysp12);

meanSE14=mean(myse14);

meanSP14=mean(mysp14);

meanSE15=mean(myse15);

meanSP15=mean(mysp15);

meanSE13=mean(myse13)

meanSP13=mean(mysp13)

meanSE16=mean(myse16)

meanSP16=mean(mysp16)

%standard dev of sensitivity 

sdSE9=std(myse9);

sdSE10=std(myse10);

sdSE11=std(myse11);

sdSE12=std(myse12);

sdSE13=std(myse13);

sdSE14=std(myse14);

sdSE15=std(myse15);

sdSE16=std(myse16);

%validation set

mSE1=mean(SE1);

mSP1=mean(SP1);

mSE2=mean(SE2);

mSP2=mean(SP2);

mSE3=mean(SE3);

mSP3=mean(SP3);

mSE4=mean(SE4);

mSP4=mean(SP4);

mSE5=mean(SE5);

mSP5=mean(SP5);

mSE6=mean(SE6);

mSP6=mean(SP6);

mSE7=mean(SE7);

mSP7=mean(SP7);

mSE8=mean(SE8);

mSP8=mean(SP8);

mSE9=mean(SE9);

mSP9=mean(SP9);

mSE10=mean(SE10);

mSP10=mean(SP10);

mSE11=mean(SE11);

mSP11=mean(SP11);

mSE12=mean(SE12);

mSP12=mean(SP12);

mSE13=mean(SE13);

mSP13=mean(SP13);

mSE14=mean(SE14);

mSP14=mean(SP14);

mSE15=mean(SE15);

mSP15=mean(SP15);

mSE16=mean(SE16);

mSP16=mean(SP16);




eb = [meanSP1 meanSP2 meanSP3 meanSP4 meanSP5 meanSP6 meanSP7 meanSP8 meanSP9 meanSP10 meanSP11 meanSP12 meanSP13 meanSP14 meanSP15 meanSP16];
Falls = [meanSE1 meanSE2 meanSE3 meanSE4 meanSE5 meanSE6 meanSE7 meanSE8 meanSE9 meanSE10 meanSE11 meanSE12 meanSE13 meanSE14 meanSE15 meanSE16];

eb1 = [mSP1 mSP2 mSP3 mSP4 mSP5 mSP6 mSP7 mSP8 mSP9 mSP10 mSP11 mSP12 mSP13 mSP14 mSP15 mSP16];
Falls1 = [mSE1 mSE2 mSE3 mSE4 mSE5 mSE6 mSE7 mSE8 mSE9 mSE10 mSE11 mSE12 mSE13 mSE14 mSE15 mSE16];


% BER = [0.0753 0.0574 0.0370 0.0222 0.0122 0.0061];

% Create a y-axis semilog plot using the semilogy function
% Plot SER data in blue and BER data in red
figure(2);
semilogx(eb, Falls, 'o-','MarkerSize',8, 'MarkerEdgeColor','red', 'MarkerFaceColor',[1 .6 .6]);
hold on
semilogx(eb1, Falls1, 'o-','MarkerSize',8, 'MarkerEdgeColor','red', 'MarkerFaceColor',[1 .6 .6]);

% Turn on the grid
grid on

% Add title and axis labels
title('Fall detection rate as a function of specificity')
ylabel('Sensitivity')
xlabel('Specificity')



figure(3);
eb=(1./(1-eb))/115;
eb1=(1./(1-eb1))/115;

plot(eb, Falls, 's-','MarkerSize',8, 'MarkerEdgeColor','red', 'MarkerFaceColor',[1 .6 .6])
hold on
plot(eb1, Falls1, 's-','MarkerSize',8, 'MarkerEdgeColor','red', 'MarkerFaceColor',[1 .6 .6])


% Turn on the grid
grid on
title('Fall detection rate as a function of false alarms per day')
ylabel('Sensitivity')
xlabel('Number of days without a false alarm')


eb = [meanSP1 meanSP2 meanSP3 meanSP4 meanSP5 meanSP6 meanSP7 meanSP8 meanSP9 meanSP10 meanSP11 meanSP12 meanSP13 meanSP14 meanSP15 meanSP16];
Falls = [meanSE1 meanSE2 meanSE3 meanSE4 meanSE5 meanSE6 meanSE7 meanSE8 meanSE9 meanSE10 meanSE11 meanSE12 meanSE13 meanSE14 meanSE15 meanSE16];

eb1 = [mSP1 mSP2 mSP3 mSP4 mSP5 mSP6 mSP7 mSP8 mSP9 mSP10 mSP11 mSP12 mSP13 mSP14 mSP15 mSP16];
Falls1 = [mSE1 mSE2 mSE3 mSE4 mSE5 mSE6 mSE7 mSE8 mSE9 mSE10 mSE11 mSE12 mSE13 mSE14 mSE15 mSE16];


figure(8);
errorbar(eb, Falls, SD_SETest, 'o-','MarkerSize',8, 'MarkerEdgeColor','blue', 'MarkerFaceColor',[  0.46  0.99   0.66]);
hold on
errorbar(eb1, Falls1, SD_SEVal,'o-','MarkerSize',8, 'MarkerEdgeColor','red', 'MarkerFaceColor',[1 .6 .6]);
title('Interpolated ROC curve with confidence intervals', 'Fontsize', 24 );
xlabel('Specificity','Fontsize', 20);
ylabel('Sensitivity','Fontsize', 20);
legend({'Real falls','Simulated falls'},'Fontsize', 16);

figure(9);
errorbar(eb, Falls, SD_SETest);
hold on
errorbar(eb1, Falls1, SD_SEVal);
title('Interpolated ROC curve with confidence intervals', 'Fontsize', 24 );
xlabel('Specificity','Fontsize', 20);
ylabel('Sensitivity','Fontsize', 20);
legend({'Real falls','Simulated falls'},'Fontsize', 16);

% close all
figure(10);
x2=[eb, fliplr(eb)];
inBetween=[Falls+(SD_SETest),fliplr(Falls-(SD_SETest))];
x22=[eb1, fliplr(eb1)];
inBetween2=[Falls1+(SD_SEVal),fliplr(Falls1-(SD_SEVal))];
patch(x22, inBetween2,[1 0.7 1]);
patch(x2, inBetween,[0.7 1 0.8]);
hold on
plot(eb, Falls, 'o-','MarkerSize',8, 'MarkerEdgeColor','blue', 'MarkerFaceColor','blue');
hold on
plot(eb1, Falls1, 'o-','MarkerSize',8, 'MarkerEdgeColor','red', 'MarkerFaceColor','red');
hold on
plot(eb1,Falls1+(SD_SEVal),'Color','black');
hold on
plot(x22,inBetween2,'Color','black');
% hold on
% plot(eb1(16),Falls1(16)-(SD_SEVal(16)):0.001:Falls1(16)+(SD_SEVal(16)),'Color','black');
title('Interpolated ROC curve with confidence intervals', 'Fontsize', 24 );
xlabel('Specificity','Fontsize', 20);
ylabel('Sensitivity','Fontsize', 20);
legend({'Simulated falls','Real falls'},'Fontsize', 16);


hyimp=[Mean_SETest;SD_SETest;Mean_SPTest;SD_SPTest]



% CI95 = bsxfun(@plus, mean(SEval,2), bsxfun(@times, [-1  1]*1.96, SEM_A));   % 95% Confidence Intervals