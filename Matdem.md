# zlw
%step1: packing the elements
clear;
fs.randSeed(1);%random model seed, 1,2,3...
B=obj_Box;%declare a box object
B.name='PoreFluidFlow';
%--------------initial model------------
B.GPUstatus='auto';%program will test the CPU and GPU speed, and choose the quicker one
B.ballR=0.1;
B.isShear=0;
B.isClump=0;%if isClump=1, particles are composed of several balls
B.distriRate=0.2;%define distribution of ball radius, 
B.sampleW=20;%width, length, height, average radius
B.sampleL=0;%when L is zero, it is a 2-dimensional model
B.sampleH=15;
B.boundaryRrate=0.999999;
B.BexpandRate=2;%boundary is 4-ball wider than 
B.PexpandRate=1;
B.isSample=1;
B.type='TriaxialCompression';
B.setType();
B.buildInitialModel();%B.show();
d=B.d;%d.breakGroup('sample');d.breakGroup('lefPlaten');
%you may change the size distribution of elements here, e.g. d.mo.aR=d.aR*0.95;
d.showB=1;
%--------------end initial model------------

d.mo.setGPU('off');
%d.aR(d.GROUP.sample)=d.aR(d.GROUP.sample)*0.88;
d.aZ(1)=d.aZ(1)*1.1;
delId=[d.GROUP.topPlaten(end-1:end);d.GROUP.botPlaten(end-1:end)];
d.delElement(delId);
d.mo.zeroBalance();
%d.show();save(['TempModel/' B.name '1.mat'],'B','d');return;

d.mo.isShear=0;

%---------- gravity sedimentation
B.gravitySediment(0.5)
B.compactSample(0.5);%input is compaction time
%------------return and save result--------------
d.status.dispEnergy();%display the energy of the model
d.show('-aR');
d.mo.bFilter(:)=1;
d.mo.zeroBalance();
d.Rrate=1;
d.mo.setGPU('off');
d.clearData(1);%clear dependent data
d.recordCalHour('Step1Finish');
save(['TempModel/' B.name '1.mat'],'B','d');
save(['TempModel/' B.name '1R' num2str(B.ballR) '-distri' num2str(B.distriRate)  'aNum' num2str(d.aNum) '.mat']);
d.calculateData();
%d.show('-Id');


%set the material of the model
clear
load('TempModel/PoreFluidFlow1.mat');
d.show('aR');
return;
B.setUIoutput();%set output of message
d=B.d;
d.calculateData();
d.mo.setGPU('off');
d.getModel();%get xyz from d.mo

%----------set material of model
matTxt=load('Mats\RockHydro.txt');
Mats{1,1}=material('RockHydro',matTxt,B.ballR);
Mats{1,1}.Id=1;
d.Mats=Mats;
%----------end set material of model

%---------assign material to layers and balance the model
B.setPlatenFixId();
d.setGroupMat('sample','RockHydro');
d.groupMat2Model({'sample'});
d.balanceBondedModel0();
d.mo.bFilter(:)=false;
d.balance('Standard',0.1);
%---------end assign material to layers and balance the model1.	

%---------save the data
d.mo.setGPU('off');
d.clearData(1);
d.recordCalHour('Step2Finish');
save(['TempModel/' B.name '2.mat'],'B','d');
save(['TempModel/' B.name '2R' num2str(B.ballR) '-distri' num2str(B.distriRate)  'aNum' num2str(d.aNum) '.mat']);
d.calculateData();


%the code is developed for the project with China Ocean Univ.
clear
fs.randSeed(2);
load('TempModel/PoreFluidFlow1.mat');
%-----------initializing the model

sampleH=12;%height of sample is 4 cm
dH=12;%water head
%---------calculate connection diameter and flow K
kRate=0.5;%change the kRate to change permeability
k=1e-9;%permeability factor
%---------end calcualte connection diameter and flow K

B.setUIoutput();
d=B.d;
d.calculateData();
d.mo.setGPU('off');
d.getModel();%get xyz from d.mo
d.showB=2;
d.deleteConnection('boundary');
d.Rrate=1;
d.resetStatus();
d.getModel();
d.mo.isCrack=1;
%-----------end initializing the model

%------------remove top elements of sample to make ocean area
sampleHcenter=mean(d.mo.aZ(d.GROUP.sample));
topModelFilter=d.mo.aZ>sampleHcenter+sampleH/2;
botModelFilter=d.mo.aZ<sampleHcenter-sampleH/2;
delId=find(topModelFilter|botModelFilter);
maxS=max(d.GROUP.sample);
delId=delId(delId<=maxS);
d.delElement(delId);
fixId=[d.GROUP.sample];
d.addFixId('X',fixId);
d.addFixId('Y',fixId);
d.addFixId('Z',fixId);
d.show('aR');
% return
%------------end remove top elements of sample to make ocean area


%------------initilizing
p=pore(d);
p.pathLimitRate=0.3;%path diameter<pathLimitRate*ballR will be connection
p.isCouple=1;%fluid-solid coupling
p.setInitialPores();
p.setPlaten('fix');%fix the coordinates of platens
p.aWaterdR=d.mo.aR*0.025;
p.setWaterdR();
%set the cracks
cdX1=d.mo.aX(p.cList(:,1)); cdX2=d.mo.aX(p.cList(:,2)); 
cdZ1=d.mo.aZ(p.cList(:,1)); cdZ2=d.mo.aZ(p.cList(:,2));
Zmax=max(d.mo.aZ((d.GROUP.sample)));
Zmin=min(d.mo.aZ((d.GROUP.sample)));
CrackFilter=find(cdZ1>Zmin&cdZ1<Zmax&cdZ2>Zmin&cdZ2<Zmax);
m = 1; %randNum
n = 5000; %crackNum
list=1:sum(size(CrackFilter,1));
% cFilter0=ones(n,m);
% rng(6);
rand('seed',1);
cFilter=CrackFilter(randperm(numel(list),n));
% sLen=length(d.GROUP.sample);
% lowKRate=0.1;
% lowKId=randperm(sLen,round(sLen*lowKRate));
% p.aWaterdR(lowKId)=p.aWaterdR(lowKId)*0.1;
% p.setWaterdR();
% d.mo.SET.aWaterdR=p.aWaterdR;
% d.show('SETaWaterdR');
% return
%-----------------end set cracks in the model
%--------------setting of the simulation
p.dT=p.d.mo.dT/kRate;
% fName=['data/step/' B.name  num2str(B.ballR) '-lowKRate' num2str(lowKRate) 'loopNum'];
topBallId=ceil(mean(d.GROUP.topPlaten));
botBallId=ceil(mean(d.GROUP.botPlaten));
% top
% pressureHigh=p.pPressure(2)+1e3*9.8*dH;
pressureHigh=p.pPressure(2)+1.2e5;
pressureLow=p.pPressure(2);
topPoreId=p.getBallConnectedPore(topBallId);
botPoreId=p.getBallConnectedPore(botBallId);
p.pPressure(topPoreId)=pressureHigh;
p.pPressure(botPoreId)=pressureLow;
p.setPressure();%set pore pressure of seawater
p.isCouple=0;%no fluid flow coupling
%--------------end setting of the simulation
% return
% save([fName '0.mat']);%return;
tic
%-----------apply high pressure

cDiameterFlow=p.cDiameter+p.cDiameterAdd;%calculate the diameter of
cDiameterFlow(cDiameterFlow<0)=0;
p.cKFlow=cDiameterFlow*k./p.cPathLength;%default K of throat is determined by diameter and path length
p.cKFlow(cFilter)=p.cKFlow(cFilter)*20;%crack K is greater
pIndex=p.getPIndex();
topPIndex=pIndex(topPoreId,:);
botPIndex=pIndex(botPoreId,:);
totalBalance=5000000;
step=1000;
% balanceRates=[];
balanceRates=zeros(totalBalance/step,1);
pMass=zeros(p.pNum,totalBalance/step);
%-------------balancing the pore pressure
for i=1:totalBalance
    p.pPressure(topPoreId)=pressureHigh;
    p.pPressure(botPoreId)=pressureLow;
    p.setPressure();%set pore pressure of seawater
    p.balance();%rate defines the balance time
    toppMass=p.poreFlowMass{topPIndex(1)}(topPIndex(2),:);
    botpMass=p.poreFlowMass{botPIndex(1)}(botPIndex(2),:);
    if mod(i,step)==0
        balanceRate=-sum(botpMass,2)/sum(toppMass,2);
        balanceRates(i/step,1)=balanceRate;
        fs.disp(['Balance rate is ' num2str(balanceRate*100) '%']);
        pMass(:,i/step)=p.pMass;
    end
end
%-------------end balancing the pore pressure
% return

stableT=p.totalT;
massI=0;
totalCircle=400;
stepNum=50;
topPoreMass=zeros(totalCircle*stepNum,1);
botPoreMass=zeros(totalCircle*stepNum,1);


for i=1:totalCircle
    for j=1:stepNum
        p.pPressure(topPoreId)=pressureHigh;
        p.pPressure(botPoreId)=pressureLow;
        p.setPressure();%set pore pressure of seawater
        p.balance();%rate defines the balance time
        massI=massI+1;
        toppMass=p.poreFlowMass{topPIndex(1)}(topPIndex(2),:);
        botpMass=p.poreFlowMass{botPIndex(1)}(botPIndex(2),:);
        topPoreMass(massI)=sum(toppMass,2);
        botPoreMass(massI)=sum(botpMass,2);
    end
    %save([fName num2str(i) '.mat']);
end
topMassAll=sum(topPoreMass);
botMassAll=sum(botPoreMass);
flowT=p.totalT-stableT;

%k = Q*L /( A*△h), calculate the hydraulic conductivity
Q=botMassAll/1e3/flowT;
L=sampleH;
A=B.sampleW*p.pThickness;
K=Q*L/(A*dH);
fs.disp(['Permeability coefficient is ' num2str(K)]);
% save([fName '0.mat']);%return;
toc
p.show('pPressure');
