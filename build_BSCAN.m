%% BUILD THE BSCAN
clearvars
% Ask for files
[file0,path0] = uigetfile('.mat','SELECT THE U0 FILE') ;
if path0==0 ; return ; end
[files,path] = uigetfile('.mat','SELECT THE BSCAN FILES',path0,'Multiselect','on') ;
if path==0 ; return ; end
files = sort(files) ;
% Import the data
DATA0 = load([path0 filesep file0]) ;
DATA = cellfun(@(f)load([path filesep f]), files) ;
% Infos
t = DATA0.t(:)*1e-6 ; 
nT = numel(t) ;
theta = cellfun(@str2double,regexprep(files,'.*_(\d+)(\.\d+)?.mat','$1$2'))*pi/180 ;
nTheta = numel(theta) ;
% Re-sort if needed
[theta,is] = sort(theta,'ascend') ;
files = files(is) ;
DATA = DATA(is) ;
% Data
s_ref = mean(DATA0.RF,2) ;
s = cat(3,DATA.RF) ;
s = reshape(mean(s,2),[nT nTheta]) ;
% Correct the gain
G = [DATA.analogGain DATA0.analogGain] ;
G = 10.^(G/20) ;
s = (G(end)./G(1:nTheta)).*s ;
% Display
clf ; axis tight
imagesc(theta*180/pi,1e6*t,s)
shading interp
clrmp = interp1([0;.5;1],[0 0 1;1 1 1;1 0 0],linspace(0,1,1000)') ;
colormap(clrmp) ;
caxis([-1 1]./500)
set(gca,'ylim',[390 405])

%% SAVE THE BSCAN
Questions={'h (mm)','rho (kg/m3)'}; %Ask measured specimen thickness and density
defaultans={'3.994','1185'} ; %Default parameters for VeroWhite
opts.Interpreter = 'tex';
Parameters=inputdlg(Questions,'Input parameters',1,defaultans,opts);

h = str2num(Parameters{1})*1e-3; % measured specimen thickness
rho = str2num(Parameters{2}); % measured specimen density

Temp = mean([DATA.T]) ; % mean T° over the experiment
tlim = 1e-6*get(gca,'ylim') ; % let the user choose a time window
indT = t>=min(tlim) & t<=max(tlim) ;

[sfile,spath] = uiputfile('.mat','SAVE THE BSCAN',[path0(1:end-1) '.mat']) ;
if spath==0 ; return ; end

t = t(indT) ;
s = s(indT,:) ;
s_ref = s_ref(indT) ;
save([spath filesep sfile],'h','t','s_ref','theta','s','rho','Temp') ;




