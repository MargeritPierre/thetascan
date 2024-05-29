% iTSCAN 
% Ultrasound properties estimation via inversion of theta-scans performed
% in double-transmission
%
% The input data (a theta-scan) must be a structure with the following
% fields (see the examples in ./data):
% - "s_ref" [nT,1] reference signal (without the specimen)
% - "s" [nT,nTheta] theta-scan
% - "t" [nT,1] time vector (sec)
% - "theta" [1,nTheta] angle vector (rad)
% - "rho" specimen density (tons/mm3)
% - "h" specimen thickness (mm)
% - "Temp" fluid temperature (°C)
% some fields can, if needed, given or overriden below
%
% Authors: 
%   Pierre Margerit: pierre.margerit@ensam.eu 
%   Anne-Sophie Poudrel: anne-sophie.poudrel@ensam.eu
% Date: 05-feb-2024
%
% Copyright: This is published under the GNU GPL v3 license ...

%% LOAD A BSCAN FILE
% Let the user choose a theta-scan
    clearvars,close all
    [file,path] = uigetfile('.mat','SELECT A TSCAN FILE') ;
    if path==0 ; return ; end
% Load the file
    load([path filesep file]) ; 
    [nT,nTheta] = size(s) ;
    if ~exist('theta','var') ; theta = (0:nTheta-1)*pi/180 ; end
% Plot the file
    clf ; plotMAP(s,theta,t) ;
    
%% INITIAL GUESS WITH PRONY METHOD

% PARAMETERS ==============================================================
% Specimen
    nu_init = 38/100 ; % initial Poisson ratio guess
%     rho = 1115 ; % (tons/mm3) specimen density
%     h = 39e-3/10 ; % (mm) specimen thickness
% Fluid
%     rho_0 = 1000 ; % (tons/mm3) fluid density
%     c0 = 1500 ; % (mm/s) wave speed in fluid
%     Temp = 20 ; % temperature in °C
    [c0,rho_0] = marzac(Temp) ;
% Signal
    omegaC = 2*pi*2.25e6 ; % central frequency of the probe
    Tmodel = 2*(max(t)-min(t)) ; % (s) signal duration for the model (avoid folding)
    zeroMean = true ; % force the input to be zero mean ?
% =========================================================================

% Time/frequency samples 
    dt = mean(diff(t)) ;
    nFmodel = ceil(Tmodel./dt/2) ;
    nTmodel = 2*nFmodel-1 ;
    fModel = (0:nFmodel-1)'/nTmodel/dt ; % /!\ no zero frequency !
    omega = 2*pi*fModel ;

% Reference signal
    s_ref = s_ref-zeroMean*mean(s_ref,1) ;

% Initialize wave speeds with Prony method at normal incidence
    cl = normalIncidenceEstim(s(:,theta==0),s_ref,t,c0,h) ;
    disp("Estimated longitudinal wave celerity: cl="+num2str(cl))
    ct = cl.*sqrt((1-2*nu_init)./(1-nu_init)./2) ;

% Deduce wavenumbers
    k0 = omega./c0 ;
    kl = omega./cl ;
    kt = omega./ct ;
    
% Compute the "initial" BSCAN
    u = TSCAN(theta,k0,rho_0,kl,kt,rho,h,s_ref) ;
    clf ; im = plotMAP(cat(3,s,u,u-s),theta,t) ;
    set(title([im.Parent],''),{'String'},{'\bf Experiment';'\bf Model';'\bf Residual'})


%% SZABO'S MODEL PARAMETERS ESTIMATION VIA RESIDUAL MINIMIZATION
% Optimization parameters
    dPmax = 1e-6 ; % maximum relative update of parameters at convergence
    maxIt = 100 ; % maximum number of update iterations
    step = .25 ; % relative updating step (damping)
    relSensitivity = true ; % use relative sensitivities for better conditionning
    variableGain = false ; % optimize the gain ?
    
% SZABO Model
    Np_dB = 20/log(10)/100 ; % unit scaling for attenuation
    v_szabo = @(v,a,y) (1./v + a/Np_dB .* tan(pi/2.*y)* ( ((omega./omegaC).^y).*(omega.^-1)-omegaC.^-1) ).^-1 ; % phase speed
    a_szabo = @(v,a,y) a.*((omega./omegaC).^y) ; % wave attenuation
    k_szabo = @(v,a,y) (omega./v_szabo(v,a,y)) - 1i*a_szabo(v,a,y)/Np_dB ; % wavenumber law

% Initial values for the parameters
args0 = struct(...
        'vlC' , real(cl)*(1 + (imag(cl)/real(cl))^2) ...
        ,'alC' , omegaC*imag(cl)/abs(cl)^2*Np_dB ...
        ,'yl' , 1.1 ...
        ,'vtC' , real(ct)*(1 + (imag(ct)/real(ct))^2) ...
        ,'atC' , omegaC*imag(ct)/abs(ct)^2*Np_dB ...
        ,'yt' , 0.7 ...
        ,'h', h ...
        ,'rho', rho ...        
        ,'theta0', 0.01 ...
            ) ;
% Parameters to optimize
optim = { ...
        'vlC' ; ...
        'alC' ; ...
        'yl' ; ...
        'vtC' ; ...
        'atC' ; ...
        'yt' ; ...
        'h' ; ...  
        ...'rho' ; ...        
        'theta0' ; ...
    } ;  

% Generator functions
bscanFun = @(arg)TSCAN(theta + arg.theta0...
                ,k0,rho_0...
                ,k_szabo(arg.vlC,arg.alC,arg.yl) ...
                ,k_szabo(arg.vtC,arg.atC,arg.yt) ...
                ,arg.rho,arg.h,s_ref) ;
sFun = @(p)bscanFun(setfields(args0,optim,p)) ;

% INITIAL BSCAN
args = args0 ; p = getfields(args0,optim) ; dP = inf ; it = 0 ;
u = sFun(p) ;
clf ; im = plotMAP(cat(3,s,u,u-s),theta,t) ;
set(title([im.Parent],''),{'String'},{'\bf Experiment';'\bf Model';'\bf Residual'})

% OPTIMIZATION
while max(abs(dP(:)./p(:)))>dPmax && it<maxIt
% Model Sensitivities
    ds_dp = dfun_dp(sFun,p) ;
    if relSensitivity ; ds_dp = ds_dp.*reshape(p,1,1,[]) ; end
% Evaluate gain from expe/model signal norms
    G_est = norm(s(:))/norm(u(:)) ;
    if variableGain ; G_it = G_est ; else G_it = 1 ; end
% Residual
    RES = u-G_it*s ;
    dRES_dp = reshape(ds_dp,[],numel(p)) ;
% Update
    jac = dRES_dp'*RES(:) ; % jacobian
    Hess = (dRES_dp'*dRES_dp) ; % Hessian
    dP = -Hess\jac ;
    if relSensitivity ; dP(:) = dP(:).*p(:) ; end
% UPDATED BSCAN
    it = it+1 ;
    p(:) = p(:) + step*dP(:) ;
    u = sFun(p) ;
    % DISPLAY
    args = setfields(args0,optim,p) ;
    disp("it:"+string(it)+" | G_it="+num2str(G_it,4)) ; disp(args) 
    plotMAP(cat(3,s,u,RES),theta,t,im) ; drawnow ;
end
    
%% FINAL VELOCITIES & UNCERTAINTIES
% Reconstruct laws with optimized parameters
    vFun = @(args,r) v_szabo(args.("v"+r+"C"),args.("a"+r+"C"),args.("y"+r)) ;
    aFun = @(args,r) a_szabo(args.("v"+r+"C"),args.("a"+r+"C"),args.("y"+r)) ;
    vl_opt =  vFun(args,'l') ; al_opt = aFun(args,'l') ;
    vt_opt =  vFun(args,'t') ; at_opt = aFun(args,'t') ;
% Evaluate uncertainties
    sigmaFun = @(dg_df,deltaf) sum((abs(dg_df.').*abs(deltaf(:))).^2,1).^(1/2) ;
    % Model
        u = sFun(p) ;
        ds_dp = dfun_dp(sFun,p) ;
    % Residual
        RES = u(:)-G_it*s(:) ;
        dRES_dp = reshape(ds_dp,[],numel(p)) ;
    % Propagate into parameters
        dp_dRES =  (dRES_dp'*dRES_dp)\dRES_dp' ; % pseudo-inverse
        sigmaP = sigmaFun(dp_dRES,abs(RES)) ; 
    % Propagate into material laws
        dvl_dp = dfun_dp(@(p)vFun(setfields(args,optim,p),'l'),p) ;
        dvt_dp = dfun_dp(@(p)vFun(setfields(args,optim,p),'t'),p) ;
        deltaVl = sigmaFun(dvl_dp(:,:),sigmaP).' ;
        deltaVt = sigmaFun(dvt_dp(:,:),sigmaP).' ;
        dal_dp = dfun_dp(@(p)aFun(setfields(args,optim,p),'l'),p) ;
        dat_dp = dfun_dp(@(p)aFun(setfields(args,optim,p),'t'),p) ;
        deltaAl = sigmaFun(dal_dp(:,:),sigmaP).' ;
        deltaAt = sigmaFun(dat_dp(:,:),sigmaP).' ;
% Display
    clf ; ax = gobjects(0) ;
    ax(end+1) = fitsubplot(1,2,1) ; hold on
        hhh = plot(omega/2e6/pi,vl_opt,'DisplayName','$v_\ell$') ;
        plot(omega/2e6/pi,vl_opt + 1.95*deltaVl,':','color',hhh.Color,'DisplayName','$v_\ell+\Delta v_\ell$')
        plot(omega/2e6/pi,vl_opt - 1.95*deltaVl,':','color',hhh.Color,'DisplayName','$v_\ell-\Delta v_\ell$')
        hhh = plot(omega/2e6/pi,vt_opt,'DisplayName','$v_t$') ;
        plot(omega/2e6/pi,vt_opt + 1.95*deltaVt,':','color',hhh.Color,'DisplayName','$v_t+\Delta v_t$')
        plot(omega/2e6/pi,vt_opt - 1.95*deltaVt,':','color',hhh.Color,'DisplayName','$v_t-\Delta v_t$')
        ax(end).YLim(1) = 0 ;
        ylabel(ax(end),'Wave Speed (m/s)') ;
        legend
    ax(end+1) = fitsubplot(1,2,2) ; hold on
        hhh = plot(omega/2e6/pi,al_opt,'DisplayName','$\alpha_\ell$') ;
        plot(omega/2e6/pi,al_opt + 1.95*deltaAl,':','color',hhh.Color,'DisplayName','$\alpha_\ell+\Delta \alpha_\ell$')
        plot(omega/2e6/pi,al_opt - 1.95*deltaAl,':','color',hhh.Color,'DisplayName','$\alpha_\ell-\Delta \alpha_\ell$')
        hhh = plot(omega/2e6/pi,at_opt,'DisplayName','$\alpha_t$') ;
        plot(omega/2e6/pi,at_opt + 1.95*deltaAt,':','color',hhh.Color,'DisplayName','$\alpha_t+\Delta \alpha_t$')
        plot(omega/2e6/pi,at_opt - 1.95*deltaAt,':','color',hhh.Color,'DisplayName','$\alpha_t-\Delta \alpha_t$')
        ylabel(ax(end),'Wave Attenuation (dB/cm)') ;
        legend
    set(ax,'Xlim',omegaC/2e6/pi*[.25 3])
    xlabel(ax,'Frequency (MHz)') ;

%% SUB-FUNCTIONS
function [u,H] = TSCAN(theta,k0,rho_0,kl,kt,rho,h,s_ref)
% Model of the theta-scan
    % Define wavevectors & modes
    % fluid waves
        chi0 = k0.*cos(theta) ;
        kappa = k0.*sin(theta) ;
    % solid waves: refraction occurs; across the frontier, k_2 is a constant==kappa
        chil = sqrt(kl.^2 - kappa.^2) ; % longitudinal
        chit = sqrt(kt.^2 - kappa.^2) ; % transverse
    % Cumpute the transfer for a simple transmission
        % Propagation operators
            Pi0 = exp(-1i*chi0*h) ;
            PiL = exp(-1i*chil*h) ;
            PiT = exp(-1i*chit*h) ;
        % Coefficients
            alpha = (chit.^2-kappa.^2).^2 ;
            beta = 4.*kappa.^2.*chil.*chit ;
            gamma = ((rho_0.*chil)./(rho.*chi0)).*((chit.^2+kappa.^2).^2) ;
        % Transfer
            HT = (-4.*gamma./Pi0).*(alpha.*PiL.*(PiT.^2-1) + beta.*PiT.*(PiL.^2-1))./(...
                        (PiL.^2-1).*(PiT.^2-1).*(alpha+beta-gamma).^2 ...
                        - 4.*beta.*gamma.*(PiL.^2-1) ...
                        - 4.*alpha.*gamma.*(PiT.^2-1) ...
                        + 4.*alpha.*beta.*(PiL-PiT).^2 ...
                    ) ;
            H = HT.^2 ; % double-transmission
            H(isnan(H)) = 1 ; % avoid problems at kappa==0
    % OUTPUT TIME SIGNAL
        nTmodel = 2*size(H,1)-1 ;
        % Input spectrum
        Sref = fft(s_ref,nTmodel) ; 
        Sref = Sref(1:size(H,1)) ; % half the spectrum
        U = H.*Sref ;
        u = ifft(U,2*size(U,1)-1,1,'symmetric') ;
        u = u(1:size(s_ref,1),:,:,:) ; % crop with the length of the input
end

function cl = normalIncidenceEstim(s,s_ref,t,c0,h)
% Estimate a complex longitudinal wave speed from normal incidence data using
% a Prony method
    % Infos
        dt = mean(diff(t)) ;
        domega = 2*pi/numel(t)/dt ;
        omega = (0:numel(t)-1)'*domega ; 
    % Signals & Transfer function
        S = fft(s,[],1) ;
        Sref = fft(s_ref,[],1) ;
        iH = Sref./S ; 
    % Water layer correction
        iH0 = iH.*exp(2i*(h./c0).*omega) ;
    % Weighting with signal bandwidth
        w = ones(size(S)) ; % initialize
        w = w.*(abs(S)>0.01*max(abs(S))) ; % thresholding
        w = w.*(abs(S)./max(abs(S))).^(2/2) ; % signal power related weight
        w = w.*(dt*omega<pi) ; % only half-spectrum is needed
    % Hankel data matrix
        hank = hankel(1:4,4:numel(iH0)).' ;
        hn = iH0(hank) ;
    % Diagonal weight matrix
        W = min(w(hank),[],2) ;
    % Prony method
        P12 = hn*[0;-1;1;0] ; % h2-h1
        P03 = hn*[-1;0;0;1] ; % h3-h0
    % Total Least-Squares
        [vec0,~] = eigs([P12 P03]'*diag(W)*[P12 P03],1,0) ;
        p2 = -vec0(1)/vec0(2) ;
    % Longitudinal complex wave speed
        cl = 2*h*domega/acos((p2-1)/2) ;
end

function df_dp = dfun_dp(fun,p,delta)
% Sensitivity estimation using finite differences
    if nargin<3 ; delta = 1e-6 ; end
    f = fun(p) ;
    df_dp = NaN([numel(f) numel(p)]) ;
    for pp = 1:numel(p)
        dp = delta*p(pp) ;
        p(pp) = p(pp) + dp ;
        df_dp(:,pp) = (reshape(fun(p),[],1)-f(:))./abs(dp) ;
        p(pp) = p(pp) - dp ;
    end
    df_dp = reshape(df_dp,[size(f) numel(p)]) ;
end

function [cw,rhow] = marzac(Temp)
% Estimation of wave celerity in water 
    % cw = 1.40285*1e3 + 5.038813*1e0*T - 5.799136*1e-2*T^2+ 3.287156*1e-4*T^3 -1.398845*1e-6*T^4;
    cw = 1.402385*1e3 + 5.038813*1e0*Temp - 5.799136*1e-2*Temp^2+ 3.287156*1e-4*Temp^3 -1.398845*1e-6*Temp^4 + 2.787860*1e-9*Temp^5;
    rho0 = 998.2071;
    T0 = 20;
    beta = 0.0002;
    rhow = rho0/(1+beta*(Temp-T0));
end

function im = plotMAP(MAP,theta,t,im,sat)
% Plot one or multiple theta-scans or equivalents
% MAP [nT nTheta nMaps] data to plot
% theta [nTheta] angle vector
% t [nT] time vector
% im (optional) vector of image graphical objects to update
% sat (optional) color saturation ratio
    if nargin<5 ; sat = 1/10 ; end % color saturation
    nAX = size(MAP,3) ;
% Image updating only
    if nargin>=4 && ~isempty(im)
        for aa = 1:nAX
            im(aa).CData = MAP(:,:,aa) ;
        end
        return ;
    end
% Build the figure
    ax = gobjects(0) ; % subplots
    im = gobjects(0) ; % image graphical handles
    for aa = 1:nAX
        ax(end+1) = fitsubplot(1,nAX,aa); 
        im(end+1) = imagesc(theta*180/pi,t*1e6,MAP(:,:,aa)) ;
    end
    axis(ax,'tight')
% Legends
    xlabel(ax,'Angle (deg)')
    ylabel(ax,'Time (us)')
% Red/blue colormap
    clrmp = interp1([0;.5;1],[0 0 1;1 1 1;1 0 0],linspace(0,1,1000)') ;
    colormap(gcf,clrmp) ;
    set(ax,'clim',max(abs(MAP(:))).*[-1 1]*sat)
% Link axes positions for easier navigation
    linkaxes(ax,'xy')
end

function values = getfields(args,fields)
% Get some fields of a structure (vectorized version of MATLAB's bultin)
    values = [] ;
    for ff = 1:numel(fields)
        values(end+1) = args.(fields{ff}) ;
    end
end

function args = setfields(args,fields,values)
% Set some fields of a structure (vectorized version of MATLAB's bultin)
    for ff = 1:numel(fields)
        args.(fields{ff}) = values(ff) ;
    end
end

function ax = fitsubplot(nR,nC,ii)
% Subplots fitted with the figure to save place !
    [cc,rr] = ind2sub([nC nR],ii) ;
    ax = axes('outerposition',[(cc-1)/nC (rr-1)/nR 1/nC 1/nR]) ;
end




