function [filt_cntX,ref_cICA] = cIOL(EEG_CNTX, REF_CNTX, fs, varargin)

props= {'idxRef'            1:size(REF_CNTX,2) 
        'mu0'               1          
        'lambda0'           1           
        'gamma'             0.0001   
        'learningRate'      0.1       
        'OverValue'         0.000001    
        'maxIter'           1000       
        'threshold'         1.75        
        'order'             20       
        'mu'                0.1      
        'gam'             	0.0001   
        'q'                 1.15     
        'Epsilon'           0.0001  
        'window_time'       2000  
        'moving_time'    	2000      
        'ref_order'         1     
        'flag_PCA'          false       
        };

varargin = reshape(varargin,2,[])';
isknown= ismember(props(:,1),varargin(:,1),'legacy');
new_props= [varargin; props(~isknown,1:2)];
opt = cell2struct(new_props(:,2), new_props(:,1));

%% PCA setting
flag_PCA = opt.flag_PCA;
%% c-ICA parameter setting
idxRef = opt.idxRef;
mu0 = opt.mu0;
lambda0 = opt.lambda0;
gamma = opt.gamma;
learningRate = opt.learningRate;
OverValue = opt.OverValue;  
maxIter = opt.maxIter;

threshold = opt.threshold;   % good ref

%% adaptive parameter setting
order = opt.order; % this number is same as the number of sample of time delay
mu = opt.mu;
gam = opt.gam;
q = opt.q;
Epsilon = opt.Epsilon;
ref_order = opt.ref_order;

%% data
X_cIOL = [];
ref_cICA = [];
ref_pca = [];

%% window setting
window_time = opt.window_time;
window_size = floor((window_time)/(1000/fs));
if window_size > size(EEG_CNTX,1)
    error('window_time should be shorter than data segment')
end
seg_term = 1:window_size+order;

moving_time = opt.moving_time;
moving_size = moving_time/(1000/fs);
if moving_time > window_time
    error('moving time should be shorter than or equal to window time')
end

X_size = size(EEG_CNTX,1);
repeatTime = (X_size - (window_size+order))/moving_size + 1;
if repeatTime <= 1
    repeatTime = 1;
    moving_size = window_size-order;
    seg_term = 1:window_size;
end
% rest_X = rem(X_size,moving_size);

%% W initialization
rCh = 2*length(idxRef);
pCh = size(EEG_CNTX,2);
% W_1 = squeeze(zeros(order*rCh,pCh));
    if ref_order == 1
        % 1차만
        W_n = squeeze(zeros(order*rCh,pCh));
    elseif ref_order == 2
        % 2차만
        W_n = squeeze(zeros((order^2)*rCh,pCh));
    elseif ref_order == 12
        % 1차 +2차
        W_n = squeeze(zeros((order + order^2)*rCh,pCh));
    else
        error('Check the reference order')
        return;
    end
    
if ndims(EEG_CNTX) == 3
    nTrials = size(EEG_CNTX,3);
else
    nTrials = 1;
end
%%
for j=1:nTrials
    % time delay manage
    X_rec = EEG_CNTX(1:order,:,j)';
    X_rec2 = EEG_CNTX(1:order,:,j)';
    % tic
    for i=1:repeatTime
        tic
        %% segment setting
        seg = seg_term + moving_size*(i-1);
        len = length(seg)/100;
        
        %% data segment
        X_cap = EEG_CNTX(seg,:,j)';
        REF = REF_CNTX(seg,idxRef,j)/1000; % /1000은 normalization

        %% Perform PCA
        if flag_PCA && size(REF,2)>1
            [coeff1,REF_PCA,~,~,~,mu1]= pca(REF,'algorithm','als');
            REF_X = REF_PCA(:,1:2)'; % 1:2]
        elseif size(REF,2)>1
            REF_X = REF(:,1:2)';
        else
            REF_X = REF';
        end
        
        %% %%%%%%%%%%%%%%%%%%%%      cICA        %%%%%%%%%%%%%%%%%%%%%%%%
        % initialization
        w_init = rand(size(X_cap,1),1);
        w_init=w_init/norm(w_init);

        [X_WH,W_WH] = whiten(X_cap);
        [REF_WH,V_ref] = whiten(REF_X);

        % Reference +,- 
         REF_all = [REF_WH; -REF_WH];

        % cICAmult
        [y_cICA, w_cICA2] = cICA_multiCh(X_WH, REF_all, threshold, w_init, learningRate, mu0, lambda0, gamma, maxIter, OverValue);
        y_ref{i} = y_cICA;


        %% %%%%%%%%%%%%%%%%%%%%  online learning  %%%%%%%%%%%%%%%%%%%%%%%%
        % parameter setting
        primary = X_WH;
        reference = y_ref{i};

        primary_size = size(primary,2);
        pCh = size(primary,1);
        rCh = size(reference,1);
        R_n = zeros(order*rCh);

        performance_curve = zeros(50000,1);

        % setting filters
        pri_wrt = primary(:, order+1:end);  %truncate primary // size=chan * time*1
        
        % 1st order
        if ref_order ~= 2
            ref_wrt1 = zeros(rCh,(primary_size - order),order);  % // size=time*30 (30개의 time delay signals)
            for di = 1: order
                ref_wrt1(:,:,di) =  reference(:,di:end-order+di-1);
            end
            ref_wrt1 = permute(ref_wrt1,[1,3,2]);
            ref_wrt1 = reshape(ref_wrt1,[rCh*order,(primary_size - order)]);
        end
        
        % 2nd order
        if ref_order ~= 1
            ref_wrt2 = zeros(rCh,(primary_size - order),order^2); % 2nd order reference
            dij = 1;
            for di = 1: order
                for dj = 1: order
                    ref_wrt2(:,:,dij) = ...
                        reference(:,di:end-order+di-1).* reference(:,dj:end-order+dj-1);
                    dij = dij+1;
                end
            end
            ref_wrt2 = permute(ref_wrt2,[1,3,2]);
            ref_wrt2 = reshape(ref_wrt2,[rCh*order^2,(primary_size - order)]);
        end
        
        if ref_order == 12
            % 1차 + 2차
            ref_wrt_U = [ref_wrt1; ref_wrt2];
        elseif ref_order == 2
            % 2차만
            ref_wrt_U = [ref_wrt2];
        elseif ref_order == 1
            % 1차만
            ref_wrt_U = [ref_wrt1];
        else
            error('Check the reference order')
            return;
        end


        %% excute AF
         [pri_out_nl,W_n] = computeAF(pri_wrt,ref_wrt_U, W_n, mu,Epsilon);

        %% result save
        Out = pri_out_nl; % 한 sample별 ref:1+2차
        Out2 = pri_wrt - (W_n' * squeeze(ref_wrt_U) ); % 한 seg 다 돈 W로

        %% inv(W_WH) * X_WH
        Y = inv(W_WH) *Out;
        X_rec = [X_rec Y(:,1:moving_size)];
        Y2 = inv(W_WH) *Out2;
        X_rec2 = [X_rec2 Y2(:,1:moving_size)];
        
        ref_pca = [ref_pca REF_X(:,1:moving_size)];
        ref_cICA = [ref_cICA y_ref{i}(:,1:moving_size)];

    end
    X_cIOL(:,:,j) = [X_rec'; EEG_CNTX(size(X_rec',1)+1:end,:,j)];
    X_cIOL2(:,:,j) = [X_rec2'; EEG_CNTX(size(X_rec2',1)+1:end,:,j) ];

end
%%
filt_cntX = X_cIOL;
filt_cntX2 = X_cIOL2;

%% Adaptive filter
function [output, w, performance_curve]=computeAF(X_s,X_r,w,mu,Epsilon)
performance_curve = 0;
output = zeros(size(X_s));
%% Full Volterra
[fea_in, L]=size(X_r);

coeflen=fea_in;
delta=0.0007;
lamda=1;
P=(1/delta)*eye(coeflen);

for n=1:L
    s_n =  X_s(:,n);
    U_n = X_r(:,n);
    
    error = s_n - w'*U_n;
    errorSquare = (error).^2;
    
    g = mu * P*U_n/(lamda+U_n'*P*U_n); % 원래 Epsilon lambda와 공유했었음.
   
    w = w + g*error';
    output(:,n) = s_n - w'*U_n;
    
    P = (1/lamda)*(P - g*U_n'*P);
    
%     figure(11)
%     plot(output')
%     ylim([-50 50])
%     MSE = sum(errorSquare)/(n);
%     performance_curve(n,1) = MSE;
end

