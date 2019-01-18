clear;clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Higham, J.E., Brevis, W & Keylock, C.J. (2018) 
% Implications of the selection of a particular modal decomposition technique for the analysis of shallow flows 
% Journal of Hydraulic Research
% Taylor and Francis
% 10.1080/00221686.2017.1419990 
% Any questions email jonny.e.higham@gmail.com

%%%%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%%%%

% When the data is input the POD modes will be calculated after which you
% are promped with a window to choose the frequency, highlighted in the
% particular POD mode, which you would like the DMD algorithm to extract. 

%%%%%%%%%%%%%%%%%% EXAMPLES %%%%%%%%%%%%%%%%%%%%

% Here we provide two examples 'cylinder' & 'groyne' change by either using
% 'cylinder' or 'groyne';
% Chose the cylinder example for one velocity component and the groyne for
% two velocity components (this can easily be extended to three components

%%%%%%%%%%%%%%%%%%  INPUT  %%%%%%%%%%%%%%%%%%%%%

freq = 5; %Hz; aquisition frequency
example = 'cylinder'; % input case (cylinder / groyne) 
mode=2; % which mode to plot and extract DMD from

%%%%%%%%%%%%%%%%%%% SCRIPT %%%%%%%%%%%%%%%%%%%%%

% Import the images from folder. The colour images are converted to
% grayscale. 
data = [];
if strcmp(example,'cylinder')
    directory = 'data_cylinder';
    count = 1;
    while exist(sprintf('%s/B_%0.4d.bmp',directory,count)) ~= 0
        data(:,:,count) = imread(sprintf('%s/B_%0.4d.bmp',directory,count));
        count = count +1;
    end
    r = size(data,1); c = size(data,2);
    % Transform the data into column vectors.
    data=reshape(data,r*c,size(data,3));
end
if strcmp(example,'groynes')
    directory = 'data_groynes';
    count = 1;
    while exist(sprintf('%s/PIV_%0.4d.mat',directory,count)) ~= 0
        load(sprintf('%s/PIV_%0.4d.mat',directory,count))
        data_u(:,:,count) = u;
        data_v(:,:,count) = v;
        count = count +1;
    end
    r = size(data_u,1); c = size(data_u,2);
    % Transform the data into column vectors.
    data=[reshape(data_u,r*c,size(data_u,3));reshape(data_v,r*c,size(data_v,3))];
end

%%%%%%%%%%%%%%%%%%% POD %%%%%%%%%%%%%%%%%%%%%%

% Perform the POD - this is mean subtracted for POD not for DMD
[Phi ,~, C]=svd(data-repmat(mean(data,2),[1 size(data,2)]),'econ');
% Plot the figures
if strcmp(example,'cylinder');figs = 3;com=1;d=1;end
if strcmp(example,'groynes');figs = 4;com=2;d=0;end
close all
figure('name','POD')
subplot(figs,1,1)
imagesc(reshape(Phi(1:r*c,mode),r,c));axis image;set(gca,'Ydir','Normal')
title(['\Phi',sprintf('_%i (u-component)',mode)])
subplot(figs,1,2-d)
if com == 2 
    imagesc(reshape(Phi(r*c+1:end,mode),r,c));axis image;set(gca,'Ydir','Normal')
    title(['\Phi',sprintf('_%i (v-component)',mode)]);
end
subplot(figs,1,3-d)
plot(C(:,mode))
title(sprintf('C_%i',mode))
subplot(figs,1,4-d)
[px, fx]=pwelch(C(:,mode),round(0.9*size(data,2)),round(0.8*size(data,2)),...
    2^12,freq);
plot(fx,px);
title(sprintf('P(C_%i)',mode))
set(figure(1),'position',[246    66   370   732])
waitfor(msgbox('Please select the peak of the POD coefficient spectra (bottom image) you would like to extract'))
[x, y]=ginput(1);
[locs, vals]=findpeaks(px); 
[v, l]=min(abs(fx(vals)-x)); 
hold on
plot(fx(vals(l)),locs(l),'ro')

%%%%%%%%%%%%%%%%%%% DMD %%%%%%%%%%%%%%%%%%%%%%

% Here we run the POD again, however this is not mean subtracted and one
% piece of data is removed from the end. 
[Phi,S,C]=svd(data(:,1:end-1),'econ');
F=(Phi'*data(:,2:end)*C)/S;
[M,Z]=eig(F) ;% Compute Eigenvalues and Eigenvectors
Q=zeros(size(F));
for loop1=1:size(data,2)-1
    Q(:,loop1)=diag(Z).^(loop1-1); %creating vandermonde with increasing powers
end
% Compute the DMD modes
Psi=data(:,1:end-1)*Q';
% Compute the frequencies
f=freq*angle(diag(Z))./(2*pi);
% Find the matching frequency from POD
[~,l]=min(abs(fx(vals(l))-f));
% Plot the DMD modes;
figure('name','DMD');
subplot(figs,1,1)
imagesc(real(reshape(Psi(1:r*c,l),r,c)));axis image;set(gca,'Ydir','Normal')
title(['\Psi',sprintf('_{f=%0.2f Hz} (u-component)',f(l))])
subplot(figs,1,2-d)
if com == 2
    imagesc(real(reshape(Psi(r*c+1:end,l),r,c)));axis image;set(gca,'Ydir','Normal')
    title(['\Psi',sprintf('_{f=%0.2f Hz} (v-component)',f(l))])
end
subplot(figs,1,3-d)
plot(real(Q(l,:)))
title(['Q',sprintf('_{f=%0.2f Hz}',f(l))])
subplot(figs,1,4-d)
[pxd, fxd]=pwelch(real(Q(l,:)),round(0.9*size(data,2)-1),round(0.8*size(data,2)-1),...
    2^12,freq);
plot(fx,mat2gray(px),'--r');
hold on
plot(fxd,mat2gray(pxd),'b');
legend('POD','DMD')
title(['P(Q',sprintf('_{f=%0.2f Hz})',f(l))])
set(figure(2),'position',[617    66   370   732])

