clc
clear
close all

readerobj=VideoReader('KO-01-300Pa.avi');

filename = 'cut_1.mat';
imgname = 'cut_1';

numFrames = get(readerobj, 'NumberOfFrames');

%% parameter initialization

para.pressure= 400     % outter pressure, unit Pa
para.phi=0.965;      % 3*phi/(2*pi), phi = 2.02, 0.965
para.a = 0.000001  % inner radius, unit meter
para.E = 1000; %elastic modulus of the cell, unit kPa

%%
for i = 1:10:numFrames
    vidFrames = read(readerobj,i);
    imshow(vidFrames(:,:,:,1));
    set(gcf,'outerposition',get(0,'screensize'));
    title('use left button on mouse to double click on target point.');
    [x,y] = getpts;
    position (i,1)=x;
    position (i,2)=y;
end

%%
% define position x and position y
position_x = position(:,1);
position_y = position(:,2);

% get rid of the zeros in x and y
position_x(position_x==0) = [];
position_y(position_y==0) = [];

% define new matrix for saving position information
position2(:,1) = position_x;
position2(:,2) = position_y;

save('testing.mat','position2');

%% curve fitting

% calculate distance
N = size(position2(:,1));
%distance = zeros(N);
distance_x = zeros(N);
distance_y = zeros(N);
distance = zeros(N);
time = 1:1:N;
for i = 1:N
    %distance(i) = (position2(i,1)-position2(1,1)).^2 + (position2(i,2)-position2(1,2).^2);
    distance_x(i) = position2(i,1)-position2(1,1);
    distance_y(i) = position2(i,2)-position2(1,2);
    distance(i) = sqrt( distance_x(i)^2+ distance_y(i)^2 ); 
end

%% calculate the inner pressure

deformation = distance(N,1) - distance(1,1); % unit in pixel
deformation = deformation*0.000002; % unit in meter
pressure = deformation*para.E/para.phi/para.a



%% save data
para.pressure = pressure;
para.deformation = deformation;

save('cell1.mat','para','position2');












%% plot the tracked data 
plot(distance)

%% fitting for the slope 

[stime,sdisplacement,ind]=manual_select_line_roi(time,distance,'select indentation roi',i,'micropipette');

%% convert the E, P_outter, diamter, to P inner
[Esample,EL,EH,cfL,gofR2]=fit_youngs_modulus_linear(sDisplacement,sForce,para,0,1)


%% save results


%%
% for i=1:numFrames
%     mov(i).cdata = vidFrames(:,:,:,i);
%     mov(i).colormap = [];
% end
% hf = figure; 
% set(hf, 'position', [150 150 readerobj.Width readerobj.Height]);
% FrameRate=25;
% movie(hf, mov, 1,FrameRate*10);
% i=1;
% for k=1:numFrames 
%     imshow(vidFrames(:,:,:,k));
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
% 
% title('use left button on mouse to double click on target point.')
% [x y] = getpts(1)
% position (i,1)=x(1);
% position (i,2)=y(1);
% 
% i=i+1;
% end

% figure(3)
% plot(position(:,1),position(:,2),'b*-')
% set(gca,'ydir','reverse');
% saveas(gcf,imgname)
% hf = figure; 
% set(hf, 'position', [150 150 readerobj.Width readerobj.Height]);
% FrameRate=50;
% movie(hf, mov, 1,FrameRate);
% imshow(vidFrames(:,:,1,i))
% hold on
% plot (position(:,1), position(:,2),'b*');
% set(gca,'ydir','reverse');
% hold off
% save(filename, 'position');

