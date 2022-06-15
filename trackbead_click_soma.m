%% initialization 
clc
clear
close all

readerobj=VideoReader('2043_shortened.avi');
numFrames = get(readerobj, 'NumberOfFrames');

%% track the soma location for calculating cell migration direction
for i = 1:1:numFrames
    vidFrames = read(readerobj,i);
    imshow(vidFrames(:,:,:,1));
    set(gcf,'outerposition',get(0,'screensize'));
    title('double click for soma location');
    [x,y] = getpts;
    position0 (i,1)=x;
    position0 (i,2)=y;
end

%% clean up zeros in tracking and save position data
% define position x and position y
position0_x = position0(:,1);
position0_y = position0(:,2);

% get rid of the zeros in x and y
position0_x(position0_x==0) = [];
position0_y(position0_y==0) = [];

% define new matrix for saving position information
position_0(:,1) = position0_x;
position_0(:,2) = position0_y;

plot(position_0(:,1), position_0(:,2))
%% calculate the migration angles
N = size(position_0(:,1));
angle0_x = zeros(N);
angle0_y = zeros(N);
angle0 = zeros(N);
time = 1:1:N;
for i = 1:N-1
    angle0_x(i) = position_0(i+1,1)-position_0(i,1);
    angle0_y(i) = position_0(i+1,2)-position_0(i,2);
    angle0(i) = angle0_y(i)/angle0_x(i);
    if angle0_x(i) > 0        
        angle0(i) = atan(angle0(i));
    elseif (angle0_x(i)<0)&(angle0_y(i)<0)
        angle0(i) = atan(angle0(i))-3.14;
    elseif (angle0_x(i)<0)&(angle0_y(i)>0)
        angle0(i) = atan(angle0(i))+3.14;
    else
        display('invalid angle')
    end
    angle0(i)=57.3*angle0(i);
end
plot(angle0)

%% save soma tracking
save('cell_1_soma.mat','position_0', 'angle0'); %remember to rename .mat file








%% initialization 
clc
clear
close all

readerobj=VideoReader('2043_shortened.avi');
numFrames = get(readerobj, 'NumberOfFrames');

%% Tracking for processes, requiring two points tracking
for i = 1:1:numFrames
    vidFrames = read(readerobj,i);
    imshow(vidFrames(:,:,:,1));
    set(gcf,'outerposition',get(0,'screensize'));
    title('use left button on mouse to double click on target point.');
    [x,y] = getpts;
    position1 (i,1)=x;
    position1 (i,2)=y;

    [x,y] = getpts;
    position2 (i,1)=x;
    position2 (i,2)=y;
end

%% clean up position tracking data
% define position x and position y
position1_x = position1(:,1);
position1_y = position1(:,2);

position2_x = position2(:,1);
position2_y = position2(:,2);

% get rid of the zeros in x and y
position1_x(position1_x==0) = [];
position1_y(position1_y==0) = [];

position2_x(position2_x==0) = [];
position2_y(position2_y==0) = [];

% define new matrix for saving position information
position_1(:,1) = position1_x;
position_1(:,2) = position1_y;

position_2(:,1) = position2_x;
position_2(:,2) = position2_y;

save('cell1_process_1_temp.mat','position_1','position_2');
%% length calculation, check whether the size of position_1 and position_2 matche

N = size(position_2(:,1));
%distance = zeros(N);
length_x = zeros(N);
length_y = zeros(N);
length = zeros(N);
time = 1:1:N;
for i = 1:N
    length_x(i) = position_2(i,1)-position_1(i,1);
    length_y(i) = position_2(i,2)-position_1(i,2);
    length(i) = sqrt( length_x(i)^2+ length_y(i)^2 ); 
end
plot(length)

%% growth cone angle calculation
N = size(position_2(:,1));
angle_x = zeros(N);
angle_y = zeros(N);
angle = zeros(N);
time = 1:1:N;
for i = 1:N
    angle_x(i) = position_2(i,1)-position_1(i,1);
    angle_y(i) = position_2(i,2)-position_1(i,2);
    angle(i) = angle_y(i)/angle_x(i);
    if angle_x(i) > 0        
        angle(i) = atan(angle(i));
    elseif (angle_x(i)<0)&(angle_y(i)<0)
        angle(i) = atan(angle(i))-3.14;
    elseif (angle_x(i)<0)&(angle_y(i)>0)
        angle(i) = atan(angle(i))+3.14;
    else
        display('invalid angle')
    end
    angle(i)=57.3*angle(i);
end
plot(angle)

%% save data
para.length = length; %unit in pixels
para.angle = angle; %absolute value
save('WT_cell1_process1.mat','para','position_1','position_2') %remember to rename .mat file
