function [] = cp1_650703()
clear all;
close all;
clc;

%% Problem 1
A = load('WA_topo.txt'); % load file into matlab

%% Problem 2
A(A<0) = 0; % convert all negative value to zero
A = A * 3.28084; % convert meter to feet

%% Problem 3
figure % create a figure
h = surf(sqrt(A)); % creat a surface plot from data
axis tight
set(h, 'edgecolor', 'none') % set the edge color
title('3D Map of Washington State Topology') % label title 
xlabel('Longitude Indices') % label x-axis
ylabel('Latitude Indices') % label y-axis
zlabel('Altitude Indices') % label z-axis
pbaspect([360, 240, 10]) %adjust the aspect ratio
view(-22,42) % set view angle

%% Problem 4
B = A(:,500:940); % keep data only for western part
P = zeros(3,3); % create an array to store result
for i = 1:3 % find the three tallest peak
    [M,I] = max(B(:)); % find the current max mountain
    [I_row, I_col] = ind2sub(size(B),I); % get the row and col index
    P(i,1:3) = [I_col+499, I_row, sqrt(M)]; % store the height and index
    B(I_row, I_col) = 0; % set the highest altitude to zero
end

%% Problem 5
hold on % hold on figure
for i = 1:3
    plot3(P(i,1),P(i,2),P(i,3),'^k') % mark top 3 peak in the same figure
end

%% Problem 7
[SeattleRow, SeattleCol] = findRC(47.6046,-122.33057); % call findRC function and input the coordinate of the city
scatter3(SeattleCol, SeattleRow, sqrt(A(SeattleRow, SeattleCol)),50,'filled','k') % plot the city on the figure
[SpokaneRow, SpokaneCol] = findRC(47.6589,-117.4250); % call findRC function and input the coordinate of the city
scatter3(SpokaneCol, SpokaneRow, sqrt(A(SpokaneRow, SpokaneCol)),75,'filled','k') % plot the city on the figure
[WenatcheeRow, WenatcheeCol] = findRC(47.4233,-120.3253); % call findRC function and input the coordinate of the city
scatter3(WenatcheeCol, WenatcheeRow, sqrt(A(WenatcheeRow, WenatcheeCol)),50,'filled','k') % plot the city on the figure
hold off % release the figure

%% Problem 15
Seattle_Row = 47.6046; % the coordinates of start point
Seattle_Col = -122.3305;
Canada_Row = 49.00; % the coordinates of end point
Canada_Col = -122.3305;
fprintf('From Seattle to Canada. \n')
bikeRides(Seattle_Row, Seattle_Col, Canada_Row, Canada_Col, A) % call function and input two sets of coordinates

%% Problem 16
Walla_Row = 46.0650; % the coordinates of start point
Walla_Col = -117.4250;
Spokane_Row = 47.6589; % the coordinates of end point
Spokane_Col = -117.4250;
fprintf('From Walla Walla to Spokane. \n')
bikeRides(Walla_Row, Walla_Col, Spokane_Row, Spokane_Col, A) % call function and input two sets of coordinates

%% Problem 17
Wenatchee_Row = 47.4233; % the coordinates of start point
Wenatchee_Col = -120.3253;
Wilderness_Row = 48.800; % the coordinates of end point
Wilderness_Col = -120.3253;
fprintf('From Wenatchee to wilderness. \n')
bikeRides(Wenatchee_Row, Wenatchee_Col, Wilderness_Row, Wilderness_Col, A) % call function and input two sets of coordinates

end

%% Problem 6 Function
function [row, column] = findRC(latitude, longitude)
southedge = 45.56666700000; % edge coordinates
northedge = 49.00000000000;
westedge = -124.71666700000;
eastedge = -116.88333300000;
gridcellwidth =  0.00833333334;
if latitude < southedge || latitude > northedge % if coordinate is outside of the state return -999
    fprintf('The coordinates is outside of the Washington state \n')
    row = -999;
    column = -999;
    return
elseif longitude < westedge || longitude > eastedge % if coordinate is outside of the state return -999
    fprintf('The coordinates is outside of the Washington state \n')
	row = -999;
	column = -999;
	return
else
    row = floor((latitude - southedge)/gridcellwidth) + 1; % calculate the index of inout location
    column = floor((longitude - westedge)/gridcellwidth) + 1;
    return
end
    
end

%% Problem 8 Function
function [] = bikeRides(lat1, lon1, lat2, lon2, A)
southedge = 45.56666700000; % edge coordinates
gridcellwidth =  0.00833333334;
[spos_Row, spos_Col] = findRC(lat1, lon1); % call findRC function to get location index
[epos_Row, epos_Col] = findRC(lat2, lon2);

%% Problem 9
z = sqrt(A(spos_Row:epos_Row, spos_Col)); % store all altitudes for the path
maxz = max(z); % find the max altitude
fprintf('The maximum altitude along the ride is %d. \n', maxz)             

%% Problem 10
avgz = sum(z)/length(z); % calculate the average altitude
fprintf('The average altitude along the ride is %d. \n \n', avgz)

%% Problem 11
a = zeros(length(z),1); % create an array with the size we needed
a(1) = z(1); % set the initial value
for i = 1:length(z) - 1
    a(i+1) = a(i) + (z(i+1)-a(i))/(i+1); % calculate the running average altitude
end

%% Problem 12
x = zeros(length(z),1); % create an array with the size we needed
%syms xlat
for i = 1:length(z)
    x(i) = (spos_Row + i - 1) * gridcellwidth + southedge; % calculate the latitude by using index
    %solve(spos_Row + i - 1 == floor((xlat - southedge)/gridcellwidth + 1), xlat);
end

%% Problem 13
figure % new figure
subplot(2,1,1); % subplot
plot(x,z,'r') % plot latitude vs. altitude
axis tight % adjust the range of plot
title('Latitude vs. Average Altitude')
xlabel('Latitude x')
ylabel('Average Altitude z')

%% Problem 14
subplot(2,1,2);
plot(x,a,':b') % plot running average altitude vs. latitude
axis tight % adjust the range of plot
title('Latitude vs. Running Average Altitude')
xlabel('Latitude x')
ylabel('Running Average Altitude a')

end


