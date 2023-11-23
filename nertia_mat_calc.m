close all; clear; clc;

% X - Runs along orthogonal to boom and Z
% Y - Direction of boom
% Z - Runs along orthogonal to boom AND panels
% Reference Point: Coordinate center of bus midpoint

%data = [m(KG), x(M), y(M), z(M)];
CDH   = [  52, -0.3, +00.50, -0.3]; 
STR   = [ 325, -0.3, -00.00, -0.3]; 
ADC   = [  78, -0.2, -00.75, -0.2]; 
CMS   = [  91, +0.2, +00.50, +0.2]; 
PWR   = [ 273, +0.3, +02.00, +0.3]; 
GNC   = [  39, +0.3, +01.50, +0.3]; 
PRP   = [ 169, 0, -05.00, 0]; 
THR   = [  78, 0, -00.75, 0]; 
SCI   = [ 195, 0, -10.00, 0]; 
PRO   = [  78, 0, -03.00, 0]; 

computeInertialMatrix(CDH, STR, ADC, CMS, PWR, GNC, PRP, THR, SCI, PRO);

function [summaryTable, CG, totalInertia] = computeInertialMatrix(varargin)
    % This function computes the inertial matrices for given sets of mass and coordinates.
    % Each input is an array [mass, x, y, z].
    % It returns a summary table with mass, coordinates, and inertial matrix elements,
    % the center of gravity (CG), and the total inertia matrix.

    n = nargin; % Number of input arguments
    summaryData = []; % Initialize summary data array
    totalMass = 0;
    weightedSumX = 0;
    weightedSumY = 0;
    weightedSumZ = 0;
    
    % Initialize total moments and products of inertia
    totalIxx = 0;
    totalIyy = 0;
    totalIzz = 0;
    totalIxy = 0;
    totalIxz = 0;
    totalIyz = 0;

    for i = 1:n
        inputData = varargin{i};
        m = inputData(1); % Mass
        x = inputData(2); % X-coordinate
        y = inputData(3); % Y-coordinate
        z = inputData(4); % Z-coordinate
        
        % Accumulate total mass and weighted position sums for CG calculation
        totalMass = totalMass + m;
        weightedSumX = weightedSumX + m * x;
        weightedSumY = weightedSumY + m * y;
        weightedSumZ = weightedSumZ + m * z;

        % Compute terms of the individual inertial matrix
        Ixx = m * (y^2 + z^2);
        Iyy = m * (x^2 + z^2);
        Izz = m * (x^2 + y^2);
        Ixy = m * x * y;
        Ixz = m * x * z;
        Iyz = m * y * z;
        
        % Accumulate total moments and products of inertia
        totalIxx = totalIxx + Ixx;
        totalIyy = totalIyy + Iyy;
        totalIzz = totalIzz + Izz;
        totalIxy = totalIxy + Ixy;
        totalIxz = totalIxz + Ixz;
        totalIyz = totalIyz + Iyz;

        % Append data to the summary array
        summaryData = [summaryData; m x y z Ixx Iyy Izz Ixy Ixz Iyz];
    end

    % Calculate CG
    CG = [weightedSumX, weightedSumY, weightedSumZ] / totalMass;

    % Calculate the inertia matrix about the CG using the parallel axis theorem
    totalInertia = [totalIxx - totalMass*(CG(2)^2 + CG(3)^2), ...
                    totalIxy - totalMass*CG(1)*CG(2), ...
                    totalIxz - totalMass*CG(1)*CG(3); ...
                    totalIxy - totalMass*CG(1)*CG(2), ...
                    totalIyy - totalMass*(CG(1)^2 + CG(3)^2), ...
                    totalIyz - totalMass*CG(2)*CG(3); ...
                    totalIxz - totalMass*CG(1)*CG(3), ...
                    totalIyz - totalMass*CG(2)*CG(3), ...
                    totalIzz - totalMass*(CG(1)^2 + CG(2)^2)];
    
    % Convert the summary data to a table
    summaryTable = array2table(summaryData, 'VariableNames', ...
                               {'Mass', 'X', 'Y', 'Z', 'Ixx', 'Iyy', 'Izz', 'Ixy', 'Ixz', 'Iyz'});
    
    % Display the CG and total inertia matrix
    disp('Center of Gravity (CG):');
    disp(CG);
    disp('Total Inertia Matrix:');
    disp(totalInertia);
    disp('Table')
    disp(summaryTable);
end


% function computeInertialMatrix(varargin)
%     % This function computes the inertial matrices for given sets of mass and coordinates.
%     % Each input is an array [mass, x, y, z].
%     % It prints a summary table with mass, coordinates, and inertial matrix elements.
% 
%     n = nargin; % Number of input arguments
%     summaryTable = []; % Initialize summary table
% 
%     for i = 1:n
%         inputData = varargin{i};
%         m = inputData(1); % Mass
%         x = inputData(2); % X-coordinate
%         y = inputData(3); % Y-coordinate
%         z = inputData(4); % Z-coordinate
% 
%         % Compute terms of the inertial matrix
%         Ixx = m*(y^2 + z^2);
%         Iyy = m*(x^2 + z^2);
%         Izz = m*(x^2 + y^2);
% 
%         % Append to the summary table
%         summaryTable = [summaryTable; m x y z Ixx Iyy Izz];
%     end
% 
%     % Convert the summary data to a table and display it
%     summaryTable = array2table(summaryTable, 'VariableNames', {'Mass', 'X', 'Y', 'Z', 'Ixx', 'Iyy', 'Izz'});
%     disp(summaryTable);
% end

% % summation as an approximation to integration
% function [inertiaTable, totalInertia] = computeInertialMatrix(varargin)
%     % This function computes the inertial matrices for given sets of mass and coordinates.
%     % Each input is an array [mass, x, y, z].
%     % It returns a summary table with mass, coordinates, and inertial matrix elements
%     % and the total inertial matrix obtained by summing the contributions of each point mass.
% 
%     n = nargin; % Number of input arguments
%     summaryTable = []; % Initialize summary table
%     totalIxx = 0;
%     totalIyy = 0;
%     totalIzz = 0;
%     totalIxy = 0;
%     totalIxz = 0;
%     totalIyz = 0;
% 
%     for i = 1:n
%         inputData = varargin{i};
%         m = inputData(1); % Mass
%         x = inputData(2); % X-coordinate
%         y = inputData(3); % Y-coordinate
%         z = inputData(4); % Z-coordinate
% 
%         % Compute terms of the inertial matrix for point masses
%         Ixx = m * (y^2 + z^2);
%         Iyy = m * (x^2 + z^2);
%         Izz = m * (x^2 + y^2);
%         Ixy = -m * x * y;
%         Ixz = -m * x * z;
%         Iyz = -m * y * z;
% 
%         % Sum the contributions for the total inertia
%         totalIxx = totalIxx + Ixx;
%         totalIyy = totalIyy + Iyy;
%         totalIzz = totalIzz + Izz;
%         totalIxy = totalIxy + Ixy;
%         totalIxz = totalIxz + Ixz;
%         totalIyz = totalIyz + Iyz;
% 
%         % Append to the summary table
%         summaryTable = [summaryTable; m x y z Ixx Iyy Izz Ixy Ixz Iyz];
%     end
% 
%     % Convert the summary data to a table and display it
%     summaryTable = array2table(summaryTable, ...
%                                'VariableNames', ...
%                                {'Mass', 'X', 'Y', 'Z', 'Ixx', 'Iyy', 'Izz', 'Ixy', 'Ixz', 'Iyz'});
%     inertiaTable = summaryTable;
% 
%     % Create the total inertial matrix
%     totalInertia = [-totalIxx, -totalIxy, -totalIxz; ...
%                     -totalIxy, -totalIyy, -totalIyz; ...
%                     -totalIxz, -totalIyz, -totalIzz];
% 
%     % Display the total inertial matrix
%     disp('Total Inertial Matrix:');
%     disp(totalInertia);
% end