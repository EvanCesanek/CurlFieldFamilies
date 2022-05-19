function R = EulerRotationMatrix(RotationAxis, Angle, AngleUnits, DebugFlag)
%% Function Description:
% This function returns a 3x3 ROTATION MATRIX (R). 
%
% AUTHOR: Sugato Ray | Created on: 11-AUG-2017 | ray.sugato[at]gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           PLEASE ACKNOWLEDGE THE AUTHOR IF YOU USE THIS CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%   rotationAxis = 'x' or, 'y' or, 'z'.                 ; Default: 'z'
%   Angle = angle of rotation (Numeric input only)      ; Default: 0
%   AngleUnits = 'D' for Degrees or, 'R' for Radians    ; Default: 'D'
%   DebugFlag = 1 (true) or, 0 (false)                  ; Default: 0
%
%   Note: DebugFlag = 1 enables printing Rotation Matrix (R) on Console.
%
% OUTPUT:
%   R = Rotation matrix about the axis of rotation
%
% EXAMPLE:
%   R = EulerRotationMatrix(rotationAxis, Angle, AngleUnits, DebugFlag);
%   R = EulerRotationMatrix('x', 60, 'D');     % Rotate by 60 degs about x-axis
%   R = EulerRotationMatrix('z', pi/3, 'R');   % Rotate by pi/3 rads about z-axis
%
%   Alternative approach: Define a function handle
%   
%       RotZ = @(ang) EulerRotationMatrix('z',ang,'D',1);
%       RotZ(60);
%
%       RotY = @(ang) EulerRotationMatrix('y',ang,'D',1);
%       RotY(60);
%
%       RotX = @(ang) EulerRotationMatrix('x',ang,'D',1);
%       RotX(60);
%
%   Test the value of RotX(60) against the following result:
%
%       Rx = [1,0,0; 0, cosd(60), -sind(60); 0, sind(60), cosd(60)]
%
%--------------------------------------------------------------------------
%% Version History
%
%   11-AUG-2017 | Version: 1.0 | First Release on Matlab FX website.
%   11-AUG-2017 | Version: 1.1 | Changed function name to EulerRotationMatrix.
%
%--------------------------------------------------------------------------
%% Code Section
    % Define default values of input parameters
        default_Axis = 'z';
        default_Angle = 0;
        default_Unit = 'D';
        default_DebugFlag = 0;
    
    % Handle exceptions in input parameters
        if nargin < 1 || isempty(RotationAxis) || strcmp(RotationAxis,'') || ~ischar(RotationAxis)
            RotationAxis = default_Axis;
            disp(['Warning: No input provided for axis of rotation: Setting axis = ',RotationAxis]);
        else
            RotationAxis = lower(strtrim(RotationAxis));
            RotationAxis = RotationAxis(1);
            if ~strcmp(RotationAxis,'x') && ~strcmp(RotationAxis,'y') && ~strcmp(RotationAxis,'z')
                RotationAxis = default_Axis;
            end    
        end 
        if nargin < 2 || isempty(Angle) || strcmp(Angle,'') || ~isnumeric(Angle)
            Angle = default_Angle;
            disp(['Warning: No input provided for angle: Setting angle = ',num2str(Angle)]);
        end
        if nargin < 3 || isempty(AngleUnits) || strcmp(AngleUnits,'') 
           AngleUnits = default_Unit; 
        else
            AngleUnits = upper(strtrim(AngleUnits)); % trim leading and trailing white spaces
            AngleUnits = AngleUnits(1);
            if ~strcmp(AngleUnits,'D') && ~strcmp(AngleUnits,'R')
                AngleUnits = default_Unit;
            end    
        end
        
        if nargin < 4 || isempty(DebugFlag) || strcmp(DebugFlag,'') || ~isnumeric(DebugFlag)
            DebugFlag = default_DebugFlag;
        else
            if (DebugFlag ~= 0) && (DebugFlag ~= 1)
                DebugFlag = default_DebugFlag;
            end    
        end
        
    % Dynamically choose Sine and Cosine functions for Units as 'D' or 'R'
        switch AngleUnits
            case 'D'
                func = {'sind','cosd'};
            case 'R'
                func = {'sin','cos'};
            otherwise
                disp('Error!!!... Angle Units are neither "D" nor "R".');
        end
    
    % Create Sine and Cosine function handles
        S = str2func(func{1}); % Sine function 
        C = str2func(func{2}); % Cosine function
    
    % Evaluate Sine and Cosine values for the given angle
        Sval = S(Angle);
        Cval = C(Angle);
    
    % Evaluate rotation matrix for the given rotationAxis
        switch RotationAxis
            case 'x'
                R = [ 1, 0, 0;...
                      0, Cval,-Sval;...
                      0, Sval, Cval ];
            case 'y'
                R = [ Cval, 0, Sval;...
                      0, 1, 0;...
                      -Sval, 0, Cval ];
            case 'z'
                R = [ Cval,-Sval, 0;...
                      Sval, Cval, 0;...
                      0, 0, 1 ];    
            otherwise
                disp('Error!!!... invalid Axis provided. Choose axis as one of (x / y / z).');
        end    
      
     % Display Rotation Matrix on Console   
     if (DebugFlag == 1)
         disp(['Rotation Matrix about axis: ',RotationAxis,' | angle: ',num2str(Angle),' | units: ',AngleUnits]);
         disp(' ');
         disp(R);
     end    
end