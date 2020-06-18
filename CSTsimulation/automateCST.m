% This script is to automate CST MWS: Starting, saving and running a new
% simulation. The new projects are not created here, that is done using CST
% UI.
% So far, this does not take into account the license issue - if CST drops
% it after a few hours, or if the VPN stops after 9 hours, so that probably
% needs to be taken care of manually.
% -------------------------------------------------------------------------
% Pragya Sharma, 11 June 2020
% ps847@cornell.edu
% -------------------------------------------------------------------------

filepath = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\TrueSize\Sim915MHz_True\';
filename = { 'RFimaging_3D_Planar_TrueSized_03_1obj1_x60y60.cst'};
nFiles = length(filename);
% Initialize CST using ActiveX framework
cst = actxserver('CSTStudio.application');

for i = 1:nFiles
    % invoke calls the units of the controlled program, CST
    % Sometimes this fails, but re-run resolves the error.
    mws = invoke(cst, 'OpenFile',[filepath,filename{i}]);
    % Calling the solver, using preset solver
    solver = invoke(mws,'Solver');
    % Defining the frequency range - most of the simulations set to 0.6-1.1 GHz
    % by default. If saved to new format, skip this.
    invoke(solver,'FrequencyRange','.8','.93'); % Only simulating in range 0.8-0.93 GHz
    % If results already exist, the next step will keep waiting for that
    % diaglog box - so take care!
    invoke(solver,'Start');

    % Selecting S-Parameters
    mws.invoke('SelectTreeItem','1D Results\S-Parameters');
    % Selecting export
    export = invoke(mws,'TOUCHSTONE');
    invoke(export,'FileName',[filepath,filename{i}(1:end-4)]);
    invoke(export,'Format','RI');
    invoke(export,'Renormalize',0);
    invoke(export,'Impedance', 50);
    invoke(export,'Write');
    invoke(mws,'save');
end
    
    
invoke(cst,'quit');