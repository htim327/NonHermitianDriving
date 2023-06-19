===============================================================
Presenting the code for the Anomalous Floquet-Anderson Insultor
===============================================================

This file presents how the MATLAB code is set up for simulating the Anomalous Floquet-Anderson Insulator and shows that the method used for calculating the pumped charge integrated over the appropriate time intervals using analytical integration approaches produces the correct results when compared to using numerical integration approaches. We start with the main file that runs the appropriate files and then plots the data.

.. code-block:: matlab

   clear;
   clc;
   % Store the error between the numerical and analytical integration
   % approaches used for calculating the charge pumped in errormat
   errormat = [];
   % Define how many noise realizations we use
   NTot = 100;
   % Define how accurate we want our numerical integration approach to be. The
   % higher the number, the more accurate.
   iter = 100;
   save('iter.mat','iter')
   for zamp = 1:NTot
       tic
       % Run TwoDimxyQ
       TwoDimxyQ
       % Load the information related to the charge pump calculated using the
       % numerical and analytical integration approaches
       load(['QTot_' num2str(iter) '.mat'])
       % Calculate the error between the analytical and numerical integration
       % approaches
       errormat = [errormat; (QTot(1,:)-QTot(2,:))];
       clearvars -except zamp errormat NTot iter
       % If you're impatient, the following lines let you become more patient.
       clc
       zamp
       toc
   end
   % Plot the data
   figure('units','normalized','outerposition',[0 0 1 1]);
   errorbar(1:length(errormat(1,:)),sum(errormat)/NTot,std(errormat)/sqrt(NTot),'LineWidth',2,'Color','m')
   title(['Error for iter=' num2str(iter)],'FontSize',40,'Interpreter','latex')
