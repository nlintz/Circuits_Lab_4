function Graphs
    close all;
    clf;
    
    set(0,'DefaultAxesFontSize', 20);
    set(0,'DefaultTextFontSize', 20);
    set(0,'defaultlinemarkersize',18);
    set(0,'defaultlinelinewidth',2);
    
    % Experiment 1 Data
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T1/Ib.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T1/Ie.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T1/Vb.mat');
    Ib_T1_exp1 = Ib;
    Ie_T1_exp1 = Ie;
    Vb_T1_exp1 = Vb;
    
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T2/Ib.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T2/Ie.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T2/Vb.mat');
    Ib_T2_exp1 = Ib;
    Ie_T2_exp1 = Ie;
    Vb_T2_exp1 = Vb;
    
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T3/Ib.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T3/Ie.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T3/Vb.mat');
    Ib_T3_exp1 = Ib;
    Ie_T3_exp1 = Ie;
    Vb_T3_exp1 = Vb;


    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T4/Ib.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T4/Ie.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_1/T4/Vb.mat');
    Ib_T4_exp1 = Ib;
    Ie_T4_exp1 = Ie;
    Vb_T4_exp1 = Vb;
    
%   Experiment 2 Data
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Sink/R_1_51k/Ix.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Sink/R_1_51k/Iz.mat');
    Ix_Sink_R1_exp2 = Ix;
    Iz_Sink_R1_exp2 = Iz;
    
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Sink/R_2_4.53k/Ix.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Sink/R_2_4.53k/Iz.mat');
    Ix_Sink_R2_exp2 = Ix;
    Iz_Sink_R2_exp2 = Iz;
    
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Sink/R_3_499k/Ix.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Sink/R_3_499k/Iz.mat');
    Ix_Sink_R3_exp2 = Ix;
    Iz_Sink_R3_exp2 = Iz;

    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Source/R_1_4.53k/Iy.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Source/R_1_4.53k/Iz.mat');
    Iy_Source_R1_exp2 = Iy;
    Iz_Source_R1_exp2 = Iz;

    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Source/R_2_51k/Iy.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Source/R_2_51k/Iz.mat');
    Iy_Source_R2_exp2 = Iy;
    Iz_Source_R2_exp2 = Iz;
    
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Source/R_3_499k/Iz.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_2/Source/R_3_499k/Iy.mat');
    Iy_Source_R3_exp2 = Iy;
    Iz_Source_R3_exp2 = Iz;
    
%   Experiment 3 Data
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Sink/R_1_51k/Ix.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Sink/R_1_51k/Iz.mat');
    Ix_Sink_R1_exp3 = Ix;
    Iz_Sink_R1_exp3 = Iz;

    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Sink/R_2_4.53k/Ix.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Sink/R_2_4.53k/Iz.mat');
    Ix_Sink_R2_exp3 = Ix;
    Iz_Sink_R2_exp3 = Iz;

    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Sink/R_3_499k/Ix.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Sink/R_3_499k/Iz.mat');
    Ix_Sink_R3_exp3 = Ix;
    Iz_Sink_R3_exp3 = Iz;

    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Source/R_1_51k/Iy.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Source/R_1_51k/Iz.mat');
    Iy_Source_R1_exp3 = Iy;
    Iz_Source_R1_exp3 = Iz;

    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Source/R_2_45.3k/Iz.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Source/R_2_45.3k/Iy.mat');
    Iy_Source_R2_exp3 = Iy;
    Iz_Source_R2_exp3 = Iz;

    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Source/R_3_499k/Iz.mat');
    load('/Users/nlintz/Documents/Olin/2013-2014/Spring/Circuits/Lab4/data/Experiment_3/Source/R_3_499k/Iy.mat');
    Iy_Source_R3_exp3 = Iy;
    Iz_Source_R3_exp3 = Iz;


    % Experiment 1
    sprintf('T1 Beta: %f, T2 Beta: %f, T3 Beta: %f, T4 Beta: %f', mean(Beta(Ib_T1_exp1, Ie_T1_exp1)), mean(Beta(Ib_T2_exp1, Ie_T2_exp1)), mean(Beta(Ib_T3_exp1, Ie_T3_exp1)), mean(Beta(Ib_T4_exp1, Ie_T4_exp1)))
    sprintf('T1 ISat: %f, T2 ISat: %f, T3 ISat: %f, T4 ISat: %f', mean(ISat(Ib_T1_exp1, Ie_T1_exp1, Vb_T1_exp1)), mean(ISat(Ib_T2_exp1, Ie_T2_exp1, Vb_T2_exp1)), mean(ISat(Ib_T3_exp1, Ie_T3_exp1, Vb_T3_exp1)), mean(ISat(Ib_T4_exp1, Ie_T4_exp1, Vb_T4_exp1)))

    figure(1);
    IV_Characteristic_Exp1([Ib_T1_exp1' Ib_T2_exp1' Ib_T3_exp1' Ib_T4_exp1'], [Ie_T1_exp1' Ie_T2_exp1' Ie_T3_exp1' Ie_T4_exp1'], [Vb_T1_exp1' Vb_T2_exp1' Vb_T3_exp1' Vb_T4_exp1']);
    
    figure(2);
    Collector_Current_Deviation_Exp1([Ib_T1_exp1' Ib_T2_exp1' Ib_T3_exp1' Ib_T4_exp1'], [Ie_T1_exp1' Ie_T2_exp1' Ie_T3_exp1' Ie_T4_exp1'], [Vb_T1_exp1' Vb_T2_exp1' Vb_T3_exp1' Vb_T4_exp1']);
    
    figure(3);
    Beta_Characteristic_Exp1([Ib_T1_exp1' Ib_T2_exp1' Ib_T3_exp1' Ib_T4_exp1'], [Ie_T1_exp1' Ie_T2_exp1' Ie_T3_exp1' Ie_T4_exp1'], [Vb_T1_exp1' Vb_T2_exp1' Vb_T3_exp1' Vb_T4_exp1']);

%     % Experiment 2
    figure(4);
    Iz_Characteristic([Ix_Sink_R1_exp2' Ix_Sink_R2_exp2' Ix_Sink_R3_exp2'], [Iz_Sink_R1_exp2' Iz_Sink_R2_exp2' Iz_Sink_R3_exp2']);
    
    figure(5);
    Iz_Characteristic([Iy_Source_R1_exp2' Iy_Source_R2_exp2' Iy_Source_R3_exp2'], [Iz_Source_R1_exp2' Iz_Source_R2_exp2' Iz_Source_R3_exp2']);

%     % Experiment 3
    figure(6);
    Iz_Characteristic([Ix_Sink_R1_exp3' Ix_Sink_R2_exp3' Ix_Sink_R3_exp3'], [Iz_Sink_R1_exp3' Iz_Sink_R2_exp3' Iz_Sink_R3_exp3']);

    figure(7);
    Iz_Characteristic([Iy_Source_R1_exp3' Iy_Source_R2_exp3' Iy_Source_R3_exp3'], [Iz_Source_R1_exp3' Iz_Source_R2_exp3' Iz_Source_R3_exp3']);
end

function res = Beta(Ib, Ie)
    Ic = Ie - Ib;
    beta = Ic./Ib;
    res = beta;
end

function res = ISat(Ib, Ie, Vb)
    Ic = Ie - Ib;
    Ut = .025;
    res = Ic./exp(Vb./Ut);
end

function IV_Characteristic_Exp1(Ib, Ie, Vb)
    Ic = Ib - Ie;
    semilogy(Vb, Ic, '.');
    hold on
    semilogy(Vb, max(Ib, 10^-10), '.');
    axis([0 .7 10e-12 10e-1]);

%     TODO - Add theoretical fits, labels, titles, legends
end

function Collector_Current_Deviation_Exp1(Ib, Ie, Vb)
    Ic = Ib - Ie;
    Ic_mean = zeros(150, 1);
    for i=1:4
        Ic_mean = Ic_mean + Ic(:, i);
    end
    Ic_Percentage_Diff = [(Ic(:,1)-Ic_mean)./Ic_mean (Ic(:,2)-Ic_mean)./Ic_mean (Ic(:,3)-Ic_mean)./Ic_mean (Ic(:,4)-Ic_mean)./Ic_mean];
    plot(Vb, Ic_Percentage_Diff, '.');
    axis([0 .7 -20 35 ]);
%     TODO - Add labels, titles, legends
end

function Beta_Characteristic_Exp1(Ib, Ie, Vb)
    Ic = Ib - Ie;
    beta_vals = Beta(Ib, Ie);
    semilogy(Ic, beta_vals, '.');
    axis([ -2e-9 6e-7 .02 1e5]);
end

function Iz_Characteristic(I, Iz)
    loglog(I, Iz, '.');
end
