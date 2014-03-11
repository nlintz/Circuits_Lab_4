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

%     figure(1)
%     hold on
    'T1 Isat:'
    ISat(Vb_T1_exp1(50:100), (Ie_T1_exp1(50:100)-Ib_T1_exp1(50:100)), Vb_T1_exp1)
    'T2 Isat:'
    ISat(Vb_T2_exp1(50:100), (Ie_T2_exp1(50:100)-Ib_T2_exp1(50:100)), Vb_T2_exp1)
    'T3 Isat:'
    ISat(Vb_T3_exp1(50:100), (Ie_T3_exp1(50:100)-Ib_T3_exp1(50:100)), Vb_T3_exp1)
    'T4 Isat:'
    ISat(Vb_T4_exp1(50:100), (Ie_T4_exp1(50:100)-Ib_T4_exp1(50:100)), Vb_T4_exp1)

    
    figure(1)
    IV_Characteristic_Exp1([Ib_T1_exp1' Ib_T2_exp1' Ib_T3_exp1' Ib_T4_exp1'], [Ie_T1_exp1' Ie_T2_exp1' Ie_T3_exp1' Ie_T4_exp1'], [Vb_T1_exp1' Vb_T2_exp1' Vb_T3_exp1' Vb_T4_exp1']);
    title('Collector/Base Current-Voltage Characteristic');
    xlabel('Voltage (Volts)');
    ylabel('Current (Amperes)');
    legend('Q_{1}', 'Q_{2}', 'Q_{3}', 'Q_{4}')
    
    figure(2);
    Collector_Current_Deviation_Exp1([Ib_T1_exp1' Ib_T2_exp1' Ib_T3_exp1' Ib_T4_exp1'], [Ie_T1_exp1' Ie_T2_exp1' Ie_T3_exp1' Ie_T4_exp1'], [Vb_T1_exp1' Vb_T2_exp1' Vb_T3_exp1' Vb_T4_exp1']);
    title('Collector Current Percentage Difference From Mean');
    xlabel('Voltage (Volts)');
    ylabel('Fractional Difference From Mean');
    legend('Q_{1}', 'Q_{2}', 'Q_{3}', 'Q_{4}');
    axis([0 .7 -15 10]);

    figure(3);
    Beta_Characteristic_Exp1([Ib_T1_exp1' Ib_T2_exp1' Ib_T3_exp1' Ib_T4_exp1'], [Ie_T1_exp1' Ie_T2_exp1' Ie_T3_exp1' Ie_T4_exp1'], [Vb_T1_exp1' Vb_T2_exp1' Vb_T3_exp1' Vb_T4_exp1']);
    title('Current Gain');
    xlabel('Voltage (Volts)');
    ylabel('Current Gain');
    legend('Q_{1}', 'Q_{2}', 'Q_{3}', 'Q_{4}')
    axis([-1e-7 7e-7 .04 1e5]);

    % Experiment 2
    figure(4);
    Iz_Characteristic([Ix_Sink_R1_exp2' Ix_Sink_R2_exp2' Ix_Sink_R3_exp2'], -1*[Iz_Sink_R1_exp2' Iz_Sink_R2_exp2' Iz_Sink_R3_exp2']);
    hold on;
    Sqrt_TheoreticalFit_Sink([Ix_Sink_R1_exp2' Ix_Sink_R2_exp2' Ix_Sink_R3_exp2'], [51e3, 4.53e3, 499e3]);
    title('Translinear Circuit 1 - Output Current (Sink)');
    xlabel('I_{x} (Amperes)');
    ylabel('I_{z} (Amperes)');
    legend('I_{y} = 4.4150e-05 (Amps)', 'I_{y} = 3.9216e-06 (Amps)', 'I_{y} = 4.0080e-07 (Amps)', 'Theoretical Fit: I_{y} = 4.4150e-05 (Amps)', 'Theoretical Fit: I_{y} = 3.9216e-06 (Amps)', 'Theoretical Fit: I_{y} = 4.0080e-07 (Amps)');

    figure(5);
    Iz_Characteristic(-1.*[Iy_Source_R1_exp2' Iy_Source_R2_exp2' Iy_Source_R3_exp2'], -1.*[Iz_Source_R1_exp2' Iz_Source_R2_exp2' Iz_Source_R3_exp2']);
    hold on
    Sqrt_TheoreticalFit_Source([Iy_Source_R1_exp2' Iy_Source_R2_exp2' Iy_Source_R3_exp2'], [4.53e3, 51e3, 499e3]);
    axis([1e-8 1e-2 1e-7 1e-2])
    title('Translinear Circuit 1 - Output Current (Source)');
    xlabel('I_{y} (Amperes)');
    ylabel('I_{z} (Amperes)');
    legend('I_{x} = 5.5188e-04 (Amps)', 'I_{x} = 4.9020e-05 (Amps)', 'I_{x} = 5.0100e-06 (Amps)', 'Theoretical Fit: I_{x} = 5.5188e-04 (Amps)', 'Theoretical Fit: I_{x} = 4.9020e-05 (Amps)', 'Theoretical Fit:  I_{x} = 5.0100e-06 (Amps)');

%     % Experiment 3
    figure(6);
    Iz_Characteristic([Ix_Sink_R1_exp3' Ix_Sink_R2_exp3' Ix_Sink_R3_exp3'], -1.*[Iz_Sink_R1_exp3' Iz_Sink_R2_exp3' Iz_Sink_R3_exp3']);
    hold on
    Squared_TheoreticalFit_Sink([Ix_Sink_R1_exp3' Ix_Sink_R2_exp3' Ix_Sink_R3_exp3'], [51e3, 4.53e3, 499e3]);
    title('Translinear Circuit 2 - Output Current (Sink)');
    xlabel('I_{x} (Amperes)');
    ylabel('I_{z} (Amperes)');
    legend('I_{y} = 4.4150e-05 (Amps)', 'I_{y} = 3.9216e-06 (Amps)', 'I_{y} = 4.0080e-07 (Amps)', 'Theoretical Fit: I_{y} = 4.4150e-05 (Amps)', 'Theoretical Fit: I_{y} = 3.9216e-06 (Amps)', 'Theoretical Fit: I_{y} = 4.0080e-07 (Amps)');

    figure(7);
    Iz_Characteristic(-1.*[Iy_Source_R1_exp3' Iy_Source_R2_exp3' Iy_Source_R3_exp3'], -1.*[Iz_Source_R1_exp3' Iz_Source_R2_exp3' Iz_Source_R3_exp3']);
    hold on;
    Inverse_TheoreticalFit_Source([Iy_Source_R1_exp2' Iy_Source_R2_exp2' Iy_Source_R3_exp2'], [51e3, 4.53e3, 499e3]);
    title('Translinear Circuit 2 - Output Current (Source)');
    xlabel('I_{y} (Amperes)');
    ylabel('I_{z} (Amperes)');
    legend('I_{x} = 5.5188e-04 (Amps)', 'I_{x} = 4.9020e-05 (Amps)', 'I_{x} = 5.0100e-06 (Amps)', 'Theoretical Fit: I_{x} = 5.5188e-04 (Amps)', 'Theoretical Fit: I_{x} = 4.9020e-05 (Amps)', 'Theoretical Fit:  I_{x} = 5.0100e-06 (Amps)');

end

function res = Beta(Ib, Ie)
    Ic = Ie - Ib;
    beta = Ic./Ib;
    res = beta;
end
    
function res = ISat(x_sample, y_sample, input) %Are we in fwd active?
    fit = polyfit(x_sample, log(y_sample), 1);
    res = exp(real(fit(2)));
end

function IV_Characteristic_Exp1(Ib, Ie, Vb)
    Ic = Ib - Ie;
    semilogy(Vb, Ic, '.');
    hold on
    semilogy(Vb, max(Ib, 10^-10), '.');
    axis([0 .7 10e-12 10e-1]);

%     TODO - Add labels, titles, legends
end

function Collector_Current_Deviation_Exp1(Ib, Ie, Vb)
    Ic = Ib - Ie;
    Ic_mean = zeros(150, 1);
    for i=1:4
        Ic_mean = Ic_mean + Ic(:, i);
    end
    Ic_mean = Ic_mean./4;
    Ic_Percentage_Diff = [(Ic(:,1)-Ic_mean)./Ic_mean (Ic(:,2)-Ic_mean)./Ic_mean (Ic(:,3)-Ic_mean)./Ic_mean (Ic(:,4)-Ic_mean)./Ic_mean];
    plot(Vb, Ic_Percentage_Diff, '.');
    axis([0 .7 -20 35 ]);
%     TODO - Add labels, titles, legends
end

function Beta_Characteristic_Exp1(Ib, Ie, Vb)
    Ic = Ib - Ie;
    beta_vals = Beta(Ib, Ie);
    semilogy(Ic, beta_vals, '.');
  % Compare against http://www.analog.com/static/imported-files/data_sheets/MAT14.pdf
%     axis([ -2e-9 6e-7 .02 1e5]);
end

function Iz_Characteristic(I, Iz)
    loglog(I, Iz, '.');
%     TODO - Add labels, titles, legends
end

function Sqrt_TheoreticalFit_Sink(Ix, R)
Vy = .2;
Iy = Vy./R;
Iz = [(Iy(:,1).*Ix(:,1)).^.5 (Iy(:,2).*Ix(:,2)).^.5 (Iy(:,3).*Ix(:,3)).^.5];
loglog(Ix, Iz);
%     TODO - Add labels, titles, legends

end

function Sqrt_TheoreticalFit_Source(Iy, R)
Vx = 2.5;
Ix = Vx./R;
Iz = [(abs(Ix(:,1).*Iy(:,1))).^.5 abs((Ix(:,2).*Iy(:,2))).^.5 abs((Ix(:,3).*Iy(:,3))).^.5];
loglog(-1.*Iy, Iz);
%     TODO - Add labels, titles, legends
end

function Squared_TheoreticalFit_Sink(Ix, R)
Vy = .2;
Iy = Vy./R;
Iz = [(Ix(:,1)).^2./(Iy(:,1)) (Ix(:,2)).^2./(Iy(:, 2)) (Ix(:,3)).^2./(Iy(:, 3))];
loglog(Ix, Iz);
%     TODO - Add labels, titles, legends
end

function Inverse_TheoreticalFit_Source(Iy, R)
Vx = 2.5;
Ix = Vx./R;
Iz = [(Ix(:,1)).^2./(Iy(:,1)) (Ix(:,2)).^2./(Iy(:, 2)) (Ix(:,3)).^2./(Iy(:, 3))];
loglog(-1.*Iy, -1.*Iz);
%     TODO - Add labels, titles, legends

end
