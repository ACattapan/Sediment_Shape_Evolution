function [varargout]=BetaFit(T,Parameter,FitTest,alpha,Lith,Orig,Location,Plotting_flag,PlotType)
% THIS FUNCTION RECEIVES INPUTS FROM A SINGLE LITHOLOGY, ORIGIN AND LOCATION
% AND PRODUCES FITS OF BETA DISTRIBUTIONS AND THEIR PLOT AGAINS DATA
% HISTOGRAM, IF REQUESTED.
%--------------------------------------------------------------------------
arguments
    T (:,:) table
    Parameter (1,1) string
    FitTest (1,1) string
    alpha (1,1) double
    Lith (:,:) string
    Orig (:,:) string
    Location (1,1) %{mustBeNumeric(LocPath)}
    Plotting_flag (1,1) string
    PlotType (1,1) string
end
%--------------------------------------------------------------------------
%                       PARAMETERS POSSIBLE VALUES
%--------------------------------------------------------------------------
% FitTest: "KS" or "AD" for Kolmogorov-Smirnov and Anderson-Darling
% respectively
% PlotFlag: "Y" or "y" AND "N" OR "n"
% PlotType: "pdf" OR "cdf"
%--------------------------------------------------------------------------
%                       EXTRACT FROM TABLE T THE DATA TO FIT
%--------------------------------------------------------------------------
Data2Fit=T.(Parameter)(T.Lithology==Lith&T.Origin==Orig&T.Location==Location);
sz=size(Data2Fit,1);
%   CHECK WHETHER DATA ARE SUITABLE FOR A FIT WITH A BETA DISTRIBUTION
if or(min(Data2Fit)<0,max(Data2Fit)>1)
    fprintf("WARNING: Input values are outside the ranges for Beta distributions!")
    return
else
%                       FIT BETA DISTRIBUTION TO THE DATA
    if sz>0
        beta=betafit(Data2Fit);
        a=beta(1);
        b=beta(2);
        % Generate a sample extracted from the Beta distribution
        % that was fit to the raw data. The sample size is the same as 
        % that of the raw data available
        Beta_randomNumbers=betarnd(a,b,sz,1);
        if FitTest=="KS"
            [hip,p]=kstest2(Data2Fit,Beta_randomNumbers,'Alpha',alpha);
        elseif FitTest=="AD"
            TestData=[Data2Fit, ones(sz,1);Beta_randomNumbers,2.*ones(sz,1)];
            [hip ADK ADKs ADKc p] = AnDarksamtestNum_Modified(TestData,alpha);
        end
        h=double(hip);
        varargout={a,b,h,p};
    end
    %--------------------------------------------------------------------------
    % PLOTTING
    %--------------------------------------------------------------------------
    if or(Plotting_flag=="Y",Plotting_flag=="y")
        f4 = figure;
        sgtitle("Distribution of "+string(Parameter)+" for "+string(Lith)+" at Location "+string(Location));
        ax=gca;
        XBeta=linspace(min(Data2Fit),max(Data2Fit),100);
        %----------------------------------------------------------------------
        %               PLOTTING THE EMPIRICAL DATA
        if PlotType=="pdf"
            hist = histogram(Data2Fit,'Normalization',PlotType);
            %%extract parameters
            counts = hist.Values;
            sum_counts = sum(counts);
            width = hist.BinWidth;
            %%area of the histogram
            area = sum_counts*width;
        elseif PlotType=="cdf"
            hist = cdfplot(Data2Fit);
            title('');
        end
        hold on;
        %----------------------------------------------------------------------
        %               PLOTTING THE BETA DISTRIBUTION FIT
        if max(Data2Fit)<=1
            YBeta = betapdf(XBeta,a,b);
            YBetacdf = betacdf(XBeta,a,b);
            if PlotType=="pdf"
                plot(XBeta,YBeta,'Color','r');
            elseif PlotType=="cdf"
                plot(XBeta,YBetacdf,'Color','r');
            end
        end
        %----------------------------------------------------------------------
        %                           COSMETICS
        xlim([min(Data2Fit)-std(Data2Fit) max(Data2Fit)+std(Data2Fit)]);
        xtickformat('%.2f');
        xlabel(Parameter);
        ylabel("Relative Frequency %");
    end
end
%--------------------------------------------------------------------------
% histogram(T.(Parameter)(T.Lithology=="W100"&T.Origin=="RIVER"&T.Location==14),'Normalization','pdf')
% hold on
% xplot=linspace(min(Data2Fit),max(Data2Fit),100);
% yplot=betapdf(xplot,a,b);
% plot(xplot,yplot);
% hold on
% histogram(T.(Parameter)(T.Lithology=="W100"&T.Origin=="RIVER"&T.Location==14))
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------