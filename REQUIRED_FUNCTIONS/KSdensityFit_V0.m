function [f,xi]=KSdensityFit(T,Parameter,alpha,Lith,Orig,Location,PlotFlag,PlotType)
% THIS FUNCTION RECEIVES INPUTS FROM A SINGLE LITHOLOGY, ORIGIN AND LOCATION
% AND PRODUCES FITS OF BETA DISTRIBUTIONS AND THEIR PLOT AGAINS DATA
% HISTOGRAM, IF REQUESTED.
%--------------------------------------------------------------------------
arguments
    T (:,:) table
    Parameter (1,1) string
    alpha (1,1) double
    Lith (:,:) string
    Orig (:,:) string
    Location (1,1) %{mustBeNumeric(LocPath)}
    PlotFlag (1,1) string
    PlotType (1,1) string
end
%--------------------------------------------------------------------------
%                       PARAMETERS POSSIBLE VALUES
%--------------------------------------------------------------------------
% PlotFlag: "Y" or "y" AND "N" OR "n"
% PlotType: "pdf" OR "cdf"
%--------------------------------------------------------------------------
%                       EXTRACT FROM TABLE T THE DATA TO FIT
%--------------------------------------------------------------------------
Data2Fit=T.(Parameter)(T.Lithology==Lith&T.Origin==Orig&T.Location==Location);
%   CHECK WHETHER DATA ARE SUITABLE FOR A FIT WITH A BETA DISTRIBUTION
if and(min(Data2Fit)>=0,max(Data2Fit)<=1)
    fprintf("WARNING: Input values are within the ranges for Beta distributions!")
    return
else
    %                       FIT KS DISTRIBUTION TO THE DATA
    if size(Data2Fit,1)>0
        [f,xi]=ksdensity(Data2Fit,'Support','positive');
        f=f';
        xi=xi';
        Ks_randomNumbers = ksdensity(Data2Fit, rand(size(Data2Fit,1),1), 'Support','positive','function', 'icdf');
        [hip,p]=kstest2(Data2Fit,Ks_randomNumbers,'Alpha',alpha);
        h=double(hip);
    end
    %--------------------------------------------------------------------------
    % PLOTTING
    %--------------------------------------------------------------------------
    if or(PlotFlag=="Y",PlotFlag=="y")
        f4 = figure;
        sgtitle("Distribution of "+string(Parameter)+" for "+string(Lith)+" at Location "+string(Location));
        ax=gca;
        Xplot=linspace(min(Data2Fit),max(Data2Fit),100);
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
        %               PLOTTING THE KERNEL SMOOTHING DISTRIBUTION FIT
        if PlotType=="pdf"
            Yplot = ksdensity(Xplot,'Support','positive','Function','pdf');
            plot(xi,f,'Color','r');
        elseif PlotType=="cdf"
            Yplot = ksdensity(Xplot,'Support','positive','Function','cdf');
            plot(Xplot,Yplot,'Color','r');
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