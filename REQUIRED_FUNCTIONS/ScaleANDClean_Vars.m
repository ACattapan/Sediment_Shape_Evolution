function T=ScaleANDClean_Vars(T0,Detect_Thresh,Plottting_flag,Fig_number)
% 
%--------------------------------------------------------------------------
% arguments
%     T0 (:,:) table
%     Detect_Thresh (1,1) double
%     Plottting_flag (1,1) string
%     Fig_number (1,1) double
% end
%--------------------------------------------------------------------------
if isempty(Detect_Thresh)
    % Threshold for particles detection according to the intermediate size
    Detect_Thresh=7;%[mm]
end
if isempty(Fig_number)
    % Label for the figure used to plot the percentage of data lost 
    Fig_number=1;
end
%--------------------------------------------------------------------------
% N.B. Data in the imoprt spreadsheet are supposed to be in pixels. These data 
% will be transpoformed in the desired unit for all output based on the average 
% scaling expressed in the spreadsheet.
% Available options for output units are: 
%                   "px" = pixels
%                   "m"  = meters
%                   "cm" = centimeters
%                   "mm" = millimiters
%                   "inch" = inches
%                   "ft" = feet
Output_data_units="mm";
%--------------------------------------------------------------------------
T=T0;
% Convert particles dimensional properties from pixel units to the unit
% expressed in the ScaleAVG variable.
% Perimeter and Area are in [mm] and [mm2] respectively in our case.
for i=1:size(T,1)
    T.Perimeter(i)=T.ScaleAVG(i)*T.Perimeter(i);
    T.Area(i)=T.ScaleAVG(i)^2*T.Area(i);
    T.ConvPerim(i)=T.ScaleAVG(i)*T.ConvPerim(i);
    T.ConvArea(i)=T.ScaleAVG(i)^2*T.ConvArea(i);
    T.R(i)=T.ScaleAVG(i)*T.R(i);
    T.a(i)=T.ScaleAVG(i)*T.a(i);
    T.b(i)=T.ScaleAVG(i)*T.b(i);
end
%--------------------------------------------------------------------------
% Check the Percentage of data lost for different threshold values, up to
% 100 mm
PercLost=zeros(100,1);
for i=1:100 %i is a threshold in [mm]
    PercLost(i,1)=size(find(T.b<=i),1)/size(T,1)*100;
end
if or(Plottting_flag=="Y",Plottting_flag=="y")
    % Plot the percentage of data lost for different thresholds in [mm]
    figure (Fig_number);
    scatter(linspace(1,i),PercLost(:,1),[],[0 0.4470 0.7410]);
    hold on;
    scatter(Detect_Thresh,PercLost(Detect_Thresh,1),[],"r","filled");
    title('Percentage of sediment outlines lost for different threshold values of intermediate size');
    xlabel(sprintf('Threshold values of intermediate size ['+Output_data_units+']'));
    ylabel('Percentage of sediment outlines lost [%]');
end
% DELETE all the particles with an intermediate size less than
% Detect_Thresh [mm]
T(T.b<Detect_Thresh,:)=[];
%    N.B. CIRCULARITY HIGHER THAN 1 SEEM TO APPEAR IN PARTICLES WITH A
%    SHAPE VERY CLOSE TO A SPHERE AND THESE ARE VERY LIMITED
% DELETE all the particles with Circularity and NormCirc bigger thasn 1
T(T.Circularity>1,:)=[];
%--------------------------------------------------------------------------
h=(T.a-T.b).^2./(T.a+T.b).^2;
EllPerim=1/2*pi*(T.a+T.b).*(1+(3.*h./(10+(4-3.*h).^(1/2))));
EllArea=1/4*pi.*T.a.*T.b;
T.TheorCirc=4*pi.*EllArea./(EllPerim.^2);
T.NormCirc=T.Circularity./T.TheorCirc;
T = movevars(T,"TheorCirc",'After',"Circularity");
T = movevars(T,"NormCirc",'After',"TheorCirc");
T(T.NormCirc>1,:)=[];