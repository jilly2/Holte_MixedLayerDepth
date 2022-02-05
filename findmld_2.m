% CODE FROM Holte et al. modified to handle a single profile salinity,
% temperature, pressure and potential density, and return MLD and
% MLD-averaged properties:
function [mld_this,mldT_this,mldS_this,mldD_this] = findmld_2(pres,temp,sal,pden)

% This subroutine calculates the mixed layer depth (MLD) of the input profile.  
% It is called by get_mld.m.

% The algorithm's parameters:          
errortol = 1*10^-10; % Error tolerance for fitting a straight line to the mixed layer -- unitless
range = 25;          % Maximum separation for searching for clusters of possible MLDs -- dbar
deltad = 100;        % Maximum separation of temperature and temperature gradient maxima for identifying
                     % intrusions at the base of the mixed layer -- dbar
tcutoffu = .5;       % Upper temperature cutoff, used to initially classify profiles as winter or summer profiles -- degrees C
tcutoffl = -.25;     % Lower temperature cutoff, used to initially classify profiles as winter or summer profiles -- degrees C
dcutoff = -.06;      % Density cutoff, used to initially classify profiles as winter or summer profiles -- kg/m^3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the MLD using a threshold method with de Boyer Montegut et al's
% criteria; a density difference of .03 kg/m^3 or a temperature difference
% of .2 degrees C.  The measurement closest to 10 dbar is used as the 
% reference value.  The threshold MLDs are interpolated to exactly match 
% the threshold criteria.

% Calculate the index of the reference value
m = length(sal); 
starti = min(find((pres-10).^2==min((pres-10).^2)));  
pres = pres(starti:m);
sal = sal(starti:m);
temp = temp(starti:m);
pden = pden(starti:m);
starti = 1;
m = length(sal); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the potential density anomaly, with a reference pressure of 0 
%pden = gsw_pden(sal,temp,pres,0)-1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Search for the first level that exceeds the potential density threshold
mldepthdens = m;
for j = starti:m
    if abs(pden(starti)-pden(j))>.03 
        mldepthdens = j;
        break;
    end
end

% Interpolate to exactly match the potential density threshold
clear pdenseg presseg presinterp pdenthreshold
presseg = [pres(mldepthdens-1) pres(mldepthdens)];
pdenseg = [pden(starti)-pden(mldepthdens-1) pden(starti) - pden(mldepthdens)];
P = polyfit(presseg,pdenseg,1);
presinterp = presseg(1):.5:presseg(2);
pdenthreshold = polyval(P,presinterp);

% The potential density threshold MLD value:
mldepthdens = presinterp(max(find(abs(pdenthreshold)<.03)));

% Search for the first level that exceeds the temperature threshold
mldepthptmp = m;
for j = starti:m
    if abs(temp(starti)-temp(j))>.2 
        mldepthptmp = j;
        break;
    end
end

% Interpolate to exactly match the temperature threshold
clear tempseg presseg presinterp tempthreshold
presseg = [pres(mldepthptmp-1) pres(mldepthptmp)];
tempseg = [temp(starti)-temp(mldepthptmp-1) temp(starti) - temp(mldepthptmp)];
P = polyfit(presseg,tempseg,1);
presinterp = presseg(1):.5:presseg(2);
tempthreshold = polyval(P,presinterp);

% The temperature threshold MLD value:
mldepthptmp = presinterp(max(find(abs(tempthreshold)<.2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the finite difference slope of the temperature, salinity and
% density profiles
clear tslope sslope dslope tslope_s sslope_s dslope_s ms
tslope = diff(temp)./diff(pres);
sslope = diff(sal)./diff(pres);
dslope = diff(pden)./diff(pres);
ms = length(tslope);

% smoothed the slope with a simple three point average using two 
% neighboring points
tslope_s = (tslope(1:ms-2) + tslope(2:ms-1) + tslope(3:ms))/3;
sslope_s = (sslope(1:ms-2) + sslope(2:ms-1) + sslope(3:ms))/3;
dslope_s = (dslope(1:ms-2) + dslope(2:ms-1) + dslope(3:ms))/3;
ms = length(tslope_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the MLD using a gradient method.  Following Dong et al., the gradient 
% criteria are .0005 kg/m^3/dbar and .005 degrees C/dbar.  If the criteria 
% are not met, the algorithm uses the temperature or density gradient extreme.
k = find( abs(dslope)>.0005 );
if any(k)
    gdmld = k(1) + 1;
else
    gdmld = min(find(abs(dslope)==max(abs(dslope))))+1;
end

l = find( abs(tslope)>.005  );
if any(l)
    gtmld = l(1) + 1;
else
    gtmld = min(find(abs(tslope)==max(abs(tslope))))+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Fit a straight line to the profile's mixed layer. Starting at the depth 
% closest to 10 dbar, use the first two points of the profile to calculate 
% a straight-line least-squares fit to the mixed layer.  Increase the depth
% and the number of points used in the fit until the bottom of the 
% profile. For each fit the error is calculated by summing the squared 
% difference between the fit and the profile over the depth of the fit.
% This step aims to accurately capture the slope of the mixed layer, and
% not its depth.  

clear errort errors errord
for j = starti+1:m
    % Fit line to temperature and calculate error
    P = polyfit(pres(starti:j),temp(starti:j),1);
    ltempfit = polyval(P,pres(starti:j));
    errort(j) = sum((temp(starti:j)-ltempfit).*(temp(starti:j)-ltempfit));

    % Fit line to salinity and calculate error
    P = polyfit(pres(starti:j),sal(starti:j),1);
    lsalfit = polyval(P,pres(starti:j));
    errors(j) = sum((sal(starti:j)-lsalfit).*(sal(starti:j)-lsalfit));

    % Fit line to potential density and calculate error
    P = polyfit(pres(starti:j),pden(starti:j),1);
    ldenfit = polyval(P,pres(starti:j));
    errord(j) = sum((pden(starti:j)-ldenfit).*(pden(starti:j)-ldenfit));

    clear ltempfit lsalfit ldenfit 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalize the errors
errort = errort/sum(errort); 
errors = errors/sum(errors);
errord = errord/sum(errord);

% Find deepest index with allowable error
upperlayert = max(find(errort < errortol));
upperlayers = max(find(errors < errortol));
upperlayerd = max(find(errord < errortol));
clear errort errors errord

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extend the mixed layer fit to the depth of the profile
P = polyfit(pres(starti:upperlayert),temp(starti:upperlayert),1);
ltempfit = polyval(P,pres(1:m));

P = polyfit(pres(starti:upperlayers),sal(starti:upperlayers),1);
lsalfit = polyval(P,pres(1:m));

P = polyfit(pres(starti:upperlayerd),pden(starti:upperlayerd),1);
ldenfit = polyval(P,pres(1:m));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit a straight line to the thermocline and extend the fit to the depth 
% of the profile.  The extreme value of each profile's smoothed gradient 
% (calculated in lines 82-84) is used to find the center of the 
% thermocline.
clear dtminfit dsmaxfit ddmaxfit

dtdzmax = max(find(abs(tslope_s) == max(abs(tslope_s))))+1;
P = polyfit(pres(dtdzmax-1:dtdzmax+1),temp(dtdzmax-1:dtdzmax+1),1);
dtminfit = polyval(P,pres(1:m));

dsdzmax = max(find(abs(sslope_s) == max(max(abs(sslope_s)))))+1;
P = polyfit(pres(dsdzmax-1:dsdzmax+1),sal(dsdzmax-1:dsdzmax+1),1);
dsmaxfit = polyval(P,pres(1:m));

dddzmax = max(find(abs(dslope_s) == max(max(abs(dslope_s)))))+1;
P = polyfit(pres(dddzmax-1:dddzmax+1),pden(dddzmax-1:dddzmax+1),1);
ddmaxfit = polyval(P,pres(1:m));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate one set of possible MLD values by finding the intersection 
% points of the mixed layer and thermocline fits.  If the fits do not
% intersect, the MLD value is set to 0.
upperdtmin = max(find(abs(dtminfit-ltempfit) == min(abs(dtminfit-ltempfit))));
if all(dtminfit-ltempfit>0)==1;
    upperdtmin = 0;
end
if all(-dtminfit+ltempfit>0)==1;
    upperdtmin = 0;
end

upperdsmax = max(find(abs(dsmaxfit-lsalfit) == min(abs(dsmaxfit-lsalfit))));
if all(-dsmaxfit+lsalfit>0)==1;
    upperdsmax = 0;
end
if all(dsmaxfit-lsalfit>0)==1;
    upperdsmax = 0;
end

upperddmax = max(find(abs(ddmaxfit-ldenfit) == min(abs(ddmaxfit-ldenfit))));
if all(ddmaxfit-ldenfit>0)==1;
    upperddmax = 0;
end
if all(-ddmaxfit+ldenfit>0)==1;
    upperddmax = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the remaining possible MLD values:  

% The maxima or minima of the temperature, salinity, and potential density
% profiles
tmax = max(find(temp == max(temp)));
smin = max(find(sal == min(sal)));
dmin = max(find(pden == min(pden)));

% The gradient MLD values
dtmax = gtmld;
dsmin = max(find(abs(sslope_s) == max(abs(sslope_s))))+1;
ddmin = gdmld;

% Sometimes subsurface temperature or salinity intrusions exist at the base 
% of the mixed layer.  For temperature, these intrusions are 
% characterized by subsurface temperature maxima located near temperature 
% gradient maxima. If the two maxima are separated by less than deltad, 
% the possible MLD value is recorded in dtandtmax.
dtmax2 = max(find(tslope_s == max(tslope_s)))+1;
if abs(pres(dtmax2)-pres(tmax)) < deltad;
    dtandtmax = min(dtmax2, tmax);
else
    dtandtmax = 0;
end
dsmin2 = max(find(sslope_s == min(sslope_s)))+1;
if abs(pres(dsmin2)-pres(smin)) < deltad;
    dsandsmin = min(dsmin2, smin);
else
    dsandsmin = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To determine if the profile resembles a typical winter or summer profile,
% the temperature change across the thermocline, tdiff, is calculated and
% compared to the temperature cutoff. tdiff is calculated as the 
% temperature change between the intersection of the mixed layer and thermocline fits and a 
% point two depth indexes deeper.  If upperdtmin is set to 0 or at the
% bottom of the profile, the points from the thermocline fit are used
% to evaluate tdiff.  
if  upperdtmin>0 && upperdtmin<(m-2) 
    tdiff = temp(upperdtmin)-temp(upperdtmin+2);
else
    tdiff = temp(dtdzmax-1)-temp(dtdzmax+1);
end
% tdiff is compared to the temperature cutoffs
if tdiff > tcutoffl && tdiff<tcutoffu
    testt = 1; % winter
else
    testt = 0; % summer
end

% For salinity and potential density profiles, the potential density 
% change across the pycnocline is calculated in a similar manner and 
% compared to a potential density cutoff.    
if upperddmax>0 && upperddmax<m-2
    ddiff = pden(upperddmax)-pden(upperddmax+2);
else
    ddiff = pden(dddzmax-1)-pden(dddzmax+1);
end
testd = testt;
if ddiff > dcutoff && tdiff > tcutoffu
    testd = 1; % winter
end
if ddiff > dcutoff && tdiff< tcutoffl
    testd = 0; % summer
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature Algorithm

% Convert the possible temperature MLDs from index to pressure
if upperdtmin > 0
    upperdtmin = pres(pres==upperdtmin);
end
tmax = pres(pres==tmax);
if dtandtmax>0
    dtandtmax = pres(pres==dtandtmax);
else
    dtandtmax = 0;
end
dtmax = pres(pres==dtmax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select the temperature MLD.  See the paper for a description of the
% steps.
mldindex=1;
if testt(mldindex) == 0 
    mixedt(mldindex) = upperdtmin(mldindex);
    analysis_t(mldindex) = 1;
    if tdiff(mldindex)<0 && mixedt(mldindex) > mldepthptmp(mldindex)
        mixedt(mldindex) = mldepthptmp(mldindex);
        analysis_t(mldindex) = 2;  
    end
    if mixedt(mldindex) > mldepthptmp(mldindex)
        if tmax(mldindex) < mldepthptmp(mldindex) && tmax(mldindex) > range 
            mixedt(mldindex) = tmax(mldindex);
            analysis_t(mldindex) = 3;
        else
            mixedt(mldindex) = mldepthptmp(mldindex);
            analysis_t(mldindex) = 4;
        end
    end       
else
    if abs(upperdtmin(mldindex)-mldepthptmp(mldindex)) < range && ...
       abs(dtandtmax(mldindex)-mldepthptmp(mldindex)) > range && ...
       upperdtmin(mldindex)<dtandtmax(mldindex)
        mixedt(mldindex) = upperdtmin(mldindex);
        analysis_t(mldindex) = 5;
    else
        if dtandtmax(mldindex) > pres(1)+range
           mixedt(mldindex) = dtandtmax(mldindex);
           analysis_t(mldindex) = 6; 
            a = [abs(dtmax(mldindex)-upperdtmin(mldindex)) ...
                 abs(dtmax(mldindex)-mldepthptmp(mldindex)) ...
                 abs(mldepthptmp(mldindex)-upperdtmin(mldindex))];
            if sum(a<range)>1
                mixedt(mldindex) = upperdtmin(mldindex);
                analysis_t(mldindex) = 7;                
            end
            if mixedt(mldindex)>mldepthptmp(mldindex)
                mixedt(mldindex) = mldepthptmp(mldindex);
                analysis_t(mldindex) = 8;
            end
        else
            if upperdtmin(mldindex)-mldepthptmp(mldindex) < range
                mixedt(mldindex) = upperdtmin(mldindex);
                analysis_t(mldindex) = 9;
            else
                mixedt(mldindex) = dtmax(mldindex);
                analysis_t(mldindex) = 10;
                if mixedt(mldindex) > mldepthptmp(mldindex)
                    mixedt(mldindex) = mldepthptmp(mldindex);
                    analysis_t(mldindex) = 11;
                end     
            end
        end
    end

    if mixedt(mldindex) == 0 && abs(mixedt(mldindex)-mldepthptmp(mldindex))>range
        mixedt(mldindex) = tmax(mldindex); 
        analysis_t(mldindex) = 12;
        if tmax(mldindex) == pres(1)
            mixedt(mldindex) = mldepthptmp(mldindex);
            analysis_t(mldindex) = 13;
        end
        if tmax(mldindex)>mldepthptmp(mldindex)
            mixedt(mldindex) = mldepthptmp(mldindex);
            analysis_t(mldindex) = 14;
        end        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Salinity Algorithm  

% Convert the possible salinity MLDs from index to pressure
if upperdsmax(mldindex)>0
    upperdsmax(mldindex) = pres(upperdsmax(mldindex));
end
dsmin(mldindex) = pres(dsmin(mldindex));
if dsandsmin(mldindex)>0
    dsandsmin(mldindex) = pres(dsandsmin(mldindex));
else
    dsandsmin(mldindex) = 0;
end

% Select the salinity MLD
if testd(mldindex) == 0
    mixeds(mldindex) = upperdsmax(mldindex);  
    analysis_s(mldindex) = 1;
    if mixeds(mldindex) - mldepthdens(mldindex) > range
        mixeds(mldindex) = mldepthdens(mldindex);
        analysis_s(mldindex) = 2;
    end
    if upperdsmax(mldindex)-dsmin(mldindex) < 0 && mldepthdens(mldindex)-dsmin(mldindex) > 0
        mixeds(mldindex) = dsmin(mldindex);
        analysis_s(mldindex) = 3;
    end
    if upperdsmax(mldindex)-dsandsmin(mldindex) < range && dsandsmin(mldindex) > range
        mixeds(mldindex) = dsandsmin(mldindex);
        analysis_s(mldindex) = 4;
    end
    if abs(mldepthdens(mldindex)-dsandsmin(mldindex)) < range && dsandsmin(mldindex) > range
        mixeds(mldindex) = dsandsmin(mldindex);
        analysis_s(mldindex) = 5;
    end  
    if mixedt(mldindex)-mldepthdens(mldindex)<0 && abs(mixedt(mldindex)-mldepthdens(mldindex))<range
        mixeds(mldindex) = mixedt(mldindex);  
        analysis_s(mldindex) = 6;
        if abs(mixedt(mldindex)-upperdsmax(mldindex))<range && upperdsmax(mldindex)-mldepthdens(mldindex)<0
            mixeds(mldindex) = upperdsmax(mldindex); 
            analysis_s(mldindex) = 7;
        end
    end
    if abs(mixedt(mldindex)-mldepthdens(mldindex))<abs(mixeds(mldindex)-mldepthdens(mldindex))
        if mixedt(mldindex)>mldepthdens(mldindex)
            mixeds(mldindex) = mldepthdens(mldindex);
            analysis_s(mldindex) = 8;
        end
    end
else
    if dsandsmin(mldindex) > range
        mixeds(mldindex) = dsandsmin(mldindex);
        analysis_s(mldindex) = 9;
        if mixeds(mldindex)>mldepthdens(mldindex)
            mixeds(mldindex) = mldepthdens(mldindex);
            analysis_s(mldindex) = 10;
        end
    else
        if dsmin(mldindex) < mldepthdens(mldindex)
            mixeds(mldindex) = dsmin(mldindex);
            analysis_s(mldindex) = 11;
            if upperdsmax(mldindex)<mixeds(mldindex) 
                mixeds(mldindex) = upperdsmax(mldindex);
                analysis_s(mldindex) = 12;
            end
        else
            mixeds(mldindex) = mldepthdens(mldindex);  
            analysis_s(mldindex) = 13;
            if upperdsmax(mldindex)<mixeds(mldindex)
                mixeds(mldindex) = upperdsmax(mldindex);
                analysis_s(mldindex) = 14;
            end
            if mixeds(mldindex) == 1 %%%%%%%%%%%%%%%%%should this be 0?
                mixeds(mldindex) = dsmin(mldindex);
                analysis_s(mldindex) = 15;
            end  
            if dsmin(mldindex) > mldepthdens(mldindex)
                mixeds(mldindex) = mldepthdens(mldindex);
                analysis_s(mldindex) = 16;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potential Density Algorithm.

% Convert the possible potential density MLDs from index to pressure
if upperddmax(mldindex)>0
    upperddmax(mldindex) = pres(upperddmax(mldindex)); 
end
dmin(mldindex) = pres(dmin(mldindex));
ddmin(mldindex) = pres(ddmin(mldindex));
 
% Select the potential density MLD
if testd(mldindex) == 0
    mixedd(mldindex) = upperddmax(mldindex);    
    analysis_d(mldindex) = 1;
    if mixedd(mldindex) > mldepthdens(mldindex)
        mixedd(mldindex) = mldepthdens(mldindex);
        analysis_d(mldindex) = 2;
    end  

    aa = [abs(mixeds(mldindex)-mixedt(mldindex)) abs(upperddmax(mldindex)-mixedt(mldindex)) abs(mixeds(mldindex)-upperddmax(mldindex))];
    if sum(aa<range)>1
        mixedd(mldindex) = upperddmax(mldindex);
        analysis_d(mldindex) = 3;                
    end
    if abs(mixeds(mldindex) - mldepthdens(mldindex)) < range && mixeds(mldindex)~=mldepthdens(mldindex)
        if mldepthdens(mldindex) < mixeds(mldindex)
            mixedd(mldindex) = mldepthdens(mldindex);
            analysis_d(mldindex) = 4;            
        else
            mixedd(mldindex) = mixeds(mldindex);
            analysis_d(mldindex) = 5;
        end
        if upperddmax(mldindex) == mldepthdens(mldindex)
            mixedd(mldindex) =  upperddmax(mldindex);
            analysis_d(mldindex) = 6;
        end
    end 
    if mixedd(mldindex)>ddmin(mldindex) && abs(ddmin(mldindex)-mixedt(mldindex))<abs(mixedd(mldindex)-mixedt(mldindex))
        mixedd(mldindex) = ddmin(mldindex);
        analysis_d(mldindex) = 7;
    end
else
    mixedd(mldindex) = mldepthdens(mldindex);
    analysis_d(mldindex) = 8;
    if mldepthptmp(mldindex)<mixedd(mldindex);
        mixedd(mldindex) = mldepthptmp(mldindex);
        analysis_d(mldindex) = 9;
    end
    if upperddmax(mldindex)<mldepthdens(mldindex) && upperddmax(mldindex)>range
        mixedd(mldindex) =  upperddmax(mldindex);
        analysis_d(mldindex) = 10;
    end
    if dtandtmax(mldindex) > range && dtandtmax(mldindex)<mldepthdens(mldindex)
        mixedd(mldindex) = dtandtmax(mldindex);
        analysis_d(mldindex) = 11;
        if abs(tmax(mldindex)-upperddmax(mldindex))<abs(dtandtmax(mldindex)-upperddmax(mldindex))
            mixedd(mldindex) = tmax(mldindex);
            analysis_d(mldindex) = 12;
        end          
        if abs(mixeds(mldindex) - mldepthdens(mldindex)) < range && mixeds(mldindex)<mldepthdens(mldindex)
            mixedd(mldindex) = min(mldepthdens(mldindex),mixeds(mldindex));
            analysis_d(mldindex) = 13;
        end 
    end
    if abs(mixedt(mldindex)-mixeds(mldindex)) < range
        if abs(min(mixedt(mldindex),mixeds(mldindex))-mixedd(mldindex)) > range
            mixedd(mldindex) = min(mixedt(mldindex),mixeds(mldindex));
            analysis_d(mldindex) = 14;
        end
    end
    if mixedd(mldindex)>ddmin(mldindex) && abs(ddmin(mldindex)-mixedt(mldindex))<abs(mixedd(mldindex)-mixedt(mldindex))
        mixedd(mldindex) = ddmin(mldindex);
        analysis_d(mldindex) = 15;
    end
    if upperddmax(mldindex)==upperdsmax(mldindex) && abs(upperdsmax(mldindex)-mldepthdens(mldindex))<range
        mixedd(mldindex) = upperddmax(mldindex);
        analysis_d(mldindex) = 16;
    end
    if mixedt(mldindex)==dmin(mldindex) 
        mixedd(mldindex) = dmin(mldindex);
        analysis_d(mldindex) = 17;
    end
end

% Return deepest mld and average T,S,density in this layer:
mld_this = max([mixedt,mixeds,mixedd,gtmld,gdmld])
ii=find(pres<=mld_this);
if isempty(ii)
    mldT_this = nan; mldS_this=nan;mldD_this=nan;
else
    mldT_this =nanmean(temp(ii));
    mldS_this = nanmean(sal(ii));
    mldD_this = nanmean(pden(ii));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output variables

% Algorithm mlds
mixedtp(mldindex) = mixedt(mldindex);
mixedsp(mldindex) = mixeds(mldindex);
mixeddp(mldindex) = mixedd(mldindex);

% Theshold method mlds
mldepthdensp(mldindex) = mldepthdens(mldindex);
mldepthptmpp(mldindex) = mldepthptmp(mldindex);

% Gradient method mlds
gtmldp(mldindex) = pres(gtmld(mldindex));
gdmldp(mldindex) = pres(gdmld(mldindex));

% Find the various methods' MLD indices for computing mixed layer average 
% properties
clear ta da tt dt

ta = find(pres<mixedt(mldindex));
da = find(pres<mixedd(mldindex));
tt = find(pres<mldepthptmp(mldindex));
dt = find(pres<mldepthdens(mldindex));
     
mixedt_ta(mldindex) = tsnanmean(temp(ta));
mixedd_ta(mldindex) = tsnanmean(temp(da));
mldepthdens_ta(mldindex) = tsnanmean(temp(dt));
mldepthptmp_ta(mldindex) = tsnanmean(temp(tt));

% Mixed layer average salinity over different MLDs
mixedt_sa(mldindex) = tsnanmean(sal(ta));
mixedd_sa(mldindex) = tsnanmean(sal(da));
mldepthdens_sa(mldindex) = tsnanmean(sal(dt));
mldepthptmp_sa(mldindex) = tsnanmean(sal(tt));

% Mixed layer average potential density over different MLDs
mixedt_da(mldindex) = tsnanmean(pden(ta));
mixedd_da(mldindex) = tsnanmean(pden(da));
mldepthdens_da(mldindex) = tsnanmean(pden(dt));
mldepthptmp_da(mldindex) = tsnanmean(pden(tt));

% Record which step selected the MLD for the temperature, salinity, and
% potential density profiles
tanalysis(mldindex) = analysis_t(mldindex);
sanalysis(mldindex) = analysis_s(mldindex);
danalysis(mldindex) = analysis_d(mldindex);

% Record which step selected the MLD for the temperature, salinity, and
% potential density profiles
tanalysis(mldindex) = analysis_t(mldindex);
sanalysis(mldindex) = analysis_s(mldindex);
danalysis(mldindex) = analysis_d(mldindex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the individual temperature, salinity, and potential density profiles, 
% as well as the mixed layer and thermocline fits and the various possible 
% MLD measures.  Turn this feature on in the 'EDIT' section of get_mld.m
yesplot=0;
if yesplot==1
        mintemp = (min(min(temp)));
        maxtemp = (max(max(temp)));
        minsal = (min(min(sal)));
        maxsal = (max(max(sal)));
        minpden = (min(min(pden)));
        maxpden = (max(max(pden)));  
        
        ta = (max(ta));
        da = (max(da));
        tt = (max(tt));
        dt = (max(dt));
        tmaxi = find(pres==tmax(mldindex));
        dtmaxi = find(pres==dtmax(mldindex));
        if upperdtmin(mldindex)>0
            upperdtmini = find(pres==upperdtmin(mldindex));
        else
            upperdtmini = length(pres);
        end
        
        
        figure
        subplot(1,3,1);
        plot(temp(:),pres(:),'ko','MarkerFaceColor','k','MarkerSize',6)
        hold on 
        tempchoice = 0:.5:length(pres);
        preschoice = mixedtp(mldindex)*ones(length(tempchoice));
        plot(tempchoice,preschoice,'k','LineWidth',4)
        %plot(ltempfit(:),pres(:),'k',dtminfit,pres,'--k','LineWidth',2)
        plot(temp(tt),mldepthptmp(mldindex),'s','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',16);
        %plot(temp(dt),mldepthdens(mldindex),'s','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',16);       
        plot(temp(gtmld(mldindex)),pres(gtmld(mldindex)),'>','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',16);              
        %plot(temp(gdmld(mldindex)),pres(gdmld(mldindex)),'>','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',16);                      
        plot(temp(tmaxi),tmax(mldindex),'o','MarkerFaceColor','c','MarkerEdgeColor','k','MarkerSize',12)
        plot(temp(dtmax2(mldindex)),pres(dtmax2(mldindex)),'s','MarkerFaceColor','c','MarkerEdgeColor','k','MarkerSize',12)
        plot(temp(dtmaxi),dtmax(mldindex),'o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',12)
        plot(temp(upperdtmini),upperdtmin(mldindex),'o','MarkerFaceColor',[1 .5 0],'MarkerEdgeColor','k','MarkerSize',12)
        set(gca, 'YDir','reverse','XAxisLocation','top')
        xlabel('Temperature (^oC)','FontSize',14)
        ylabel('Pressure (dbar)','FontSize',14')
        axis([mintemp maxtemp 0 max(pres)])
        grid on
        hold off
        
        
        if upperdsmax(mldindex)>0
            upperdsmaxi = find(pres==upperdsmax(mldindex));
        else
            upperdsmaxi = length(pres);
        end
        dsmini = find(pres==dsmin(mldindex));
        
        subplot(1,3,2);
        plot(sal(:),pres(:),'ko','MarkerFaceColor','k','MarkerSize',6)
        hold on
        tempchoice = 0:.5:length(pres);
        preschoice = mixedsp(mldindex)*ones(length(tempchoice));
        plot(tempchoice,preschoice,'k','LineWidth',4)
        %plot(lsalfit,pres,'k',dsmaxfit,pres,'--k','LineWidth',2)
        plot(sal(tt),mldepthptmp(mldindex),'s','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',16);
        %plot(sal(dt),mldepthdens(mldindex),'s','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',16);       
        plot(sal(gtmld(mldindex)),pres(gtmld(mldindex)),'>','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',16);              
        plot(sal(gdmld(mldindex)),pres(gdmld(mldindex)),'>','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',16);                 
        plot(sal(upperdsmaxi),upperdsmax(mldindex),'o','MarkerFaceColor',[1 .5 0],'MarkerEdgeColor','k','MarkerSize',12)
        plot(sal(smin(mldindex)),pres(smin(mldindex)),'o','MarkerFaceColor','c','MarkerEdgeColor','k','MarkerSize',12)
        plot(sal(dsmini),dsmin(mldindex),'o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',12);
        plot(temp(dsmin2(mldindex)),pres(dsmin2(mldindex)),'s','MarkerFaceColor','c','MarkerEdgeColor','k','MarkerSize',12)    
        set(gca, 'YDir','reverse','XAxisLocation','top')
        grid on 
        %titlename = [' profile ' int2str(mldindex) ' of float ' floatnumber];
        %title(titlename,'FontSize',18)
        xlabel('Salinity (PSU)','FontSize',14')
        axis([minsal maxsal 0 max(pres)])
        hold off

        if upperddmax(mldindex)>0
            upperddmaxi = find(pres==upperddmax(mldindex));
        else
            upperddmaxi = length(pres);
        end        
        dmini = find(pres==dmin(mldindex));
        ddmini = find(pres==ddmin(mldindex));
        
        subplot(1,3,3);
        plot(pden(:),pres(:),'ko','MarkerFaceColor','k','MarkerSize',6)
        hold on
        tempchoice = 0:.5:length(pres);
        preschoice = mixeddp(mldindex)*ones(length(tempchoice));
        plot(tempchoice,preschoice,'k','LineWidth',4)
        %plot(ldenfit,pres,'k',ddmaxfit,pres,'--k','LineWidth',2)
        plot(pden(tt),mldepthptmp(mldindex),'s','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',16);
        %plot(pden(dt),mldepthdens(mldindex),'s','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',16);       
        plot(pden(gtmld(mldindex)),pres(gtmld(mldindex)),'>','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',16);              
        plot(pden(gdmld(mldindex)),pres(gdmld(mldindex)),'>','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',16);               
        plot(pden(upperddmaxi),upperddmax(mldindex),'o','MarkerFaceColor',[1 .5 0],'MarkerEdgeColor','k','MarkerSize',12)
        plot(pden(dmini),dmin(mldindex),'o','MarkerFaceColor','c','MarkerEdgeColor','k','MarkerSize',12)
        plot(pden(ddmini),ddmin(mldindex),'o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',12);
        set(gca, 'YDir','reverse','XAxisLocation','top')
        grid on 
        xlabel(texlabel('Density (sigma_theta)'),'FontSize',14')
        axis([minpden maxpden 0 max(pres)])
        hold off
end
